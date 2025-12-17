# MCP-MD: Molecular Dynamics Input File Generation Agent

Amber系に特化したMD入力ファイル生成AIエージェントシステム。LangGraph + FastMCPで構築された3フェーズワークフロー（Clarification → Setup → Validation）。

## 特徴

- **LangGraph統合**: ステートフルなワークフロー、永続化、人間フィードバック
  - LangChain 1.0準拠のStateGraphベースの実装
  - langchain-mcp-adaptersで公式MCP統合
  - チェックポイント機能で中断・再開可能
- **ReActパターン**: Phase 1でPDB構造を事前検査してから適切な質問を生成
  - `fetch_molecules`/`inspect_molecules`ツールで構造を分析
  - マルチチェーン構造やリガンドの有無を自動検出
  - シンプルな単一チェーンタンパク質は自動で処理進行
- **Boltz-2統合**: FASTAやSMILESから高精度な構造予測と結合親和性予測
- **AmberTools完結**: 配位子パラメータ化に外部QMソフト不要（AM1-BCC電荷計算）
- **FastMCP統合**: モジュラーな5つの独立サーバー、型安全な自動スキーマ生成
- **OpenMM専用**: Pythonプログラマブルなプロダクションレディなスクリプト生成

## 📚 ドキュメント

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - プロジェクト全体のアーキテクチャ・実装プラン・技術仕様
- **[CLAUDE.md](CLAUDE.md)** - Claude Code用ガイダンス・開発パターン
- **[AGENTS.md](AGENTS.md)** - Cursor AI Agent設定とガイドライン
- **[.cursor/rules/](.cursor/rules/)** - プロジェクトルールと開発ワークフロー

## インストール

### 前提条件

- Python 3.11以上
- [conda](https://docs.conda.io/en/latest/) または [mamba](https://mamba.readthedocs.io/)
- GPU推奨（Boltz-2、OpenMM高速化）

### 手順

#### 1. conda環境のセットアップ

```bash
# conda環境作成
conda create -n mcp-md python=3.11
conda activate mcp-md

# 科学計算パッケージをインストール
conda install -c conda-forge openmm rdkit mdanalysis biopython pandas numpy scipy openblas pdbfixer

# MD準備ツール
conda install -c conda-forge ambertools packmol smina
```

#### 2. Pythonパッケージのインストール

```bash
# プロジェクトのクローン
git clone https://github.com/matsunagalab/mcp-md.git
cd mcp-md

# パッケージをインストール（editable mode）
pip install -e .
```

#### 3. Boltz-2のインストール（オプション）

Boltz-2は Phase 2-3（Setup/Validation）で使用します。必要になったときにインストールしてください：

```bash
# CUDA対応GPUがある場合
pip install 'boltz[cuda]' --no-deps

# その後、不足している依存関係を個別にインストール
pip install torch hydra-core pytorch-lightning einops einx mashumaro modelcif wandb

# または、scipyをダウングレードしてから通常インストール
conda install -c conda-forge scipy=1.13.1
pip install 'boltz[cuda]'
```

> **注意**: Boltz-2の依存関係の一つ（fairscale）がscipy==1.13.1を厳密に要求するため、condaで既にインストールされているscipyと競合する場合があります。`--no-deps`オプションを使用することで、既存のパッケージを保持したまま、不足しているものだけを追加できます。

#### 4. Ollamaのインストール（オプション）
OllamaはLocal LLMのローカル実行環境です。デフォルトではOllamaの`gpt-oss:20b`モデルを使用します。

```bash
# Macの場合
brew install ollama
brew pull gpt-oss:20b
brew services start ollama
```

## 使用方法

### CLI (main.py)

```bash
# Interactive mode - エージェントと対話しながらセットアップ（推奨）
python main.py interactive
python main.py interactive "Setup MD for PDB 1AKE"

# Batch mode - 完全自動でワークフロー実行
python main.py batch "Setup MD for PDB 1AKE in explicit water, 1 ns at 300K"

# JSON出力付きバッチ処理
python main.py batch "Setup MD for 1AKE" --output-json results.json

# 中断したセッションを再開
python main.py resume --thread-id md_session_xxxxx

# Phase 1のみ（SimulationBrief生成）
python main.py clarify "Setup MD for PDB 1AKE"

# MCPサーバー一覧
python main.py list-servers

# ヘルプ
python main.py --help
python main.py info
```

### Notebook開発

```bash
jupyter notebook notebooks/md_agent_v2.ipynb
```

### MCPサーバーのテスト

各FastMCPサーバーを単独でテスト可能：

```bash
# MCP Inspector起動（Structure Serverを例に）
mcp dev servers/structure_server.py

# 別のサーバーをテストする場合
mcp dev servers/genesis_server.py
mcp dev servers/solvation_server.py
mcp dev servers/amber_server.py
mcp dev servers/md_simulation_server.py
```

### MCPサーバー一覧

| サーバー | 説明 |
|---------|------|
| `structure_server` | PDB/AlphaFold/PDB-REDOからの構造取得、チェーン分離、構造修復、リガンドGAFF2パラメータ化 |
| `genesis_server` | Boltz-2によるFASTA配列からの構造予測（単量体・多量体対応） |
| `solvation_server` | packmol-memgenによる溶媒和（水ボックス）・脂質膜埋め込み |
| `amber_server` | tleapによるAmberトポロジー（parm7）・座標（rst7）ファイル生成 |
| `md_simulation_server` | OpenMMによるMD実行、MDTrajによるトラジェクトリ解析 |

## ディレクトリ構造

```
mcp-md/
├── main.py               # CLI エントリポイント
│
├── src/mcp_md/           # ソースコード（直接編集）
│   ├── prompts.py                  # プロンプトテンプレート
│   ├── utils.py                    # ユーティリティ
│   ├── state_scope.py              # Phase 1状態定義
│   ├── state_setup.py              # Phase 2状態定義
│   ├── state_validation.py         # Phase 3状態定義
│   ├── state_full.py               # 統合状態定義
│   ├── clarification_agent.py      # Phase 1: ReAct Agent（構造検査→質問）
│   ├── setup_agent.py              # Phase 2: ReAct Setup Agent
│   ├── validation_agent.py         # Phase 3: 検証・レポート
│   ├── mcp_integration.py          # MCP統合
│   └── full_agent.py               # 3フェーズ統合
│
├── notebooks/            # テスト・デモ用
│   ├── 1_clarification.ipynb       # Phase 1 テスト
│   ├── md_agent_v2.ipynb           # 統合テスト
│   └── test_*.ipynb                # MCPサーバーテスト
│
├── servers/              # FastMCPサーバー（5サーバー）
│   ├── structure_server.py         # PDB取得・構造修復・リガンドGAFF2パラメータ化
│   ├── genesis_server.py           # Boltz-2構造生成（FASTA→PDB）
│   ├── solvation_server.py         # 溶媒和・膜埋め込み（packmol-memgen）
│   ├── amber_server.py             # Amberトポロジー・座標生成（tleap）
│   └── md_simulation_server.py     # MD実行・解析（OpenMM/MDTraj）
│
├── common/               # 共通ライブラリ
│   ├── base.py                     # BaseToolWrapper
│   └── utils.py                    # 共通ユーティリティ
│
├── checkpoints/          # LangGraphチェックポイント
├── ARCHITECTURE.md       # 詳細アーキテクチャ
├── CLAUDE.md             # Claude Code ガイダンス
├── AGENTS.md             # Cursor AI Agent設定
└── README.md             # このファイル
```

## 開発ワークフロー

### Direct Python Files

このプロジェクトは **Direct Python Files** パターンを採用しています：

```
✅ src/mcp_md/ を直接編集
✅ notebooks/ でテスト・デモ
✅ ruff check src/mcp_md/ でフォーマットチェック

🚫 %%writefile でのコード生成は非推奨
```

### コードフォーマット

```bash
# フォーマットチェック
ruff check src/mcp_md/

# 自動修正
ruff check src/mcp_md/ --fix
```

## ライセンス

MIT License

## 引用

このツールを使用する場合、以下を引用してください：

### Boltz-2

```
S. Passaro et al., Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.
bioRxiv (2025). doi:10.1101/2025.06.14.659707
```

### AmberTools

```
D. A. Case et al., AmberTools, J. Chem. Inf. Model. 63, 6183 (2023).
```

### OpenMM

```
P. Eastman et al., OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials,
J. Phys. Chem. B 128, 109 (2024).
```
