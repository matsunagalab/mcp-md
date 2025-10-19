# MCP-MD アーキテクチャ・実装プラン

## 1. プロジェクト概要

### 目標

CHARMM-GUIに代わる、お手軽でフレクシブルな様々なMD入力システムを作成できるエージェント。

### 主要技術スタック

- **Boltz-2**: AlphaFold3超える精度の構造予測 + FEP法匹敵の親和性予測（1000倍高速）
- **AmberTools**: 配位子パラメータ化（GAFF2 + AM1-BCC、外部QMソフト不要）
- **OpenMM**: Pythonプログラマブルなプロダクション対応MD
- **LM Studio**: ローカルLLMによる自然言語ワークフロー生成
- **MCP (Model Context Protocol)**: モジュラーなツール統合フレームワーク

### 主要機能

- FASTAやSMILESから高精度な構造・親和性予測（Boltz-2）
- バーチャルスクリーニング（複数リガンド同時評価）
- 配位子の完全自動パラメータ化（SMILES → GAFF2）
- pH指定プロトネーション（PDB2PQR+PROPKA）
- ジスルフィド結合・金属サイト自動検出
- 膜タンパク質系構築（Packmol-Memgen）
- OpenMM MDスクリプト自動生成
- 複数MD形式エクスポート（Amber、GROMACS、OpenMM）

---

## 2. 全体アーキテクチャ

### MCPサーバー構成

```
[ユーザ]
    ↓ (自然言語/コマンド)
[LLMエージェント (未実装)]
    ├─ Planner (プロンプト→DAG化)
    ├─ Validator (QC/再試行)
    └─ WorkflowEngine (実行オーケストレーション)
        ↓ (MCP ツール呼び出し)
    ┌──────────────────────────┐
    │  MCPサーバー群 (実装済み)   │
    ├──────────────────────────┤
    │ 1. Structure Server       │ ← Boltz-2, PDBFixer, PDB2PQR
    │ 2. Ligand Server          │ ← RDKit, AmberTools
    │ 3. Docking Server         │ ← smina
    │ 4. Assembly Server        │ ← tleap, Packmol-Memgen
    │ 5. Protocol Server        │ ← OpenMM
    │ 6. Export Server          │ ← ParmEd
    └──────────────────────────┘
        ↓
    [MD入力ファイル (prmtop/inpcrd/script)]
```

### ディレクトリ構造

```
mcp-md/
├── servers/                  # MCPサーバー実装
│   ├── base_server.py       # ベースクラス
│   ├── structure_server.py  # Phase 1
│   ├── ligand_server.py     # Phase 2
│   ├── docking_server.py    # Phase 3
│   ├── assembly_server.py   # Phase 4
│   ├── protocol_server.py   # Phase 5
│   └── export_server.py     # Phase 6
├── core/                     # コアエンジン
│   ├── llm_client.py        # LM Studio統合
│   ├── planner.py           # プランニング (未完成)
│   ├── validator.py         # QC/バリデーション (未完成)
│   ├── workflow.py          # ワークフロー実行 (未完成)
│   ├── models.py            # データモデル
│   └── utils.py             # ユーティリティ
├── tools/                    # 外部ツールラッパー
│   ├── base_wrapper.py      # ベースクラス
│   ├── boltz2_wrapper.py
│   ├── pdbfixer_wrapper.py
│   ├── pdb2pqr_wrapper.py
│   ├── rdkit_wrapper.py
│   ├── ambertools_wrapper.py
│   ├── smina_wrapper.py
│   ├── packmol_wrapper.py
│   └── openmm_wrapper.py
├── tests/                    # テストコード
├── examples/                 # 使用例
├── docs/                     # ドキュメント
├── main.py                   # エントリーポイント
└── pyproject.toml           # プロジェクト設定
```

### データフロー例

**典型的なワークフロー**: PDB + SMILES → 複合体MD系

```
1. [Structure] fetch_pdb("1ABC") 
   → PDB構造取得

2. [Structure] clean_structure() 
   → 欠損原子補完、水除去

3. [Structure] protonate_structure(ph=7.4)
   → pH指定プロトネーション

4. [Ligand] smiles_to_3d("CC(=O)...")
   → SMILES → 3D構造

5. [Ligand] generate_gaff_params()
   → GAFF2パラメータ生成（AM1-BCC）

6. [Ligand] create_ligand_lib()
   → tleapライブラリ作成

7. [Docking] dock_ligand_smina()
   → リガンドドッキング（オプション）

8. [Assembly] build_complex_system()
   → タンパク質-リガンド複合体構築

9. [Assembly] solvate_box()
   → 溶媒化（TIP3P/OPC）

10. [Assembly] add_salt(0.15M)
    → イオン付与

11. [Protocol] create_openmm_workflow()
    → OpenMM MDスクリプト生成

12. [Export] package_system()
    → ZIP化（README・スクリプト含む）
```

---

## 3. 実装状況サマリー

### Phase別実装状況

| Phase | サーバー | 状態 | 主要機能 | ファイル |
|-------|---------|------|---------|---------|
| 1 | Structure | ✅ 完全実装 | Boltz-2予測、PDBFixer、PDB2PQR、修飾検出 | `servers/structure_server.py` |
| 2 | Ligand | ✅ 完全実装 | RDKit 3D生成、AmberTools GAFF2 (AM1-BCC) | `servers/ligand_server.py` |
| 3 | Docking | ✅ 完全実装 | smina分子ドッキング | `servers/docking_server.py` |
| 4 | Assembly | ✅ 完全実装 | tleap系構築、Packmol-Memgen膜系 | `servers/assembly_server.py` |
| 5 | Protocol | ✅ 完全実装 | OpenMM MDスクリプト生成 | `servers/protocol_server.py` |
| 6 | Export | ✅ 完全実装 | ParmEd形式変換、パッケージング | `servers/export_server.py` |
| - | Planner | ❌ 未実装 | LLMベースプラン生成 | `core/planner.py` (骨格のみ) |
| - | Validator | ❌ 未実装 | QC/バリデーション | `core/validator.py` (骨格のみ) |
| - | WorkflowEngine | ❌ 未実装 | 実行オーケストレーション | `core/workflow.py` (骨格のみ) |

### 実装済み機能

#### ✅ MCPサーバー（全6つ）
- 各サーバーは独立して起動・実行可能
- MCPプロトコルに準拠したツール登録
- エラーハンドリング・ロギング実装

#### ✅ 外部ツールラッパー（全9つ）
- Boltz-2、PDBFixer、PDB2PQR、RDKit、AmberTools、smina、Packmol、OpenMM、ParmEd
- 統一されたインターフェース（BaseToolWrapper）
- conda環境統合

#### ✅ コアユーティリティ
- LM Studioクライアント（OpenAI互換API）
- データモデル定義（Pydantic）
- 共通ユーティリティ関数

### 未実装機能

#### ❌ 統合ワークフローシステム
- **Planner**: 自然言語→DAG変換、LLMベースプラン生成
- **Validator**: 各ステップの出力検証、QCチェック、エラー再試行
- **WorkflowEngine**: DAG実行、MCPサーバー呼び出しオーケストレーション、進捗管理

#### ❌ CLIコマンド
- `mcp-md run workflow_plan.md`: ワークフロー実行
- `mcp-md chat`: LLMインタラクティブモード

---

## 4. Phase 1: Structure Server（実装済み）

### 目的

構造の取得・予測・前処理を一元化。Boltz-2による最先端の構造・親和性予測を統合。

### 実装ファイル

- `servers/structure_server.py` (494行) - MCPサーバー本体
- `tools/boltz2_wrapper.py` (210行) - Boltz-2ラッパー
- `tools/pdbfixer_wrapper.py` (127行) - PDBFixerラッパー
- `tools/pdb2pqr_wrapper.py` (129行) - PDB2PQR+PROPKAラッパー

### 主要機能

#### 1. 構造取得
- **`fetch_pdb`**: PDB、AlphaFold DB、PDB-REDOから構造取得
  - 対応ソース: `pdb`, `alphafold`, `pdb-redo`
  - 自動ダウンロード・検証

#### 2. Boltz-2構造予測

**Boltz-2統合**: [GitHub](https://github.com/jwohlwend/boltz)
- AlphaFold3超える精度
- FEP法匹敵の親和性予測（1000倍高速）
- 結合親和性: `affinity_probability_binary`（ヒット探索）、`affinity_pred_value`（リガンド最適化）

**ツール**:
- **`predict_structure_boltz2`**: FASTAから構造予測
  - MSA使用オプション
  - 複数モデル生成
  - 信頼度評価（pLDDT、pAE）
  
- **`predict_complex_with_affinity`**: FASTA + SMILES → 複合体 + 親和性
  - タンパク質-リガンド複合体同時予測
  - 結合親和性計算（IC50、バインダー確率）
  - 例: `IC50 = 0.63 μM`, `binder_prob = 0.85`
  
- **`screen_ligands_boltz2`**: バーチャルスクリーニング
  - 複数リガンド同時評価（100+化合物）
  - 親和性ランキング
  
- **`complete_missing_residues_boltz2`**: 欠損残基補完
  - 高精度補完（テンプレートベース）

#### 3. 構造クリーニング
- **`clean_structure`**: PDBFixerによるクリーニング
  - 欠損原子・残基補完
  - 異常altloc処理
  - 水分子除去

- **`add_hydrogens`**: 水素付与（pH指定）

#### 4. 高度なプロトネーション
- **`protonate_structure`**: PDB2PQR + PROPKA
  - pH指定プロトネーション（デフォルト: 7.4）
  - pKa計算統合
  - 力場指定（AMBER、CHARMM、PARSE）

#### 5. 修飾検出
- **`detect_modifications`**: 構造修飾の自動検出
  - ジスルフィド結合（SSBOND）
  - 修飾残基（MODRES）
  - 金属サイト（Zn、Mg、Ca、Fe、Cu、Mn、Na、K）

#### 6. 構造検証
- **`validate_structure`**: PDB構造の健全性チェック

### 技術仕様

#### Boltz-2 YAML入力

```yaml
- name: target_protein
  sequences:
    - protein:
        id: ["A"]
        sequence: "MKTAYIAKQRQISFVKSHFSRQ..."
- name: ligand
  smiles: "CC(=O)Oc1ccccc1C(=O)O"
```

#### PDB2PQR使用例

```bash
pdb2pqr --ff amber --with-ph 7.4 --titration-state-method propka \
  input.pdb output.pqr
```

### サポート力場

- タンパク質: ff19SB, ff14SB
- 核酸: OL15, OL3
- 脂質: lipid17, CHARMM36
- 糖鎖: GLYCAM06

---

## 5. Phase 2: Ligand Server（実装済み）

### 目的

配位子の3D構造生成からGAFF2パラメータ化まで完全自動化。**外部QMソフト不要**（AmberToolsのみで完結）。

### 実装ファイル

- `servers/ligand_server.py` (187行) - MCPサーバー本体
- `tools/rdkit_wrapper.py` (88行) - RDKitラッパー
- `tools/ambertools_wrapper.py` (223行) - AmberToolsラッパー

### 主要機能

#### 1. 3D構造生成
- **`smiles_to_3d`**: SMILES → 3D構造
  - RDKit ETKDG法
  - MMFF94コンフォーマ最適化

#### 2. AmberTools GAFF2パラメータ化

**電荷計算方法**（外部QMソフト不要）:
- **AM1-BCC**: antechamber内蔵（推奨、バランス型）
- **Gasteiger**: 高速・低精度
- **RESP**: sqmでHF/6-31G* ESP計算（高精度）

**ツール**:
- **`generate_gaff_params`**: antechamber + parmchk2
  - GAFF2力場パラメータ
  - AM1-BCC電荷計算
  - 欠損パラメータ自動補完

- **`create_ligand_lib`**: tleap用ライブラリ作成
  - .lib + .frcmod生成

- **`parameterize_ligand_complete`**: SMILES → tleapライブラリ（一括）
  - 3D生成 → パラメータ化 → ライブラリ作成（ワンコマンド）

### 技術仕様

#### AmberTools電荷計算

```bash
# AM1-BCC（推奨）
antechamber -i input.pdb -fi pdb -o output.mol2 -fo mol2 \
  -c bcc -at gaff2 -rn LIG -nc 0

# RESP（高精度）
antechamber -i input.pdb -fi pdb -o output.mol2 -fo mol2 \
  -c resp -at gaff2 -rn LIG -nc 0

# Gasteiger（高速）
antechamber -i input.pdb -fi pdb -o output.mol2 -fo mol2 \
  -c gas -at gaff2 -rn LIG -nc 0
```

#### パラメータ補完

```bash
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -a Y
```

---

## 6. Phase 3: Docking Server（実装済み）

### 目的

smina専用の分子ドッキング（AutoDock Vina不要）。

### 実装ファイル

- `servers/docking_server.py` (135行) - MCPサーバー本体
- `tools/smina_wrapper.py` (162行) - sminaラッパー

### 主要機能

- **`prepare_receptor`**: PDBQT変換（smina用）
- **`prepare_ligand`**: PDBQT変換
- **`dock_ligand_smina`**: sminaドッキング
  - 結合サイト定義（center + size）
  - 複数ポーズ生成
  - スコア関数: Vina, Vinardo等
- **`align_to_reference`**: 既知リガンド整列
- **`cluster_poses`**: ポーズクラスタリング

### sminaの利点

- Vina互換 + 拡張スコア関数
- 高速計算
- カスタムスコア関数対応

---

## 7. Phase 4: Assembly Server（実装済み）

### 目的

tleapスクリプト自動生成による系組立。溶媒化、イオン付与、膜系構築を統合。

### 実装ファイル

- `servers/assembly_server.py` (156行) - MCPサーバー本体
- `tools/ambertools_wrapper.py` - tleap統合
- `tools/packmol_wrapper.py` (144行) - Packmol-Memgenラッパー

### 主要機能

#### 1. tleap系構築
- **`build_system_tleap`**: 完全なMD系構築
  - タンパク質読み込み
  - 配位子統合
  - 力場適用（ff19SB、ff14SB等）
  - 溶媒化（TIP3P、OPC）
  - イオン付与（中性化 + 塩濃度指定）

#### 2. Packmol-Memgen膜系構築
- **`build_membrane_system`**: 脂質二重層構築
  - 脂質組成指定（POPC、POPE、CHOL等）
  - タンパク質の膜埋め込み
  - 自動配向最適化
  - ミセル系対応

### 技術仕様

#### tleapスクリプト例

```tcl
# 力場読み込み
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p

# 構造読み込み
protein = loadpdb protein.pdb
loadamberparams ligand.frcmod
loadoff ligand.lib
ligand = loadpdb ligand.pdb

# 複合体作成
complex = combine {protein ligand}

# 溶媒化（12Åパディング）
solvateBox complex TIP3PBOX 12.0

# イオン付与（0.15M NaCl）
addIonsRand complex Na+ 0
addIonsRand complex Cl- 0
addIonsRand complex Na+ 24
addIonsRand complex Cl- 24

# 出力
saveamberparm complex system.prmtop system.inpcrd
savepdb complex system.pdb

quit
```

#### Packmol-Memgen使用

```bash
packmol-memgen \
  --pdb protein.pdb \
  --lipids POPC:0.7:POPE:0.3 \
  --membrane-type bilayer \
  --dist 15.0 \
  --output membrane_system
```

---

## 8. Phase 5: Protocol Server（実装済み）

### 目的

OpenMM専用MDスクリプト生成（Amber/GROMACS/NAMD不要）。

### 実装ファイル

- `servers/protocol_server.py` (220行) - MCPサーバー本体
- `tools/openmm_wrapper.py` (125行) - OpenMMラッパー

### 主要機能

- **`generate_openmm_minimization`**: 最小化スクリプト
- **`generate_openmm_equilibration`**: NPT平衡化
- **`generate_openmm_production`**: プロダクションMD
- **`create_openmm_workflow`**: 統合ワークフロー（最小化→平衡化→プロダクション）
- **`add_position_restraints`**: 位置拘束設定

### OpenMMの利点

- Pythonプログラマブル
- GPU/CPU自動最適化
- tleap生成prmtop/inpcrd直接読込
- リアルタイム監視・制御

### 典型的プロトコル

1. Minimization (1000-5000 steps)
2. NVT Equilibration (100 ps, Langevin 300K)
3. NPT Equilibration (500 ps, MonteCarloBarostat)
4. NPT Production (10-100+ ns, DCD出力)

---

## 9. Phase 6: Export Server（実装済み）

### 目的

形式変換・解析準備・パッケージング。

### 実装ファイル

- `servers/export_server.py` (178行) - MCPサーバー本体
- ParmEd統合（依存関係）

### 主要機能

- **`export_amber`**: prmtop/inpcrd出力
- **`export_gromacs`**: ParmEd → top/gro変換
- **`export_openmm`**: XML形式
- **`export_charmm`**: PSF/CRD変換
- **`setup_analysis`**: cpptraj解析スクリプト雛形
- **`package_system`**: ZIP化（README・実行スクリプト含む）

### サポートMDエンジン

- **OpenMM**: フルサポート（推奨）
- **Amber**: prmtop/inpcrd出力
- **GROMACS**: ParmEdで変換
- **CHARMM**: PSF/CRD変換
- **NAMD**: PSF/PDB + 設定ファイル

---

## 10. 統合ワークフローエンジン（未実装）

### 目的

自然言語クエリから完全自動でMD入力を生成するエンドツーエンドシステム。

### アーキテクチャ

```
[ユーザ: 自然言語クエリ]
    ↓
[Planner: LLMベースプラン生成]
    ↓
[WorkflowEngine: DAG実行]
    ├─ 各ステップでMCPツール呼び出し
    ├─ [Validator: 出力検証]
    ├─ エラー時: 再試行/代替パス
    └─ 進捗管理・ログ
    ↓
[MD入力ファイル + レポート]
```

### Planner設計（未実装）

**ファイル**: `core/planner.py`（骨格のみ存在）

**機能**:
1. 自然言語クエリをパース
2. LLM（LM Studio）でワークフロー生成
3. DAG（Directed Acyclic Graph）構築
4. 各ノードにMCPツール割り当て
5. Markdown形式でプラン保存

**LM Studio統合**:
- OpenAI互換API（`http://localhost:1234/v1`）
- プライバシー保護（データがローカルに留まる）
- 推奨モデル: `gpt-oss-20b`（設定済み）

**実装優先度**: 高

**例**:
```python
planner = MDWorkflowPlanner()
plan = planner.plan_from_query(
    "PDB 1ABCからタンパク質を取得して、Aspirinをドッキング後、MD系を構築"
)
planner.save_plan(plan, "workflow_plan.md")
```

### Validator設計（未実装）

**ファイル**: `core/validator.py`（骨格のみ存在）

**機能**:
1. 各ステップの出力検証
2. QCチェック（原子数、エネルギー、構造健全性）
3. エラー検出・分類
4. 再試行ロジック
5. 品質レポート生成

**実装優先度**: 高

**例**:
```python
validator = MDWorkflowValidator()
result = await validator.validate_step(
    step_name="clean_structure",
    output={"pdb_file": "cleaned.pdb", "num_atoms": 1234}
)
if not result.is_valid:
    # 再試行または代替パス
```

### WorkflowEngine設計（未実装）

**ファイル**: `core/workflow.py`（骨格のみ存在）

**機能**:
1. DAGトポロジカルソート
2. 各ステップを順次実行
3. MCPサーバー呼び出しオーケストレーション
4. ステップ間データ受け渡し
5. 進捗管理・ロギング
6. エラーハンドリング・リカバリ

**実装優先度**: 高

**例**:
```python
engine = WorkflowEngine()
await engine.run_workflow("workflow_plan.md")
```

### 未実装CLIコマンド

```bash
# ワークフロー実行
mcp-md run workflow_plan.md

# インタラクティブモード
mcp-md chat "Boltz-2でFASTA配列から構造予測して"
```

---

## 11. 外部ツール依存関係

### conda経由（推奨）

```bash
conda create -n mcp-md python=3.11
conda activate mcp-md
conda install -c conda-forge ambertools packmol smina pdbfixer
```

### pip経由

```bash
pip install -e .  # 以下を自動インストール:
# boltz>=2.0.0
# pdb2pqr>=3.1.0
# propka>=3.5.0
# rdkit>=2023.9.1
# openmm>=8.3.1
# parmed>=4.3.0
# mdanalysis>=2.10.0
# openai>=1.0.0 (LM Studio用)
# pydantic>=2.12.3
# typer>=0.19.2
# rich>=14.2.0
# mcp>=1.18.0
```

### 環境構築

詳細は`README.md`の「インストール手順」を参照。

---

## 12. 開発ロードマップ

### 次のステップ（優先度順）

#### 1. Planner/Validator/WorkflowEngine実装（最優先）

**目標**: 自然言語からMD入力まで完全自動化

**タスク**:
- [ ] Planner: LM Studio統合、DAG生成、Markdown保存
- [ ] Validator: 出力検証ロジック、QCチェック
- [ ] WorkflowEngine: DAG実行、MCPサーバー呼び出し
- [ ] `mcp-md run`コマンド実装
- [ ] 統合テスト

**期間**: 2-3週間

#### 2. ドキュメント・テスト充実

**タスク**:
- [ ] 各MCPサーバーのユニットテスト
- [ ] エンドツーエンド統合テスト
- [ ] ユーザーガイド拡充
- [ ] APIリファレンス生成

**期間**: 1-2週間

#### 3. 追加機能

**低優先度**:
- [ ] GUI（Webベース）
- [ ] 他MDエンジン対応（GROMACS、NAMD等のネイティブサポート）
- [ ] クラウドデプロイ対応
- [ ] ワークフロー可視化

---

## 13. 使用例

### 例1: Boltz-2でFASTAから構造予測 → MD入力

```bash
conda activate mcp-md

# 1. Structure Server起動
python -m servers.structure_server &

# 2. Ligand Server起動
python -m servers.ligand_server &

# 3. Assembly Server起動
python -m servers.assembly_server &

# 4. Protocol Server起動
python -m servers.protocol_server &

# 5. ワークフロー実行（未実装）
# mcp-md run examples/boltz2_workflow.md
```

### 例2: PDB + SMILESから複合体MD系

詳細なPythonコード例は`examples/phase_124_workflow.md`を参照。

---

## 14. 参考資料

### 外部ツールドキュメント

- **Boltz-2**: https://github.com/jwohlwend/boltz
- **AmberTools**: https://ambermd.org/AmberTools.php
- **OpenMM**: https://openmm.org/
- **PDB2PQR**: https://pdb2pqr.readthedocs.io/
- **smina**: https://sourceforge.net/projects/smina/
- **Packmol-Memgen**: https://github.com/BIOS/packmol-memgen

### プロジェクトドキュメント

- **README.md**: クイックスタート・インストール
- **examples/phase_124_workflow.md**: 実践的なワークフロー例
- **ARCHITECTURE.md**: このファイル（マスタープラン）

### 論文引用

#### Boltz-2
```bibtex
@article{passaro2025boltz2,
  author = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and ...},
  title = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  year = {2025},
  journal = {bioRxiv}
}
```

#### AmberTools
```bibtex
@article{case2023ambertools,
  title={AmberTools},
  author={Case, D.A. and ...},
  journal={Journal of Chemical Information and Modeling},
  year={2023}
}
```

#### OpenMM
```bibtex
@article{eastman2017openmm,
  title={OpenMM 7: Rapid development of high performance algorithms for molecular dynamics},
  author={Eastman, Peter and ...},
  journal={PLOS Computational Biology},
  year={2017}
}
```

---

## 15. まとめ

### 実装完了部分

- ✅ 全6つのMCPサーバー（Structure, Ligand, Docking, Assembly, Protocol, Export）
- ✅ 全9つの外部ツールラッパー
- ✅ LM Studioクライアント
- ✅ データモデル・ユーティリティ

### 未完成部分

- ❌ Planner/Validator/WorkflowEngine統合
- ❌ `mcp-md run`コマンド
- ❌ エンドツーエンドテスト

### 現在の使用方法

各MCPサーバーは独立して起動可能。統合ワークフローシステムが未実装のため、現在は各サーバーを個別に呼び出す必要がある。

**次のマイルストーン**: Planner/Validator/WorkflowEngine実装により、完全自動化を達成。

