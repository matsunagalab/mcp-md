# MCP-MD アーキテクチャ・実装プラン

## 1. プロジェクト概要

### 目的とポジショニング

**Amber系に最適化したAIエージェント＋MCPツール群**

- **主軸**: Amber/GAFF/OpenFF/ParmEd/OpenMM エコシステムに特化
- **非競合**: CHARMM-GUIとは棲み分け（CHARMM系は変換経由で二次対応、将来拡張）
- **永続化**: MCP標準でツール接続を維持可能（将来のLLM/実行基盤の更新に強い）
- **ホスト/クライアント**: [LangChain](https://github.com/langchain-ai/langchain) + [LangGraph](https://github.com/langchain-ai/langgraph)に統一（MCPツール統合）

### 主要技術スタック

- **LangChain 1.0+**: LLM統合、ツール抽象化、プロンプト管理
  - LangChain 1.0では全てのchainsとagentsがLangGraph上に統一
  - `langchain-core`, `langchain-openai` (or `langchain-anthropic`)を使用
- **LangGraph**: ステートフルなグラフベースのワークフロー構築
  - チェックポイント機能で永続化とtime-travel可能
  - 複雑な制御フローと人間フィードバック統合をネイティブサポート
- **FastMCP**: MCPサーバーの実装とツール提供
- **Boltz-2**: 構造予測・複合体生成ツール
- **AmberTools**: 完全OSS、配位子パラメータ化（GAFF2 + AM1-BCC）
- **OpenMM**: Pythonプログラマブル、GPU最適化、プロダクション対応MD
- **MCP (Model Context Protocol)**: 標準化されたツール統合（ツールの永続性・相互運用性）

### 主要機能

1. **Hybrid Planning**: 固定スケルトン + 自律サブルーチン（意思決定ログ記録）
2. **品質保証**: 自作MolProbity等による物理化学的一貫性チェック
3. **再現性**: Plan/決定/生成物をJSON保存

---

## 2. 全体アーキテクチャ（ハイブリッド設計）

### システム構成

```
[Chat UI / Jupyter / CLI]
    ↓
[LangGraph Agent] ────────┐
  │                       │ (State Graph)
  ├─ StateGraph           │  - workflow_state (current step, params)
  │   ├─ Planner Node     │  - user_preferences (pH, 塩, box)
  │   ├─ Tool Nodes       │  - execution_history (過去の実行)
  │   └─ Decision Node    │  - decision_log (決定根拠)
  │                       │
  ├─ Checkpointer         │ (永続化)
  │   └─ SQLite/Postgres  │  - グラフステートの保存
  │                       │  - リトライ・巻き戻し可能
  │                       │
  └─ MCP Tools ───────────├─── [FastMCP Servers] 🆕
                          │
                          ├─ Structure Server
                          │   ├─ fetch_pdb
                          │   ├─ clean_structure
                          │   ├─ add_hydrogens
                          │   ├─ protonate_structure
                          │   ├─ detect_modifications
                          │   └─ validate_structure
                          │
                          ├─ Genesis Server 🆕
                          │   ├─ boltz2_protein_from_seq
                          │   ├─ boltz2_protein_from_fasta
                          │   └─ boltz2_multimer
                          │
                          ├─ Complex Server 🆕
                          │   ├─ boltz2_complex
                          │   ├─ boltz2_screen_ligands
                          │   ├─ smina_dock
                          │   └─ refine_poses
                          │
                          ├─ Ligand Server
                          │   ├─ smiles_to_3d
                          │   ├─ generate_gaff_params
                          │   ├─ create_ligand_lib
                          │   └─ parameterize_ligand_complete
                          │
                          ├─ Assembly Server
                          │   ├─ build_system_tleap
                          │   ├─ build_membrane_system
                          │   └─ build_mixed_solvent
                          │
                          ├─ Export Server
                          │   ├─ export_amber
                          │   ├─ export_gromacs
                          │   ├─ export_openmm
                          │   ├─ package_system
                          │   └─ convert_format
                          │
                          └─ QC/Min Server 🆕
                              ├─ openmm_minimize
                              ├─ clash_check
                              ├─ bond_check
                              ├─ chirality_check
                              ├─ run_full_qc
                              └─ posebusters_check

    ↓
[Persistent Storage]
  ├─ checkpoints/         (LangGraph state snapshots)
  │   └─ <thread_id>/     (会話セッション単位)
  └─ runs/<timestamp>/
      ├─ plan.json        (固定スケルトン + 決定ログ)
      ├─ outputs/         (PDB, prmtop, inpcrd, etc.)
      ├─ qc_report.json   (PoseBusters, 最小化指標)
      └─ metadata.json    (seed, hash, 再現用)
```

### FastMCP統合の特徴

- **モジュラー設計**: 各サーバーファイルが完全に独立して動作可能
- **自動スキーマ生成**: 型ヒントとdocstringから自動的にMCPツールスキーマを生成
- **標準準拠**: MCP標準プロトコルに完全準拠
- **開発効率**: デコレータベースのシンプルなAPI（`@mcp.tool`）
- **独立実行**: 各サーバーが `python -m servers.{server_name}` で単独起動可能
- **共通ライブラリ**: `common/` モジュールで外部ツール実行とユーティリティを共有
- **LangChain統合**: MCPツールをLangChain `Tool`オブジェクトにラップして利用

### ハイブリッド設計の核心

#### 1. 固定スケルトン（Plan）
```
fetch/generate → repair/protonate → ligand_param → 
complex → assemble → solvate/ions → export → 
minimize → package
```

**特徴**:
- 工程は固定（Amber特化パイプライン）
- 順序は決定論的
- 学術的再現性を担保

#### 2. 自律サブルーチン（LangGraph Node）
各工程内でツール・パラメータを動的選択：

```python
# 例: 複合体生成の意思決定（LangGraphノード内）
def complex_generation_node(state: WorkflowState):
    """複合体生成ノード"""
    
    if state.pdb_exists:
        tool_result = fetch_pdb_tool.invoke({"pdb_id": state.pdb_id})
    elif state.fasta_provided:
        tool_result = boltz2_protein_from_seq.invoke({
            "sequence": state.fasta_sequence
        })
        # 決定ログを状態に追加
        state.decision_log.append({
            "step": "protein_generation",
            "tool": "boltz2_protein_from_seq",
            "reason": "No PDB available, using Boltz-2"
        })
    
    # 複合体ポーズ生成
    if state.use_ai_model:
        poses = boltz2_complex.invoke({
            "protein": tool_result["output"],
            "ligand": state.ligand_smiles,
            "top_k": 5
        })
        state.decision_log.append({
            "step": "complex_generation",
            "tool": "boltz2_complex",
            "affinity": poses[0]["affinity"]
        })
        
        if state.refine_poses:
            poses = smina_dock.invoke({
                "receptor": poses[0]["protein"],
                "ligands": poses,
                "local_search": True
            })
            state.decision_log.append({
                "step": "pose_refinement",
                "tool": "smina_dock",
                "reason": "Improve local geometry"
            })
    
    return {"poses": poses, "decision_log": state.decision_log}
```

**特徴**:
- 失敗時は自動で1回再試行（LangGraphの条件付きエッジ）
- それでもNGなら人間にフィードバック（`interrupt_before/after`）
- すべての決定をStateに記録

---

## 3. MCPツール統合

### Genesis MCP: 構造生成

FASTA配列からPDB構造を生成：

```python
# FASTA → PDB
protein_pdb = boltz2_protein_from_seq(
    sequence="MKTAYIAKQRQISFVKSHFSRQ...",
    num_models=5
)
```

### Complex MCP: 複合体生成

タンパク質-配位子複合体の姿勢予測：

```python
# 受容体 + SMILES → 複合体候補
complexes = boltz2_complex(
    protein_pdb="receptor.pdb",
    ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
    top_k=10
)

# Sminaで局所精密化（オプション）
refined_poses = smina_dock(
    receptor="receptor.pdb",
    ligands=complexes[:5],
    local_search=True
)
```

### QC/Min MCP: 品質保証

物理化学的一貫性チェック：

```python
# PoseBustersチェック
qc_report = posebusters_check(pdb_file="complex.pdb")

# OpenMM最小化
minimized = openmm_minimize(
    prmtop="system.prmtop",
    inpcrd="system.inpcrd",
    max_iterations=5000
)
```

---

## 4. 永続エージェント（LangGraph × MCP）

### LangChain 1.0とLangGraphの関係

- **LangChain 1.0の変更**: 従来の`chains`と`agents`を廃止、全てLangGraph上に統一
- **推奨アプローチ**: 
  - シンプルなReActエージェント → `create_react_agent()` (高レベル抽象化)
  - 複雑なワークフロー → LangGraphのStateGraphを直接使用（推奨）
- **本プロジェクトの選択**: 固定スケルトン + 条件分岐のため、StateGraphを直接使用

### LangGraphの特徴

- **公式**: https://github.com/langchain-ai/langgraph
- **特徴**: 
  - ステートフルグラフベースのワークフロー
  - チェックポイント機能（永続化、time-travel、分岐実行）
  - 条件分岐とサイクル（ループ）のネイティブサポート
  - 人間フィードバック統合（`interrupt_before/after`）
  - サブグラフによるモジュラー設計
- **LangChain統合**: LangChain `Tool`をノード内で直接利用可能
- **MCP統合**: MCPサーバーをLangChain `Tool`にラップして利用

### 運用のキーポイント

#### 1. StateGraph定義
```python
from typing import TypedDict, Annotated, Sequence
from langgraph.graph import StateGraph, END
from langgraph.checkpoint.sqlite import SqliteSaver
from langchain_core.tools import Tool

# ワークフローステート定義
class WorkflowState(TypedDict):
    # 入力パラメータ
    query: str
    pdb_id: str | None
    ligand_smiles: str | None
    
    # 実行状態
    current_step: str
    outputs: dict
    
    # ユーザー設定
    user_preferences: dict  # {ph: 7.4, salt: 0.15, ...}
    
    # 決定ログ
    decision_log: Annotated[Sequence[dict], "append"]
    
    # エラーハンドリング
    retry_count: int
    error: str | None

# MCPツールをLangChain Toolにラップ
def create_mcp_tool(server_module: str, tool_name: str) -> Tool:
    """MCPサーバーからツールを作成
    
    Args:
        server_module: サーバーモジュール名 (e.g., "servers.structure_server")
        tool_name: ツール名 (e.g., "fetch_pdb")
    
    Returns:
        LangChain Tool object
    """
    async def run_mcp_tool(**kwargs):
        # MCPサーバーをStdioTransportで起動してツール呼び出し
        import subprocess
        import json
        
        cmd = ["python", "-m", server_module]
        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        # MCP JSONRPCリクエスト
        request = {
            "jsonrpc": "2.0",
            "method": f"tools/{tool_name}",
            "params": kwargs,
            "id": 1
        }
        
        stdout, _ = proc.communicate(json.dumps(request).encode())
        response = json.loads(stdout.decode())
        return response.get("result")
    
    return Tool(
        name=tool_name,
        description=f"MCP tool: {tool_name}",
        coroutine=run_mcp_tool
    )

# グラフ構築
graph = StateGraph(WorkflowState)
```

#### 2. ノード定義とグラフ構築
```python
# MCPツールの作成
fetch_pdb_tool = create_mcp_tool("servers.structure_server", "fetch_pdb")
clean_structure_tool = create_mcp_tool("servers.structure_server", "clean_structure")
boltz2_complex_tool = create_mcp_tool("servers.complex_server", "boltz2_complex")
# ... 他のツール

# ノード定義
def planner_node(state: WorkflowState):
    """プランニングノード: 固定スケルトンを準備"""
    return {
        "current_step": "fetch",
        "user_preferences": state.get("user_preferences", {
            "ph": 7.4,
            "salt_concentration": 0.15,
            "water_model": "TIP3P",
            "force_field": "ff19SB"
        })
    }

async def structure_fetch_node(state: WorkflowState):
    """構造取得ノード"""
    result = await fetch_pdb_tool.ainvoke({"pdb_id": state["pdb_id"]})
    
    return {
        "outputs": {**state["outputs"], "structure": result},
        "current_step": "repair",
        "decision_log": [{
            "step": "fetch",
            "tool": "fetch_pdb",
            "params": {"pdb_id": state["pdb_id"]}
        }]
    }

# ... 他のノード定義

# グラフ構築
graph.add_node("planner", planner_node)
graph.add_node("fetch", structure_fetch_node)
graph.add_node("repair", repair_node)
graph.add_node("ligand_param", ligand_param_node)
graph.add_node("complex", complex_node)
graph.add_node("assemble", assemble_node)
graph.add_node("qc", qc_node)

# エッジ定義（固定スケルトン）
graph.set_entry_point("planner")
graph.add_edge("planner", "fetch")
graph.add_edge("fetch", "repair")
graph.add_edge("repair", "ligand_param")
graph.add_edge("ligand_param", "complex")
graph.add_edge("complex", "assemble")
graph.add_edge("assemble", "qc")
graph.add_edge("qc", END)

# 条件付きエッジ（エラーハンドリング）
def should_retry(state: WorkflowState):
    if state.get("error") and state["retry_count"] < 1:
        return "retry"
    elif state.get("error"):
        return "human_feedback"
    return "continue"

graph.add_conditional_edges(
    "complex",
    should_retry,
    {
        "retry": "complex",
        "human_feedback": END,
        "continue": "assemble"
    }
)
```

#### 3. チェックポイント機能（永続化）
```python
from langgraph.checkpoint.sqlite import SqliteSaver

# チェックポイント保存用（永続化）
memory = SqliteSaver.from_conn_string("checkpoints/workflow.db")

# グラフのコンパイル
app = graph.compile(checkpointer=memory)

# 実行（スレッドIDで会話セッションを管理）
config = {"configurable": {"thread_id": "session_123"}}
result = app.invoke(
    {
        "query": "Fetch PDB 1ABC and dock Aspirin",
        "pdb_id": "1ABC",
        "ligand_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "user_preferences": {"ph": 7.4, "salt_concentration": 0.15},
        "decision_log": [],
        "retry_count": 0
    },
    config=config
)

# 中断からの再開
# グラフは自動的に前回の状態から再開される
result = app.invoke(None, config=config)

# 過去の実行履歴取得
for state in app.get_state_history(config):
    print(f"Step: {state.values['current_step']}")
    print(f"Decisions: {state.values['decision_log']}")
```

#### 4. 人間フィードバック統合
```python
# 人間の確認が必要なノードを指定
app = graph.compile(
    checkpointer=memory,
    interrupt_before=["complex"]  # 複合体生成前に確認
)

# 実行（複合体生成前で停止）
result = app.invoke(initial_state, config=config)

# ユーザーに確認を求める
print(f"Current state: {result}")
user_approval = input("Continue with complex generation? (y/n): ")

if user_approval == "y":
    # 実行再開
    result = app.invoke(None, config=config)
```

#### 5. 安全性・認証
- **MCP統合**: MCPサーバーは標準プロトコルで動作（認証は別管理）
- **実装方針**:
  - 環境変数で機密値（API key等）を管理
  - 一時クレデンシャル運用
  - ローカル実行前提（外部API最小化）

---

## 5. ユーザ体験（UX）設計

### 典型的な対話フロー

#### 1. 最初の返答
```
User: "PDB 1ABCにAspirinをドッキングしてMD系を作って"

Agent:
了解しました。以下のプランで進めます：

【固定スケルトン】
1. PDB取得（1ABC）
2. 構造修復・プロトネーション
3. リガンド処理（Aspirin, GAFF2/AM1-BCC）
4. 複合体生成（Boltz-2 + Smina）
5. 系構築（tleap）
6. 溶媒化・イオン付与
7. Amber出力（prmtop/inpcrd）
8. 最小化・QC
9. パッケージング

【質問: 以下を教えてください】
- pH: [default: 7.4]
- 塩濃度: [default: 0.15 M]
- Box padding: [default: 12.0 Å]
- 出力形式: [Amber / GROMACS / OpenMM]
- 既知結合部位: [あれば指定、なければBoltz-2で推定]
```

#### 2. 実行中の可視化
```
[Step 1/9] PDB取得 ✅ (1ABC.pdb, 1234 atoms)
[Step 2/9] 構造修復 ✅ (欠損残基 3箇所補完)
[Step 3/9] リガンド処理 ⏳
  - SMILES → 3D (RDKit) ✅
  - GAFF2パラメータ化 (AM1-BCC) ⏳
    Decision: AM1-BCC選択（バランス重視、計算時間 < 1min）

[中間プレビュー]
🔗 http://localhost:8080/view/intermediate.pdb (NGLビューワ)
```

#### 3. 失敗時の誘導
```
[Step 4/9] 複合体生成 ⚠️ エラー
  - Boltz-2予測: 親和性が極めて低い（binder_prob = 0.12）
  
【自動リトライ】
  - Smina局所サーチで代替候補を探索中...
  - 結果: 候補なし

【提案: 以下から選択してください】
a) 結合サイトを手動指定（推奨残基: SER195, HIS57, ASP102）
b) リガンドの3Dコンフォーマを変更（ETKDG → UFF）
c) Boltz-2の設定変更（MSA使用、top_k=10）
```

---

## 6. ロードマップ（Amber特化 → 拡張）

### Phase 1: MVP（最小実装）🎯

**目標**: Amber特化の基本パイプライン

**機能**:
- PDB/FASTA → Boltz-2 or fetch → 修復 → GAFF2/AM1-BCC
- Boltz-2複合体 ± Smina → tleap → Amber出力（prmtop/inpcrd）
- 基本的なQC（最小化、電荷整合）

**MCPサーバー**:
- Structure MCP
- Genesis MCP
- Complex MCP
- Ligand MCP
- Assembly MCP
- Export MCP（Amberのみ）

**期間**: 4-6週間

### Phase 2: QC強化

**目標**: 学術発表レベルの品質保証

**追加機能**:
- 自作MolProbityによる物理化学的一貫性チェック
- QCレポート生成（JSON + Markdown）

**期間**: 2-3週間

### Phase 3: 出力拡張

**目標**: GROMACS/OpenMM対応、膜系

**追加機能**:
- GROMACS出力（ParmEd）
- OpenMM XML出力
- 膜系構築（Packmol-Memgen）
- 混合溶媒対応（Packmol）

**期間**: 3-4週間

### Phase 4: HPC/永続運用

**目標**: 大規模スクリーニング、長期運用

**追加機能**:
- Strandsエージェントのキュー実行
- 結果キャッシュ（同一入力の再利用）
- 結果索引（検索可能なDB）
- HPC連携（Slurmジョブ投入）

**期間**: 4-6週間

### Phase 5: 将来拡張（低優先度）

**CHARMM系対応**:
- CHARMMパラメータ化（CGenFF）
- CHARMM-GUI出力への変換

**特殊系対応**:
- 糖鎖（GLYCAM）
- 金属中心（MCPB.py）
- RNA特化（OL3）

**プラグインアーキテクチャ**:
- 外部開発者がMCPサーバーを追加可能
- コミュニティ貢献の促進

---

## 7. 現在の実装状況（FastMCP統合完了）

### 実装済み（7 FastMCP Servers）✅

| Component | Status | 主要機能 |
|-----------|--------|---------|
| Structure Server | ✅ | PDB取得、PDBFixer、PDB2PQR、構造検証 |
| Genesis Server | ✅ 🆕 | Boltz-2タンパク質生成（FASTA→PDB、マルチマー） |
| Complex Server | ✅ 🆕 | Boltz-2複合体予測、Sminaドッキング、ポーズ精密化 |
| Ligand Server | ✅ | RDKit 3D生成、AmberTools GAFF2/AM1-BCC |
| Assembly Server | ✅ | tleap系構築、Packmol-Memgen膜系 |
| Export Server | ✅ | Amber/GROMACS/OpenMM形式変換、パッケージング |
| QC/Min Server | ✅ 🆕 | OpenMM最小化、衝突検出、結合長・キラリティチェック |

### FastMCP統合アーキテクチャ

#### 1. ディレクトリ構造
```
mcp-md/
├── common/              # 共通ライブラリ
│   ├── __init__.py
│   ├── base.py         # BaseToolWrapper
│   └── utils.py        # 共通ユーティリティ
├── servers/            # FastMCPサーバー（7ファイル）
│   ├── __init__.py
│   ├── structure_server.py
│   ├── genesis_server.py
│   ├── complex_server.py
│   ├── ligand_server.py
│   ├── assembly_server.py
│   ├── export_server.py
│   └── qc_min_server.py
├── core/               # LangGraphエージェント実装
│   ├── __init__.py
│   ├── workflow_graph.py      # StateGraph定義
│   ├── workflow_nodes.py      # ノード実装
│   ├── workflow_state.py      # WorkflowState定義
│   ├── mcp_integration.py     # MCPツール統合
│   ├── decision_logger.py     # 決定ログユーティリティ
│   └── models.py             # データモデル
├── checkpoints/        # LangGraphチェックポイント
│   └── workflow.db     # SQLiteストレージ
├── examples/
│   ├── phase_124_workflow.md
│   └── langgraph_example.py  # 実行例
└── pyproject.toml      # fastmcp, langchain, langgraph依存
```

#### 2. FastMCP統合の実装パターン

各サーバーは以下の標準パターンで実装：

```python
from fastmcp import FastMCP
from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)
mcp = FastMCP("Server Name")

# 外部ツールラッパー初期化
tool_wrapper = BaseToolWrapper("tool_name", conda_env="mcp-md")

@mcp.tool
def tool_name(param1: str, param2: int = 0) -> dict:
    """Tool description
    
    Args:
        param1: Parameter description
        param2: Optional parameter
    
    Returns:
        Result dictionary
    """
    # 実装コード
    return result

if __name__ == "__main__":
    mcp.run()  # STDIO transport (default)
```

#### 3. LangGraph × MCP統合パターン

```python
# core/mcp_integration.py
from langchain_core.tools import Tool
import subprocess
import json

def create_mcp_tool(server_module: str, tool_name: str, description: str = "") -> Tool:
    """MCPサーバーをLangChain Toolにラップ
    
    Args:
        server_module: サーバーモジュール名 (e.g., "servers.structure_server")
        tool_name: ツール名 (e.g., "fetch_pdb")
        description: ツールの説明文
    
    Returns:
        LangChain Tool object
    """
    async def run_mcp_tool(**kwargs):
        cmd = ["python", "-m", server_module]
        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # MCP JSONRPCリクエスト
        request = {
            "jsonrpc": "2.0",
            "method": f"tools/{tool_name}",
            "params": kwargs,
            "id": 1
        }
        
        stdout, stderr = proc.communicate(json.dumps(request))
        
        if proc.returncode != 0:
            raise RuntimeError(f"MCP tool error: {stderr}")
        
        response = json.loads(stdout)
        return response.get("result")
    
    return Tool(
        name=tool_name,
        description=description or f"MCP tool: {tool_name}",
        coroutine=run_mcp_tool
    )

def load_all_mcp_tools() -> dict[str, Tool]:
    """全MCPツールを読み込み"""
    tools = {}
    
    # Structure Server
    tools["fetch_pdb"] = create_mcp_tool(
        "servers.structure_server",
        "fetch_pdb",
        "Fetch PDB structure from RCSB PDB"
    )
    tools["clean_structure"] = create_mcp_tool(
        "servers.structure_server",
        "clean_structure",
        "Clean and repair PDB structure using PDBFixer"
    )
    
    # Genesis Server
    tools["boltz2_protein_from_seq"] = create_mcp_tool(
        "servers.genesis_server",
        "boltz2_protein_from_seq",
        "Generate protein structure from FASTA sequence using Boltz-2"
    )
    
    # Complex Server
    tools["boltz2_complex"] = create_mcp_tool(
        "servers.complex_server",
        "boltz2_complex",
        "Generate protein-ligand complex using Boltz-2"
    )
    
    # ... 他のツール定義
    
    return tools

# core/workflow_graph.py
from langgraph.graph import StateGraph, END
from langgraph.checkpoint.sqlite import SqliteSaver
from .workflow_state import WorkflowState
from .workflow_nodes import (
    planner_node,
    create_structure_fetch_node,
    create_repair_node,
    # ...
)
from .mcp_integration import load_all_mcp_tools

def create_workflow_graph():
    """ワークフローグラフを構築"""
    # MCPツール読み込み
    mcp_tools = load_all_mcp_tools()
    
    # グラフ構築
    graph = StateGraph(WorkflowState)
    
    # ノード追加（ツールを渡す）
    graph.add_node("planner", planner_node)
    graph.add_node("fetch", create_structure_fetch_node(mcp_tools))
    graph.add_node("repair", create_repair_node(mcp_tools))
    # ...
    
    # エッジ定義
    graph.set_entry_point("planner")
    graph.add_edge("planner", "fetch")
    # ...
    graph.add_edge("qc", END)
    
    # チェックポイント設定
    memory = SqliteSaver.from_conn_string("checkpoints/workflow.db")
    
    return graph.compile(checkpointer=memory)
```

### 削除されたファイル（FastMCPに置き換え）

- ~~`tools/`~~ ディレクトリ全体（10ファイル） → `common/`に統合
- ~~`servers/base_server.py`~~ → FastMCP標準機能で代替
- ~~`servers/archive/`~~ → 旧実装削除

---

## 8. 技術詳細（Phase別）

以下は既存実装の技術詳細です。新設計への移行時に参照してください。

### Phase 1: Structure Server

**実装ファイル**:
- `servers/structure_server.py` (494行)
- `tools/boltz2_wrapper.py` (210行)
- `tools/pdbfixer_wrapper.py` (127行)
- `tools/pdb2pqr_wrapper.py` (129行)

**主要ツール**:
1. `fetch_pdb`: PDB/AlphaFold/PDB-REDO取得
2. `predict_structure_boltz2`: FASTA→構造（→ Genesis MCPに移行）
3. `predict_complex_with_affinity`: FASTA+SMILES→複合体（→ Complex MCPに移行）
4. `clean_structure`: PDBFixerクリーニング
5. `protonate_structure`: PDB2PQR+PROPKA
6. `detect_modifications`: ジスルフィド・金属検出

**リファクタリング**:
- Boltz-2関連ツールを Genesis/Complex MCP に分離
- Structure MCP は構造取得・修復・プロトネーションに専念

### Phase 2: Ligand Server

**実装ファイル**:
- `servers/ligand_server.py` (187行)
- `tools/rdkit_wrapper.py` (88行)
- `tools/ambertools_wrapper.py` (223行)

**主要ツール**:
1. `smiles_to_3d`: SMILES → 3D（RDKit ETKDG）
2. `generate_gaff_params`: antechamber + parmchk2（GAFF2/AM1-BCC）
3. `create_ligand_lib`: tleap用ライブラリ

**新設計での位置づけ**:
- そのまま維持（Amber特化の核心部分）
- OpenFF統合は Phase 3以降

### Phase 3: Docking Server

**実装ファイル**:
- `servers/docking_server.py` (135行)
- `tools/smina_wrapper.py` (162行)

**主要ツール**:
1. `prepare_receptor/ligand`: PDBQT変換
2. `dock_ligand_smina`: Sminaドッキング
3. `align_to_reference`: 既知リガンド整列

**新設計での位置づけ**:
- Complex MCP に統合
- Boltz-2複合体予測の補助ツールとして位置づけ

### Phase 4: Assembly Server

**実装ファイル**:
- `servers/assembly_server.py` (156行)
- `tools/ambertools_wrapper.py` - tleap統合
- `tools/packmol_wrapper.py` (144行)

**主要ツール**:
1. `build_system_tleap`: 完全MD系構築
2. `build_membrane_system`: Packmol-Memgen膜系

**新設計での位置づけ**:
- そのまま維持（Amber特化の核心）

### Phase 5: Protocol Server

**実装ファイル**:
- `servers/protocol_server.py` (220行)
- `tools/openmm_wrapper.py` (125行)

**主要ツール**:
1. `generate_openmm_minimization`: 最小化
2. `generate_openmm_equilibration`: 平衡化
3. `generate_openmm_production`: プロダクションMD

**新設計での位置づけ**:
- そのまま維持
- 最小化機能は QC/Min MCP にも複製

### Phase 6: Export Server

**実装ファイル**:
- `servers/export_server.py` (178行)
- ParmEd統合

**主要ツール**:
1. `export_amber`: prmtop/inpcrd
2. `export_gromacs`: ParmEd変換
3. `export_openmm`: XML
4. `package_system`: ZIP化

**新設計での位置づけ**:
- そのまま維持
- Phase 1はAmberのみ、Phase 3でGROMACS/OpenMM追加

---

## 9. 外部ツール依存関係

### 環境セットアップ（推奨: conda + uv）

#### 1. conda環境作成と科学計算ツールのインストール

```bash
# conda環境作成
conda create -n mcp-md python=3.11
conda activate mcp-md

# 科学計算ツール（conda-forge推奨）
conda install -c conda-forge ambertools packmol smina pdbfixer

# uvのインストール（conda環境内）
pip install uv
```

#### 2. conda環境内でuvを使ってPythonパッケージをインストール

```bash
# conda環境がアクティブな状態で実行
conda activate mcp-md

# 基本依存関係のインストール（conda環境に直接インストール）
uv pip install -e .

# または、pyproject.tomlから直接インストール
uv pip install --project pyproject.toml

# 特定のLLMプロバイダーも含める場合
uv pip install -e ".[openai]"      # OpenAI/LM Studio
uv pip install -e ".[anthropic]"   # Claude
uv pip install -e ".[google]"      # Gemini

# 全てのオプション依存関係
uv pip install -e ".[openai,anthropic,google,dev]"
```

#### 3. 実行方法

```bash
# conda環境がアクティブな状態で
conda activate mcp-md

# uv runを使って実行（高速起動）
uv run python main.py

# またはMCPサーバーの起動
uv run python -m servers.structure_server

# LangGraphワークフローの実行
uv run python -m core.workflow_graph

# 通常のpythonコマンドも使用可能
python main.py
```

#### 4. pyproject.toml設定例

```toml
[project]
name = "mcp-md"
version = "0.1.0"
description = "Amber-focused MD setup with LangGraph + MCP"
requires-python = ">=3.11"
dependencies = [
    "boltz>=2.0.0",
    "pdb2pqr>=3.1.0",
    "propka>=3.5.0",
    "rdkit>=2023.9.1",
    "openmm>=8.3.1",
    "parmed>=4.3.0",
    "fastmcp>=0.1.0",
    "langchain-core>=1.0.0",
    "langgraph>=0.2.0",
]

[project.optional-dependencies]
openai = ["langchain-openai>=0.2.0"]
anthropic = ["langchain-anthropic>=0.3.0"]
google = ["langchain-google-genai>=0.1.0"]
dev = ["pytest>=7.0", "black>=24.0", "ruff>=0.1.0"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

### 主要パッケージ一覧

#### 科学計算ツール（conda経由）
- **AmberTools**: 完全OSSのAmberツール群（tleap, antechamber, parmchk2等）
- **Packmol**: 溶媒・膜系の構築
- **Smina**: ドッキングツール（AutoDock Vina fork）
- **PDBFixer**: PDB構造修復

#### Pythonパッケージ（uv経由）
- **Boltz-2**: 構造予測・複合体生成
- **PDB2PQR + PROPKA**: プロトネーション
- **RDKit**: ケモインフォマティクス
- **OpenMM + ParmEd**: MD計算とトポロジー変換
- **FastMCP**: MCPサーバー実装
- **LangChain Core + LangGraph**: ワークフロー構築

### 注意事項

1. **conda + uv併用の方針**: 
   - **conda環境**: 科学計算ツール（C/C++バイナリ）+ Python本体
   - **uv pip**: conda環境内でPythonパッケージをインストール（高速）
   - **uv run**: conda環境内でスクリプト実行（キャッシュ活用で高速起動）
   
2. **uv独自の仮想環境は使わない**: 
   - `uv sync` は実行しない（`.venv`を作成してしまう）
   - `uv pip install` を使ってconda環境に直接インストール
   - `uv run` はconda環境のPythonを使用
   
3. **依存関係のロック**: 
   - conda環境では `conda env export > environment.yml` でロック
   - Pythonパッケージは `uv pip compile pyproject.toml -o requirements.txt` でロック可能
   - または `pip freeze > requirements.txt`

4. **MCP統合**: 
   - `langchain-mcp`パッケージは不要
   - MCPツールを`langchain_core.tools.Tool`として直接実装

---

## 10. 参考資料

### 主要論文

#### Boltz-2
```bibtex
@article{passaro2025boltz2,
  author = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and others},
  title = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  journal = {bioRxiv},
  year = {2025}
}
```

### 外部リンク

#### 主要フレームワーク
- **LangChain**: https://github.com/langchain-ai/langchain
  - **LangChain 1.0 ドキュメント**: https://python.langchain.com/docs/
  - **LangChain 1.0 移行ガイド**: https://python.langchain.com/docs/versions/v0_3/migrating_chains/
- **LangGraph**: https://github.com/langchain-ai/langgraph
  - **LangGraph ドキュメント**: https://langchain-ai.github.io/langgraph/
  - **チェックポイント機能**: https://langchain-ai.github.io/langgraph/concepts/persistence/
- **FastMCP**: https://github.com/jlowin/fastmcp
- **MCP Protocol**: https://modelcontextprotocol.io/
- **uv**: https://github.com/astral-sh/uv
  - **uvドキュメント**: https://docs.astral.sh/uv/

#### 科学計算ツール
- **Boltz-2**: https://github.com/jwohlwend/boltz
- **AmberTools**: https://ambermd.org/AmberTools.php
- **OpenMM**: https://openmm.org/
- **PoseBusters**: https://github.com/maabuu/posebusters

---

## 11. まとめ

### 現在地

- ✅ 7つのFastMCPサーバー実装済み（基本機能）
- ✅ 共通ライブラリ（`common/`）完成
- ❌ LangGraph統合未実装（最優先）
- ❌ ワークフローノード未実装
- ❌ チェックポイント機能未実装

### 次のステップ

1. **LangGraph Agent統合**（最優先、2週間）
   - StateGraph定義（`core/workflow_graph.py`）
   - WorkflowState定義（`core/workflow_state.py`）
   - MCP統合（`core/mcp_integration.py`）
   
2. **ワークフローノード実装**（2週間）
   - 各工程のノード実装（`core/workflow_nodes.py`）
   - 条件付きエッジとエラーハンドリング
   - 決定ログ統合
   
3. **チェックポイント機能**（1週間）
   - SQLiteベースの永続化
   - 中断・再開機能
   - 実行履歴表示
   
4. **MVP完成**（Phase 1完了、4週間）

### プロジェクト特性

- **非競合**: CHARMM-GUIと棲み分け（Amber特化）
- **将来性**: MCP標準でツール永続化、LLM/実行基盤の更新に強い
- **拡張性**: LangGraphのモジュラー設計により、ノード/エッジの追加が容易
- **標準準拠**: LangChain 1.0の設計思想に準拠
  - 従来のchains/agentsは使用せず、LangGraphのStateGraphを直接利用
  - 固定スケルトンワークフローに最適化
  - MCPツールをLangChain `Tool`として統合（`langchain-core.tools.Tool`）
