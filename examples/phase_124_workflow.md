# Phase 1, 2, 4統合ワークフロー例

このドキュメントでは、Phase 1（Structure Server）、Phase 2（Ligand Server）、Phase 4（Assembly Server）の主要機能を組み合わせた実践的なワークフロー例を示します。

## ワークフロー1: PDB + SMILES → 完全なMD系構築

### 概要
既存のPDB構造にリガンドを追加し、完全なMD入力系を構築します。

### ステップ

#### 1. PDB構造の取得とクリーニング（Phase 1）

```python
# Structure Serverを使用
from servers.structure_server import StructureServer

server = StructureServer()

# PDBから構造取得
pdb_result = await server.fetch_pdb(pdb_id="1ABC", source="pdb")
# => {"pdb_id": "1ABC", "file_path": "output/1ABC.pdb", "num_atoms": 1234, ...}

# 構造のクリーニング
cleaned = await server.clean_structure(
    pdb_file=pdb_result["file_path"],
    remove_water=True,
    fix_missing=True
)
# => {"output": "output/cleaned.pdb", "operations": [...]}

# pH指定プロトネーション（PDB2PQR+PROPKA）
protonated = await server.protonate_structure(
    pdb_file=cleaned["output"],
    ph=7.4,
    forcefield="AMBER"
)
# => {"output": "output/protonated_pdb2pqr.pdb", "ph": 7.4, ...}

# 修飾の検出
mods = await server.detect_modifications(pdb_file=protonated["output"])
# => {"disulfide_bonds": [...], "modified_residues": [...], "metal_sites": [...]}
```

#### 2. 配位子の準備とパラメータ化（Phase 2）

```python
# Ligand Serverを使用
from servers.ligand_server import LigandServer

ligand_server = LigandServer()

# SMILESから3D構造生成 → GAFF2パラメータ化 → ライブラリ作成（一括実行）
ligand_params = await ligand_server.parameterize_ligand_complete(
    smiles="CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    net_charge=0,
    residue_name="ASP",
    charge_method="bcc"  # AM1-BCC（推奨）
)
# => {
#     "gaff_mol2": "output/gaff_params/ASP_gaff.mol2",
#     "frcmod": "output/gaff_params/ASP.frcmod",
#     "library": "output/ligand_lib/ASP.lib",
#     "charges": [-0.5, 0.1, ...],
#     "total_charge": 0.0012,
#     ...
# }
```

#### 3. 系の組立（Phase 4）

```python
# Assembly Serverを使用
from servers.assembly_server import AssemblyServer

assembly_server = AssemblyServer()

# タンパク質-リガンド複合体の構築（溶媒化+イオン）
system = await assembly_server.build_system_tleap(
    protein_pdb=protonated["output"],
    ligand_lib=ligand_params["library"],
    forcefield="leaprc.protein.ff19SB",
    water_model="tip3p",
    box_padding=12.0,  # Angstroms
    salt_conc=0.15     # M (生理的塩濃度)
)
# => {
#     "prmtop": "output/system/system.prmtop",
#     "inpcrd": "output/system/system.inpcrd",
#     "leap_in": "output/system/tleap.in"
# }
```

---

## ワークフロー2: Boltz-2予測 → MD系構築

### 概要
Boltz-2でFASTAから構造予測し、MD系を構築します。

### ステップ

#### 1. Boltz-2構造予測（Phase 1）

```python
from servers.structure_server import StructureServer

server = StructureServer()

# FASTAから構造予測
prediction = await server.predict_structure_boltz2(
    fasta="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK...",
    use_msa=True,
    num_models=5
)
# => {
#     "structures": ["boltz2_prediction/model_0.pdb", ...],
#     "confidence": {"plddt": [85.2, 83.1, ...], "pae": [...]},
#     "yaml_input": "boltz2_prediction/boltz_input.yaml"
# }

# 最高信頼度モデルを選択
best_model = prediction["structures"][0]

# 構造検証
validation = await server.validate_structure(pdb_file=best_model)
# => {"valid": True, "num_atoms": 1450, "chains": ["A"]}
```

#### 2. MD系構築（Phase 4）

```python
from servers.assembly_server import AssemblyServer

assembly_server = AssemblyServer()

# タンパク質のみ系（溶媒化+イオン）
system = await assembly_server.build_system_tleap(
    protein_pdb=best_model,
    forcefield="leaprc.protein.ff19SB",
    water_model="opc",  # OPC水モデル（高精度）
    box_padding=10.0,
    salt_conc=0.15
)
```

---

## ワークフロー3: Boltz-2複合体予測 + 親和性評価

### 概要
FASTA + SMILESからBoltz-2で複合体構造と結合親和性を同時予測します。

### ステップ

#### 1. Boltz-2複合体+親和性予測（Phase 1）

```python
from servers.structure_server import StructureServer

server = StructureServer()

# タンパク質-リガンド複合体予測 + 親和性計算
complex_pred = await server.predict_complex_with_affinity(
    protein_fasta="MKTAYIAK...",
    ligand_smiles=["CC(=O)Oc1ccccc1C(=O)O"],  # Aspirin
    use_msa=True,
    num_models=3
)
# => {
#     "structures": ["boltz2_complex/model_0.pdb", ...],
#     "confidence": {...},
#     "affinity": {
#         "probability_binary": 0.85,  # バインダー確率（0-1）
#         "pred_value": -6.2,          # log10(IC50) in μM
#         "ic50_um": 0.63              # IC50 = 0.63 μM
#     },
#     ...
# }

print(f"結合親和性: IC50 = {complex_pred['affinity']['ic50_um']:.2f} μM")
print(f"バインダー確率: {complex_pred['affinity']['probability_binary']:.2%}")
```

#### 2. 配位子の再パラメータ化（Phase 2）

Boltz-2で予測された複合体構造から配位子を抽出し、AmberToolsでパラメータ化：

```python
from servers.ligand_server import LigandServer

ligand_server = LigandServer()

# 複合体から配位子を抽出（手動またはスクリプト）
# ... ligand.pdb作成 ...

# GAFF2パラメータ生成
ligand_params = await ligand_server.generate_gaff_params(
    ligand_file="ligand.pdb",
    net_charge=0,
    charge_method="bcc",
    residue_name="LIG"
)

# ライブラリ作成
lib = await ligand_server.create_ligand_lib(
    mol2_file=ligand_params["mol2"],
    frcmod_file=ligand_params["frcmod"],
    residue_name="LIG"
)
```

#### 3. MD系構築（Phase 4）

```python
from servers.assembly_server import AssemblyServer

assembly_server = AssemblyServer()

# 複合体系の構築
system = await assembly_server.build_system_tleap(
    protein_pdb=complex_pred["structures"][0],
    ligand_lib=lib["lib"],
    forcefield="leaprc.protein.ff19SB",
    water_model="tip3p",
    box_padding=12.0,
    salt_conc=0.15
)
```

---

## ワークフロー4: 膜タンパク質系の構築

### 概要
膜タンパク質をPackmol-Memgenで脂質二重層に埋め込みます。

### ステップ

#### 1. 構造準備（Phase 1）

```python
from servers.structure_server import StructureServer

server = StructureServer()

# 膜タンパク質のPDB取得
membrane_protein = await server.fetch_pdb(pdb_id="2RH1", source="pdb")

# クリーニング
cleaned = await server.clean_structure(
    pdb_file=membrane_protein["file_path"],
    remove_water=True
)
```

#### 2. 膜系構築（Phase 4）

```python
from servers.assembly_server import AssemblyServer

assembly_server = AssemblyServer()

# 脂質二重層膜系の構築
membrane_system = await assembly_server.build_membrane_system(
    protein_pdb=cleaned["output"],
    lipid_composition={
        "POPC": 0.7,  # 70% POPC
        "POPE": 0.2,  # 20% POPE
        "CHOL": 0.1   # 10% コレステロール
    },
    membrane_type="bilayer",
    dist_to_bilayer=15.0  # タンパク質から膜までの距離（Å）
)
# => {
#     "output_pdb": "output/membrane_system/membrane_system.pdb",
#     "membrane_type": "bilayer",
#     "lipid_composition": {...},
#     ...
# }
```

---

## 高度な使用例: バーチャルスクリーニング

### 複数リガンドの同時スクリーニング（Phase 1 + 2）

```python
from servers.structure_server import StructureServer

server = StructureServer()

# 候補化合物のSMILESリスト
candidate_smiles = [
    "CC(=O)Oc1ccccc1C(=O)O",     # Aspirin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",   # Ibuprofen
    # ... 100個の化合物 ...
]

# Boltz-2で一括スクリーニング
screening_results = await server.screen_ligands_boltz2(
    protein_fasta="MKTAYIAK...",
    ligand_smiles_list=candidate_smiles,
    screening_mode="binary"  # ヒット探索モード
)
# => {
#     "results": [
#         {
#             "smiles": "...",
#             "affinity_score": 0.92,
#             "structure": "complex_0.pdb",
#             "rank": 1
#         },
#         ...
#     ],
#     "ranked_by": "affinity_probability_binary",
#     "num_ligands": 100
# }

# トップ10ヒットを抽出
top_hits = screening_results["results"][:10]

# 各ヒットに対してMD系を構築
for hit in top_hits:
    # Phase 2: パラメータ化
    params = await ligand_server.parameterize_ligand_complete(
        smiles=hit["smiles"],
        residue_name=f"HIT{hit['rank']}"
    )
    
    # Phase 4: MD系構築
    md_system = await assembly_server.build_system_tleap(
        protein_pdb=hit["structure"],
        ligand_lib=params["library"]
    )
```

---

## まとめ

### Phase 1: Structure Server
- PDB取得（PDB、AlphaFold、PDB-REDO）
- Boltz-2構造予測（FASTA → 構造）
- Boltz-2複合体予測（FASTA + SMILES → 複合体 + 親和性）
- バーチャルスクリーニング
- PDBFixer構造クリーニング
- PDB2PQR+PROPKAプロトネーション
- 修飾検出（ジスルフィド結合、金属サイト）

### Phase 2: Ligand Server
- RDKit 3D構造生成
- AmberTools GAFF2パラメータ化（AM1-BCC電荷）
- tleapライブラリ作成
- 完全自動ワークフロー

### Phase 4: Assembly Server
- tleap系構築（タンパク質、複合体）
- 溶媒化（TIP3P、OPC）
- イオン付与（中性化、塩濃度指定）
- Packmol-Memgen膜系構築

これらを組み合わせることで、様々なMD入力系を効率的に構築できます。

