# Phase 1, 2, 4 実装詳細

## 実装概要

Phase 1（Structure Server）、Phase 2（Ligand Server）、Phase 4（Assembly Server）の主要機能を完全実装しました。

---

## Phase 1: Structure Server

### 実装ファイル
- `servers/structure_server.py` - MCPサーバー本体
- `tools/boltz2_wrapper.py` - Boltz-2ラッパー
- `tools/pdbfixer_wrapper.py` - PDBFixerラッパー
- `tools/pdb2pqr_wrapper.py` - PDB2PQR+PROPKAラッパー **[新規追加]**

### 実装機能

#### 1. 構造取得
- **`fetch_pdb`**: PDB、AlphaFold、PDB-REDOから構造取得
  - 対応ソース: `pdb`, `alphafold`, `pdb-redo`
  - 自動ダウンロード・検証

#### 2. Boltz-2構造予測
- **`predict_structure_boltz2`**: FASTAから構造予測
  - MSA使用オプション
  - 複数モデル生成
  - 信頼度評価（pLDDT、pAE）
  
- **`predict_complex_with_affinity`**: FASTA + SMILES → 複合体 + 親和性
  - タンパク質-リガンド複合体同時予測
  - 結合親和性計算（IC50、バインダー確率）
  
- **`screen_ligands_boltz2`**: バーチャルスクリーニング
  - 複数リガンド同時評価
  - ランキング機能

- **`complete_missing_residues_boltz2`**: 欠損残基補完
  - 高精度補完
  - テンプレートベース

#### 3. 構造クリーニング
- **`clean_structure`**: PDBFixerによるクリーニング
  - 欠損原子・残基補完
  - 異常altloc処理
  - 水分子除去

- **`add_hydrogens`**: 水素付与（pH指定）
  - pH依存プロトネーション
  - PDBFixer統合

#### 4. 高度なプロトネーション **[新規追加]**
- **`protonate_structure`**: PDB2PQR + PROPKA
  - pH指定プロトネーション
  - pKa計算統合
  - 力場指定（AMBER、CHARMM、PARSE）

#### 5. 修飾検出 **[新規追加]**
- **`detect_modifications`**: 構造修飾の自動検出
  - ジスルフィド結合
  - 修飾残基（MODRES）
  - 金属サイト（Zn、Mg、Ca、Fe等）

#### 6. 構造検証
- **`validate_structure`**: PDB構造の健全性チェック
  - 原子数カウント
  - チェーン情報抽出

### 技術仕様

#### Boltz-2統合
```python
# YAML入力生成
- name: target_protein
  sequences:
    - protein:
        id: ["A"]
        sequence: "MKTAYIAK..."
- name: ligand
  smiles: "CC(=O)Oc1ccccc1C(=O)O"
```

#### PDB2PQR使用例
```bash
pdb2pqr --ff amber --with-ph 7.4 --titration-state-method propka input.pdb output.pqr
```

---

## Phase 2: Ligand Server

### 実装ファイル
- `servers/ligand_server.py` - MCPサーバー本体
- `tools/rdkit_wrapper.py` - RDKitラッパー
- `tools/ambertools_wrapper.py` - AmberToolsラッパー

### 実装機能

#### 1. 3D構造生成
- **`smiles_to_3d`**: SMILES → 3D構造
  - RDKit ETKDG法
  - コンフォーマ最適化（MMFF94）

#### 2. AmberTools GAFF2パラメータ化
- **`generate_gaff_params`**: antechamber + parmchk2
  - 電荷計算: AM1-BCC（推奨）、Gasteiger、RESP
  - GAFF2力場パラメータ
  - 欠損パラメータ自動補完

- **`create_ligand_lib`**: tleap用ライブラリ作成
  - .lib + .frcmodファイル生成
  - tleap統合準備

#### 3. 完全自動ワークフロー
- **`parameterize_ligand_complete`**: SMILES → tleapライブラリ（一括実行）
  - 3D生成 → パラメータ化 → ライブラリ作成
  - ワンコマンド実行

### 技術仕様

#### AmberTools電荷計算
```bash
# AM1-BCC（推奨、最もバランスが良い）
antechamber -i input.pdb -fi pdb -o output.mol2 -fo mol2 \
  -c bcc -at gaff2 -rn LIG -nc 0

# RESP（高精度、計算コスト高）
antechamber -i input.pdb -fi pdb -o output.mol2 -fo mol2 \
  -c resp -at gaff2 -rn LIG -nc 0

# Gasteiger（高速、低精度）
antechamber -i input.pdb -fi pdb -o output.mol2 -fo mol2 \
  -c gas -at gaff2 -rn LIG -nc 0
```

#### パラメータ補完
```bash
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -a Y
```

---

## Phase 4: Assembly Server

### 実装ファイル
- `servers/assembly_server.py` - MCPサーバー本体
- `tools/ambertools_wrapper.py` - AmberTools tleapラッパー
- `tools/packmol_wrapper.py` - Packmol-Memgenラッパー **[新規追加]**

### 実装機能

#### 1. tleap系構築
- **`build_system_tleap`**: 完全なMD系構築
  - タンパク質読み込み
  - 配位子統合
  - 力場適用（ff19SB、ff14SB等）
  - 溶媒化（TIP3P、OPC）
  - イオン付与（中性化 + 塩濃度）

#### 2. 膜タンパク質系 **[新規追加]**
- **`build_membrane_system`**: Packmol-Memgen統合
  - 脂質二重層構築
  - 脂質組成指定（POPC、POPE、CHOL等）
  - 膜タンパク質埋め込み
  - 自動配向最適化

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

## 統合ワークフロー例

### 1. PDB + SMILES → MD系

```python
# Phase 1: 構造取得・クリーニング
pdb = await structure_server.fetch_pdb(pdb_id="1ABC")
cleaned = await structure_server.clean_structure(pdb["file_path"])
protonated = await structure_server.protonate_structure(cleaned["output"], ph=7.4)

# Phase 2: 配位子パラメータ化
ligand = await ligand_server.parameterize_ligand_complete(
    smiles="CC(=O)Oc1ccccc1C(=O)O",
    residue_name="ASP"
)

# Phase 4: 系構築
system = await assembly_server.build_system_tleap(
    protein_pdb=protonated["output"],
    ligand_lib=ligand["library"],
    salt_conc=0.15
)
```

### 2. Boltz-2予測 → MD系

```python
# Phase 1: Boltz-2構造予測
prediction = await structure_server.predict_structure_boltz2(
    fasta="MKTAYIAK...",
    use_msa=True
)

# Phase 4: MD系構築
system = await assembly_server.build_system_tleap(
    protein_pdb=prediction["structures"][0],
    water_model="opc"
)
```

### 3. 膜タンパク質系

```python
# Phase 1: 構造準備
membrane_prot = await structure_server.fetch_pdb(pdb_id="2RH1")
cleaned = await structure_server.clean_structure(membrane_prot["file_path"])

# Phase 4: 膜系構築
membrane_sys = await assembly_server.build_membrane_system(
    protein_pdb=cleaned["output"],
    lipid_composition={"POPC": 0.7, "POPE": 0.3}
)
```

---

## 出力ファイル

### Phase 1出力
- `*.pdb` - クリーニング済み構造
- `boltz2_prediction/model_*.pdb` - Boltz-2予測構造
- `boltz2_prediction/boltz_input.yaml` - Boltz-2入力
- `boltz2_prediction/predictions.json` - 親和性データ

### Phase 2出力
- `gaff_params/*.mol2` - GAFF2パラメータ付きMOL2
- `gaff_params/*.frcmod` - 欠損パラメータ
- `ligand_lib/*.lib` - tleapライブラリ

### Phase 4出力
- `system/system.prmtop` - Amber topology
- `system/system.inpcrd` - Amber coordinates
- `system/tleap.in` - tleapスクリプト
- `membrane_system/membrane_system.pdb` - 膜系構造

---

## 依存関係

### Python パッケージ
```
boltz>=2.0.0          # Boltz-2構造予測
rdkit>=2023.9.1       # 配位子3D生成
openmm>=8.3.1         # OpenMM統合
pdbfixer>=1.9         # PDB修復
```

### 外部ツール（conda経由）
```bash
conda install -c conda-forge ambertools packmol smina
pip install pdb2pqr propka
```

---

## 次のステップ

Phase 1、2、4の実装が完了しました。次は：

1. **Phase 3: Docking Server** - smina統合ドッキング
2. **Phase 5: Protocol Server** - OpenMM MDスクリプト生成
3. **Phase 6: Export Server** - 形式変換・パッケージング
4. **Planner/Validator** - ワークフロー自動化

詳細は`examples/phase_124_workflow.md`を参照してください。

