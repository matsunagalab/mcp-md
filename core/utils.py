"""
Common utility functions for MCP-MD.
"""

import os
import logging
import subprocess
from pathlib import Path
from typing import Optional, Union
from datetime import datetime

logger = logging.getLogger(__name__)


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """ロガーをセットアップ
    
    Args:
        name: ロガー名
        level: ログレベル
    
    Returns:
        設定済みロガー
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger


def ensure_directory(path: Union[str, Path]) -> Path:
    """ディレクトリが存在することを確認、なければ作成
    
    Args:
        path: ディレクトリパス
    
    Returns:
        Pathオブジェクト
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def check_file_exists(path: Union[str, Path]) -> bool:
    """ファイルの存在確認
    
    Args:
        path: ファイルパス
    
    Returns:
        True if file exists
    """
    return Path(path).is_file()


def run_command(
    cmd: list[str],
    cwd: Optional[Union[str, Path]] = None,
    timeout: Optional[int] = None,
    capture_output: bool = True
) -> subprocess.CompletedProcess:
    """外部コマンドを実行
    
    Args:
        cmd: コマンドと引数のリスト
        cwd: 作業ディレクトリ
        timeout: タイムアウト（秒）
        capture_output: 出力をキャプチャするか
    
    Returns:
        CompletedProcessオブジェクト
    
    Raises:
        subprocess.CalledProcessError: コマンド実行失敗
        subprocess.TimeoutExpired: タイムアウト
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=capture_output,
            text=True,
            timeout=timeout,
            check=True
        )
        logger.debug(f"Command completed successfully")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.stderr}")
        raise
    except subprocess.TimeoutExpired as e:
        logger.error(f"Command timed out after {timeout}s")
        raise


def generate_timestamp() -> str:
    """現在時刻のタイムスタンプ文字列を生成
    
    Returns:
        ISO 8601形式のタイムスタンプ
    """
    return datetime.now().isoformat()


def generate_unique_id(prefix: str = "") -> str:
    """ユニークIDを生成
    
    Args:
        prefix: IDのプレフィックス
    
    Returns:
        ユニークID文字列
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if prefix:
        return f"{prefix}_{timestamp}"
    return timestamp


def read_fasta(fasta_path: Union[str, Path]) -> dict[str, str]:
    """FASTAファイルを読み込み
    
    Args:
        fasta_path: FASTAファイルパス
    
    Returns:
        {sequence_id: sequence}の辞書
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def write_fasta(sequences: dict[str, str], output_path: Union[str, Path]):
    """FASTAファイルを書き込み
    
    Args:
        sequences: {sequence_id: sequence}の辞書
        output_path: 出力FASTAファイルパス
    """
    with open(output_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            # 80文字ごとに改行
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def count_atoms_in_pdb(pdb_path: Union[str, Path]) -> int:
    """PDBファイル内の原子数をカウント
    
    Args:
        pdb_path: PDBファイルパス
    
    Returns:
        原子数
    """
    count = 0
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                count += 1
    return count


def get_pdb_chains(pdb_path: Union[str, Path]) -> list[str]:
    """PDBファイルのチェインIDを取得
    
    Args:
        pdb_path: PDBファイルパス
    
    Returns:
        チェインIDのリスト
    """
    chains = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21:22].strip()
                if chain_id:
                    chains.add(chain_id)
    return sorted(chains)


def check_external_tool(tool_name: str) -> bool:
    """外部ツールが利用可能かチェック
    
    Args:
        tool_name: ツール名（コマンド名）
    
    Returns:
        True if tool is available in PATH
    """
    try:
        result = subprocess.run(
            ['which', tool_name],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except Exception:
        return False


def get_conda_env_path(env_name: str) -> Optional[str]:
    """conda環境のパスを取得
    
    Args:
        env_name: conda環境名
    
    Returns:
        環境パス（存在しない場合None）
    """
    try:
        result = subprocess.run(
            ['conda', 'env', 'list'],
            capture_output=True,
            text=True,
            check=True
        )
        for line in result.stdout.split('\n'):
            if line.strip().startswith(env_name):
                parts = line.split()
                if len(parts) >= 2:
                    return parts[-1]
    except Exception as e:
        logger.warning(f"Failed to get conda env path: {e}")
    
    return None


def validate_smiles(smiles: str) -> bool:
    """SMILES文字列の妥当性を簡易チェック
    
    Args:
        smiles: SMILES文字列
    
    Returns:
        True if valid (basic check)
    """
    # 基本的な文字チェック
    if not smiles or not smiles.strip():
        return False
    
    # 許可される文字
    allowed_chars = set("CNOPSFClBrI[]()=#-+\\/@0123456789")
    return all(c in allowed_chars for c in smiles)


class WorkingDirectory:
    """一時的に作業ディレクトリを変更するコンテキストマネージャー"""
    
    def __init__(self, path: Union[str, Path]):
        self.path = Path(path)
        self.original_dir = Path.cwd()
    
    def __enter__(self):
        os.chdir(self.path)
        return self.path
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.original_dir)

