#!/usr/bin/env python3
"""Clarification Agent Demo Script.

このスクリプトは、notebooks/1_clarification.ipynb で実装された
clarification & brief作成フローを対話的に体験できます。
"""

import asyncio
from dotenv import load_dotenv
from langchain_core.messages import HumanMessage

# 環境変数をロード
load_dotenv()

from mcp_md.clarification_agent import clarification_graph


def print_separator():
    """セパレーターを表示"""
    print("\n" + "="*60 + "\n")


async def run_clarification_demo():
    """Clarificationフローのデモを実行"""
    print("🚀 MD Setup Clarification Agent Demo")
    print_separator()
    print("このデモでは、MDシミュレーションのセットアップに必要な情報を")
    print("対話的に収集し、最終的に構造化されたSimulationBriefを生成します。")
    print_separator()
    
    # 初期メッセージの収集
    messages = []
    
    print("📝 MDシミュレーションの要件を入力してください。")
    print("   (例: PDB ID 1ABC のタンパク質でシミュレーションしたい)")
    print("   (例: アスピリンとP450の複合体をシミュレートしたい)")
    
    user_input = input("\n> ")
    messages.append(HumanMessage(content=user_input))
    
    # Clarificationフローを実行
    iteration = 0
    max_iterations = 5  # 無限ループ防止
    
    while iteration < max_iterations:
        print_separator()
        print(f"🤖 エージェント実行中... (反復 {iteration + 1})")
        
        # グラフを実行
        result = await clarification_graph.ainvoke({"messages": messages})
        
        # 最後のAIメッセージを取得
        ai_messages = [msg for msg in result["messages"] if msg.type == "ai"]
        if ai_messages:
            last_ai_message = ai_messages[-1]
            print(f"\n💬 エージェント: {last_ai_message.content}")
        
        # Simulation Briefが生成されたかチェック
        if result.get("simulation_brief"):
            print_separator()
            print("✅ 情報収集完了！Simulation Briefを生成しました。")
            print_separator()
            
            brief = result["simulation_brief"]
            print("📋 Simulation Brief:")
            print(f"  - PDB ID: {brief.pdb_id}")
            print(f"  - FASTA Sequence: {brief.fasta_sequence}")
            print(f"  - Ligand SMILES: {brief.ligand_smiles}")
            print(f"  - pH: {brief.ph}")
            print(f"  - Salt Concentration: {brief.salt_concentration} M")
            print(f"  - Water Model: {brief.water_model}")
            print(f"  - Box Padding: {brief.box_padding} Å")
            print(f"  - Force Field: {brief.force_field}")
            print(f"  - Use Boltz-2 Docking: {brief.use_boltz2_docking}")
            print(f"  - Refine with Smina: {brief.refine_with_smina}")
            print(f"  - Output Formats: {brief.output_formats}")
            
            print_separator()
            print("🎉 デモ完了！次は notebooks/2_setup_agent.ipynb でSetup Agentを実装します。")
            break
        
        # 追加の質問がある場合、ユーザーに回答を促す
        print("\n📝 回答を入力してください:")
        user_input = input("> ")
        
        # メッセージを更新
        messages = result["messages"] + [HumanMessage(content=user_input)]
        iteration += 1
    
    if iteration >= max_iterations:
        print_separator()
        print("⚠️  最大反復回数に達しました。")


def main():
    """メイン関数"""
    try:
        asyncio.run(run_clarification_demo())
    except KeyboardInterrupt:
        print("\n\n👋 デモを中断しました。")
    except Exception as e:
        print(f"\n❌ エラーが発生しました: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

