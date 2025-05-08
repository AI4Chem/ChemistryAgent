from tools import ToolPool, ToolSet

# 数据集生成
from case.toolcase_single import SingleToolCaseGeneration
from case.toolcase_multiple import MultiToolCaseGeneration

from model import GPT, InternLM, LLaMA2, LLaMA3, FMModel

# 评测
from eval.process import construct_dataset_with_retriever_result
from eval import SingleToolEvaluation, MultiToolEvaluation, ResponseEvaluation

# Agent
from chat import Agent

model_info_list = [
    # {
    #     "model_name": "Llama3_1Instruct",
    #     "model_type": "Transformers_CausalLM",
    #     "model_path": "/mnt/hwfile/ai4chem/wumengsong/PTM/Meta-Llama-3_1-8B-Instruct"
    # },
    # {
    #     "model_name": "Llama2Chat",
    #     "model_type": "Transformers_CausalLM",
    #     "model_path": "/mnt/hwfile/ai4chem/wumengsong/PTM/Meta-Llama-2-7B-Chat"
    # },
    # {
    #     "model_name": "ChemLLMLoRAv1",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "ChemLLM",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/sft-v1"
    # },
    # {
    #     "model_name": "ChemLLMLoRAv2",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "ChemLLM",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/sft-v2"
    # },
    # {
    #     "model_name": "ChemLLMLoRAv2_1",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "ChemLLM",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/sft-v2_1"
    # },
    # {
    #     "model_name": "Llama3LoRAv1",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_1Instruct",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/sft_v1"
    # },
    # {
    #     "model_name": "Llama3LoRAv2",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_1Instruct",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/sft_v2"
    # },
    # {
    #     "model_name": "Llama3LoRAv2_1",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_1Instruct",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/sft_v2_1"
    # },
    {
        "model_name": "Qwen2_5-7B-Instruct",
        "model_type": "Transformers_CausalLM",
        "model_path": "/mnt/hwfile/ai4chem/wumengsong/PTM/Qwen_2_5-7B-Instruct"
    },
    # {
    #     "model_name": "Llama3DPOv1",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v1"
    # },
    # {
    #     "model_name": "Llama3DPOv2",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v2"
    # },
    # {
    #     "model_name": "Llama3DPOv3",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v3"
    # },
    # {
    #     "model_name": "Llama3DPOv4",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v4"
    # },
    # {
    #     "model_name": "Llama3DPOv4CP3",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v4/checkpoint-30"
    # },
    # {
    #     "model_name": "Llama3DPOv4CP6",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v4/checkpoint-60"
    # },
    # {
    #     "model_name": "Llama3DPOv4CP9",
    #     "model_type": "Peft_CausalLM",
    #     "model_foundation": "Llama3_LoRA_Merge",
    #     "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/DPO_v4/checkpoint-90"
    # },
    # {
    #     "model_name": "Llama3_LoRA_Merge",
    #     "model_type": "Transformers_CausalLM",
    #     "model_path": "/mnt/hwfile/ai4chem/wumengsong/PTM/llama3_lora_v2_1"
    # },
    {
        "model_name": "ChemLLM",
        "model_type": "Transformers_CausalLM",
        "model_path": "/mnt/hwfile/ai4chem/CKPT/ChemLLM-20B-Chat-DPO"
    },
]

# * = = =  加载 工具集 
toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])

def main():
    # tools = ToolPool("chem_lib")
    # print(tools.call_tool("acidity_calculation",[9.47]))
    # print(tools.call_tool("get_element_properties",["Zn"]))
    # print(tools.call_tool("analyze_combustion",[44,18]))

    # tools = ToolPool("chemcrow")
    # print(tools.call_tool("Query2CAS",["CC"]))

    # tools = ToolPool("chem_lib")
    # print(tools.toolinfo["reaction_stoichiometry_amounts"])

    # tools = ToolPool("chem_lib")
    # print(tools.call_tool("judge_the_balance_of_reaction",["2Al + 3H2SO4 --> Al2(SO4)3 + 3H2"]))


    # * = = = 生成调用结果 
    # for tooldir in toolset.tool_pool_dict:
    #     toolset.tool_pool_dict[tooldir].generate_tool_calling_result()

    # * = = = 生成single用例 

    # case_genera = SingleToolCaseGeneration(toolset)
    # case_genera.generate_single_case_with_calling()

    # * = = =  生成multi用例 

    # case_genera = MultiToolCaseGeneration(toolset)

    # - step:1
    # # exception_list = ["chem_lib/limiting_reagent_of_reaction"]
    # exception_list = ["chem_lib/limiting_reagent_of_reaction", "chemistrytools/get_compound_CID"]
    # case_genera.mix_all_toolpools(exception_list)
    # # case_genera.mix_all_toolpools()
    # case_genera.choose_candidate_tools()

    # # - step:2
    # case_genera.generate_tool_calling_chain()

    # # - step:3
    # case_genera.generate_query()

    # * = = = 

    # input_data_path_list = [
    #     "data/single/single_train.jsonl",
    #     "data/single/single_dev.jsonl",
    #     "data/single/single_test.jsonl",
    #     "data/multiple/method_SealTools/multiple_train.jsonl",
    #     "data/multiple/method_SealTools/multiple_dev.jsonl",
    #     "data/multiple/method_SealTools/multiple_test.jsonl"]
    # output_data_path_list = [
    #     "data/retrieved/single_train_NVEmbed.jsonl",
    #     "data/retrieved/single_dev_NVEmbed.jsonl",
    #     "data/retrieved/single_test_NVEmbed.jsonl",
    #     "data/retrieved/multiple_train_NVEmbed.jsonl",
    #     "data/retrieved/multiple_dev_NVEmbed.jsonl",
    #     "data/retrieved/multiple_test_NVEmbed.jsonl"]
    # construct_dataset_with_retriever_result(toolset, input_data_path_list, output_data_path_list)

    # checkpoint_path = "/home/bingxing2/ailab/group/ai4phys/EXPORT/20b_FT_4_28_zhangdi_DPO_5_3"
    # lora_path = "/home/bingxing2/ailab/wumengsong/LLaMA-Factory/saves/chemllm-20b-DPO/lora/sft-1e5"
    # model = InternLM(checkpoint_path, lora_path)

    # checkpoint_path = "/home/bingxing2/ailab/group/ai4phys/EXPORT/20b_FT_4_28_zhangdi_DPO_5_3"
    # checkpoint_path = "/home/bingxing2/ailab/group/ai4phys/CKPT/internlm2_5-20b-chat"
    # model = InternLM(checkpoint_path)

    # checkpoint_path = "/home/bingxing2/ailab/group/ai4phys/CKPT/Meta-Llama-2-7B-Chat"
    # model = LLaMA2(checkpoint_path)

    # checkpoint_path = "/home/bingxing2/ailab/group/ai4phys/wxz/Meta-Llama-3.1-8B-Instruct"
    # model = LLaMA3(checkpoint_path)

    # model = GPT()

    import fastmindapi as FM

    # * 评测

    # 模型：FM
    client = FM.Client(IP="10.140.24.115",PORT=10003)
    client.add_model_info_list(model_info_list)
    model = FMModel(client)
    for model_info in model_info_list[:-1]: # *
        for unload_model in model_info_list:
            client.unload_model(unload_model["model_name"])
        print('【Model】:',model_info)
        model_name = model_info["model_name"] # *
        if "Peft" in model_info["model_type"]:
            client.load_model(model_info["model_foundation"])
        client.load_model(model_name)
        model.model_name = model_name

    # 模型：GPT API
    # model_name_list = ["gpt-4o-mini", "claude-3-5-sonnet-20240620"] # *
    # for model_name in model_name_list: # *
    #     model = GPT(model_name)

        # 单工具评测
        input_data_path = "data/single/single_test_DPR.jsonl"
        output_data_path = f"result/single_test_DPR_{model_name}.jsonl"
        singleeval = SingleToolEvaluation(toolset, model)
        singleeval.eval_0_shot(input_data_path, output_data_path)
        SingleToolEvaluation.calculate_score(output_data_path)
        responseeval = ResponseEvaluation(toolset, model)
        responseeval.eval(output_data_path)

        # 多工具评测
        input_data_path = "data/multiple/method_SealTools/multiple_test_DPR.jsonl"
        output_data_path = f"result/multiple_test_DPR_{model_name}.jsonl"
        multieval = MultiToolEvaluation(toolset, model)
        multieval.eval_0_shot(input_data_path, output_data_path)
        MultiToolEvaluation.calculate_score(output_data_path)
        responseeval = ResponseEvaluation(toolset, model)
        responseeval.eval(output_data_path)

    # * 调用工具与否对比测试
    # model = GPT("gpt-4o-mini")
    # singleeval = SingleToolEvaluation(toolset, model)
    # singleeval.generate_without_tool_calling("data/single/single_test_comparison.jsonl", "result/single_test_gpt-4o-mini_WoT_response.jsonl")
    # responseeval = ResponseEvaluation(toolset, model)
    # responseeval.score("result/single_test_gpt-4o-mini_WoT_response.jsonl", "data/single/single_test_response.jsonl")

    # * GPT打分
    model_name_list = ["Llama3_1Instruct", 
                    #    "Llama2Chat", 
                       "Qwen2_5-7B-Instruct",
                    #    "ChemLLMLoRAv1", 
                    #    "ChemLLMLoRAv2", 
                    #    "ChemLLMLoRAv2_1",
                    #    "Llama3LoRAv1",
                    #    "Llama3LoRAv2",
                    #    "Llama3LoRAv2_1",
                    #    "Llama3DPOv1"
                    #    "ChemLLM",
                       "gpt-4o-mini", 
                       "claude-3-5-sonnet-20240620"]
    for model_name in model_name_list:
        model = GPT("gpt-4o-mini")
        responseeval = ResponseEvaluation(toolset, model)

        gold_data_path = "data/single/single_test_response.jsonl"
        input_data_path = f"result/single_test_DPR_{model_name}_response.jsonl"
        responseeval.score(input_data_path, gold_data_path)

        gold_data_path = "data/multiple/method_SealTools/multiple_test_response.jsonl"
        input_data_path = f"result/multiple_test_DPR_{model_name}_response.jsonl"
        responseeval.score(input_data_path, gold_data_path)

    # * 打印最终得分
    import glob
    for file_path in glob.glob("result/*"):
        if file_path.endswith("score.jsonl"):
            print(file_path)
            score = ResponseEvaluation.count(file_path)
            print(score)
            print(" = = = = = = ")

    # * = = =
    # import os
    # os.environ["TOKENIZERS_PARALLELISM"] = "false"
    # retriever_checkpoint = "/mnt/hwfile/ai4chem/wumengsong/PTM/all-MiniLM-L6-v2"

    # agent = Agent(GPT(),toolset, retriever_checkpoint)
    # # agent = Agent(GPT(),toolset)
    # # query = "What are all the properties of the element with the atomic symbol \"Au\"?"
    # query = "Can you tell me about the properties of Au?"
    # # query = "告诉我一些关于金元素的相关属性信息?"
    # agent(query)


if __name__ == "__main__":
    main()
