import tools
from tools import ToolPool, ToolSet
from tools.toolinfo import matgen_tool_package

# 数据集生成
from case.toolcase_single import SingleToolCaseGeneration
from case.toolcase_multiple import MultiToolCaseGeneration

from model import GPT, FMModel

# 评测
from eval.process import construct_dataset_with_retriever_result
from eval import SingleToolEvaluation, MultiToolEvaluation, ResponseEvaluation

# Agent
from chat import Agent

model_info_list = [
    {
        "model_name": "Llama3_1Instruct",
        "model_type": "Transformers_CausalLM",
        "model_path": "/mnt/hwfile/ai4chem/wumengsong/PTM/Meta-Llama-3_1-8B-Instruct"
    },
    {
        "model_name": "Llama3LoRAv2_1",
        "model_type": "Peft_CausalLM",
        "model_foundation": "Llama3_1Instruct",
        "model_path": "/mnt/petrelfs/wumengsong/checkpoint/llama3_lora/sft_v2_1"
    },
    {
        "model_name": "Qwen2_5-7B-Instruct",
        "model_type": "Transformers_CausalLM",
        "model_path": "/mnt/hwfile/ai4chem/wumengsong/PTM/Qwen_2_5-7B-Instruct"
    },
]

# * = = =  加载 工具集 
toolset = ToolSet(matgen_tool_package, 
                  toolpackage_path="src/toolpool/pymatgen_common_PROCESSED")

def main():
    # tool_name_list = toolset.get_all_tool_names_list()
    # print(len(tool_name_list))

    # * = = = 生成调用结果 
    # for tooldir in toolset.tool_pool_dict:
    #     toolset.tool_pool_dict[tooldir].generate_tool_calling_result(dir="pymatgen/tool_calling_result", execution=False)

    # * = = = 生成single用例 

    # case_genera = SingleToolCaseGeneration(toolset)
    # case_genera.generate_single_case_with_calling(calling_dir="pymatgen/tool_calling_result", case_dir="pymatgen/single/single_case")

    # * = = =  生成multi用例 

    # case_genera = MultiToolCaseGeneration(toolset)

    # path_1 = "data/pymatgen/multiple/method_SealTools/candidate_tool.jsonl"
    # path_2 = "data/pymatgen/multiple/method_SealTools/tool_calling_chain.jsonl"
    # path_2_1 = "data/pymatgen/multiple/method_SealTools/tool_calling_chain_p.json"
    # path_3 = "data/pymatgen/multiple/method_SealTools/multiple_case.jsonl"

    # # - step:1
    # case_genera.mix_all_toolpools()
    # case_genera.choose_candidate_tools(candidate_tool_path=path_1)

    # - step:2
    # case_genera.generate_tool_calling_chain(calling_dir="pymatgen/tool_calling_result", candidate_tool_path=path_1, tool_calling_chain_path=path_2, execution=False)

    # case_genera.sp_process(path_2)
    # case_genera.postprocess_calling_chain_format(path_2, path_2_1)

    # # - step:3
    # case_genera.generate_query(tool_calling_chain_path=path_2_1, multiple_case_path=path_3)

    # *
    # input_data_path_list = [
    #     "data/pymatgen/single/single_train.jsonl",
    #     "data/pymatgen/single/single_dev.jsonl",
    #     "data/pymatgen/single/single_test.jsonl",
    #     "data/pymatgen/multiple/method_SealTools/multiple_test.jsonl",
    #     "data/pymatgen/multiple/method_SealTools/multiple_notest.jsonl"]
    # output_data_path_list = [
    #     "data/retrieved/single_mat_train_NVEmbed.jsonl",
    #     "data/retrieved/single_mat_dev_NVEmbed.jsonl",
    #     "data/retrieved/single_mat_test_NVEmbed.jsonl",
    #     "data/retrieved/multiple_mat_test_NVEmbed.jsonl",
    #     "data/retrieved/multiple_mat_notest_NVEmbed.jsonl"]
    # construct_dataset_with_retriever_result(toolset, input_data_path_list, output_data_path_list)

    import fastmindapi as FM

    # * 评测
    # 模型：FM
    # client = FM.Client(IP="10.140.24.28",PORT=10003)
    # client.add_model_info_list(model_info_list)
    # model = FMModel(client)
    # for model_info in model_info_list[:]: # *
    #     for unload_model in model_info_list:
    #         client.unload_model(unload_model["model_name"])
    #     print('【Model】:',model_info)
    #     model_name = model_info["model_name"] # *
    #     if "Peft" in model_info["model_type"]:
    #         client.load_model(model_info["model_foundation"])
    #     client.load_model(model_name)
    #     model.model_name = model_name

    # 模型：GPT API
    model_name_list = ["gpt-4o-mini", "deepseek-r1"] # * , "claude-3-5-sonnet-20240620"
    # model_name_list = ["Llama3_1Instruct"]
    for model_name in model_name_list: # *
        model = GPT(model_name)

        # 单工具评测
        input_data_path = "data/retrieved/single_mat_test_NVEmbed.jsonl"
        output_data_path = f"result_NV/mat_single_test_NV_{model_name}.jsonl"
        singleeval = SingleToolEvaluation(toolset, model)
        singleeval.eval_0_shot(input_data_path, output_data_path,stop_strings=[]) # only for deepseek-r1
        SingleToolEvaluation.calculate_score(output_data_path)
        # responseeval = ResponseEvaluation(toolset, model)
        # responseeval.eval(output_data_path)

        # 多工具评测
        input_data_path = "data/retrieved/multiple_mat_test_NVEmbed.jsonl"
        output_data_path = f"result_NV/mat_multiple_test_NV_{model_name}.jsonl"
        multieval = MultiToolEvaluation(toolset, model)
        multieval.eval_0_shot(input_data_path, output_data_path, stop_strings=[])
        MultiToolEvaluation.calculate_score(output_data_path)
        # responseeval = ResponseEvaluation(toolset, model)
        # responseeval.eval(output_data_path)

if __name__ == "__main__":
    main()
