import os
import random

from tools import ToolSet
from utils.file_io import read_JSON, write_JSON
from retrieval import arrange_tool_list_and_corpus, select_candidate_tool_list
from retrieval.model import MiniLM, NV_Embed
from retrieval.retriever import DenseRetriever
from eval import process_calling_list_text
from prompt import prompt_eval_0_shot

from model.GPT import call_openai_api

toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])

prompt_instruction=prompt_eval_0_shot["instruction"]
prompt_input=prompt_eval_0_shot["input"]

def process_raw_calling_dict(toolset, tool_calling_dict):
    tooldir, toolname = tool_calling_dict["tool"].split("/")[0],tool_calling_dict["tool"].split("/")[1]
    toolinfo = toolset(tooldir).toolinfo[toolname]
    processed_param_list = []
    tool_idx = 0
    for param_name in toolinfo["parameters"]["properties"]:
        processed_param_list.append(str(param_name)+"="+str(tool_calling_dict["params"][tool_idx]) if not isinstance(tool_calling_dict["params"][tool_idx], str) else str(param_name)+"='"+tool_calling_dict["params"][tool_idx]+"'")
        tool_idx += 1
    tool_calling_command = str(tool_calling_dict["tool"])+"("+", ".join(processed_param_list)+")"
    result = tool_calling_dict["return"]
    return [tool_calling_command, result]


def process_dataset_in_Alpaca_format(input_dataset_path_list, output_dataset_path):
    # retriever_checkpoint = "/home/bingxing2/ailab/group/ai4phys/wumengsong/all-MiniLM-L6-v2"
    # retrieval_model = MiniLM(retriever_checkpoint)
    retriever_checkpoint = "/home/bingxing2/ailab/group/ai4phys/wumengsong/NV-Embed-v2"
    retrieval_model = NV_Embed(retriever_checkpoint)
    retriever = DenseRetriever(*arrange_tool_list_and_corpus(toolset), retrieval_model)

    output_dataset = []
    for input_dataset_path in input_dataset_path_list:
        print(input_dataset_path)
        input_dataset = read_JSON(input_dataset_path)
        retriever_metric_counter = {
            "num": len(input_dataset),
            "top_1": 0,
            "top_5": 0,
            "top_10": 0
        }
        for input_data in input_dataset:
            gold_tool_list = [data["tool"] for data in input_data["calling_chain"]]
            retrieved_tool_list = retriever.retrieve(input_data["query"],top_k=10)

            if set(gold_tool_list).issubset(retrieved_tool_list[:10]):
                retriever_metric_counter["top_10"] += 1
                if set(gold_tool_list).issubset(retrieved_tool_list[:5]):
                    retriever_metric_counter["top_5"] += 1
                    if set(gold_tool_list).issubset(retrieved_tool_list[0]):
                        retriever_metric_counter["top_1"] += 1

            candidate_tool_list = select_candidate_tool_list(10, gold_tool_list, retrieved_tool_list)
            candidate_tool_information_list = toolset.get_tool_information_text_list(candidate_tool_list)

            tool_calling_list = [process_raw_calling_dict(toolset, data) for data in input_data["calling_chain"]]
            for tool_idx in range(len(tool_calling_list)):
                random.shuffle(candidate_tool_information_list)
                output_data = {}
                output_data["instruction"] = prompt_instruction
                output_data["input"] = prompt_input.format(input_data["query"], process_calling_list_text(tool_calling_list[:tool_idx]), "\n".join(candidate_tool_information_list))
                output_data["output"] = "Tool Calling: "+tool_calling_list[tool_idx][0]
                output_dataset.append(output_data)
                write_JSON(output_dataset_path, output_dataset, indent=4)
            # 构造负例，将完成工具调用后的对话状态加入，根据Prompt需返回None
            if len(tool_calling_list)>1:
                random.shuffle(candidate_tool_information_list)
                output_data = {}
                output_data["instruction"] = prompt_instruction
                output_data["input"] = prompt_input.format(input_data["query"], process_calling_list_text(tool_calling_list[:]), "\n".join(candidate_tool_information_list))
                output_data["output"] = "None"
                output_dataset.append(output_data)
                write_JSON(output_dataset_path, output_dataset, indent=4)
        print("Top 1: " + str(retriever_metric_counter["top_1"]/retriever_metric_counter["num"]))
        print("Top 5: " + str(retriever_metric_counter["top_5"]/retriever_metric_counter["num"]))
        print("Top 10: " + str(retriever_metric_counter["top_10"]/retriever_metric_counter["num"]))
        
    write_JSON(output_dataset_path, output_dataset, indent=4)

# Please generate the wrongly formatted function calling which looks right. You can add or modify the parameter name, remove the quotes for strings, change '/' for '.' and so on.

DPO_format_prompt = '''Please generate a function call that appears correct at first glance but contains subtle formatting errors. These errors may include:
- Modifying or misspelling parameter names.
- Removing or misplacing quotation marks around string values.
- Using incorrect delimiters, such as replacing / with . or vice versa.
- Altering the function's syntax or structure slightly while keeping it seemingly plausible.
- ...
You only need to output the wrong calling.

Right calling:
chemistrytools/calculate_compound_molar_mass(compound='C₂₃H₃₅O₁₅')
Wrong calling:
chemistrytools/calculate_compound_molar_mass(compound=C₂₃H₃₅O₁₅)

Right calling:
chem_lib/galvanic_cell_properties(element_electrode1='Cu', element_electrode2='H')
Wrong calling:
chem_lib/galvanic_cell_properties(element_electrode1='Cu', element_electrode2='H', temperature=300)

Right calling:
{}
Wrong calling:
'''

DPO_content_prompt = '''Please modify the parameter value to make the function calling wrong according to the given information.
You only need to output the wrong calling.

{}

Right calling:
{}
Wrong calling:
'''

def construct_DPO_dataset(input_dataset_path, output_dataset_path, mode="content"):
    input_dataset = read_JSON(input_dataset_path)
    generated_id_list = []
    if os.path.exists(output_dataset_path):
        output_dataset = read_JSON(output_dataset_path)
        for output_data in output_dataset:
            generated_id_list.append(output_data["id"])
    else:
        output_dataset = []
        write_JSON(output_dataset_path, output_dataset)
    for _ in range(1000):
        idx = random.randint(0,len(input_dataset)-1)
        while idx in generated_id_list:
            idx = random.randint(0,len(input_dataset)-1)

        if "None" == input_dataset[idx]["output"]:
            continue

        if mode == "format":
            input = DPO_format_prompt.format(input_dataset[idx]["output"][14:])
            output = call_openai_api(input)
        else:
            assert mode == "content"
            input = DPO_content_prompt.format(input_dataset[idx]["input"], input_dataset[idx]["output"][14:])
            output = call_openai_api(input)

        print("GPT input: ", input_dataset[idx]["output"][14:])
        print("GPT output: ", output)
        print()

        wrong_calling = output
        if "Wrong calling:" in wrong_calling:
            wrong_calling = wrong_calling.split("Wrong calling:")[-1].strip()

        output_data = {}
        output_data["id"] = idx
        output_data["instruction"] = input_dataset[idx]["instruction"]
        output_data["input"] = input_dataset[idx]["input"]
        output_data["chosen"] = input_dataset[idx]["output"]
        output_data["rejected"] = "Tool Calling: "+ wrong_calling
        output_dataset.append(output_data)
        generated_id_list.append(idx)
        write_JSON(output_dataset_path, output_data, mode='a')

def mix_DPO_dataset():
    def process(raw_dataset):
        dataset = []
        for raw_data in raw_dataset:
            dataset.append({
                "instruction": raw_data["instruction"],
                "input": raw_data["input"],
                "chosen": raw_data["chosen"],
                "rejected": raw_data["rejected"]
            })
        return dataset

    output_dataset_path = "data/chem_tools_nvembed_DPO_v1.json"
    content_data_amount = 2000
    format_data_amount = 1000

    dataset = []
    dataset.extend(random.sample(process(read_JSON("data/chem_tools_nvembed_DPO_content_raw.jsonl")), content_data_amount))
    dataset.extend(random.sample(process(read_JSON("data/chem_tools_nvembed_DPO_format_raw.jsonl")), format_data_amount))
    write_JSON(output_dataset_path, dataset, indent=4)

if __name__ == "__main__":
    # process_dataset_in_Alpaca_format(["data/single/single_train.jsonl","data/multiple/method_SealTools/multiple_train.jsonl"],"data/chem_tools_nvembed_v2_1.json")
    # construct_DPO_dataset("data/chem_tools_nvembed_v2_1.json","data/chem_tools_nvembed_DPO_format_raw.jsonl", mode="format")
    # construct_DPO_dataset("data/chem_tools_nvembed_v2_1.json","data/chem_tools_nvembed_DPO_content_raw.jsonl")
    mix_DPO_dataset()