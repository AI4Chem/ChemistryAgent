import sys
sys.path.append('./src/')

import os
from tools.toolset import ToolSet

from utils.file_io import read_JSON, write_JSON
from model.GPT import call_openai_api
from prompt import prompt_response_generation

class ResponseGeneration:
    def __init__(self, toolset):
        self.toolset = toolset

    def generate_response(self, input_dataset_path):
        output_dataset_path = input_dataset_path.replace(".jsonl", "_response.jsonl")
        input_dataset = read_JSON(input_dataset_path)
        if os.path.exists(output_dataset_path):
            output_dataset = read_JSON(output_dataset_path)
        else:
            output_dataset = []
        finished_data_id = [data["id"] for data in output_dataset]
        for input_data in input_dataset:
            if input_data["id"] in finished_data_id:
                continue
            output_data = {}

            calling_text_list = []
            for calling in input_data["calling_chain"]:
                tool_name = calling["tool"]
                tool_info = self.toolset(tool_name.split("/")[0]).toolinfo[tool_name.split("/")[1]]
                tool_description = tool_info["description"]
                tool_param_list = [param for param in tool_info["parameters"]["properties"]]

                calling_text_list.append(f"Tool Name: {tool_name}, Tool Description: {tool_description}")
                calling_text_list.append("Input: "+tool_name+"()"+", ".join([tool_param_list[idx]+"="+(str(calling["params"][idx]) if calling["params"][idx] is not str else '"'+calling["params"][idx]+'"') for idx in range(len(tool_param_list)) ])+")")
                calling_text_list.append("Output: "+str(calling["return"]))

            tool_calling_text = "\n".join(calling_text_list)

            prompt = prompt_response_generation.format(input_data["query"], tool_calling_text)
            output = call_openai_api(prompt)
            print("【GPT input 】:", prompt)
            print("【GPT output】:", output)
            response = output.split("Question:")[0].strip()

            output_data["id"] = input_data["id"]
            output_data["query"] = input_data["query"]
            output_data["response"] = response
            output_data["calling_chain"] = input_data["calling_chain"]
            output_dataset.append(output_data)
            write_JSON(output_dataset_path, output_data, mode="a")
        # write_JSON(output_dataset_path, output_dataset)


if __name__ == "__main__":
    toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])
    case_genera = ResponseGeneration(toolset)
    for input_dataset_path in ["data/single/single_test.jsonl", "data/multiple/method_SealTools/multiple_test.jsonl"]:
        case_genera.generate_response(input_dataset_path)