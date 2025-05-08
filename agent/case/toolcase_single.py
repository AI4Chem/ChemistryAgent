import os
import json

from utils.file_io import read_JSON, write_JSON
from model.GPT import call_openai_api
from prompt import single_case_prompt

class SingleToolCaseGeneration:
    def __init__(self, toolset):
        self.toolset = toolset

    def generate_single_case_with_calling(self,calling_dir: str="tool_calling_result", case_dir: str="single/single_case"):
        def generate_query(toolinfo, input_param, output_result):
            tool_information = json.dumps(toolinfo, indent=4)
            tool_calling = toolinfo["name"]+"("+", ".join([str(param) if type(param) is not str else '"'+param+'"' for param in input_param])+")"
            prompt = single_case_prompt().format(tool_information, tool_calling)
            output = call_openai_api(prompt)
            print("GPT input :", prompt)
            print("GPT output:", output)
            return output
        for tool_dir in self.toolset.tool_package_list:
            toolpool = self.toolset(tool_dir)
            for tool_name in toolpool.toolinfo:
                tool_calling_dataset = read_JSON("data/"+calling_dir+"/"+tool_dir+"/"+tool_name+".jsonl")
                single_case_dir = "data/"+case_dir+"/"+tool_dir
                os.makedirs(single_case_dir, exist_ok=True)
                single_case_path = os.path.join(single_case_dir, tool_name+".jsonl")
                print(single_case_path)
                if os.path.exists(single_case_path):
                    case_dataset = read_JSON(single_case_path)
                else:
                    case_dataset = []
                input_set = list()
                for case in case_dataset:
                    input_set.append(case[1])
                for [input_param, output_result] in tool_calling_dataset:
                        if input_param not in input_set:
                            query = generate_query(toolpool.toolinfo[tool_name], input_param, output_result)
                            if type(query) is not Exception:
                                try:
                                    user_query = query.split("User Query:")[1].strip()
                                except:
                                    user_query = query
                                case_dataset.append([user_query, input_param, output_result])
                                write_JSON(single_case_path, case_dataset)


if __name__ == "__main__":
    from tools.toolpool import ToolPool
    tooldir_list = ["chem_lib", "cactus", "chemcrow", "chemistrytools"]
    toolpool_list = [ToolPool(tooldir) for tooldir in tooldir_list]
    case_genera = SingleToolCaseGeneration(toolpool_list)
    case_genera.generate_single_case_with_calling()
