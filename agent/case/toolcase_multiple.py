import os
import json
import random

from utils import match_given_pattern
from utils.file_io import read_JSON, write_JSON
from tools.toolpool import ToolPool
from tools.toolparam import generate_random_string
from model.GPT import call_openai_api
from prompt import prompt_candidate_tool_selection, prompt_tool_calling_chain, prompt_chain_to_query

# candidate_tool_path = "data/multiple/method_SealTools/candidate_tool.jsonl"
# tool_calling_chain_path = "data/multiple/method_SealTools/tool_calling_chain.json"
# multiple_case_path = "data/multiple/method_SealTools/multiple_case.jsonl"

class MultiToolCaseGeneration:
    def __init__(self, toolset):
        self.toolset = toolset
        dataset_path = "data/multiple/method_SealTools/"
        os.makedirs(dataset_path,exist_ok=True)

    # 将工具池混合
    def mix_all_toolpools(self, exception_list = []):
        all_tools_list = []
        for tooldir in self.toolset.tool_package_list:
            for toolname in self.toolset(tooldir).toolinfo:
                if f"{tooldir}/{toolname}" not in exception_list:
                    tooldescription = self.toolset(tooldir).toolinfo[toolname]["description"]
                    all_tools_list.append([tooldir, toolname, tooldescription])
        self.all_tools_list = all_tools_list
    
    # STEP 1：GPT挑选工具组合
    def sample_tools(self, max_num=200):
        sample_num = len(self.all_tools_list) if max_num > len(self.all_tools_list) else max_num
        sampled_list = random.sample(self.all_tools_list, sample_num)
        random.shuffle(sampled_list)
        sampled_dict = {}
        for sampled_tool in sampled_list:
            sampled_dict[sampled_tool[0]+"/"+sampled_tool[1]] = sampled_tool[2]
        return sampled_dict
    def choose_candidate_tools(self, attempt_times=1000, candidate_tool_path: str="data/multiple/method_SealTools/candidate_tool.jsonl"):
        if os.path.exists(candidate_tool_path):
            candidate_tool_list = read_JSON(candidate_tool_path)
        else:
            candidate_tool_list = []
        for _ in range(attempt_times):
            sampled_tools_dict = self.sample_tools()
            input_text = prompt_candidate_tool_selection.format(str(sampled_tools_dict))
            output_text = call_openai_api(input_text)
            print("GPT input:", input_text)
            print("GPT output:", output_text)
            try:
                check_flag = 1
                # selected_tools = eval(match_given_pattern(output_text.split("Tool calling chain")[-1].split("Query")[0].replace("\n",""),r'\[.+\]'))
                # task_description = eval(match_given_pattern(output_text.split("Query")[1].replace("\n",""),r'\[.+\]'))[0]
                selected_tools = eval(match_given_pattern(output_text.split("reordered_calling_chain =")[-1].replace("\n",""),r'\[.+\]'))
                task_description = eval(match_given_pattern(output_text.split("query =")[1].split("reordered_calling_chain =")[0].replace("\n",""),r'\[.+\]'))[0]
                for tool in selected_tools:
                    if tool not in sampled_tools_dict:
                        check_flag = 0
                        break
                candidate_tool_set = [set(tool[0]) for tool in candidate_tool_list]
                selected_tools_no_duplicate = []
                for tool in selected_tools:
                    if tool not in selected_tools_no_duplicate:
                        selected_tools_no_duplicate.append(tool)
                selected_tools = selected_tools_no_duplicate

                if check_flag:
                    if set(selected_tools) in candidate_tool_set:
                        candidate_tool_list[candidate_tool_set.index(set(selected_tools))][1].append(task_description)
                    else:
                        candidate_tool_list.append([selected_tools,[task_description]])

                    print("** selected tools: ", selected_tools)
                    print("** task description: ", task_description)

                    write_JSON(candidate_tool_path, sorted(candidate_tool_list, reverse=True, key=lambda x: len(x[1])))
            except Exception as e:
                print("choose candidate tools 后处理失败：", e)

    # STEP 2：由工具组合生成工具调用链
    def generate_tool_calling_chain(self,
                                    calling_dir = "tool_calling_result",
                                    candidate_tool_path = "data/multiple/method_SealTools/candidate_tool.jsonl",
                                    tool_calling_chain_path = "data/multiple/method_SealTools/tool_calling_chain.jsonl",
                                    execution = True):
        expected_num = 3000
        max_single_num = 1 # 5 # 40
        def get_tool_examples(tool_full_name, calling_dir = "tool_calling_chain"):
            tooldir, toolname = tool_full_name.split('/')[0], tool_full_name.split('/')[1]
            tool_calling_dataset = read_JSON(f"data/{calling_dir}/{tooldir}/{toolname}.jsonl")
            tool_examples = random.sample(tool_calling_dataset, min(5,len(tool_calling_dataset)))
            return [example[0] for example in tool_examples]
        def get_detailed_tool_information(tool_information_dict):
            tool_information_text = ""
            for tool_full_name in tool_information_dict:
                tool_information_text += f"Name: {tool_full_name}\n"
                tool_information_text += "Description: "+tool_information_dict[tool_full_name]["description"]+"\n"
                tool_information_text += "Calling examples: "+str(tool_information_dict[tool_full_name]["tool_examples"])[1:-1]+"\n"
            return tool_information_text[:-1]
        def get_tool_calling_chain_text(calling_chain_attempt, execution = True):
            tool_calling_chain_text = ""
            for tool_calling in calling_chain_attempt:
                tool_calling_chain_text += f"Calling: {tool_calling[0]}[{str(tool_calling[1])[1:-1]}]\n"
                if execution:
                    tool_calling_chain_text += "Return: "+(str(tool_calling[2]) if not isinstance(tool_calling[2], str) else '"'+str(tool_calling[2])+'"')+'\n'
                else:
                    tool_calling_chain_text += "Return: 【Not Calling Yet】\n"
            return tool_calling_chain_text[:-1]

        candidate_tool_dataset = read_JSON(candidate_tool_path)
        strict_counter = sum([len(data[1]) if len(data[1])>1 else 0 for data in candidate_tool_dataset])

        calling_chain_dataset = []
        calling_chain_counter = {}
        if os.path.exists(tool_calling_chain_path):
            calling_chain_dataset = read_JSON(tool_calling_chain_path)
            for data in calling_chain_dataset:
                if data[0] not in calling_chain_counter:
                    calling_chain_counter[data[0]] = 1
                else:
                    calling_chain_counter[data[0]] += 1

        # 遍历上一步生成的候选工具集
        for candidate_tool in candidate_tool_dataset:
            if len(candidate_tool[1]) >= 1: # > 1 确保上一步生成candidate_tool的合理性(多次生成)
                # 初始化Prompt内相关信息
                tool_information_dict = {}
                for tool in candidate_tool[0]:
                    tool_information_dict[tool] = {}
                    tool_information_dict[tool]["tooldir"] = tool.split('/')[0]
                    tool_information_dict[tool]["toolname"] = tool.split('/')[1]
                    tool_information_dict[tool]["description"] = self.toolset(tool.split('/')[0]).toolinfo[tool.split('/')[1]]["description"]

                failing_flag = 0
                if str(candidate_tool[0]) not in calling_chain_counter:
                    calling_chain_counter[str(candidate_tool[0])] = 0
                while (calling_chain_counter[str(candidate_tool[0])] < min(int(expected_num*len(candidate_tool[1])/strict_counter)+1, max_single_num)) and (failing_flag<=3):
                    calling_chain_attempt = []
                    attempt_times = 0
                    while ((len(calling_chain_attempt)<len(candidate_tool[0])) and (attempt_times < 10)):
                        attempt_times += 1
                        query = random.choice(candidate_tool[1])
                        # 每次尝试重新给tool选择调用例子
                        for tool in tool_information_dict:
                            tool_information_dict[tool]["tool_examples"] = get_tool_examples(tool, calling_dir = calling_dir)
                        tool_information_text = get_detailed_tool_information(tool_information_dict)
                        tool_calling_chain_text = get_tool_calling_chain_text(calling_chain_attempt, execution=execution)

                        try:
                            input_text = prompt_tool_calling_chain.format(query, str(set(candidate_tool[0])), tool_information_text, tool_calling_chain_text)
                            output_text = call_openai_api(input_text)
                            print("* * * * * *")
                            print("GPT Input : ", input_text)
                            print("GPT Output: ", output_text)
                            print("* * * * * *")
                            generated_calling = output_text.split("Calling:")[1].split("\n")[0].strip().strip("`")
                            generated_tool = generated_calling.split("[")[0].strip()
                            generated_parameter_list = eval(match_given_pattern(generated_calling, r'\[.+\]'))
                            print("Generated calling: ", generated_calling)
                            print("Generated tool: ", generated_tool)
                            print("Generated parameters: ", generated_parameter_list)
                            if generated_tool not in candidate_tool[0]: # 判断工具名是否有误
                                failing_flag += 1
                                continue
                            if generated_tool not in [tool[0] for tool in calling_chain_attempt]: # 判断该工具是否已经在链中调用过
                                if execution:
                                    result = self.toolset(generated_tool.split('/')[0]).call_tool(generated_tool.split('/')[1], generated_parameter_list)
                                    print("Generated result: ", result)
                                else:
                                    result = None
                                    if "example" in str(generated_parameter_list):
                                        failing_flag += 1
                                        continue
                                calling_chain_attempt.append([generated_tool, generated_parameter_list, result])
                        except Exception as e:
                            failing_flag += 1
                            print("生成工具调用链报错：",e)
                            print("尝试计数：", attempt_times)
                            print("错误计数：",failing_flag)

                        if len(calling_chain_attempt) == len(candidate_tool[0]):
                            output_flag = 1
                            for calling_idx in range(len(calling_chain_attempt)):
                                if calling_chain_attempt[calling_idx][0] != candidate_tool[0][calling_idx]:
                                    failing_flag += 1
                                    output_flag = 0
                                    print("工具调用顺序有误！")
                                    break
                                package_name = candidate_tool[0][calling_idx].split('/')[0]
                                tool_name = candidate_tool[0][calling_idx].split('/')[1]
                                if candidate_tool[0][calling_idx] in ["IO_Operations/get_structure_by_material_id", "Basic_Functionality/read_structure_from_file"]:
                                    pass
                                elif len(calling_chain_attempt[calling_idx][1]) != len(self.toolset(package_name).toolinfo[tool_name]["parameters"]["properties"]):
                                    failing_flag += 1
                                    output_flag = 0
                                    print("工具参数填充有误！")
                                    break
                            if not output_flag:
                                break

                            if [str(candidate_tool[0]),calling_chain_attempt] not in calling_chain_dataset:
                                calling_chain_dataset.append([str(candidate_tool[0]),calling_chain_attempt])
                                calling_chain_counter[str(candidate_tool[0])] += 1
                                write_JSON(tool_calling_chain_path, [str(candidate_tool[0]),calling_chain_attempt], mode="a")
                                failing_flag = 0
                                break

    def sp_process(self, data_path):
        dataset = read_JSON(data_path)
        for data_idx in range(len(dataset)):
            key = generate_random_string(length=12)
            for calling_idx in range(len(dataset[data_idx][1])):
                if dataset[data_idx][1][calling_idx][0] == "IO_Operations/get_structure_by_material_id":
                    dataset[data_idx][1][calling_idx][1]  = ["sk-"+key, str(random.randint(1,1000))]
                elif dataset[data_idx][1][calling_idx][0] == "Basic_Functionality/read_structure_from_file":
                    if len(dataset[data_idx][1][calling_idx][1]) == 1:
                        dataset[data_idx][1][calling_idx][1] = [dataset[data_idx][1][calling_idx][1][0], dataset[data_idx][1][calling_idx][1][0].split(".")[0]+'.pkl']
        write_JSON(data_path, dataset)

    def postprocess_calling_chain_format(self, raw_tool_calling_chain_path, new_tool_calling_chain_path):
        raw_dataset = read_JSON(raw_tool_calling_chain_path)
        new_dataset = {}

        for raw_data in raw_dataset:
            if raw_data[0] not in new_dataset:
                new_dataset[raw_data[0]] = []
            new_dataset[raw_data[0]].append(raw_data[1])
        write_JSON(new_tool_calling_chain_path, new_dataset, indent=4)

    # STEP 3：由工具调用链生成具体问题
    def generate_query(self,
                       tool_calling_chain_path = "data/multiple/method_SealTools/tool_calling_chain.json",
                       multiple_case_path = "data/multiple/method_SealTools/multiple_case.jsonl"):
        def get_tool_information_text(tool_info_list):
            tool_information_text = ""
            for tool_info in tool_info_list:
                tool_information_text += str(tool_info)+"\n"
            return tool_information_text[:-1]
        def get_tool_calling_chain_text(tool_chain_case):
            tool_calling_chain_text = ""
            for tool_calling in tool_chain_case:
                tool_calling_chain_text += f"Calling: {tool_calling[0]}[{str(tool_calling[1])[1:-1]}]\n"
                tool_calling_chain_text += "Return: "+(str(tool_calling[2]) if not isinstance(tool_calling[2], str) else '"'+str(tool_calling[2])+'"')+'\n'
            return tool_calling_chain_text[:-1]

        calling_chain_dataset = read_JSON(tool_calling_chain_path)
        generated_calling_chain = []
        if os.path.exists(multiple_case_path):
            multiple_case_dataset = read_JSON(multiple_case_path)
            for multiple_case in multiple_case_dataset:
                generated_calling_chain.append(multiple_case[1])
        else:
            multiple_case_dataset = []
        for tool_chain_template in calling_chain_dataset:
            for tool_chain_case in calling_chain_dataset[tool_chain_template]:
                if tool_chain_case not in generated_calling_chain:
                    print("Processing tool chain:", tool_chain_case)
                    new_case = [None, tool_chain_case]

                    try:
                        tool_information_text = get_tool_information_text([self.toolset(tool_calling[0].split("/")[0]).toolinfo[tool_calling[0].split("/")[1]] for tool_calling in tool_chain_case])
                        tool_calling_chain_text = get_tool_calling_chain_text(tool_chain_case)

                        input_text = prompt_chain_to_query().format(tool_information_text, tool_calling_chain_text)
                        output_text = call_openai_api(input_text)
                        print("* * * * * *")
                        print("GPT Input : ", input_text)
                        print("GPT Output: ", output_text)
                        print("* * * * * *")
                        generated_query = output_text.split("Query:")[1].split("\n")[0].strip()
                        if generated_query != "":
                            print("Generated query: ", generated_query)
                            new_case[0] = generated_query
                            multiple_case_dataset.append(new_case)
                            generated_calling_chain.append(tool_chain_case)
                            write_JSON(multiple_case_path, multiple_case_dataset)

                    except Exception as e:
                        print("生成多工具用户Query报错：", e)
        pass

if __name__ == "__main__":
    tooldir_list = ["chem_lib", "cactus", "chemcrow", "chemistrytools"]
    toolpool_list = [ToolPool(tooldir) for tooldir in tooldir_list]
    case_genera = MultiToolCaseGeneration(toolpool_list)
    case_genera.choose_candidate_tools()