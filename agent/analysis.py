import os

# from tools import ToolSet
# from tools.toolinfo import basic_tool_package, matgen_tool_package

from utils.file_io import read_JSON, write_JSON

def count_all_single_case(tooldir_list):
    counter = 0
    for tooldir in tooldir_list:
        case_dir = os.path.join("data/single/single_case/", tooldir)
        case_files = os.listdir(case_dir)
        for case_file in case_files:
            case_path = os.path.join(case_dir, case_file)
            counter += len(read_JSON(case_path))
    return counter

def count_all_multi_tool_calling_chain(dataset_path):
    dataset = read_JSON(dataset_path)
    return sum([len(dataset[template]) for template in dataset])

def count_all_multi_case(dataset_path):
    dataset = read_JSON(dataset_path)
    return len(dataset), sum([len(data[1]) for data in dataset])

def count_per_dataset(dataset_path):
    dataset = read_JSON(dataset_path)
    return len(dataset), sum([len(data["calling_chain"]) for data in dataset]), max([len(data["calling_chain"]) for data in dataset])

def count_candidate_tools(dataset_path):
    dataset = read_JSON(dataset_path)
    tools_counter = {}
    for data in dataset:
        for tool in data[0]:
            if tool not in tools_counter:
                tools_counter[tool] = 1
            else:
                tools_counter[tool] += 1
    print(len(tools_counter))
    return dict(sorted(tools_counter.items(), key=lambda item: item[1], reverse=True))

# ------------------------

from collections import Counter, defaultdict
import json

def process_tool_file_optimized(input_file, output_file, max_tool_count):
    # 统计所有工具的全局出现频率
    tool_counter = Counter()
    line_tools = []

    # 读取文件并解析每行的工具列表
    lines = read_JSON(input_file)

    # 解析 JSON 数据
    for line in lines:
        tools = line[0]  # 工具列表
        line_tools.append((tools, line))  # 保存工具和原始数据
        tool_counter.update(tools)

    # 行优先级队列：按稀有工具的最小频率排序
    def line_priority(line):
        tools, _ = line
        return min(tool_counter[tool] for tool in tools)

    # 筛选结果
    filtered_lines = []
    current_counter = Counter()

    while line_tools:
        # 按优先级筛选当前行
        line_tools.sort(key=line_priority, reverse=True)
        tools, data = line_tools.pop(0)

        # 检查是否可以加入结果集
        if all(current_counter[tool] < max_tool_count for tool in tools):
            filtered_lines.append(data)
            current_counter.update(tools)
        else:
            # 更新全局统计，降低稀有工具误删的风险
            for tool in tools:
                tool_counter[tool] -= 1

    # 写回结果文件
    write_JSON(output_file, filtered_lines)

def print_func_info(toolset, output_path):
    dataset = {}
    for package_name in toolset.tool_pool_dict:
        for tool_name in toolset.tool_pool_dict[package_name].toolinfo:
            full_tool_name = package_name+"/"+tool_name
            dataset[full_tool_name] = toolset.tool_pool_dict[package_name].toolinfo[tool_name]["description"]
            # dataset[full_tool_name]["Description"] = toolset.tool_pool_dict[package_name].toolinfo[tool_name]["description"]
    write_JSON(output_path, dataset, indent=4)

def tool_calling_per_data(input_file):
    input_dataset = read_JSON(input_file)
    counter = {}
    for data in input_dataset:
        number = len(data[1])
        if number not in counter:
            counter[number] = 1
        else:
            counter[number] += 1
    print(counter)

if __name__ == "__main__":
    # tooldir_list = ["chem_lib", "chemcrow", "cactus", "chemistrytools"]
    # print("单工具用例数量：", count_all_single_case(tooldir_list))

    # multi_case_path = "data/multiple/method_SealTools/tool_calling_chain.json"
    # print("多工具调用链数量：", count_all_multi_tool_calling_chain(multi_case_path))

    # multi_case_path = "data/multiple/method_SealTools/multiple_case.jsonl"
    # print("多工具用例数量/调用次数：", count_all_multi_case(multi_case_path))

    # dataset = "data/single/single_train.jsonl"
    # print(count_per_dataset(dataset))
    # dataset = "data/multiple/method_SealTools/multiple_train.jsonl"
    # print(count_per_dataset(dataset))

    # dataset = "data/pymatgen/multiple/method_SealTools/candidate_tool.jsonl"
    # print(count_candidate_tools(dataset))

    # input_file = 'data/pymatgen/multiple/method_SealTools/candidate_tool.jsonl'
    # output_file = 'data/pymatgen/multiple/method_SealTools/candidate_tool_p.jsonl'
    # max_tool_count = 50
    # process_tool_file_optimized(input_file, output_file, max_tool_count)

    # raw_dataset = read_JSON("data/pymatgen/multiple/method_SealTools/backup/tool_calling_chain_2.json")
    # processed_dataset = []
    # for template in raw_dataset:
    #     for case in raw_dataset[template]:
    #         processed_dataset.append([template, case])
    # write_JSON("data/pymatgen/multiple/method_SealTools/tool_calling_chain.jsonl", processed_dataset)

    # toolset_list = []
    # toolset_list.append(ToolSet(basic_tool_package))
    # toolset_list.append(ToolSet(matgen_tool_package, toolpackage_path="src/toolpool/pymatgen_common_PROCESSED"))
    # for toolset in toolset_list:
    #     for tool_package in toolset.tool_pool_dict:
    #         toolpool = toolset(tool_package)
    #         for tool_name in toolpool.toolinfo:
    #             print("Name: ",tool_name)
    #             print("Description: ", toolpool.toolinfo[tool_name]["description"])

    # toolset = ToolSet([
    #     "chem_lib", 
    #     "cactus", 
    #     "chemcrow", 
    #     "chemistrytools",
    # ])
    # # toolset = ToolSet(matgen_tool_package, toolpackage_path="src/toolpool/pymatgen_common_PROCESSED")
    # tool_name_list = toolset.get_all_tool_names_list()
    # print(len(tool_name_list))
    
    # toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])
    # print_func_info(toolset, "./chemistry.json")
    # toolset = ToolSet(matgen_tool_package, 
    #               toolpackage_path="src/toolpool/pymatgen_common_PROCESSED")
    # print_func_info(toolset, "./materials.json")

    chem_multi_case_file = "/Users/wumengsong/Code/PJLab/ChemAgent/ChemToolLearning/data/multiple/method_SealTools/multiple_case.jsonl"
    mat_multi_case_file = "/Users/wumengsong/Code/PJLab/ChemAgent/ChemToolLearning/data/pymatgen/multiple/method_SealTools/multiple_case.jsonl"
    tool_calling_per_data(chem_multi_case_file)
    tool_calling_per_data(mat_multi_case_file)
