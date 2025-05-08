from tools import ToolPool, ToolSet
from tools.toolparam import generate_random_string
from model import GPT

from utils.file_io import read_JSON

# 评测
from eval import SingleToolEvaluation, MultiToolEvaluation

import requests
import json

def main():
    # toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])
    # model = GPT()

    # # output_data_path = "result/single_test_DPR_chemllm_raw.jsonl"
    # # output_data_path = "result/single_test_DPR_chemllm_lora.jsonl"
    # output_data_path = "result/single_test_DPR_internlm.jsonl"
    # singleeval = SingleToolEvaluation(toolset, model)
    # singleeval.calculate_score(output_data_path)


    # # output_data_path = "result/multiple_test_DPR_chemllm_raw.jsonl"
    # # # output_data_path = "result/multiple_test_DPR_chemllm_lora.jsonl"
    # # multieval = MultiToolEvaluation(toolset, model)
    # # multieval.calculate_score(output_data_path)

    # 设置URL
    url = "http://172.16.4.2:10001/call_tool/"
    # 设置请求头
    headers = {"Content-Type": "application/json"}
    # 请求体数据
    data = {
        "command": "chem_lib/analyze_combustion(CO2=1408.0, H2O=1036.8)"
        }

    # 发送POST请求
    response = requests.post(url, headers=headers, data=json.dumps(data))

    # 打印响应
    print(response.status_code)
    print(response.json())  # 如果响应是JSON格式的数据

if __name__ == "__main__":
    # main()

    # dataset = read_JSON("./result/multiple_test_DPR_gpt-4o-mini.jsonl")
    # counter = 0
    # for data in dataset:
    #     # counter += len(data["pred_exec_result"])
    #     counter += len(data["gold"])
    # print(counter)

    print(generate_random_string())