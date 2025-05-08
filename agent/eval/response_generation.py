import os
import random

from utils.file_io import read_JSON, write_JSON
from prompt import prompt_response_generation, prompt_response_evaluation

class ResponseEvaluation:
    def __init__(self, toolset, model):
        self.model = model
        self.toolset = toolset
        self.response_generation_prompt = prompt_response_generation
        self.response_evaluation_prompt = prompt_response_evaluation

    def eval(self, input_dataset_path):
        output_dataset_path = input_dataset_path.replace(".jsonl", "_response.jsonl")
        eval_dataset = read_JSON(input_dataset_path)
        if os.path.exists(output_dataset_path):
            output_dataset = read_JSON(output_dataset_path)
        else:
            output_dataset = []
        for eval_data in eval_dataset:
            if eval_data["query"] not in [output_data["query"] for output_data in output_dataset]:

                output_data = eval_data

                calling_text_list = []
                for calling in eval_data["pred_exec_result"]:
                    tool_name = calling[0].split("(")[0]
                    try:
                        tool_info = self.toolset(tool_name.split("/")[0]).toolinfo[tool_name.split("/")[1]]
                        tool_description = tool_info["description"]
                        tool_param_list = [param for param in tool_info["parameters"]["properties"]]

                        # calling_text_list.append(f"Tool Name: {tool_name}, Tool Description: {tool_description}")
                        calling_text_list.append("Calling: "+calling[0])
                        calling_text_list.append("Return: "+str(calling[1]))
                    except:
                        pass

                tool_calling_text = "\n".join(calling_text_list) if calling_text_list != [] else "None"

                prompt = self.response_generation_prompt.format(eval_data["query"], tool_calling_text)
                print("【Input 】:", prompt)
                output = self.model.answer(prompt, stop_strings=["\n"])
                print("【Output】:", output)
                response = output.split("Question:")[0].strip()

                output_data["response"] = response
                output_dataset.append(output_data)
                write_JSON(output_dataset_path, output_data, mode='a')
        write_JSON(output_dataset_path, output_dataset)

    def score(self, input_dataset_path, gold_dataset_path):
        output_dataset_path = input_dataset_path.replace("response.jsonl", "score.jsonl")
        input_dataset = read_JSON(input_dataset_path)
        gold_dataset = read_JSON(gold_dataset_path)
        gold_hash_dict = {data["query"]:idx for idx, data in enumerate(gold_dataset)}

        if os.path.exists(output_dataset_path):
            output_dataset = read_JSON(output_dataset_path)
        else:
            output_dataset = []
        for input_data in input_dataset:
            if input_data["query"] not in [output_data["query"] for output_data in output_dataset]:
                output_data = input_data

                query = input_data["query"]
                gold_answer = gold_dataset[gold_hash_dict[query]]["response"]
                pred_answer = input_data["response"]

                prompt = self.response_evaluation_prompt.format(query, gold_answer, pred_answer)
                print("【Input 】:", prompt)
                output = self.model.answer(prompt, stop_strings=["\n"])
                print("【Output】:", output)
                response = output.split("Question:")[0].strip()

                score = None
                if "Yes" in response:
                    score = "Yes"
                elif "No" in response:
                    score = "No"

                output_data["score"] = score
                output_dataset.append(output_data)
                write_JSON(output_dataset_path, output_data, mode='a')
        write_JSON(output_dataset_path, output_dataset)

    @staticmethod
    def count(input_dataset_path):
        dataset = read_JSON(input_dataset_path)
        counter = {"Yes":0, "No":0, "None":0}
        for data in dataset:
            if data["score"] == "No":
                counter["No"] += 1
            elif data["score"] == "Yes":
                counter["Yes"] += 1
            elif data["score"] == None:
                counter["None"] += 1
        counter["Total"] = len(dataset)
        if counter["Total"] != counter["Yes"] + counter["No"] + counter["None"]:
            print("Error: Total != Yes + No + None")
        counter["Pass"] = counter["Yes"] / counter["Total"]
        return counter
