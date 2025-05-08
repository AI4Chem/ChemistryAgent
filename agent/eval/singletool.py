import os
import random

from utils.file_io import read_JSON, write_JSON, is_json_serializable
from prompt import prompt_eval_0_shot, prompt_answer
from tools.parse_tool_calling import parse_function_calling

class SingleToolEvaluation:
    def __init__(self, toolset, model):
        self.model = model
        self.toolset = toolset
        self.eval_0_shot_prompt = prompt_eval_0_shot["instruction"]+prompt_eval_0_shot["input"]
        self.answer_prompt = prompt_answer

    def eval_0_shot(self, input_path, output_path, execution=True, stop_strings = ["\n"]):
        eval_dataset = read_JSON(input_path)
        if os.path.exists(output_path):
            output_dataset = read_JSON(output_path)
        else:
            output_dataset = []
        for eval_data in eval_dataset:
            if eval_data["query"] not in [output_data["query"] for output_data in output_dataset]:
                query = eval_data["query"]
                candidate_tools = eval_data["candidate_tools"]

                tool_information = self.toolset.get_tool_information_text_list(candidate_tools)

                input_text = self.eval_0_shot_prompt.format(query,"","\n".join(tool_information))
                print("INPUT:\n"+input_text)
                output_text = self.model.answer(input_text, stop_strings=stop_strings)
                print("OUTPUT:\n"+output_text)
                try:
                    calling_command = output_text.split("Tool Calling:")[1] if "Tool Calling:" in output_text else output_text
                    calling_command = calling_command.split("\n")[0].strip()
                    if execution:
                        execution_result = self.toolset.exec_tool_calling(calling_command)
                    else:
                        execution_result = "【Not Calling】"
                except Exception as e:
                    execution_result = "【Error: "+str(e)+"】"
                print("RESULT:\n"+str(execution_result))
                output_data = {}
                output_data["query"] = query
                # output_data["input_text"] = input_text
                output_data["pred_text"] = output_text
                if is_json_serializable(execution_result):
                    output_data["pred_exec_result"] = [[calling_command, execution_result]]
                else:
                    output_data["pred_exec_result"] = [[calling_command, str(execution_result)]]
                output_data["gold"] = eval_data["calling_chain"]
                output_dataset.append(output_data)
                write_JSON(output_path, output_data, mode='a')
        write_JSON(output_path, output_dataset)

    @staticmethod
    def calculate_score(eval_data_path, execution=True):
        eval_dataset = read_JSON(eval_data_path)

        # * Tool
        correct_counter = 0
        all_counter = 0
        parse_error_counter = 0
        for eval_data in eval_dataset:
            pred_tool_name_list = []
            for pred_calling in eval_data["pred_exec_result"]:
                if pred_calling[0] == 'None': # 跳过未调用的轮次
                    continue
                all_counter += 1
                try:
                    func_name, args, kwargs = parse_function_calling(pred_calling[0])
                    pred_tool_name_list.append(func_name)
                except:
                    parse_error_counter += 1
                    # print("【Failed Parsing】:",pred_calling)
            gold_tool_name_list = []
            for gold_calling in eval_data["gold"]:
                gold_tool_name_list.append(gold_calling["tool"])
            
            if pred_tool_name_list != []:
                if pred_tool_name_list[0] == gold_tool_name_list[0]:
                    correct_counter += 1

        print("【file】:",eval_data_path)
        print(f"【Parse】:{str(all_counter-parse_error_counter)}/{str(all_counter)}={str((all_counter-parse_error_counter)/all_counter)}")
        print("【Tool ACC】: ",str(correct_counter/len(eval_dataset)))

        # * Parameter
        correct_counter = 0
        pred_counter = 0
        gold_counter = 0
        parse_error_counter = 0
        for eval_data in eval_dataset:
            pred_tool_name_list = []
            pred_tool_param_list = []
            for pred_calling in eval_data["pred_exec_result"]:
                if pred_calling[0] == 'None': # 跳过未调用的轮次
                    continue
                try:
                    func_name, args, kwargs = parse_function_calling(pred_calling[0])
                    pred_tool_name_list.append(func_name)
                    param_list = []
                    param_list.extend(args)
                    param_list.extend([kwargs[key] for key in kwargs])
                    pred_tool_param_list.append(param_list)
                except:
                    parse_error_counter += 1
            gold_tool_name_list = []
            gold_tool_param_list = []
            for gold_calling in eval_data["gold"]:
                gold_tool_name_list.append(gold_calling["tool"])
                gold_tool_param_list.append(gold_calling["params"])

            gold_counter += len(gold_tool_param_list[0])
            if pred_tool_param_list != []:
                pred_counter += len(pred_tool_param_list[0])
                correct_counter += sum(1 for a, b in zip(pred_tool_param_list[0], gold_tool_param_list[0]) if a == b)
            # print(pred_tool_param_list[pred_index], gold_tool_param_list[gold_index])


        print("【file】:",eval_data_path)
        print("【Param P】: ",str(correct_counter/pred_counter))
        print("【Param R】: ",str(correct_counter/gold_counter))
        print("【Param F1】:",str(2*correct_counter/(pred_counter+gold_counter)))

        # * Return
        if execution:
            counter = 0
            for eval_data in eval_dataset:
                pred = eval_data["pred_exec_result"][0][1]
                gold = eval_data["gold"][0]["return"]
                if type(pred) is float and type(gold) is float:
                    if round(pred, 4) == round(gold, 4):
                        counter += 1
                elif pred == gold:
                    counter += 1
            print("【file】:",eval_data_path)
            print("【Return ACC】:",str(counter/len(eval_dataset)))

    def generate_without_tool_calling(self, input_path, output_path):
        eval_dataset = read_JSON(input_path)
        if os.path.exists(output_path):
            output_dataset = read_JSON(output_path)
        else:
            output_dataset = []
        for eval_data in eval_dataset:
            if eval_data["query"] not in [output_data["query"] for output_data in output_dataset]:
                output_data = eval_data

                query = eval_data["query"]
                input_text = self.answer_prompt.format(query)
                print("【INPUT】:"+input_text)
                output_text = self.model.answer(input_text, stop_strings=["\n"])
                print("【OUTPUT】:"+output_text)
                response = output_text.split("Question:")[0].strip()
                
                output_data["response"] = response
                output_dataset.append(output_data)
                write_JSON(output_path, output_data, mode='a')
        write_JSON(output_path, output_dataset)