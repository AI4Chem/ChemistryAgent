import os
import random

from utils.file_io import read_JSON, write_JSON, is_json_serializable
from prompt import prompt_eval_0_shot
from tools.parse_tool_calling import parse_function_calling
from . import process_calling_list_text

class MultiToolEvaluation:
    def __init__(self, toolset, model):
        self.model = model
        self.toolset = toolset
        self.eval_0_shot_prompt = prompt_eval_0_shot["instruction"]+prompt_eval_0_shot["input"]

    def eval_0_shot(self, data_set_path, output_path, execution=True, stop_strings = ["\n"]):
        eval_dataset = read_JSON(data_set_path)
        if os.path.exists(output_path):
            output_dataset = read_JSON(output_path)
        else:
            output_dataset = []
        for eval_data in eval_dataset:
            if eval_data["query"] not in [output_data["query"] for output_data in output_dataset]:
                query = eval_data["query"]
                candidate_tools = eval_data["candidate_tools"]

                tool_information = self.toolset.get_tool_information_text_list(candidate_tools)
                
                already_calling_list = []
                raw_output_text_list = []
                continue_flag = 1
                retry_times = 0
                while (continue_flag and (retry_times<=10)):
                    retry_times += 1
                    input_text = self.eval_0_shot_prompt.format(query,"\n"+process_calling_list_text(already_calling_list),"\n".join(tool_information))
                    print("INPUT:\n"+input_text)
                    output_text = self.model.answer(input_text, stop_strings=stop_strings)
                    raw_output_text_list.append(output_text)
                    print("OUTPUT:\n"+output_text)
                    calling_flag = 0
                    try:
                        calling_command = output_text.split("Tool Calling:")[1] if "Tool Calling:" in output_text else output_text
                        calling_command = calling_command.split("\n")[0].strip()
                        calling_flag = 1
                        if already_calling_list != []:
                            if calling_command == already_calling_list[-1][0]: # * 避免重复调用 
                                continue_flag = 0
                                calling_flag = 0
                        if execution:
                            execution_result = self.toolset.exec_tool_calling(calling_command)
                        else:
                            execution_result = "【Not Calling】"
                    except Exception as e:
                        execution_result = "【Error: "+str(e)+"】"
                    print("RESULT:\n"+str(execution_result))
                    if calling_flag:
                        if is_json_serializable(execution_result):
                            already_calling_list.append([calling_command,execution_result])
                        else:
                            already_calling_list.append([calling_command,str(execution_result)])
                output_data = {}
                output_data["query"] = query
                # output_data["input_text"] = input_text
                output_data["pred_text"] = raw_output_text_list
                output_data["pred_exec_result"] = already_calling_list
                output_data["gold"] = eval_data["calling_chain"]
                output_dataset.append(output_data)
                write_JSON(output_path, output_data, mode='a')
        write_JSON(output_path, output_dataset)

    @staticmethod
    def calculate_score(eval_data_path, execution=True):
        eval_dataset = read_JSON(eval_data_path)

        # * Tool
        correct_counter = 0
        pred_counter = 0
        gold_counter = 0
        parse_error_counter = 0
        for eval_data in eval_dataset:
            pred_tool_name_list = []
            for pred_calling in eval_data["pred_exec_result"]:
                if pred_calling[0] == 'None': # 跳过未调用的轮次
                    continue
                try:
                    func_name, args, kwargs = parse_function_calling(pred_calling[0])
                    pred_tool_name_list.append(func_name)
                except:
                    parse_error_counter += 1
                    # print("【Failed Parsing】:",pred_calling)
            gold_tool_name_list = []
            for gold_calling in eval_data["gold"]:
                gold_tool_name_list.append(gold_calling["tool"])
            
            pred_counter += len(pred_tool_name_list)
            gold_counter += len(gold_tool_name_list)
            for pred_tool_name in pred_tool_name_list:
                if pred_tool_name in gold_tool_name_list:
                    correct_counter += 1
                    del gold_tool_name_list[gold_tool_name_list.index(pred_tool_name)]

        print("【file】:",eval_data_path)
        print(f"【Parse】:{str(pred_counter)}/{str(pred_counter+parse_error_counter)}={str(pred_counter/(pred_counter+parse_error_counter))}")
        print("【Tool P】: ",str(correct_counter/pred_counter))
        print("【Tool R】: ",str(correct_counter/gold_counter))
        print("【Tool F1】:",str(2*correct_counter/(pred_counter+gold_counter)))

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

            pred_counter += sum([len(pred_param) for pred_param in pred_tool_param_list])
            gold_counter += sum([len(gold_param) for gold_param in gold_tool_param_list])
            for pred_index in range(len(pred_tool_name_list)):
                pred_tool_name = pred_tool_name_list[pred_index]
                if pred_tool_name in gold_tool_name_list:
                    gold_index = gold_tool_name_list.index(pred_tool_name)
                    correct_counter += sum(1 for a, b in zip(pred_tool_param_list[pred_index], gold_tool_param_list[gold_index]) if a == b)
                    # print(pred_tool_param_list[pred_index], gold_tool_param_list[gold_index])
                    del gold_tool_name_list[gold_index]
                    del gold_tool_param_list[gold_index]

        print("【file】:",eval_data_path)
        print("【Param P】: ",str(correct_counter/pred_counter))
        print("【Param R】: ",str(correct_counter/gold_counter))
        print("【Param F1】:",str(2*correct_counter/(pred_counter+gold_counter)))

        # * Return
        if execution:
            correct_counter = 0
            pred_counter = 0
            gold_counter = 0
            for eval_data in eval_dataset:
                pred = []
                for pred_tool in eval_data["pred_exec_result"]:
                    if pred_tool[0] == "None":
                        continue
                    pred.append([pred_tool[0].split("(")[0],pred_tool[1]])
                gold = [ [gold_tool["tool"], gold_tool["return"]] for gold_tool in eval_data["gold"] ]

                pred_counter += len(pred)
                gold_counter += len(gold)
                for pred_result in pred:
                    if pred_result[0] in [gold_tool[0] for gold_tool in gold]:
                        t_index = [gold_tool[0] for gold_tool in gold].index(pred_result[0])
                        if type(pred_result[1]) is float and type(gold[t_index][1]) is float:
                            if round(pred_result[1], 4) == round(gold[t_index][1], 4):
                                correct_counter += 1
                                gold.pop(t_index)
                        elif pred_result[1] == gold[t_index][1]:
                            correct_counter += 1
                            gold.pop(t_index)

            print("【file】:",eval_data_path)
            print("【Result P】: ",str(correct_counter/pred_counter))
            print("【Result R】: ",str(correct_counter/gold_counter))
            print("【Result F1】:",str(2*correct_counter/(pred_counter+gold_counter)))
