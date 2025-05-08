import os
import csv
import json
import pickle

def is_json_serializable(obj):
    try:
        json.dumps(obj)
        return True
    except (TypeError, OverflowError):
        return False

def read_json(data_path):
    with open(data_path,'r', encoding='UTF-8') as f:
        dataset = json.load(f)
    return dataset


def read_jsonl(data_path):
    dataset=[]
    with open(data_path,'r', encoding='UTF-8') as f:
        for line in f:
            dataset.append(json.loads(line))
    return dataset


def write_json(data_path, dataset, indent=0, mode='w'):
    with open(data_path, mode, encoding='UTF-8') as f:
        if indent != 0:
            json.dump(dataset, f, ensure_ascii=False, indent=indent)
        else:
            json.dump(dataset, f, ensure_ascii=False)


def write_jsonl(data_path, dataset, indent=0, mode='w'):
    with open(data_path, mode, encoding='UTF-8') as f:
        if mode == 'w':
            for data in dataset:
                if indent != 0:
                    f.writelines(json.dumps(data, ensure_ascii=False, indent=indent))
                else:
                    f.writelines(json.dumps(data, ensure_ascii=False))
                f.write('\n')
        else:
            assert mode == 'a'
            if indent != 0:
                f.writelines(json.dumps(dataset, ensure_ascii=False, indent=indent))
            else:
                f.writelines(json.dumps(dataset, ensure_ascii=False))
            f.write('\n')


def read_JSON(data_path):
    if data_path.split('.')[-1].lower() == "json":
        try:
            return read_json(data_path)
        except:
            return read_jsonl(data_path)
    elif data_path.split('.')[-1].lower() == "jsonl":
        try:
            return read_jsonl(data_path)
        except:
            return read_json(data_path)
    else:
        print("data_path error !!!")


def write_JSON(data_path, dataset, indent=0, mode='w'):
    if data_path.split('.')[-1].lower() == "json":
        try:
            return write_json(data_path, dataset, indent, mode)
        except:
            return write_jsonl(data_path, dataset, indent, mode)
    elif data_path.split('.')[-1].lower() == "jsonl":
        try:
            return write_jsonl(data_path, dataset, indent, mode)
        except:
            return write_json(data_path, dataset, indent, mode)
    else:
        print("data_path error !!!")


def str_to_JSON(input_text_string):
    def replace_in_JSON(data, old, new):
        if isinstance(data, dict):
            return {k: replace_in_JSON(v, old, new) for k, v in data.items()}
        elif isinstance(data, list):
            return [replace_in_JSON(item, old, new) for item in data]
        elif isinstance(data, str):
            return data.replace(old, new)
        else:
            return data
    output_dict = json.loads(json.dumps(eval(input_text_string.replace("\n","【换行符】"))))
    output_dict = replace_in_JSON(output_dict, "【换行符】", "\n")
    return output_dict


def read_pickle(data_path):
    with open(data_path,'rb') as f:
        dataset = pickle.load(f)
    return dataset


def write_pickle(data_path, dataset):
    with open(data_path,'wb') as f:
        pickle.dump(dataset, f)


def read_csv(data_path, delimiter=","): # "\t" for tsv
    with open(data_path, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter=delimiter)
        data_list = [row for row in reader]
    return data_list


def write_csv(data_path, dataset, delimiter=","): # "\t" for tsv
    with open(data_path, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file, delimiter=delimiter)
        writer.writerows(dataset)


