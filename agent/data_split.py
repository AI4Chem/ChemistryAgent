import os
import random
import copy

from utils.file_io import read_JSON, write_JSON


# * single case split 
def sample_single_case(raw_data_dir, output_dataset_path_list):
    train_dataset = []
    dev_dataset = []
    test_dataset = []

    for _,dirs,_ in os.walk(raw_data_dir):
        for dir_name in dirs:
            folder_path = os.path.join(raw_data_dir, dir_name)
            for _,_,files in os.walk(folder_path):
                for tool_file in files:
                    tool_name = tool_file.split(".")[0]
                    tool_full_name = f"{dir_name}/{tool_name}"
                    tool_dataset = read_JSON(os.path.join(folder_path,tool_file))

                    sampled_dev_dataset = random.sample(tool_dataset, max(int(0.1*len(tool_dataset)),1))
                    remaining_list = [item for item in tool_dataset if item not in sampled_dev_dataset]
                    sampled_test_dataset = random.sample(remaining_list, max(int(0.1*len(tool_dataset)),1))
                    sampled_train_dataset = [item for item in remaining_list if item not in sampled_test_dataset]

                    for sampled_data in sampled_train_dataset:
                        train_dataset.append({
                            "id": None,
                            "query": sampled_data[0],
                            "calling_chain": [{"tool":tool_full_name, "params":sampled_data[1], "return":sampled_data[2]}]
                        })
                    for sampled_data in sampled_dev_dataset:
                        dev_dataset.append({
                            "id": None,
                            "query": sampled_data[0],
                            "calling_chain": [{"tool":tool_full_name, "params":sampled_data[1], "return":sampled_data[2]}]
                        })            
                    for sampled_data in sampled_test_dataset:
                        test_dataset.append({
                            "id": None,
                            "query": sampled_data[0],
                            "calling_chain": [{"tool":tool_full_name, "params":sampled_data[1], "return":sampled_data[2]}]
                        })
    data_id = 0
    for idx in range(len(train_dataset)):
        train_dataset[idx]["id"] = "single_"+str(data_id)
        data_id+=1
    for idx in range(len(dev_dataset)):
        dev_dataset[idx]["id"] = "single_"+str(data_id)
        data_id+=1
    for idx in range(len(test_dataset)):
        test_dataset[idx]["id"] = "single_"+str(data_id)
        data_id+=1
    write_JSON(output_dataset_path_list[0], train_dataset)
    write_JSON(output_dataset_path_list[1], dev_dataset)
    write_JSON(output_dataset_path_list[2], test_dataset)

# * multiple case split 
def sample_multi_case(raw_data_path, output_dataset_path_list):
    train_dataset = []
    dev_dataset = []
    test_dataset = []

    tool_dataset = read_JSON(raw_data_path)

    sampled_dev_dataset = random.sample(tool_dataset, max(int(0.1*len(tool_dataset)),1))
    remaining_list = [item for item in tool_dataset if item not in sampled_dev_dataset]
    sampled_test_dataset = random.sample(tool_dataset, max(int(0.1*len(tool_dataset)),1))
    sampled_train_dataset = [item for item in remaining_list if item not in sampled_test_dataset]

    data_id = 0
    for sampled_data in sampled_train_dataset:
        train_dataset.append({
            "id": "multiple_"+str(data_id),
            "query": sampled_data[0],
            "calling_chain": [{"tool":calling[0], "params":calling[1], "return":calling[2]} for calling in sampled_data[1]]
        })
        data_id += 1
    for sampled_data in sampled_dev_dataset:
        dev_dataset.append({
            "id": "multiple_"+str(data_id),
            "query": sampled_data[0],
            "calling_chain": [{"tool":calling[0], "params":calling[1], "return":calling[2]} for calling in sampled_data[1]]
        })
        data_id += 1              
    for sampled_data in sampled_test_dataset:
        test_dataset.append({
            "id": "multiple_"+str(data_id),
            "query": sampled_data[0],
            "calling_chain": [{"tool":calling[0], "params":calling[1], "return":calling[2]} for calling in sampled_data[1]]
        })
        data_id += 1 
    write_JSON(output_dataset_path_list[0], train_dataset)
    write_JSON(output_dataset_path_list[1], dev_dataset)
    write_JSON(output_dataset_path_list[2], test_dataset)

# * pymatgen data split

def sample_single_mat_case(raw_data_dir, output_dataset_path_list):
    train_dataset = []
    dev_dataset = []
    test_dataset = []

    for _,dirs,_ in os.walk(raw_data_dir):
        for dir_name in dirs:
            folder_path = os.path.join(raw_data_dir, dir_name)
            for _,_,files in os.walk(folder_path):
                for tool_file in files:
                    tool_name = tool_file.split(".")[0]
                    tool_full_name = f"{dir_name}/{tool_name}"
                    tool_dataset = read_JSON(os.path.join(folder_path,tool_file))

                    sampled_dev_dataset = random.sample(tool_dataset, min(int(0.1*len(tool_dataset)),10))
                    remaining_list = [item for item in tool_dataset if item not in sampled_dev_dataset]
                    sampled_test_dataset = random.sample(remaining_list, min(int(0.1*len(tool_dataset)),10))
                    sampled_train_dataset = [item for item in remaining_list if item not in sampled_test_dataset]

                    for sampled_data in sampled_train_dataset:
                        train_dataset.append({
                            "id": None,
                            "query": sampled_data[0],
                            "calling_chain": [{"tool":tool_full_name, "params":sampled_data[1], "return":sampled_data[2]}]
                        })
                    for sampled_data in sampled_dev_dataset:
                        dev_dataset.append({
                            "id": None,
                            "query": sampled_data[0],
                            "calling_chain": [{"tool":tool_full_name, "params":sampled_data[1], "return":sampled_data[2]}]
                        })            
                    for sampled_data in sampled_test_dataset:
                        test_dataset.append({
                            "id": None,
                            "query": sampled_data[0],
                            "calling_chain": [{"tool":tool_full_name, "params":sampled_data[1], "return":sampled_data[2]}]
                        })
    data_id = 0
    for idx in range(len(train_dataset)):
        train_dataset[idx]["id"] = "single_mat_"+str(data_id)
        data_id+=1
    for idx in range(len(dev_dataset)):
        dev_dataset[idx]["id"] = "single_mat_"+str(data_id)
        data_id+=1
    for idx in range(len(test_dataset)):
        test_dataset[idx]["id"] = "single_mat_"+str(data_id)
        data_id+=1
    write_JSON(output_dataset_path_list[0], train_dataset)
    write_JSON(output_dataset_path_list[1], dev_dataset)
    write_JSON(output_dataset_path_list[2], test_dataset)
    
def sample_multiple_mat_case(raw_data_path, output_dataset_path_list):
    sampled_test_dataset = []
    sampled_notest_dataset = []

    test_dataset = []
    notest_dataset = []

    classified_dataset = {}

    raw_dataset = read_JSON(raw_data_path)
    for data in raw_dataset:
        tool_chain = str([calling[0] for calling in data[1]])
        if tool_chain not in classified_dataset:
            classified_dataset[tool_chain] = [data]
        else:
            classified_dataset[tool_chain].append(data)
    
    for tool_chain in classified_dataset:
        data_list = copy.deepcopy(classified_dataset[tool_chain])
        sampled_test_dataset.append(data_list.pop(random.randint(0, len(data_list) - 1)))
        sampled_notest_dataset.extend(data_list)

    data_id = 0
    for sampled_data in sampled_test_dataset:
        test_dataset.append({
            "id": "multiple_mat_"+str(data_id),
            "query": sampled_data[0],
            "calling_chain": [{"tool":calling[0], "params":calling[1], "return":calling[2]} for calling in sampled_data[1]]
        })
        data_id += 1
    for sampled_data in sampled_notest_dataset:
        notest_dataset.append({
            "id": "multiple_mat_"+str(data_id),
            "query": sampled_data[0],
            "calling_chain": [{"tool":calling[0], "params":calling[1], "return":calling[2]} for calling in sampled_data[1]]
        })
        data_id += 1
    
    write_JSON(output_dataset_path_list[0], test_dataset)
    write_JSON(output_dataset_path_list[1], notest_dataset)

if __name__ == "__main__":
    # data_dir = "data/single/single_case"
    # output_path_list = [
    #     "data/single/single_train.jsonl",
    #     "data/single/single_dev.jsonl",
    #     "data/single/single_test.jsonl"
    # ]
    # sample_single_case(data_dir, output_path_list)

    # data_path = "data/multiple/method_SealTools/multiple_case.jsonl"
    # output_path_list = [
    #     "data/multiple/method_SealTools/multiple_train.jsonl",
    #     "data/multiple/method_SealTools/multiple_dev.jsonl",
    #     "data/multiple/method_SealTools/multiple_test.jsonl"
    # ]
    # sample_multi_case(data_path, output_path_list)

    data_dir = "data/pymatgen/single/single_case"
    output_path_list = [
        "data/pymatgen/single/single_train.jsonl",
        "data/pymatgen/single/single_dev.jsonl",
        "data/pymatgen/single/single_test.jsonl"
    ]
    sample_single_mat_case(data_dir, output_path_list)

    data_path = "data/pymatgen/multiple/method_SealTools/multiple_case.jsonl"
    output_path_list = [
        "data/pymatgen/multiple/method_SealTools/multiple_test.jsonl",
        "data/pymatgen/multiple/method_SealTools/multiple_notest.jsonl",
    ]
    sample_multiple_mat_case(data_path, output_path_list)