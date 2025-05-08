import sys
sys.path.append('./src/')

import random

from retrieval import arrange_tool_list_and_corpus, select_candidate_tool_list
from retrieval.model import NV_Embed
from retrieval.retriever import BM25_retriever, DenseRetriever
from utils.file_io import read_JSON, write_JSON

def construct_dataset_with_retriever_result(toolset,input_data_file_path_list,output_file_path_list):
     # retriever = BM25_retriever(*arrange_tool_list_and_corpus(toolset))
     # retriever_checkpoint = "/Users/wumengsong/Resource/all-MiniLM-L6-v2"
     retriever_checkpoint = "/home/bingxing2/ailab/group/ai4phys/wumengsong/NV-Embed-v2"

     retriever = DenseRetriever(*arrange_tool_list_and_corpus(toolset), NV_Embed(retriever_checkpoint))
     for file_idx in range(len(input_data_file_path_list)):
          input_dataset = read_JSON(input_data_file_path_list[file_idx])
          output_dataset = []
          for input_data in input_dataset:
               print("Retrieving: ", input_data["query"])
               retrieved_tools = retriever.retrieve(input_data["query"], top_k=10)
               gold_tools = [tool["tool"] for tool in input_data["calling_chain"]]
               candidate_tools = select_candidate_tool_list(10, gold_tools, retrieved_tools)
               random.shuffle(candidate_tools)

               output_dataset.append({
                    "id": input_data["id"],
                    "query": input_data["query"],
                    "calling_chain": input_data["calling_chain"],
                    "candidate_tools": candidate_tools
               })
          write_JSON(output_file_path_list[file_idx], output_dataset)

# 从单工具评测集中筛选出适合用于调用工具与否效果对比的子集
def filter_dataset_for_tool_calling_effect(input_data_file_path, output_file_path):
     exception_list = [
          "chem_lib/get_element_properties",
          "chem_lib/galvanic_cell_properties",
          "chem_lib/perform_electrolysis",
          "chem_lib/convert_compound_stoichiometry_amount",
          "chemistrytools/get_element_information"
     ]
     input_dataset = read_JSON(input_data_file_path)
     output_dataset = []
     for data in input_dataset:
          if data["calling_chain"][0]["tool"] not in exception_list:
               output_dataset.append(data)
     write_JSON(output_file_path, output_dataset)

def filter_result(filtered_dataset_file_path,input_data_file_path, output_file_path):
     filtered_dataset = read_JSON(filtered_dataset_file_path)
     filtered_query_list = [data["query"] for data in filtered_dataset]
     input_dataset = read_JSON(input_data_file_path)
     output_dataset = []
     for data in input_dataset:
          if data["query"] in filtered_query_list:
               output_dataset.append(data)
     write_JSON(output_file_path, output_dataset)

if __name__ == "__main__":
     # filter_dataset_for_tool_calling_effect("data/single/single_test.jsonl", "data/single/single_test_comparison.jsonl")
     filter_result("data/single/single_test_comparison.jsonl", "result/single_test_DPR_gpt-4o-mini_score.jsonl", "result/single_test_DPR_gpt-4o-mini_WT_score.jsonl")