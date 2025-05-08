

def arrange_tool_list_and_corpus(toolset):
    tool_corpus = []
    tool_list = []
    for tool_dir in toolset.tool_package_list:
        for tool_name in toolset(tool_dir).toolinfo:
            tool_description = toolset(tool_dir).toolinfo[tool_name]["description"]
            tool_corpus.append("Tool Name: {}\nTool Description: {}".format(tool_dir+"/"+tool_name, tool_description))
            tool_list.append(tool_dir+"/"+tool_name)
    return tool_list, tool_corpus

def select_candidate_tool_list(top_n, gold_tool_list, retrieved_tool_list):
    candidate_tool_list = gold_tool_list
    while len(candidate_tool_list) < top_n:
            for tool in retrieved_tool_list:
                if tool not in candidate_tool_list:
                    candidate_tool_list.append(tool)
    return candidate_tool_list