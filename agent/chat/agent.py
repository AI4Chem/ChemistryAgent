from retrieval import arrange_tool_list_and_corpus
from retrieval.retriever import BM25_retriever, DenseRetriever
from retrieval.model import MiniLM
from prompt import prompt_eval_0_shot

from eval import process_calling_list_text

class Agent:
    def __init__(self, model, toolset, retriever=None):
        self.model = model
        self.toolset = toolset

        # * 检索模块 初始化 
        if retriever is None:
            self.retriever = BM25_retriever(*arrange_tool_list_and_corpus(self.toolset))
        else:
            assert isinstance(retriever, str)
            self.retriever = DenseRetriever(*arrange_tool_list_and_corpus(self.toolset), MiniLM(retriever))

        self.tool_calling_prompt = prompt_eval_0_shot["instruction"]+prompt_eval_0_shot["input"]

    def __call__(self, input_text):
        calling_status, result, calling_command = self.call_tool(input_text)
        # return calling_status, result, calling_command
        if calling_status:
            print("Need to call the tool:",calling_command)
            print(result)
        else:
            print("No need to call a tool.")
            print(result)


    def call_tool(self, raw_input_text, retrieved_tools=[], tool_calling=[], prompt=None):
        if retrieved_tools == []:
            retrieved_tools = self.retriever.retrieve(raw_input_text, top_k=10)
        tool_information = self.toolset.get_tool_information_text_list(retrieved_tools)
        if prompt is not None:
            input_text = (prompt["instruction"]+prompt["input"]).format(raw_input_text,
                                                                        process_calling_list_text([[data[1],data[2]] for data in tool_calling]),
                                                                        "\n".join(tool_information))
        else:
            input_text = self.tool_calling_prompt.format(raw_input_text,process_calling_list_text([[data[1],data[2]] for data in tool_calling]),"\n".join(tool_information))
        output_text = self.model(input_text)
        try:
            calling_command = output_text.split("Tool Calling:")[1].strip() if "Tool Calling:" in output_text else output_text.strip()
            tool_name = calling_command.split("(")[0]
            if tool_name not in retrieved_tools:
                tool_calling.append(["None", calling_command, ""])
            else:
                execution_result = self.toolset.exec_tool_calling(calling_command)
                tool_calling.append([tool_name, calling_command, execution_result])
        except Exception as e:
            execution_result = "【Error: "+str(e)+"】"
            tool_calling.append(["Error", output_text, execution_result])
        return {"query": raw_input_text, "retrieved_tools": retrieved_tools, "tool_calling": str(tool_calling)}
    
    def chat(self, input_text):
        return self.model(input_text)
