import sys
sys.path.append('.')

from fastapi import FastAPI
from pydantic import BaseModel

from .model import GPT
from .tools import ToolSet
from .chat import Agent

toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])

import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"
retriever_checkpoint = "/mnt/hwfile/ai4chem/wumengsong/PTM/all-MiniLM-L6-v2"

agent = Agent(GPT(),toolset, retriever_checkpoint)

# * = = =

class Tool_List(BaseModel):
    tool_list: list

class Chat(BaseModel):
    query: str

class Chat_Context(BaseModel):
    query: str
    retrieved_tools: list
    tool_calling: list
    prompt: dict = None

app = FastAPI()


@app.get("/")
async def index():
    return {"msg": "The server has been deployed successfully!"}

@app.get("/get_all_tools/")
async def get_all_available_tools():
    return agent.toolset.get_all_tool_names_list()

@app.post("/get_tool_information/")
async def get_tool_information(item: Tool_List):
    return agent.toolset.get_tool_information_text_list(item.tool_list)

@app.post("/chat/")
async def chat(item: Chat):
    return agent.chat(item.query)

@app.post("/call_tool/")
async def chat_with_the_tool(item: Chat_Context):
    return agent.call_tool(raw_input_text=item.query, retrieved_tools=item.retrieved_tools, tool_calling=item.tool_calling, prompt=item.prompt)
