import sys
sys.path.append('./src/')

from fastapi import FastAPI
from pydantic import BaseModel

from .tools import ToolSet

toolset = ToolSet(["chem_lib", "cactus", "chemcrow", "chemistrytools"])

import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"

# * = = =

class Calling(BaseModel):
    command: str

app = FastAPI()

@app.get("/")
async def index():
    return {"msg": "The server has been deployed successfully!"}

@app.post("/call_tool/")
async def chat_with_the_tool(item: Calling):
    return toolset.exec_tool_calling(item.command)

