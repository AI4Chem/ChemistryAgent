import os
import re

from utils.file_io import read_JSON, write_JSON
from . import post_process_execution_result
from .parse_func_file import get_function_info
from .toolparam import get_tool_param, get_matgen_tool_param
from .toolinfo import basic_tool_package, matgen_tool_package


def parse_docstring(input_string: str) -> dict | str:
    # 使用正则表达式解析字符串
    pattern = re.compile(
        r'Name:\s*(?P<name>.*?)\n'
        r'Description:\s*(?P<description>.*?)\n'
        r'Parameters:\n(?P<parameters>.*?)\n'
        r'Returns:\n(?P<returns>.*?)(?=\n[A-Za-z]|$)', re.DOTALL)
    match = pattern.match(input_string)

    if match:
        # 提取匹配的内容
        name = match.group('name')
        description = match.group('description')
        parameters = match.group('parameters').strip().split('\n')
        returns = match.group('returns').strip().split('\n')

        return {
            "Name": name,
            "Description": description,
            "Parameters": parameters,
            "Returns": returns
            }
    else:
        return "No match found."


class ToolPool:
    def __init__(self, tool_package_name, tool_package_path = "src/toolpool"):
        self.tooldir = tool_package_name
        self.packagepath = tool_package_path
        self.toolpath = read_JSON(os.path.join(self.packagepath,self.tooldir,"tools.json"))
        self.toolrawinfo = {}
        self.get_toolrawinfo_from_toolpath()
        self.toolinfo = {}
        self.process_toolinfo()

    def get_toolrawinfo_from_toolpath(self):
        for tool_name in self.toolpath:
            tool_relative_path = self.toolpath[tool_name]["path"]
            self.toolrawinfo[tool_name] = get_function_info(os.path.join(self.packagepath,self.tooldir,tool_relative_path))[tool_name]

    def process_toolinfo(self):
        for tool_name in self.toolrawinfo:
            print("Getting Tool Info: ", tool_name)
            raw_info = self.toolrawinfo[tool_name]
            doc_string = parse_docstring(raw_info['description'])
            tool_info = {}
            tool_info["name"] = raw_info['name']
            tool_info["description"] = doc_string['Description']
            tool_info["parameters"] = {"type": "object", "properties":{}}
            param_description_dict = {}
            for param_text in doc_string["Parameters"]:
                param_name = param_text.split(":")[0]
                param_description = param_text[len(param_name)+1:].strip()
                param_description_dict[param_name.strip()] = param_description
            for raw_param_dict in raw_info['parameters']:
                param_name = raw_param_dict["name"]
                tool_info["parameters"]["properties"][param_name] = {}
                tool_info["parameters"]["properties"][param_name]["type"] = raw_param_dict["type"]
                if param_name in param_description_dict:
                    tool_info["parameters"]["properties"][param_name]["description"] = param_description_dict[param_name][len(raw_param_dict["type"])+1:].strip()
            self.toolinfo[tool_name] = tool_info


    def call_tool(self, tool_name: str, parameters_list: list):
        tool_path = self.toolpath[tool_name]["path"]
        import_cmd = "from toolpool."+self.tooldir+"."+tool_path.split(".py")[0].replace("/",".")+" import "+tool_name
        exec(import_cmd)
        calling_cmd = tool_name+"("+",".join([str(param) if type(param) is not str else '"'+param+'"' for param in parameters_list])+")"
        calling_result = eval(calling_cmd)
        calling_result = post_process_execution_result(calling_result)
        return calling_result
        # try:
        #     return eval(calling_cmd)
        # except:  # noqa: E722
        #     return "Calling Tool Failed."

    def generate_tool_calling_result(self, num: int = 200, dir: str="tool_calling_result", execution: bool=True):
        dir_path = "data/"+dir+"/"+self.tooldir+"/"
        os.makedirs(dir_path, exist_ok=True)
        for tool_name in self.toolinfo:
            print("Generating calling result with Tool: ", tool_name)
            file_path = os.path.join(dir_path, tool_name+".jsonl")
            if os.path.exists(file_path):
                calling_dataset = read_JSON(file_path)
            else:
                calling_dataset = []
            input_set = set([str(calling_data[0]) for calling_data in calling_dataset])
            RETRY_MAX = 100 #5
            retry = RETRY_MAX
            while retry > 0 and (len(calling_dataset) < num):

                if execution:
                    input_param = get_tool_param(tool_name, self.toolinfo[tool_name]["parameters"]["properties"])
                    print(f"Generating  Tool {tool_name} with Input {input_param}")
                    result = None
                    try:
                        result = self.call_tool(tool_name, input_param)
                    except Exception as e:  # noqa: E722
                        print("Calling Error: ", e)
                        retry -= 1
                    # if self.tooldir=="chem_lib" and tool_name == "galvanic_cell_potential":
                    #     if isinstance(result, str):
                    #         retry += 1
                    #         continue
                    # elif self.tooldir=="chem_lib" and tool_name == "galvanic_cell_properties":
                    #     if isinstance(result, str):
                    #         retry += 1
                    #         continue
                    # if self.tooldir=="chemcrow" and tool_name == "ExplosiveCheck":
                    #     if result.startswith("Explosive Check Error."):
                    #         retry += 1
                    #         continue
                    # if self.tooldir=="chemistrytools" and tool_name == "get_compound_charge_by_CID":
                    #     if result==0:
                    #         retry += 1
                    #         continue

                else:
                    input_param = get_matgen_tool_param(tool_name, self.toolinfo[tool_name]["parameters"]["properties"])
                    print(f"Generating  Tool {tool_name} with Input {input_param}")
                    result = None # - 不再尝试调用工具

                if str(input_param) not in input_set and (result is not None if execution else True) :
                    retry = RETRY_MAX
                    input_set.add(str(input_param))
                    calling_dataset.append([input_param, result])
                    write_JSON(file_path, calling_dataset)
                    print(f"Tool: {tool_name}, Params: {input_param}, Returns: {result}")
                else:
                    retry -= 1


# {
#     "name": "get_current_temperature",
#     "description": "Get the current temperature for a specific location",
#     "parameters": 
#     {
#           "type": "object",
#           "properties": 
#            {
#             "location": {
#               "type": "string",
#               "description": "The city and state, e.g., San Francisco, CA"
#             },
#             "unit": {
#               "type": "string",
#               "enum": ["Celsius", "Fahrenheit"],
#               "description": "The temperature unit to use. Infer this from the user's location."
#             }
#      },
#      "required": ["location", "unit"]
# }


if __name__ == "__main__":
    # tooldir_list = ["chem_lib", "cactus", "chemcrow", "chemistrytools"]
    # for tooldir in tooldir_list:
    #     tools = ToolPool(tooldir)
    #     tools.generate_tool_calling_result()

    tools = ToolPool("chem_lib")
    print(tools.call_tool("get_element_properties",["Zn"]))
    print(tools.call_tool("analyze_combustion",[44,18]))
