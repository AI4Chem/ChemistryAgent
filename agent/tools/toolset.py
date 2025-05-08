from . import post_process_execution_result
from .toolpool import ToolPool

class ToolSet:
    def __init__(self,toolpackage_list, toolpackage_path = "src/toolpool"):
        self.tool_package_list = toolpackage_list
        self.tool_pool_dict = {tool_package: ToolPool(tool_package, tool_package_path=toolpackage_path) for tool_package in self.tool_package_list}

    def __call__(self, tool_package_name):
        return self.tool_pool_dict[tool_package_name]
    
    def exec_tool_calling(self, tool_calling: str): # 不将参数列表化，直接调用
        index = tool_calling.find('/')
        tool_dir = tool_calling[:index]
        calling_command = tool_calling[index+1:]
        tool_name = calling_command.split("(")[0]
        tool_path = self.tool_pool_dict[tool_dir].toolpath[tool_name]["path"]
        import_cmd = "from toolpool."+self.tool_pool_dict[tool_dir].tooldir+"."+tool_path.split(".py")[0].replace("/",".")+" import "+tool_name
        exec(import_cmd)
        calling_cmd = calling_command
        calling_result = eval(calling_cmd)
        calling_result = post_process_execution_result(calling_result)
        if calling_result == calling_cmd[1:-1]:
            return "【Error: Output Format】"
        else:
            return calling_result

    def get_all_tool_names_list(self):
        full_tool_names_list = []
        for package_name in self.tool_pool_dict:
            for tool_name in self.tool_pool_dict[package_name].toolinfo:
                full_tool_names_list.append(package_name+"/"+tool_name)
        return full_tool_names_list
        
    def get_tool_information_text_list(self, full_tool_name_list):
        tool_information_list = []
        for candidate_tool in full_tool_name_list:
            tool_information_list.append(self.tool_pool_dict[candidate_tool.split("/")[0]].toolinfo[candidate_tool.split("/")[1]])
            tool_information_list[-1]["name"] = candidate_tool
            tool_information_list[-1] = str(tool_information_list[-1])
        return tool_information_list