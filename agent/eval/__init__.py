def process_calling_list_text(already_calling_list):
    if already_calling_list != []:
        calling_text_list = ["Existing Tool Calling Chain:"]
        for tool_call in already_calling_list:
            calling_text_list.append("Tool Calling: "+tool_call[0])
            calling_text_list.append("Result: "+str(tool_call[1]))
        return "\n".join(calling_text_list)
    else:
        return ""

from .singletool import SingleToolEvaluation
from .multitool import MultiToolEvaluation
from .response_generation import ResponseEvaluation