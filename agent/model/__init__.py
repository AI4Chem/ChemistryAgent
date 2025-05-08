from transformers import StoppingCriteria

class StopOnNewLine(StoppingCriteria):
    def __init__(self, tokenizer, input_length):
        self.tokenizer = tokenizer
        self.input_length = input_length

    def __call__(self, input_ids, scores, **kwargs):
        # 只检查新生成的部分，忽略输入的内容
        decoded_text = self.tokenizer.decode(input_ids[0][self.input_length:], skip_special_tokens=True)
        return '\n' in decoded_text

from .GPT import GPT
from .InternLM import InternLM
from .LLaMA import LLaMA2, LLaMA3
from .FastMindAPI import FMModel

def AutoModel(model_name, *args, **kwargs):
    match model_name:
        case "GPT":
            return GPT(*args, **kwargs)
        case "InternLM":
            return InternLM(*args, **kwargs)
        case "LLaMA2":
            return LLaMA2(*args, **kwargs)
        case "LLaMA3":
            return LLaMA3(*args, **kwargs)
        case "FaseMindAPI":
            return FMModel(*args, **kwargs)
