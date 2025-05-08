import torch
from transformers import AutoTokenizer, AutoModelForCausalLM

class InternLM:
    def __init__(self, checkpoint_path, lora_path = None):
        self.tokenizer = AutoTokenizer.from_pretrained(checkpoint_path, trust_remote_code=True)
        # Set `torch_dtype=torch.float16` to load model in float16, otherwise it will be loaded as float32 and cause OOM Error.
        if lora_path is not None:
            # from peft import PeftModelForCausalLM
            # self.raw_model = AutoModelForCausalLM.from_pretrained(checkpoint_path, torch_dtype=torch.float16, trust_remote_code=True, device_map="auto")
            # self.model = PeftModelForCausalLM.from_pretrained(model=self.raw_model, model_id=lora_path, trust_remote_code=True, device_map="auto")
            from peft import AutoPeftModelForCausalLM
            self.model = AutoPeftModelForCausalLM.from_pretrained(lora_path, trust_remote_code=True, device_map="auto") 
        else:
            self.model = AutoModelForCausalLM.from_pretrained(checkpoint_path, torch_dtype=torch.float16, trust_remote_code=True, device_map="auto")
        self.model = self.model.eval()

    def answer(self,input_text, stop_strings=None):
        inputs = self.tokenizer(input_text, return_tensors="pt").to("cuda")
        input_length = inputs.input_ids.shape[1]
        generate_ids = self.model.generate(inputs.input_ids, max_new_tokens=1024, stop_strings=stop_strings, tokenizer=self.tokenizer)
        response = self.tokenizer.batch_decode(generate_ids, skip_special_tokens=True)[0]
        re_inputs = self.tokenizer.batch_decode(inputs.input_ids, skip_special_tokens=True)[0]
        return response[len(re_inputs):]

    def __call__(self, input_text):
        return self.answer(input_text)
