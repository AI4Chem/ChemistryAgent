from transformers import AutoTokenizer, AutoModelForCausalLM

class LLaMA2:
    def __init__(self,checkpoint_path=''):
        self.checkpoint_path = checkpoint_path
        self.tokenizer = AutoTokenizer.from_pretrained(self.checkpoint_path)
        self.model = AutoModelForCausalLM.from_pretrained(self.checkpoint_path, device_map="auto")

    def answer(self,input_text, stop_strings=None):
        inputs = self.tokenizer(input_text, return_tensors="pt").to("cuda")
        input_length = inputs.input_ids.shape[1]

        outputs = self.model.generate(**inputs, max_new_tokens=1024, stop_string=stop_strings, tokenizer=self.tokenizer)#2048
        output_text = self.tokenizer.batch_decode(outputs, skip_special_tokens=True, clean_up_tokenization_spaces=False)[0]
        return output_text[len(input_text):]

    def __call__(self, input_text):
        return self.answer(input_text)

class LLaMA3:
    def __init__(self,checkpoint_path=''):
        self.checkpoint_path = checkpoint_path
        self.tokenizer = AutoTokenizer.from_pretrained(self.checkpoint_path)
        self.model = AutoModelForCausalLM.from_pretrained(self.checkpoint_path, device_map="auto")

    def answer(self,input_text, stop_strings=None):
        inputs = self.tokenizer(input_text, return_tensors="pt").to("cuda")
        input_length = inputs.input_ids.shape[1]

        outputs = self.model.generate(**inputs, max_new_tokens=1024, stop_strings=stop_strings, tokenizer=self.tokenizer)#2048
        output_text = self.tokenizer.batch_decode(outputs, skip_special_tokens=True, clean_up_tokenization_spaces=False)[0]
        return output_text[len(input_text):]

    def __call__(self, input_text):
        return self.answer(input_text)