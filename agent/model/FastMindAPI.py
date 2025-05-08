
class FMModel:
    def __init__(self, client):
        self.client = client
        self.model_name = None

    def answer(self, input_text, stop_strings=None):
        data = {'input_text':input_text, 'stop_strings':stop_strings, 'max_new_tokens':1024}
        response = self.client.generate(self.model_name, data=data)
        print("【response】:"+str(response))
        return response["output_text"]

    def __call__(self, input_text):
        return self.answer(input_text)
