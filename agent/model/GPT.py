import openai
from openai import OpenAI

from .api_key import OPENAI_KEY, BASE_URL

# 设置OpenAI API密钥
client = OpenAI(
    api_key= OPENAI_KEY,
    base_url= BASE_URL)

# 定义一个函数来调用ChatGPT API
def call_openai_api(prompt):
    try:
        completion = client.chat.completions.create(
        model= 'gpt-4o-mini', # "gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": prompt}
        ]
        )
        return completion.choices[0].message.content
    except Exception as e:
        return e

class GPT:
    def __init__(self, model_name):
        self.model_name = model_name

    def answer(self,input_text, stop_strings=None):
        try:
            completion = client.chat.completions.create(
            model = self.model_name,
            messages = [
                {"role": "system", "content": "You are a helpful assistant."},
                {"role": "user", "content": input_text}
            ],
            stop = stop_strings
            )
            return completion.choices[0].message.content
        except Exception as e:
            return str(e)
        
    def __call__(self, input_text):
        return self.answer(input_text)