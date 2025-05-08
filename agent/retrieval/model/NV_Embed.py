import torch
import torch.nn.functional as F
from transformers import AutoModel

class NV_Embed():
    def __init__(self, checkpoint_path: str):
        self.checkpoint_path = checkpoint_path
        self.model = AutoModel.from_pretrained(self.checkpoint_path, trust_remote_code=True, device_map="auto")
        
    def encode(self, text: str | list[str]):
        if isinstance(text, str):
            text = [text]
        text_embeddings = self.model.encode(text)
        text_embeddings = F.normalize(text_embeddings, p=2, dim=1)
        return text_embeddings if text_embeddings.size()[0]!= 1 else text_embeddings[0]
