import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, AutoModel

# Mean Pooling - Take attention mask into account for correct averaging
def mean_pooling(model_output, attention_mask):
    token_embeddings = model_output[0] # First element of model_output contains all token embeddings
    input_mask_expanded = attention_mask.unsqueeze(-1).expand(token_embeddings.size()).float()
    return torch.sum(token_embeddings * input_mask_expanded, 1) / torch.clamp(input_mask_expanded.sum(1), min=1e-9)

class MiniLM:
    def __init__(self, checkpoint_path):
        self.checkpoint_path = checkpoint_path
        # Load model from HuggingFace Hub
        self.tokenizer = AutoTokenizer.from_pretrained(self.checkpoint_path)
        self.model = AutoModel.from_pretrained(self.checkpoint_path)

    def encode(self, sentence):
        # Tokenize sentences
        encoded_input = self.tokenizer([sentence], padding=True, truncation=True, return_tensors='pt')
        # Compute token embedding
        with torch.no_grad():
            model_output = self.model(**encoded_input)
        # Perform pooling
        sentence_embedding = mean_pooling(model_output, encoded_input['attention_mask'])
        # Normalize embedding
        sentence_embedding = F.normalize(sentence_embedding, p=2, dim=1)
        return sentence_embedding[0]
