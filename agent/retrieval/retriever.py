import numpy as np
import torch.nn.functional as F

from .rank_bm25 import BM25Okapi


class BM25_retriever:
    def __init__(self, tool_list, corpus):
        self.corpus = [self.tokenize(doc) for doc in corpus]
        self.tool_list = tool_list
        self.retriever = BM25Okapi(self.corpus)

    def tokenize(self, text):
        tokenized_text = text.replace(":","").replace("_"," ").split(" ")
        return tokenized_text

    def retrieve(self, query, n=5):
        tokenized_query = self.tokenize(query)
        doc_scores = self.retriever.get_scores(tokenized_query)
        top_n = np.argsort(doc_scores)[::-1][:n]
        return [self.tool_list[i] for i in top_n]

class DenseRetriever:
    def __init__(self, tool_list, corpus, retriever):
        self.corpus = corpus
        self.corpus_vectors = []
        self.tool_list = tool_list
        self.retriever = retriever

        print("计算 工具库 检索向量中~")
        self.calculate_corpus()
        print("工具库 检索向量 计算完毕")

    def calculate_corpus(self):
        self.corpus_vectors = []
        for sentence in self.corpus:
            self.corpus_vectors.append(self.retriever.encode(sentence))

    def retrieve(self, query, top_k=10):
        if self.corpus_vectors is []:
            raise AssertionError("Corpus of the retriever is none.")
        
        query_vector = self.retriever.encode(query)
        similarities = []
        for i, corpus_vector in enumerate(self.corpus_vectors):
            similarity = F.cosine_similarity(query_vector, corpus_vector, dim=0)
            similarities.append((i, similarity.item()))
        # 按相似度排序并获取 top k
        top_k_object = sorted(similarities, key=lambda x: x[1], reverse=True)[:top_k]
        # 返回 top k 的序号
        top_k_indices = [index for index, _ in top_k_object]
        return [self.tool_list[i] for i in top_k_indices]

    # def calculate_score(self, query):
    #     query_vector = self.retriever.encode(query)
    #     similarities = []
    #     for i, corpus_vector in enumerate(self.corpus_vectors):
    #         similarity = F.cosine_similarity(query_vector, corpus_vector, dim=0)
    #         similarities.append((i, similarity.item()))
    #     return similarities

    # def retrieve(self, query, n=10):
    #     similarities = self.retriever.calculate_score(query)
    #     # 按相似度排序并获取 top k
    #     top_k = sorted(similarities, key=lambda x: x[1], reverse=True)[:n]
    #     # 返回 top k 的序号
    #     top_k_indices = [index for index, _ in top_k]
    #     return [self.tool_list[i] for i in top_k_indices]
