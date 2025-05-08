import os
import re
import random
import string

def match_given_pattern(text,pattern):
    location = re.search(pattern,text).span()
    matched_text = text[location[0]:location[1]]
    return matched_text

def generate_random_string(length=10):
    # 生成一个包含字母和数字的字符集
    characters = string.ascii_letters + string.digits
    # 随机选择字符并生成编码
    code = ''.join(random.choice(characters) for _ in range(length))
    return code

def generate_random_filepath(base_dir="/tmp", depth=3, filename="file.txt"):
    """
    生成随机文件路径
    :param base_dir: 基础目录
    :param depth: 文件夹深度
    :param filename: 文件名
    :return: 随机文件路径
    """
    folders = [generate_random_string() for _ in range(depth)]
    file_path = os.path.join(base_dir, *folders, filename)
    return file_path