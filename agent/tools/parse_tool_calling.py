import ast

def parse_function_calling(call_str):
    # 替换路径分隔符，使其可被解析为模块和函数调用
    if '/' in call_str:
        call_str = call_str.replace('/', '.')
    try:
        # 解析字符串为抽象语法树
        tree = ast.parse(call_str, mode='eval')
        # 检查是否为函数调用
        if isinstance(tree.body, ast.Call):
            func = tree.body
            # 获取函数名
            if isinstance(func.func, ast.Attribute):
                func_name = f"{func.func.value.id}/{func.func.attr}"
            elif isinstance(func.func, ast.Name):
                func_name = func.func.id
            else:
                raise ValueError("Unsupported function format.")
            
            # 获取参数
            args = [ast.literal_eval(arg) for arg in func.args]
            kwargs = {kw.arg: ast.literal_eval(kw.value) for kw in func.keywords}

            return func_name, args, kwargs
        else:
            raise ValueError("Not a function call.")
    except Exception as e:
        raise ValueError(f"Invalid function call string: {e}")