import ast

# def get_function_info(file_path):
#     with open(file_path, "r", encoding="utf-8") as file:
#         tree = ast.parse(file.read(), filename=file_path)

#     functions_dict = {}

#     for node in ast.walk(tree):
#         if isinstance(node, ast.FunctionDef):
#             function_info = {
#                 "name": node.name,
#                 "description": ast.get_docstring(node),
#                 "parameters": [arg.arg for arg in node.args.args],
#                 "return_type": None  # We will handle this later if possible
#             }

#             # Check if the function has a type annotation for the return type
#             if node.returns and isinstance(node.returns, ast.Name):
#                 function_info["return_type"] = node.returns.id
#             elif node.returns and isinstance(node.returns, ast.Subscript):
#                 function_info["return_type"] = ast.unparse(node.returns)

#             functions_dict[function_info["name"]] = function_info

#     return functions_dict


def get_function_info(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        tree = ast.parse(file.read(), filename=file_path)

    functions_dict = {}

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            function_info = {
                "name": node.name,
                "description": ast.get_docstring(node),
                "parameters": [],
                "return_type": None
            }

            # Check the function parameters for type annotations
            for arg in node.args.args:
                param_info = {
                    "name": arg.arg,
                    "type": None
                }
                if arg.annotation:
                    if isinstance(arg.annotation, ast.Name):
                        param_info["type"] = arg.annotation.id
                    elif isinstance(arg.annotation, ast.Subscript):
                        param_info["type"] = ast.unparse(arg.annotation)
                function_info["parameters"].append(param_info)

            # Check if the function has a type annotation for the return type
            if node.returns:
                if isinstance(node.returns, ast.Name):
                    function_info["return_type"] = node.returns.id
                elif isinstance(node.returns, ast.Subscript):
                    function_info["return_type"] = ast.unparse(node.returns)

            functions_dict[function_info["name"]] = function_info

    return functions_dict


if __name__ == "__main__":
    # Example usage
    file_path = "src/cactus/tools/calculate_sa.py"  # Replace with the path to your Python file
    function_info = get_function_info(file_path)
    print(function_info)

    exec("from cactus.tools.calculate_sa import CalculateSA")
    print(eval('''CalculateSA("C(=O)(O)[O-].[Na+]")'''))

