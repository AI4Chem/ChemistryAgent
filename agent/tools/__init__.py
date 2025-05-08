from decimal import Decimal, getcontext

def post_process_execution_result(execution_result):
    if isinstance(execution_result, float):
        # 设置上下文精度
        getcontext().prec = 10
        number = Decimal(str(execution_result))
        approx = number.normalize()
        return float(approx)
    else:
        return execution_result

from .toolpool import ToolPool
from .toolset import ToolSet