import threading
from queue import Queue

# 需要执行的函数
def worker(t_func, task):
    t_func(task)

# 线程工作函数
def thread_worker(task_function, queue):
    t_func = task_function
    while not queue.empty():
        task = queue.get()
        worker(t_func, task)
        queue.task_done()

def work_with_multithreading(task_list, task_function, num_threads):
    # 任务列表
    tasks = task_list
    
    # 创建任务队列
    task_queue = Queue()
    for task in tasks:
        task_queue.put(task)
    
    # 控制并行线程数量
    num_threads = num_threads
    threads = []

    # 创建并启动线程
    for _ in range(num_threads):
        thread = threading.Thread(target=thread_worker, args=(task_function, task_queue,))
        thread.start()
        threads.append(thread)
    
    # 等待所有任务完成
    task_queue.join()
    
    # 确保所有线程都结束
    for thread in threads:
        thread.join()
    
    print("All tasks completed.")

if __name__ == "__main__":
    pass