import timeit
import wrapper

def time(func,number=5000,*args,**kwargs):
    func=wrapper.wrap(func,*args,**kwargs)
    time=timeit.timeit(func,number=number)
    return time
