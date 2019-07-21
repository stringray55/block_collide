def wrap(func,*args,**kwargs):
    def wrapped():
        return func(*args,**kwargs)
    return wrapped
