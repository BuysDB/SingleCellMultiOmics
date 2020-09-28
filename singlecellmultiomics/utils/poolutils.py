def pool_wrapper(args):
    f, kwargs =  args
    return f(**kwargs)
