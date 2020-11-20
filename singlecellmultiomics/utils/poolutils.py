def pool_wrapper(args):
    """
    Wraps a function, so it can be natively called from a Pool

    Args:
        f, kwargs : f is the function to call, kwargs a dictionary of arguments to pass to f

    Returns:
        return value of f called with **kwargs
    """
    f, kwargs =  args
    return f(**kwargs)
