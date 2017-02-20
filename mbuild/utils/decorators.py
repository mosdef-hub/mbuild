"""Some helpful decorators"""


from warnings import warn


def deprecated(warning_string=None):
    def old_function(fcn):
        def wrapper(*args, **kwargs):
            printed_message = '{0} is deprecated. '.format(fcn.__name__)
            if warning_string:
                printed_message += warning_string
            warn(printed_message, DeprecationWarning)
            fcn(*args, **kwargs)
        return wrapper
    return old_function
