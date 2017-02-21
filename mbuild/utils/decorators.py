"""Some helpful decorators"""


from warnings import warn


def deprecated(warning_string=''):
    def old_function(fcn):
        def wrapper(*args, **kwargs):
            printed_message = '{0} is deprecated. {1}'.format(fcn.__name__, warning_string)
            warn(printed_message, DeprecationWarning)
            fcn(*args, **kwargs)
        return wrapper
    return old_function
