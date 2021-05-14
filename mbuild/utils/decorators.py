"""Some helpful decorators."""


from warnings import warn


def deprecated(warning_string=""):
    """Decorate deprecated functions."""

    def old_function(fcn):
        def wrapper(*args, **kwargs):
            printed_message = "{0} is deprecated. {1}".format(
                fcn.__name__, warning_string
            )
            warn(printed_message, DeprecationWarning)
            fcn(*args, **kwargs)

        return wrapper

    return old_function


def breaking_change(warning_string=""):
    """Decorate functions with breaking changes."""

    def old_function(fcn):
        def wrapper(*args, **kwargs):
            printed_message = "{0} has breaking change. {1}".format(
                fcn.__name__, warning_string
            )
            warn(printed_message, Warning)
            fcn(*args, **kwargs)

        return wrapper

    return old_function
