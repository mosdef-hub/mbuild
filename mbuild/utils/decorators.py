"""Some helpful decorators."""

from functools import wraps
from warnings import warn


def make_callable(decorator):  # pragma: no cover
    """Make function callable.

    See https://stackoverflow.com/questions/3888158/making-decorators-with-
    optional-arguments
    """

    @wraps(decorator)
    def inner(*args, **kw):  # pragma: no cover
        if len(args) == 1 and not kw and callable(args[0]):
            return decorator()(args[0])
        else:
            return decorator(*args, **kw)

    return inner


def deprecated(warning_string=""):  # pragma: no cover
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


@make_callable
def deprecated_property(
    warning_msg="Property deprecated", always=False
):  # pragma: no cover
    """Alert users to deprecated properties of an object.

    Deprecation messages are shown only once per runtime by default

    Parameters
    ----------
    warning_msg : str, Property deprecated
        The deprecation warning to show
    always: bool, default=False
        If true, always show deprecation warning, default behavior shows
        deprecation warning once per runtime
    """

    def old_func(fcn):
        class Property(object):
            def __init__(self, fget=None, fset=None, fdel=None, doc=None):
                self.fget = fget
                self.fset = fset
                self.fdel = fdel
                if doc is None and fget is not None:
                    doc = fget.__doc__
                self.__doc__ = doc
                self.warning_msg = warning_msg
                self.always_show = always

            def __get__(self, obj, objtype=None):
                self.show_warning_msg()
                if obj is None:
                    return self
                if self.fget is None:
                    raise AttributeError("unreadable attribute")
                return self.fget(obj)

            def __set__(self, obj, value):
                self.show_warning_msg()
                if self.fset is None:
                    raise AttributeError("can't set attribute")
                self.fset(obj, value)

            def __delete__(self, obj):
                self.show_warning_msg()
                if self.fdel is None:
                    raise AttributeError("can't delete attribute")
                self.fdel(obj)

            def getter(self, fget):
                return type(self)(fget, self.fset, self.fdel, self.__doc__)

            def setter(self, fset):
                return type(self)(self.fget, fset, self.fdel, self.__doc__)

            def deleter(self, fdel):
                return type(self)(self.fget, self.fset, fdel, self.__doc__)

            def show_warning_msg(self):
                if self.warning_msg is not None:
                    warn(
                        f"Property {fcn.__name__} is deprecated. "
                        f"{self.warning_msg}",
                        DeprecationWarning,
                    )
                    if not self.always_show:
                        self.warning_msg = None

        prop = Property(fget=fcn)
        return prop

    return old_func


def breaking_change(warning_string=""):
    """Decorate functions with breaking changes."""

    def old_function(fcn):
        def wrapper(*args, **kwargs):
            warn(
                f"{fcn.__name__} has breaking change. {warning_string}", Warning
            )
            fcn(*args, **kwargs)

        return wrapper

    return old_function


def experimental_feature(warning_string=""):  # pragma no cover
    """Decorate experimental methods."""

    def experimental_function(fcn):
        @wraps(fcn)
        def wrapper(*args, **kwargs):
            warn(
                f"{fcn.__name__} is an experimental feature and is not subject to follow standard deprecation cycles. Use at your own risk! {warning_string}",
                Warning,
            )
            fcn(*args, **kwargs)

        return wrapper

    return experimental_function
