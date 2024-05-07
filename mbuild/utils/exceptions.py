"""Exceptions for mBuild."""


class RemovedFuncError(Exception):
    """Exception for mBuild functions that have been deprecated and removed."""

    def __init__(
        self, deprecated_func, new_func, version_deprecated, version_removed
    ):
        self.deprecated_func = deprecated_func
        self.new_func = new_func
        self.version_deprecated = version_deprecated
        self.version_removed = version_removed

    def __str__(self):
        """Return the string representation of the error."""
        msg = f"""
        {self.deprecated_func} has been removed. Please use {self.new_func}.
        Deprecated since mBuild version {self.version_deprecated} and removed
        in version {self.version_removed}."""

        return msg
