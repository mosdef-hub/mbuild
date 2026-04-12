"""mBuild specific exceptions."""


class MBuildError(Exception):
    """Base class for all non-trivial errors raised by mBuild."""

class PathConvergenceError(MBuildError):
    """Class for providing key information to path builders that valid points were inaccessible."""
    # TODO: Passable path or state that provides relevant information on error.

