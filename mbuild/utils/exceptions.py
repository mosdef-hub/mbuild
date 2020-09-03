class DeprecationError(Exception):
    """ Error for Deprecated mBuild functions """

    def __init__(self, depr_func, new_func, version):
        self.depr_func = depr_func
        self.new_func = new_func
        self.version = version

    def __str__(self):
        msg = f"{self.depr_func} has been deprecated. Please use {self.new_func}. Deprecated since mBuild version {self.version}."

        return msg
