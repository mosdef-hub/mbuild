"""mBuild sorting utilities."""

import re


def _atoi(text):
    return int(text) if text.isdigit() else text


def natural_sort(text):
    """Break apart a string containing letters and digits."""
    return [_atoi(a) for a in re.split(r"(\d+)", text)]
