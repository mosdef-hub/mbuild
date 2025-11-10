# ruff: noqa: F401
# ruff: noqa: F403
"""mBuild: a hierarchical, component based molecule builder."""

import logging
import sys
from logging.handlers import RotatingFileHandler

from mbuild.box import Box
from mbuild.coarse_graining import coarse_grain
from mbuild.compound import *
from mbuild.conversion import load
from mbuild.coordinate_transform import *
from mbuild.lattice import Lattice
from mbuild.packing import *
from mbuild.pattern import *
from mbuild.port import Port
from mbuild.recipes import recipes

__version__ = '1.3.1'
__date__ = '2025-11-10'


class DeduplicationFilter(logging.Filter):
    """A logging filter that suppresses duplicate messages."""

    def __init__(self):
        super().__init__()
        self.logged_messages = set()

    def filter(self, record):
        log_entry = (record.name, record.levelno, record.msg)
        if log_entry not in self.logged_messages:
            self.logged_messages.add(log_entry)
            return True
        return False


class HeaderRotatingFileHandler(RotatingFileHandler):
    def __init__(
        self,
        filename,
        mode="w",
        maxBytes=0,
        backupCount=0,
        encoding=None,
        delay=False,
        header="",
    ):
        self.header = header
        super().__init__(filename, mode, maxBytes, backupCount, encoding, delay)

    def _open(self):
        """
        Open the current base log file, with the header written.
        """
        stream = super()._open()
        if stream.tell() == 0 and self.header:  # Only write header if file is empty
            stream.write(self.header + "\n")
        return stream


class mBuildLogger:
    def __init__(self):
        self.library_logger = logging.getLogger("mbuild")
        self.library_logger.setLevel(logging.DEBUG)

        # Create handlers
        self.console_handler = logging.StreamHandler(sys.stdout)
        self.console_handler.setLevel(logging.WARNING)

        # Create a formatter
        self.formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

        # Add formatter to handlers
        self.console_handler.setFormatter(self.formatter)

        # Initialize and add the deduplication filter
        self.dedup_filter = DeduplicationFilter()
        self.console_handler.addFilter(self.dedup_filter)

        # Clear any previous handlers to avoid duplicates in Jupyter
        self._clear_handlers()

        # Add handlers to the library logger
        self.library_logger.addHandler(self.console_handler)

    def _clear_handlers(self):
        handlers = self.library_logger.handlers[:]
        for handler in handlers:
            self.library_logger.removeHandler(handler)

    def debug_file(self, filename: str):
        """Print logging Debug messages to file `filename`."""
        # Get the path to the Python interpreter
        python_executable = sys.executable

        # Get the list of command-line arguments
        command_arguments = sys.argv

        # Construct the full command
        full_command = [python_executable] + command_arguments
        header = f"Log details for mBuild {__version__} from running \n{full_command}"
        self.file_handler = HeaderRotatingFileHandler(
            filename, mode="a", maxBytes=10**6, backupCount=2, header=header
        )
        self.file_handler.setLevel(logging.DEBUG)

        self.file_handler.addFilter(DeduplicationFilter())  # fresh duplication handler
        self.file_handler.setFormatter(self.formatter)
        self.library_logger.addHandler(self.file_handler)

    def print_level(self, level: str):
        """Print sys.stdout screen based on the logging `level` passed."""
        levelDict = {
            "notset": logging.NOTSET,
            "debug": logging.DEBUG,
            "info": logging.INFO,
            "warning": logging.WARNING,
            "error": logging.ERROR,
            "critical": logging.CRITICAL,
        }
        logLevel = levelDict.get(level.lower())
        if logLevel:
            self.console_handler.setLevel(logLevel)  # sets stdout
        else:
            raise ValueError(
                f"INCORRECT {level=}. Please set level of {levelDict.keys()}"
            )


# Example usage in __init__.py
mbuild_logger = mBuildLogger()
mbuild_logger.library_logger.setLevel(logging.INFO)
