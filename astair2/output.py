"""Helpers shared by the commands for naming, opening and writing output files."""

import gzip
import logging
import os
import sys

logs = logging.getLogger(__name__)


def output_directory(directory):
    """The absolute path of an existing output directory."""
    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError("The output directory does not exist.")
    return directory


def input_name(input_file):
    """The input file name without its directory and last extension, used to name outputs."""
    return os.path.splitext(os.path.basename(input_file))[0]


def exit_if_exists(file_names, message):
    """Stops with an error instead of overwriting existing output."""
    if any(os.path.isfile(file_name) for file_name in file_names):
        logs.error(message)
        sys.exit(1)


def no_information_value(no_information):
    """The placeholder for values that cannot be computed; '0' means the number 0."""
    return 0 if no_information == "0" else no_information


def open_text(file_name, compress=False):
    """A text handle for writing, gzip-compressed if requested."""
    if compress:
        return gzip.open(
            file_name, "wt", compresslevel=9, encoding="utf8", newline="\n"
        )
    return open(file_name, "w", newline="\n")


def tab_line(values):
    return "\t".join(map(str, values)) + "\n"
