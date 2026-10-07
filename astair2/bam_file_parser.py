import logging
import os

import pysam

logs = logging.getLogger(__name__)


def open_alignments(input_file, threads=1):
    """Opens a coordinate-sorted BAM or CRAM file, building its index first if it has none."""
    if not os.path.isfile(input_file):
        raise FileNotFoundError("The input file {} does not exist.".format(input_file))
    index = {".bam": ".bai", ".cram": ".crai"}.get(os.path.splitext(input_file)[1])
    if index and not os.path.isfile(input_file + index):
        logs.info("Building index for the input file.")
        try:
            pysam.index(input_file)
        except pysam.SamtoolsError as error:
            raise ValueError(
                "The input file {} is not sorted by coordinates. Please sort it before running again.".format(
                    input_file
                )
            ) from error
    return pysam.AlignmentFile(input_file, "rb", threads=threads)


def iterate_reads(input_file, threads=1, region=None):
    """Yields the reads of a region (chromosome, start, end), or all reads including unmapped ones."""
    with open_alignments(input_file, threads) as alignments:
        yield from (
            alignments.fetch(*region) if region else alignments.fetch(until_eof=True)
        )
