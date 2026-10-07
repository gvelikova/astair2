import logging
import re
import sys
from os import path

import click
import pysam

from astair2.bam_file_parser import iterate_reads, open_alignments
from astair2.cigar_search import correct_positions_for_cigar
from astair2.output import input_name, output_directory
from astair2.simple_fasta_parser import read_reference

logs = logging.getLogger(__name__)

# read flag -> (reference base, converted base, offset of the cytosine within a CG)
STRANDS = {
    99: ("C", "T", 0),
    147: ("C", "T", 0),
    0: ("C", "T", 0),
    83: ("G", "A", 1),
    163: ("G", "A", 1),
    16: ("G", "A", 1),
}


@click.command()
@click.option(
    "reference",
    "--reference",
    "-f",
    required=True,
    help="Reference DNA sequence in FASTA format.",
)
@click.option(
    "input_file",
    "--input_file",
    "-i",
    required=True,
    help="BAM|CRAM format file containing sequencing reads.",
)
@click.option(
    "method",
    "--method",
    "-m",
    required=False,
    default="mCtoT",
    type=click.Choice(["CtoT", "mCtoT"]),
    help="Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).",
)
@click.option(
    "bases_noncpg",
    "--bases_noncpg",
    default=3,
    type=int,
    help="The number of cytosines conversion events in CpH content to consider the read for removal. Default value is 3.",
)
@click.option(
    "per_chromosome",
    "--per_chromosome",
    "-chr",
    default=None,
    type=str,
    help="When used, it calculates the modification rates only per the chromosome given. (Default None).",
)
@click.option(
    "N_threads",
    "--N_threads",
    "-t",
    default=1,
    required=True,
    help="The number of threads to spawn (Default 1).",
)
@click.option(
    "single_end",
    "--se",
    "-se",
    default=False,
    is_flag=True,
    required=False,
    help="Indicates single-end sequencing reads (Default False).",
)
@click.option(
    "directory",
    "--directory",
    "-d",
    required=True,
    type=str,
    help="Output directory to save files.",
)
def filter(
    reference,
    input_file,
    method,
    bases_noncpg,
    per_chromosome,
    N_threads,
    single_end,
    directory,
):
    """Look for sequencing reads with more than N CpH modifications."""
    removing_mod_err(
        reference,
        input_file,
        method,
        bases_noncpg,
        per_chromosome,
        N_threads,
        single_end,
        directory,
    )


def base_positions(base, sequence):
    """Indices of a base in a sequence, ignoring case."""
    return [match.start() for match in re.finditer(base, sequence, re.IGNORECASE)]


def has_indels_or_clips(read):
    # TODO check: by operator precedence, insertions and deletions only count for reads with tags.
    return (
        (
            len(read.get_tags()) != 0
            and ("I" in read.cigarstring or "D" in read.cigarstring)
        )
        or "S" in read.cigarstring
        or "H" in read.cigarstring
    )


def non_cpg_conversions(read, reference_sequence, method):
    """Read positions of non-CpG cytosines that look modified (mCtoT) or unmodified (CtoT)."""
    reference, converted, cpg_offset = STRANDS[read.flag]
    window = reference_sequence[
        read.reference_start : read.reference_start + read.query_alignment_length
    ]
    cpg = [index + cpg_offset for index in base_positions("CG", window)]
    references = base_positions(reference, window)
    if has_indels_or_clips(read):
        references = correct_positions_for_cigar(read.cigarstring, references)
        cpg = correct_positions_for_cigar(read.cigarstring, cpg)
    conversions = set(base_positions(converted, read.query_sequence))
    if method == "mCtoT":
        return (set(references) & conversions) - set(cpg)
    return set(references) - set(cpg) - conversions


def is_high_cph(read, sequences, method, bases_noncpg):
    """Whether a read has at least bases_noncpg non-CpG conversion events."""
    reference = STRANDS[read.flag][0]
    # mismatches against a reference C (or G) are a quick necessary condition
    if read.get_tag("MD").upper().count(reference) < bases_noncpg:
        return False
    return (
        len(non_cpg_conversions(read, sequences[read.reference_name], method))
        >= bases_noncpg
    )


def removing_mod_err(
    reference,
    input_file,
    method,
    bases_noncpg,
    per_chromosome,
    N_threads,
    single_end,
    directory,
):
    """Splits the reads into those with fewer and those with at least bases_noncpg non-CpG
    conversions. Reads that are not primary, properly paired alignments are left out."""
    logs.info(
        "asTair's excessive non-CpG modification read removal function started running."
    )
    directory = output_directory(directory)
    name = path.join(directory, input_name(input_file))
    # TODO check: per_chromosome has no effect.
    sequences = read_reference(reference)
    read_flags = (0, 16) if single_end else (83, 99, 147, 163)
    with (
        open_alignments(input_file, N_threads) as template,
        pysam.AlignmentFile(
            name + "_high_CpH_filtered.bam", "wb", template=template
        ) as kept,
        pysam.AlignmentFile(
            name + "_high_CpH_removed.bam", "wb", template=template
        ) as removed,
    ):
        for read in iterate_reads(input_file, N_threads):
            if read.flag not in read_flags:
                continue
            try:
                high_cph = is_high_cph(read, sequences, method, bases_noncpg)
            except KeyError:
                logs.error(
                    "The input file does not contain a MD tag column.", exc_info=True
                )
                sys.exit(1)
            (removed if high_cph else kept).write(read)
    logs.info(
        "asTair's excessive non-CpG modification read removal function finished running."
    )
