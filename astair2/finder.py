import logging
from os import path

import click

from astair2.context_search import context_keys, find_cytosine_contexts
from astair2.output import (
    exit_if_exists,
    input_name,
    open_text,
    output_directory,
    tab_line,
)
from astair2.simple_fasta_parser import read_reference

logs = logging.getLogger(__name__)


@click.command()
@click.option(
    "reference",
    "--reference",
    "-f",
    required=True,
    help="Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.",
)
@click.option(
    "context",
    "--context",
    "-co",
    required=False,
    default="all",
    type=click.Choice(["all", "CpG", "CHG", "CHH"]),
    help="Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis. (Default all).",
)
@click.option(
    "user_defined_context",
    "--user_defined_context",
    "-uc",
    default=None,
    required=False,
    type=str,
    help="At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.",
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
    "compress",
    "--gz",
    "-z",
    default=False,
    is_flag=True,
    required=False,
    help="Indicates whether the mods file output will be compressed with gzip (Default False).",
)
@click.option(
    "directory",
    "--directory",
    "-d",
    required=True,
    type=str,
    help="Output directory to save files.",
)
def find(reference, context, user_defined_context, per_chromosome, compress, directory):
    """Output positions of Cs from fasta file per context."""
    find_contexts(
        reference, context, user_defined_context, per_chromosome, compress, directory
    )


def find_contexts(
    reference, context, user_defined_context, per_chromosome, compress, directory
):
    """Writes the coordinates of the cytosines of the requested contexts in a bed-like format."""
    logs.info("asTair cytosine contexts positions finder started running.")
    directory = output_directory(directory)
    chromosome_part = "" if per_chromosome is None else per_chromosome + "_"
    file_name = path.join(
        directory,
        input_name(reference)
        + "_"
        + chromosome_part
        + context
        + ".bed"
        + (".gz" if compress else ""),
    )
    exit_if_exists(
        [file_name, file_name + ".gz"],
        "Bed file with this name exists. Please rename before rerunning.",
    )
    keys = context_keys(context, user_defined_context)
    # TODO known bug: per_chromosome only changes the file name, all sequences are written.
    sequences = read_reference(reference)
    with open_text(file_name, compress) as output:
        header = ["#CHROM", "START", "END", "STRAND", "SPECIFIC_CONTEXT", "CONTEXT"]
        output.write(tab_line(header[:5] if compress else header))
        for name, sequence in sequences.items():
            logs.info(
                "Looking for cytosine positions on {} chromosome (sequence).".format(
                    name
                )
            )
            positions, _ = find_cytosine_contexts(
                sequence, name, keys, user_defined_context, stranded=True
            )
            for position in sorted(positions):
                row = position + positions[position][:2]
                # TODO known bug: the compressed output loses the CONTEXT column.
                output.write(tab_line(row[:5] if compress else row))
    logs.info("asTair cytosine contexts positions finder finished running.")
