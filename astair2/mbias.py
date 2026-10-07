import logging
from math import ceil
from os import path

import click
import numpy

from astair2.bam_file_parser import iterate_reads
from astair2.output import input_name, no_information_value, output_directory, tab_line
from astair2.plotting import plot_colors, pyplot
from astair2.read_context import read_cytosine_calls
from astair2.safe_division import non_zero_division, safe_rounder
from astair2.simple_fasta_parser import read_reference

logs = logging.getLogger(__name__)

CONTEXTS = ("CpG", "CHG", "CHH")
DEFAULT_COLORS = ["teal", "gray", "maroon"]


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
    "directory",
    "--directory",
    "-d",
    required=True,
    help="Output directory to save files.",
)
@click.option(
    "read_length",
    "--read_length",
    "-l",
    type=int,
    required=True,
    help="The read length is needed to calculate the Mbias.",
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
    "per_chromosome",
    "--per_chromosome",
    "-chr",
    type=str,
    help="When used, it calculates the modification rates only per the chromosome given.",
)
@click.option(
    "single_end",
    "--single_end",
    "-se",
    default=False,
    is_flag=True,
    required=False,
    help="Indicates single-end sequencing reads (Default False).",
)
@click.option(
    "plot",
    "--plot",
    "-p",
    required=False,
    is_flag=True,
    help="Phred scores will be visualised and output as a pdf file. Requires installed matplotlib.",
)
@click.option(
    "colors",
    "--colors",
    "-c",
    default=["teal", "gray", "maroon"],
    type=list,
    required=False,
    help="List of color values used for visualistion of CpG, CHG and CHH modification levels per read, which are given as color1,color2,color3. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default 'teal','gray','maroon').",
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
    "no_information",
    "--no_information",
    "-ni",
    default="0",
    type=click.Choice([".", "0", "*", "NA"]),
    required=False,
    help="What symbol should be used for a value where no enough quantative information is used. (Default *).",
)
def mbias(
    reference,
    input_file,
    directory,
    read_length,
    method,
    single_end,
    plot,
    colors,
    N_threads,
    per_chromosome,
    no_information,
):
    """Generate modification per read length information (Mbias). This is a quality-control measure."""
    Mbias_plotting(
        reference,
        input_file,
        directory,
        read_length,
        method,
        single_end,
        plot,
        colors,
        N_threads,
        per_chromosome,
        no_information,
    )


def read_orientations(single_end):
    """(flags of first reads or top-strand reads, flags of second reads or bottom-strand reads)."""
    return ((0,), (16,)) if single_end else ((99, 83), (147, 163))


def empty_counts(read_length):
    """Per read position counts of modified and unmodified cytosines, per context."""
    return {
        (context, modified): [0] * read_length
        for context in CONTEXTS
        for modified in (True, False)
    }


def mbias_counts(reads, sequences, read_length, method, single_end, per_chromosome):
    """Per-position modified and unmodified cytosine counts for both read orientations."""
    first, second = read_orientations(single_end)
    top_flags = (0,) if single_end else (99, 147)
    counts = (empty_counts(read_length), empty_counts(read_length))
    for read in reads:
        if read.reference_length == 0 or read.flag not in first + second:
            continue
        if per_chromosome is not None and read.reference_name != per_chromosome:
            continue
        calls = read_cytosine_calls(
            read, sequences[read.reference_name], read.flag in top_flags, method
        )
        if len(read.query_sequence) <= read_length:
            orientation_counts = counts[0 if read.flag in first else 1]
            for read_index, call in calls.items():
                orientation_counts[call][read_index] += 1
    return counts


def combined_single_end(counts):
    """Single-end reads of both strands are reported together as the first orientation."""
    top, bottom = counts
    return {key: [a + b for a, b in zip(top[key], bottom[key])] for key in top}, bottom


def mbias_columns(counts, no_information):
    """(modification levels, unmodified counts, modified counts) per context, as lists over read positions.

    TODO known bug: the counts are multiplied by 100 like the levels."""

    def columns(context):
        modified, unmodified = counts[(context, True)], counts[(context, False)]
        levels = [
            safe_rounder(non_zero_division(m, u + m, no_information), 3, True)
            for m, u in zip(modified, unmodified)
        ]
        return (
            levels,
            [safe_rounder(u, 3, True) for u in unmodified],
            [safe_rounder(m, 3, True) for m in modified],
        )

    return {context: columns(context) for context in CONTEXTS}


def mbias_header(single_end):
    first, second = ("OT", "OB") if single_end else ("1", "2")
    return ["#POSITION_(bp)"] + [
        "{}_{}_READ_{}".format(column, context, read)
        for context in CONTEXTS
        for read in (first, second)
        for column in ("MOD_LVL", "UNMOD_COUNT", "MOD_COUNT")
    ]


def mbias_rows(columns, read_length):
    first, second = columns
    for position in range(read_length):
        yield [position + 1] + [
            values[position]
            for context in CONTEXTS
            for orientation in (first, second)
            for values in orientation[context]
        ]


def plot_mbias(file_name, levels, read_length, single_end, colors):
    """Line plot of the modification level per read position and context."""
    plt = pyplot()
    if plt is None:
        return
    x_axis = list(range(1, read_length + 1))
    plt.figure()
    fig, axes = plt.subplots(1 if single_end else 2, 1)
    fig.suptitle("Sequencing M-bias", fontsize=14)
    panels = (
        [(axes, "Base positions", levels[0])]
        if single_end
        else [
            (axes[0], "First in pair base positions", levels[0]),
            (axes[1], "Second in pair base positions", levels[1]),
        ]
    )
    if not single_end:
        plt.subplots_adjust(hspace=0.4)
        plt.subplots_adjust(right=1)
    for axis, label, orientation_levels in panels:
        axis.set_ylabel("Modification level, %", fontsize=12)
        axis.set_xlabel(label, fontsize=12)
        for context, color in zip(CONTEXTS, colors):
            axis.plot(
                x_axis,
                orientation_levels[context][0],
                linewidth=1.0,
                linestyle="-",
                color=color,
            )
        axis.xaxis.set_ticks(
            numpy.arange(0, read_length + 1, step=ceil(read_length / 10))
        )
        axis.yaxis.set_ticks(numpy.arange(0, 101, step=10))
        axis.grid(color="lightgray", linestyle="solid", linewidth=1)
    plt.figlegend(
        list(CONTEXTS),
        loc="center left",
        bbox_to_anchor=(0.9, 0.5) if single_end else (1, 0.5),
    )
    plt.savefig(file_name, dpi=330, bbox_inches="tight", pad_inches=0.15)
    plt.close()


def Mbias_plotting(
    reference,
    input_file,
    directory,
    read_length,
    method,
    single_end,
    plot,
    colors,
    N_threads,
    per_chromosome,
    no_information,
):
    """Writes the modification level per read position, read orientation and context (M-bias) and
    optionally plots it."""
    logs.info("asTair's M-bias summary function started running.")
    name = path.join(output_directory(directory), input_name(input_file))
    sequences = read_reference(reference)
    counts = mbias_counts(
        iterate_reads(input_file, N_threads),
        sequences,
        read_length,
        method,
        single_end,
        per_chromosome,
    )
    if single_end:
        counts = combined_single_end(counts)
    columns = [
        mbias_columns(orientation, no_information_value(no_information))
        for orientation in counts
    ]
    with open(name + "_Mbias.txt", "w") as stats:
        stats.write(tab_line(mbias_header(single_end)))
        stats.writelines(tab_line(row) for row in mbias_rows(columns, read_length))
    if plot:
        try:
            plot_mbias(
                name + "_M-bias_plot.pdf",
                columns,
                read_length,
                single_end,
                plot_colors(colors, DEFAULT_COLORS),
            )
        except Exception:
            logs.error("asTair cannot output the Mbias plot.", exc_info=True)
    logs.info("asTair's M-bias summary function finished running.")
