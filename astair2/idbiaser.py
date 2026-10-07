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

STATS_CONTEXTS = ("CHH", "CHG", "CpG")
PLOT_CONTEXTS = ("CpG", "CHG", "CHH")
KINDS = ("insert", "deletion")
DEFAULT_COLORS = ["teal", "deepskyblue", "mediumblue", "orange", "gold", "sienna"]


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
    help="The read length is needed to calculate the IDbias.",
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
    default=["teal", "deepskyblue", "mediumblue", "orange", "gold", "sienna"],
    type=list,
    required=False,
    help="List of color values used for visualistion of CpG, CHG and CHH modification levels per read, which are given as color1,color2,color3. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default 'teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna').",
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
    help="What symbol should be used for a value where no enough quantative information is used. (Default 0).",
)
def idbias(
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
    """Generate indel count per read length information (IDbias). This is a quality-control measure."""
    IDbias_plotting(
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


def empty_profile(read_length):
    profile = {
        "reads": 0,
        "insert_reads": 0,
        "deletion_reads": 0,
        "indel_reads": 0,
        "coverage": [0] * read_length,
        "insert": [0] * read_length,
        "deletion": [0] * read_length,
    }
    profile.update(
        {
            (kind, context, modified): [0] * read_length
            for kind in KINDS
            for context in STATS_CONTEXTS
            for modified in (True, False)
        }
    )
    return profile


def read_indels(read):
    """(read indices of inserted or clipped bases, read-relative reference indices of deleted bases)."""
    if "D" not in str(read.cigarstring) and "I" not in str(read.cigarstring):
        return set(), set()
    pairs = read.get_aligned_pairs()
    return (
        {
            read_index
            for read_index, reference_index in pairs
            if reference_index is None
        },
        {
            reference_index - read.reference_start
            for read_index, reference_index in pairs
            if read_index is None
        },
    )


def add_read(profile, read, calls):
    """Adds the indels of one read to the profile of its orientation."""
    inserts, deletions = read_indels(read)
    length = len(read.query_sequence)
    profile["reads"] += 1
    profile["insert_reads"] += bool(inserts) and not deletions
    profile["deletion_reads"] += bool(deletions) and not inserts
    profile["indel_reads"] += bool(inserts) and bool(deletions)
    # TODO known bug: meant to look for modified cytosines within 5 bp of each indel, the original
    # condition `any(range(i - 5, i + 5)) in calls` tests whether read position 1 is such a cytosine.
    category = calls.get(1)
    for index in range(length):
        profile["coverage"][index] += 1
    for index in range(length):
        kind = (
            "insert" if index in inserts else "deletion" if index in deletions else None
        )
        if kind:
            profile[kind][index] += 1
            if category:
                profile[(kind,) + category][index] += 1


def idbias_profiles(reads, sequences, read_length, method, single_end, per_chromosome):
    """Indel counts per read position for both read orientations."""
    first, second = read_orientations(single_end)
    top_flags = (0,) if single_end else (99, 147)
    profiles = (empty_profile(read_length), empty_profile(read_length))
    for read in reads:
        if read.reference_length == 0 or read.flag not in first + second:
            continue
        if per_chromosome is not None and read.reference_name != per_chromosome:
            continue
        calls = read_cytosine_calls(
            read, sequences[read.reference_name], read.flag in top_flags, method
        )
        add_read(profiles[0 if read.flag in first else 1], read, calls)
    return profiles


def combined_single_end(profiles):
    """Single-end reads of both strands are reported together as the first orientation."""
    top, bottom = profiles
    return {
        key: (
            [a + b for a, b in zip(value, bottom[key])]
            if isinstance(value, list)
            else value + bottom[key]
        )
        for key, value in top.items()
    }, bottom


def unmodified_counts(profile, orientation, kind, context):
    # TODO known bug: the second reads' CHH insertion totals use the CHG unmodified counts.
    if orientation == 1 and kind == "insert" and context == "CHH":
        return profile[("insert", "CHG", False)]
    return profile[(kind, context, False)]


def totals(profile, orientation):
    """{(kind, context): (cytosines next to indels, modified ones)} over all read positions."""
    return {
        (kind, context): (
            sum(unmodified_counts(profile, orientation, kind, context))
            + sum(profile[(kind, context, True)]),
            sum(profile[(kind, context, True)]),
        )
        for kind in KINDS
        for context in STATS_CONTEXTS
    }


def stats_header():
    columns = ["total_reads", "total_insertions", "total_deletions", "total_both"] + [
        "{}_{}_{}".format(context, kind, measure)
        for kind in KINDS
        for measure in ("total", "modified")
        for context in STATS_CONTEXTS
    ]
    return ["#POS"] + [
        "{}_{}".format(read, column) for read in ("R1", "R2") for column in columns
    ]


def stats_rows(profiles, read_length, no_information):
    summary = [no_information]
    for orientation, profile in enumerate(profiles):
        orientation_totals = totals(profile, orientation)
        summary += [
            profile["reads"],
            profile["insert_reads"],
            profile["deletion_reads"],
            profile["indel_reads"],
        ]
        summary += [
            orientation_totals[(kind, context)][measure]
            for kind in KINDS
            for measure in (0, 1)
            for context in STATS_CONTEXTS
        ]
    yield summary
    for index in range(read_length):
        row = [index + 1]
        for profile in profiles:
            row += [
                profile["coverage"][index],
                profile["insert"][index],
                profile["deletion"][index],
                no_information,
            ]
            row += [
                profile[(kind, context, modified)][index]
                for kind in KINDS
                for modified in (False, True)
                for context in STATS_CONTEXTS
            ]
        yield row


def modification_rates(profile, orientation, kind, context, no_information):
    """Percentage of modified cytosines next to indels per read position."""
    modified, unmodified = profile[(kind, context, True)], unmodified_counts(
        profile, orientation, kind, context
    )
    return [
        safe_rounder(non_zero_division(m, m + u, no_information), 12, True)
        for m, u in zip(modified, unmodified)
    ]


def _panels(plt, single_end, sharey, layout):
    fig, axes = plt.subplots(
        *((1, 1) if single_end else layout), sharey=sharey and not single_end
    )
    return fig, [axes] if single_end else list(axes)


def plot_idbias(name, profiles, read_length, single_end, colors, no_information):
    """The four ID-bias plots."""
    plt = pyplot()
    if plt is None:
        return
    x_axis = list(range(1, read_length + 1))
    ticks = numpy.arange(0, read_length + 1, step=ceil(read_length / 10))
    orientations = (0,) if single_end else (0, 1)
    pair_labels = ("First in pair", "Second in pair")

    fig, axes = _panels(plt, single_end, True, (1, 2))
    fig.suptitle("Sequencing ID-bias: relative indel abundance", fontsize=14)
    if not single_end:
        plt.subplots_adjust(wspace=0.4)
        plt.subplots_adjust(right=1)
    for axis, orientation in zip(axes, orientations):
        profile = profiles[orientation]
        axis.set_ylabel("Relative abundance, %", fontsize=12)
        axis.set_xlabel(
            "Indel type" if single_end else pair_labels[orientation] + " indel type",
            fontsize=12,
        )
        axis.bar(
            [0, 1, 2],
            [
                non_zero_division(profile[key], profile["reads"], 0) * 100
                for key in ("insert_reads", "deletion_reads", "indel_reads")
            ],
            color=["lightgray", "deepskyblue", "mediumblue"],
        )
        axis.set_xticklabels(
            (["insert"] if single_end else ["", "insert"]) + ["deletion", "both"],
            fontsize=12,
        )
        axis.grid(color="lightgray", linestyle="solid", linewidth=1)
    plt.savefig(
        name + "_ID-bias_abundance_plot.pdf",
        dpi=330,
        bbox_inches="tight",
        pad_inches=0.25,
    )
    plt.close()

    bar_labels = [
        "{} {}".format(context, kind) for kind in KINDS for context in PLOT_CONTEXTS
    ]
    fig, axes = _panels(plt, single_end, True, (1, 2))
    fig.suptitle(
        "Sequencing ID-bias: relative indel abundance in 10bp from a modified cytosine",
        fontsize=14,
    )
    if not single_end:
        plt.subplots_adjust(wspace=0.4)
        plt.subplots_adjust(right=1)
    for axis, orientation in zip(axes, orientations):
        orientation_totals = totals(profiles[orientation], orientation)
        axis.set_ylabel("Relative abundance, %", fontsize=12)
        axis.set_xlabel(
            "Indel type" if single_end else pair_labels[orientation], fontsize=12
        )
        axis.bar(
            range(6),
            [
                non_zero_division(
                    orientation_totals[(kind, context)][1],
                    orientation_totals[(kind, context)][0],
                    0,
                )
                * 100
                for kind in KINDS
                for context in PLOT_CONTEXTS
            ],
            color=colors,
        )
        axis.set_xticks(range(6))
        axis.set_xticklabels(bar_labels, fontsize=12, rotation=90)
        axis.grid(color="lightgray", linestyle="solid", linewidth=1)
    plt.savefig(
        name + "_ID-bias_abundance_10bp_mod_site_plot.pdf",
        dpi=330,
        bbox_inches="tight",
        pad_inches=0.25,
    )
    plt.close()

    fig, axes = _panels(plt, single_end, True, (2, 1))
    fig.suptitle(
        "Sequencing ID-bias: indel co-localisation at 10bp from a modified position",
        fontsize=14,
    )
    if not single_end:
        plt.subplots_adjust(hspace=0.4)
        plt.subplots_adjust(right=1)
    for axis, orientation in zip(axes, orientations):
        axis.set_ylabel("Indel rate, %", fontsize=12)
        axis.set_xlabel(
            (
                "Base positions"
                if single_end
                else pair_labels[orientation] + " base positions"
            ),
            fontsize=12,
        )
        for (kind, context), color in zip(
            [(kind, context) for kind in KINDS for context in PLOT_CONTEXTS], colors
        ):
            # TODO known bug: the first-in-pair CHH lines show the second reads' rates.
            source = 1 if context == "CHH" else orientation
            axis.plot(
                x_axis,
                modification_rates(
                    profiles[source], source, kind, context, no_information
                ),
                linewidth=1.5,
                linestyle="-",
                color=color,
            )
        axis.xaxis.set_ticks(ticks)
        axis.grid(color="lightgray", linestyle="solid", linewidth=1)
    plt.figlegend(
        bar_labels,
        loc="center left",
        bbox_to_anchor=(1, 0.5) if single_end else (1.1, 0.5),
    )
    plt.savefig(
        name + "_ID-bias_modification_colocalisation_plot.pdf",
        dpi=330,
        bbox_inches="tight",
        pad_inches=0.25,
    )
    plt.close()

    fig, axes = _panels(plt, single_end, True, (2, 1))
    fig.suptitle("Sequencing ID-bias: indel rate per read length", fontsize=14)
    if not single_end:
        plt.subplots_adjust(hspace=0.4)
        plt.subplots_adjust(right=1)
    for axis, orientation in zip(axes, orientations):
        profile = profiles[orientation]
        axis.set_ylabel("Indel rate, %", fontsize=12)
        axis.set_xlabel(
            (
                "Base positions"
                if single_end
                else pair_labels[orientation] + " base positions"
            ),
            fontsize=12,
        )
        for kind, color in (("insert", colors[1]), ("deletion", colors[3])):
            axis.plot(
                x_axis,
                [
                    non_zero_division(count, coverage, no_information) * 100
                    for count, coverage in zip(profile[kind], profile["coverage"])
                ],
                linewidth=1.5,
                linestyle="-",
                color=color,
            )
        axis.xaxis.set_ticks(ticks)
        axis.grid(color="lightgray", linestyle="solid", linewidth=1)
    plt.figlegend(
        ["insert", "deletion"],
        loc="center left",
        bbox_to_anchor=(0.9, 0.5) if single_end else (1.05, 0.5),
    )
    plt.savefig(
        name + "_ID-bias_indel_rate_plot.pdf",
        dpi=330,
        bbox_inches="tight",
        pad_inches=0.25,
    )
    plt.close()


def IDbias_plotting(
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
    """Writes indel counts per read position and read orientation (ID-bias), and optionally plots them."""
    logs.info("asTair's ID-bias summary function started running.")
    name = path.join(output_directory(directory), input_name(input_file))
    no_information = no_information_value(no_information)
    sequences = read_reference(reference)
    profiles = idbias_profiles(
        iterate_reads(input_file, N_threads),
        sequences,
        read_length,
        method,
        single_end,
        per_chromosome,
    )
    if single_end:
        profiles = combined_single_end(profiles)
    chromosome_part = "" if per_chromosome is None else "_" + per_chromosome
    with open(name + chromosome_part + "_ID-bias.stats", "w") as stats:
        stats.write(tab_line(stats_header()))
        stats.writelines(
            tab_line(row) for row in stats_rows(profiles, read_length, no_information)
        )
    if plot:
        try:
            plot_idbias(
                name,
                profiles,
                read_length,
                single_end,
                plot_colors(colors, DEFAULT_COLORS),
                no_information,
            )
        except Exception:
            logs.error("asTair cannot output the IDbias plot.", exc_info=True)
    logs.info("asTair's ID-bias summary function finished running.")
