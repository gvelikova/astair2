import gzip
import logging
import random
import re
from concurrent.futures import ProcessPoolExecutor
from itertools import chain, islice
from os import path

import click
import pysam

from astair2.output import output_directory, tab_line
from astair2.plotting import plot_colors, pyplot
from astair2.safe_division import non_zero_division_NA
from astair2.statistics_summary import general_statistics_summary

logs = logging.getLogger(__name__)

BASES = ("T", "C", "A", "G")
DEFAULT_COLORS = ["skyblue", "mediumaquamarine", "khaki", "lightcoral"]
DEFAULT_SAMPLE_SIZE = 10000000


@click.command()
@click.option(
    "fq1",
    "--fq1",
    "-1",
    required=True,
    help="First in pair (R1) sequencing reads file in fastq.gz format.",
)
@click.option(
    "fq2",
    "--fq2",
    "-2",
    required=False,
    help="Second in pair (R2) sequencing reads file in fastq.gz format.",
)
@click.option(
    "calculation_mode",
    "--calculation_mode",
    "-cm",
    required=False,
    default="means",
    type=click.Choice(["means", "absolute"]),
    help="Gives the mode of computation used for the Phred scores summary, where means runs faster. (Default means)",
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
@click.option(
    "sample_size",
    "--sample_size",
    "-s",
    default=10000000,
    type=int,
    required=False,
    help="The number of reads to sample for the analysis. (Default 10 000 000).",
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
    "minimum_score",
    "--minimum_score",
    "-q",
    required=False,
    default=15,
    type=int,
    help="Minimum Phred score used for visualisation only. (Default 15).",
)
@click.option(
    "colors",
    "--colors",
    "-c",
    default=["skyblue", "mediumaquamarine", "khaki", "lightcoral"],
    type=list,
    required=False,
    help="List of color values used for visualistion of A, C, G, T, they are given as color1,color2,color3,color4. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default skyblue,mediumaquamarine,khaki,lightcoral).",
)
def phred(
    fq1,
    fq2,
    calculation_mode,
    directory,
    sample_size,
    minimum_score,
    colors,
    plot,
    single_end,
):
    """Calculate per base (A, C, T, G) Phred scores for each strand."""
    Phred_scores_plotting(
        fq1,
        fq2,
        calculation_mode,
        directory,
        sample_size,
        minimum_score,
        colors,
        plot,
        single_end,
    )


def numeric_Phred_score(score):
    """Converts ASCII fastq sequencing scores to numeric scores."""
    return list(pysam.qualitystring_to_array(score))


def fastq_records(lines):
    """Yields (sequence, quality string) for each four-line FASTQ record."""
    lines = iter(lines)
    while True:
        record = list(islice(lines, 4))
        if len(record) < 4:
            return
        yield record[1].rstrip("\r\n"), record[3].rstrip("\r\n")


def read_base_qualities(sequence, qualities):
    """The Phred scores of the T, C, A and G bases of one read."""
    scores = numeric_Phred_score(qualities)
    return tuple(
        [score for score, base in zip(scores, sequence) if base == wanted]
        for wanted in BASES
    )


def read_summary(base_scores, calculation_mode):
    """Per read: the mean score per base ('means') or all scores per base ('absolute')."""
    if calculation_mode == "means":
        return tuple(
            non_zero_division_NA(sum(scores), len(scores)) for scores in base_scores
        )
    return base_scores


def sample(values, sample_size):
    """Keeps the first sample_size values, then lets each of the next 9 * sample_size values replace a
    random kept one, and ignores the rest."""
    kept = []
    for index, value in enumerate(values):
        if index < sample_size:
            kept.append(value)
        elif index < 10 * sample_size:
            kept[random.randint(0, len(kept) - 1)] = value
        else:
            break
    return kept


def Phred_score_main_body(file_to_load, calculation_mode, sample_size):
    """Per-read Phred score summaries of a FASTQ file, sampled down to sample_size reads."""
    cutoff = DEFAULT_SAMPLE_SIZE if sample_size is None else int(sample_size)
    summaries = (
        read_summary(read_base_qualities(sequence, qualities), calculation_mode)
        for sequence, qualities in fastq_records(file_to_load)
    )
    return sample(summaries, cutoff)


def fastq_read_values(fastq_file, sample_size, calculation_mode):
    with gzip.open(fastq_file, "rt") as handle:
        return Phred_score_main_body(handle, calculation_mode, sample_size)


def base_score_lists(read_values, calculation_mode):
    """All scores per base (T, C, A, G) over the sampled reads, leaving out reads without that base."""
    per_base = [
        [values[index] for values in read_values if values[index] != "NA"]
        for index in range(len(BASES))
    ]
    if calculation_mode == "absolute":
        per_base = [list(chain.from_iterable(scores)) for scores in per_base]
    return dict(zip(BASES, per_base))


def summary_rows(scores, title):
    """The rows of one section of the _total_Phred.txt file."""
    statistics = {
        base: general_statistics_summary(scores[base]) for base in ("A", "C", "T", "G")
    }
    names = {"A": "adenines", "C": "cytosines", "T": "thymines", "G": "guanines"}
    yield [
        "____________________________________{}____________________________________".format(
            title
        )
    ]
    yield [
        "______________________________________________________________________________________"
    ]
    for label, index in (
        ("mean ", 0),
        ("median ", 1),
        ("q25 ", 3),
        ("q75 ", 4),
        ("sd ", 2),
        ("min ", 5),
        ("max ", 6),
    ):
        yield [label] + [
            "{}: {}".format(names[base], statistics[base][index])
            for base in ("A", "C", "T", "G")
        ]
    yield [
        "______________________________________________________________________________________"
    ]


def output_name(fq1, paired):
    name = path.splitext(path.basename(fq1))[0]
    return re.sub("_(R1|1).fq", "", name) if paired else re.sub(".fq", "", name)


def plot_phred(file_name, scores, minimum_score, colors):
    """Box plots of the Phred scores per base for each read file."""
    plt = pyplot()
    if plt is None:
        return
    import matplotlib.ticker as ticker

    fig, axes = plt.subplots(1, len(scores))
    axes = list(axes) if len(scores) > 1 else [axes]
    fig.suptitle("Sequencing base quality", fontsize=14)
    if len(scores) > 1:
        plt.subplots_adjust(wspace=0.4)
    maxy = max(35, max(max(chain.from_iterable(s.values())) for s in scores) + 1)
    labels = (
        ["Single-end read"] if len(scores) == 1 else ["First in pair", "Second in pair"]
    )
    for index, (axis, file_scores, label) in enumerate(zip(axes, scores, labels)):
        box = axis.boxplot(
            [file_scores[base] for base in ("A", "C", "G", "T")],
            tick_labels=["A", "C", "G", "T"],
            patch_artist=True,
        )
        if index == 0:
            axis.set_ylabel("Phred score", fontsize=12)
        axis.set_xlabel(label, fontsize=12)
        axis.axis([0, 5, minimum_score, maxy])
        axis.yaxis.set_major_locator(ticker.MultipleLocator(5))
        axis.grid(color="lightgray", linestyle="solid", linewidth=1)
        for patch, color in zip(box["boxes"], colors):
            patch.set(color="black", linewidth=1)
            patch.set_facecolor(color)
        for line in box["whiskers"] + box["medians"]:
            line.set(color="black", linewidth=1)
        for flier in box["fliers"]:
            flier.set(
                marker="", markersize=1, markerfacecolor=None, markeredgecolor=None
            )
    if len(scores) > 1:
        plt.vlines(
            -1,
            minimum_score,
            maxy,
            alpha=0.3,
            linewidth=1,
            linestyle="--",
            color="gray",
            clip_on=False,
        )
    plt.savefig(file_name, dpi=330, bbox_inches="tight")
    plt.close()


def Phred_scores_plotting(
    fq1,
    fq2,
    calculation_mode,
    directory,
    sample_size,
    minimum_score,
    colors,
    plot,
    single_end,
):
    """Writes summary statistics of the Phred scores per base for one or two FASTQ files, and optionally plots them."""
    logs.info("asTair's Phred scores statistics summary function started running.")
    paired = not single_end or fq2 is not None
    fastq_files = [fq1, fq2] if paired else [fq1]
    for fastq_file in fastq_files:
        if fastq_file is None or not path.isfile(fastq_file):
            raise FileNotFoundError(
                "The input fastq file {} does not exist.".format(fastq_file)
            )
    name = path.join(output_directory(directory), output_name(fq1, paired))
    with ProcessPoolExecutor(max_workers=len(fastq_files)) as executor:
        read_values = list(
            executor.map(
                fastq_read_values,
                fastq_files,
                [sample_size] * len(fastq_files),
                [calculation_mode] * len(fastq_files),
            )
        )
    scores = [base_score_lists(values, calculation_mode) for values in read_values]
    titles = ["First in pair_", "Second in pair"] if paired else ["__Single-end__"]
    with open(name + "_total_Phred.txt", "w") as summary:
        for file_scores, title in zip(scores, titles):
            summary.writelines(
                tab_line(row) for row in summary_rows(file_scores, title)
            )
    if plot:
        try:
            plot_phred(
                name + "_phred_scores_plot.pdf",
                scores,
                minimum_score,
                plot_colors(colors, DEFAULT_COLORS),
            )
        except Exception:
            logs.error("asTair cannot output the Phred scores plot.", exc_info=True)
    logs.info("asTair's Phred scores statistics summary function finished running.")
