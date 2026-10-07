import logging
from functools import lru_cache
from collections import Counter
from itertools import zip_longest
from os import path
from typing import TypeAlias

import click

from astair2.bam_file_parser import open_alignments
from astair2.context_search import USER_CONTEXT, context_keys, find_cytosine_contexts
from astair2.output import (
    exit_if_exists,
    input_name,
    no_information_value,
    open_text,
    output_directory,
    tab_line,
)
from astair2.safe_division import non_zero_division, safe_rounder
from astair2.simple_fasta_parser import (
    read_reference,
    reference_names,
    write_reference_with_underscores,
)
from astair2.vcf_reader import read_vcf

logs = logging.getLogger(__name__)

MODS_HEADER = (
    "#CHROM",
    "START",
    "END",
    "MOD_LEVEL",
    "MOD",
    "UNMOD",
    "REF",
    "ALT",
    "SPECIFIC_CONTEXT",
    "CONTEXT",
    "SNV",
    "TOTAL_DEPTH",
)
STATS_HEADER = (
    "#CONTEXT",
    "SPECIFIC_CONTEXT",
    "MEAN_MODIFICATION_RATE_PERCENT",
    "TOTAL_POSITIONS",
    "COVERED_POSITIONS",
    "MODIFIED",
    "UNMODIFIED",
)

SPECIFIC_CONTEXTS = (
    "CAG",
    "CCG",
    "CTG",
    "CTT",
    "CCT",
    "CAT",
    "CTA",
    "CTC",
    "CAC",
    "CAA",
    "CCA",
    "CCC",
    "CGA",
    "CGT",
    "CGC",
    "CGG",
)

# (row label, motif key of the total counts, specific contexts) per --context choice
StatisticsGroupRow: TypeAlias = tuple[str, str, tuple[str, ...]]
StatisticsGroup: TypeAlias = tuple[StatisticsGroupRow, ...]

STATISTICS_GROUPS: dict[str, StatisticsGroup] = {
    "CpG": (("CpG", "CG", ("CGA", "CGC", "CGG", "CGT")),),
    "CHG": (("CHG", "CHG", ("CAG", "CCG", "CTG")),),
    "CHH": (
        ("CHH", "CHH", ("CTT", "CAT", "CCT", "CTA", "CAA", "CCA", "CTC", "CAC", "CCC")),
    ),
}
STATISTICS_GROUPS["all"] = (
    STATISTICS_GROUPS["CpG"]
    + STATISTICS_GROUPS["CHG"]
    + STATISTICS_GROUPS["CHH"]
    + (("CNN", "CN", ()),)
)

# reference base -> (expected bases in the order unmodified, T, opposite base, A; converted base)
STRAND_BASES = {"C": (("C", "T", "G", "A"), "T"), "G": (("G", "T", "C", "A"), "A")}


class PositionDropped(Exception):
    """A covered position that is left out of the output (see the TODOs where it is raised)."""


@click.command()
@click.option(
    "input_file",
    "--input_file",
    "-i",
    required=True,
    help="BAM|CRAM format file containing sequencing reads.",
)
# @click.option('control_file', '--control_file', '-c', required=False, help='BAM|CRAM format file containing sequencing reads used as a matched control.')
@click.option(
    "known_snp",
    "--known_snp",
    "-ks",
    required=False,
    help="VCF format file containing genotyped WGS high quality variants or known common variants in VCF format (dbSNP, 1000 genomes, etc.).",
)
@click.option(
    "model",
    "--model",
    "-mo",
    default="none",
    type=click.Choice(["none"]),
    required=False,
    help="Decide on model for class estimation.",
)
@click.option(
    "reference",
    "--reference",
    "-f",
    required=True,
    help="Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.",
)
# @click.option('exclude_variants', '--exclude_variants', '-ev', default=False, is_flag=True, help='When set to true does not output variants.)
@click.option(
    "zero_coverage",
    "--zero_coverage",
    "-zc",
    default=False,
    is_flag=True,
    help="When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.",
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
    required=False,
    type=str,
    help="At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.",
)
@click.option(
    "library",
    "--library",
    "-li",
    required=False,
    default="directional",
    type=click.Choice(["directional", "reverse"]),
    help="Provides information for the library preparation protocol (Default directional).",
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
    "skip_clip_overlap",
    "--skip_clip_overlap",
    "-sc",
    required=False,
    default=True,
    is_flag=True,
    help="Random removal of overlapping bases between pair-end reads. Skipping is recommended for pair-end libraries, unless the overlaps are removed prior to calling. (Default True)",
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
    "minimum_base_quality",
    "--minimum_base_quality",
    "-bq",
    required=False,
    type=int,
    default=20,
    help="Set the minimum base quality for a read base to be used in the pileup (Default 20).",
)
@click.option(
    "minimum_mapping_quality",
    "--minimum_mapping_quality",
    "-mq",
    required=False,
    type=int,
    default=0,
    help="Set the minimum mapping quality for a read to be used in the pileup (Default 0).",
)
@click.option(
    "adjust_acapq_threshold",
    "--adjust_capq_threshold",
    "-amq",
    required=False,
    type=int,
    default=0,
    help="Used to adjust the mapping quality with default 0 for no adjustment and a recommended value for adjustment 50. (Default 0).",
)
@click.option(
    "add_indels",
    "--add_indels",
    "-ai",
    required=False,
    default=True,
    type=bool,
    help="Adds inserted bases and Ns for base skipped from the reference (Default True).",
)
@click.option(
    "redo_baq",
    "--redo_baq",
    "-rbq",
    required=False,
    default=False,
    type=bool,
    help="Re-calculates per-Base Alignment Qualities ignoring existing base qualities (Default False).",
)
@click.option(
    "compute_baq",
    "--compute_baq",
    "-cbq",
    required=False,
    default=True,
    type=bool,
    help="Performs re-alignment computing of per-Base Alignment Qualities (Default True).",
)
@click.option(
    "ignore_orphans",
    "--ignore_orphans",
    "-io",
    required=False,
    default=True,
    type=bool,
    help="Ignore reads not in proper pairs (Default True).",
)
@click.option(
    "max_depth",
    "--max_depth",
    "-md",
    required=False,
    type=int,
    default=250,
    help="Set the maximum read depth for the pileup. Please increase the maximum value for spike-ins and other highly-covered sequences. (Default 250).",
)
@click.option(
    "per_chromosome",
    "--per_chromosome",
    "-chr",
    type=str,
    help="When used, it calculates the modification rates only per the chromosome given.",
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
@click.option(
    "add_underscores",
    "--add_underscores",
    "-au",
    default=False,
    is_flag=True,
    required=False,
    help="Indicates outputting a new reference fasta file with added underscores in the sequence names that is afterwards used for calling. (Default False).",
)
@click.option(
    "no_information",
    "--no_information",
    "-ni",
    default="*",
    type=click.Choice([".", "0", "*", "NA"]),
    required=False,
    help="What symbol should be used for a value where no enough quantative information is used. (Default *).",
)
@click.option(
    "start_clip",
    "--start_clip",
    "-scl",
    required=False,
    type=int,
    default=0,
    help="Set the length of the bases in the start (5 prime) of the reads that will not be used for calling. (Default 0).",
)
@click.option(
    "end_clip",
    "--end_clip",
    "-ecl",
    required=False,
    type=int,
    default=0,
    help="Set the length of the bases in the end (3 prime) of the reads that will not be used for calling. (Default 0).",
)
def call(
    input_file,
    known_snp,
    model,
    reference,
    context,
    zero_coverage,
    skip_clip_overlap,
    minimum_base_quality,
    user_defined_context,
    library,
    method,
    minimum_mapping_quality,
    adjust_acapq_threshold,
    add_indels,
    redo_baq,
    compute_baq,
    ignore_orphans,
    max_depth,
    per_chromosome,
    N_threads,
    directory,
    compress,
    single_end,
    add_underscores,
    no_information,
    start_clip,
    end_clip,
):
    """Call modified cytosines from a bam or cram file. The output consists of two files, one containing modification counts per nucleotide, the other providing genome-wide statistics per context."""
    cytosine_modification_finder(
        input_file,
        known_snp,
        model,
        reference,
        context,
        zero_coverage,
        skip_clip_overlap,
        minimum_base_quality,
        user_defined_context,
        library,
        method,
        minimum_mapping_quality,
        adjust_acapq_threshold,
        add_indels,
        redo_baq,
        compute_baq,
        ignore_orphans,
        max_depth,
        per_chromosome,
        N_threads,
        directory,
        compress,
        single_end,
        add_underscores,
        no_information,
        start_clip,
        end_clip,
    )


def strand_flags(reference, single_end, library):
    """(flags of reads informative for the modification, flags of reads from the opposite strand)."""
    top, bottom = ((0,), (16,)) if single_end else ((99, 147), (83, 163))
    informative, opposite = (top, bottom) if reference == "C" else (bottom, top)
    return (
        (informative, opposite) if library == "directional" else (opposite, informative)
    )


def base_counts(read_counts, flags, bases):
    """Number of reads with the given flags that show each of the bases."""
    return [sum(read_counts[(flag, base)] for flag in flags) for base in bases]


def most_common_base(read_counts):
    """The base of the first (flag, base) pair with the highest count."""
    return max(read_counts, key=read_counts.get)[1]


def alternative_allele(read_counts, reference):
    """The most common base if it may indicate a variant, else None."""
    if not read_counts:
        return None
    alternative = most_common_base(read_counts)
    # TODO check: T is not accepted as an alternative allele here, and A is listed twice.
    return (
        alternative
        if alternative in ("G", "A", "C", "A", "N") and alternative != reference
        else None
    )


def variant_heuristic(read_counts, opposite_flags, reference, no_information):
    """(snv label or None, alternative allele or None) from the reads of the strand opposite to the
    modification. A position is called homozygous when at least 80% of the informative opposite-strand
    reads show the alternative base."""
    bases = STRAND_BASES[reference][0]
    if not read_counts or most_common_base(read_counts) == reference:
        return None, None
    unmodified, thymines, _, adenines = base_counts(read_counts, opposite_flags, bases)
    if (
        non_zero_division(thymines, unmodified + thymines, 0) >= 0.8
        and reference == "C"
    ):
        return "homozygous", alternative_allele(read_counts, reference)
    adenine_fraction = non_zero_division(
        adenines, unmodified + adenines, no_information
    )
    if isinstance(adenine_fraction, str):
        # TODO known bug: comparing the no-information symbol with 0.8 raised a TypeError that dropped
        # the position whenever no opposite-strand read showed the reference base or A.
        raise PositionDropped
    snv = "homozygous" if adenine_fraction >= 0.8 and reference == "G" else None
    return snv, alternative_allele(read_counts, reference)


def modification_counts(
    read_counts, informative_flags, reference, method, no_information
):
    """(modification level, modified count, unmodified count) from the informative reads."""
    bases, converted_base = STRAND_BASES[reference]
    counts = dict(zip(bases, base_counts(read_counts, informative_flags, bases)))
    original, converted = counts[reference], counts[converted_base]
    if method == "mCtoT":
        modified, unmodified = converted, original
    else:
        modified, unmodified = original, converted
    level = non_zero_division(modified, modified + unmodified, no_information)
    return safe_rounder(level, 3, False), modified, unmodified


def counted_read(flag, query_position, read_length, start_clip, end_clip):
    """Whether a read base lies outside the clipped read ends."""
    if start_clip == 0 and end_clip == 0:
        return True
    if flag in (83, 99):
        return start_clip < query_position < read_length - end_clip
    if flag in (147, 163):
        return end_clip < query_position < read_length - start_clip
    return False


def pileup_read_counts(column, add_indels, start_clip, end_clip):
    """Counter of (read flag, base) over the reads of a pileup column, in order of first appearance."""
    sequences = column.get_query_sequences(
        mark_matches=False, mark_ends=False, add_indels=add_indels
    )
    reads = zip_longest(column.pileups, sequences, column.get_query_positions())
    return Counter(
        (read.alignment.flag, base.upper())
        for read, base, query_position in reads
        if counted_read(
            read.alignment.flag,
            query_position,
            read.alignment.query_length,
            start_clip,
            end_clip,
        )
    )


def call_position(column, position, information, known_variant, settings):
    """The .mods record of one covered cytosine position."""
    specific_context, context, reference = information
    if not settings["single_end"] and not settings["ignore_orphans"]:
        # TODO known bug: including orphan reads failed with a TypeError that dropped every position.
        raise PositionDropped
    informative_flags, opposite_flags = strand_flags(
        reference, settings["single_end"], settings["library"]
    )
    read_counts = pileup_read_counts(
        column, settings["add_indels"], settings["start_clip"], settings["end_clip"]
    )
    if known_variant:
        snv, alternative = "WGS_known", None
    else:
        snv, alternative = variant_heuristic(
            read_counts, opposite_flags, reference, settings["no_information"]
        )
    level, modified, unmodified = modification_counts(
        read_counts,
        informative_flags,
        reference,
        settings["method"],
        settings["no_information"],
    )
    return (
        position[0],
        position[1],
        position[1] + 1,
        level,
        modified,
        unmodified,
        reference,
        alternative or STRAND_BASES[reference][1],
        specific_context,
        context,
        snv or "No",
        column.get_num_aligned(),
    )


def call_chromosome(pileups, positions, true_variants, settings):
    """Yields the .mods records of the covered positions of one chromosome."""
    for column in pileups:
        position = (
            column.reference_name,
            column.reference_pos,
            column.reference_pos + 1,
        )
        if position not in positions:
            continue
        try:
            yield call_position(
                column,
                position,
                positions[position],
                position in true_variants,
                settings,
            )
        except PositionDropped:
            continue
        except AssertionError:
            logs.exception(
                "Failed getting query sequences (AssertionError, pysam). Please decrease the max_depth parameter."
            )
        except TypeError:
            # TODO known bug: reads with deletions at a position used with read end clipping raise a TypeError.
            continue


@lru_cache(maxsize=None)
def statistics_keys(specific_context, context, user_defined_context):
    """The summary statistics a call contributes its modified and unmodified counts to."""
    keys = [key for key in ("CpG", "CHH", "CHG") if context.startswith(key)]
    keys += [key for key in SPECIFIC_CONTEXTS if specific_context.startswith(key)]
    if context.startswith("CN"):
        keys.append("CNN")
    if user_defined_context and context.startswith(USER_CONTEXT):
        keys.append(USER_CONTEXT)
    return tuple(keys)


def add_to_statistics(statistics, record, user_defined_context):
    """Adds one .mods record to the (covered, modified, unmodified) Counters."""
    covered, modified, unmodified = statistics
    specific_context, context, snv = record[8], record[9], record[10]
    covered[context] += 1
    covered[specific_context] += 1
    if snv == "No":
        for key in statistics_keys(specific_context, context, user_defined_context):
            modified[key] += record[4]
            unmodified[key] += record[5]


def statistics_rows(
    statistics, total_counts, context, user_defined_context, no_information
):
    """The rows of the .stats file: modification rates per context and specific context."""
    covered, modified, unmodified = statistics

    def rate(key):
        return safe_rounder(
            non_zero_division(
                modified[key], modified[key] + unmodified[key], no_information
            ),
            3,
            True,
        )

    yield STATS_HEADER
    if user_defined_context:
        yield (
            user_defined_context,
            "*",
            rate(USER_CONTEXT),
            total_counts[USER_CONTEXT],
            covered[USER_CONTEXT],
            modified[USER_CONTEXT],
            unmodified[USER_CONTEXT],
        )
    for label, motif, specific_contexts in STATISTICS_GROUPS[context]:
        yield (
            label,
            "*",
            rate(label),
            total_counts[motif] + total_counts[motif + "b"],
            covered[label],
            modified[label],
            unmodified[label],
        )
        for specific in specific_contexts:
            yield (
                "*",
                specific,
                rate(specific),
                total_counts[specific],
                covered[specific],
                modified[specific],
                unmodified[specific],
            )


def cytosine_modification_finder(
    input_file,
    known_snp,
    model,
    reference,
    context,
    zero_coverage,
    skip_clip_overlap,
    minimum_base_quality,
    user_defined_context,
    library,
    method,
    minimum_mapping_quality,
    adjust_acapq_threshold,
    add_indels,
    redo_baq,
    compute_baq,
    ignore_orphans,
    max_depth,
    per_chromosome,
    N_threads,
    directory,
    compress,
    single_end,
    add_underscores,
    no_information,
    start_clip,
    end_clip,
):
    """Calls the modification level of every covered cytosine of the requested contexts and writes
    them to a .mods file, with summary statistics per context in a .stats file."""
    logs.info("asTair modification finder started running.")
    directory = output_directory(directory)
    chromosome_part = "" if per_chromosome is None else per_chromosome + "_"
    base_name = path.join(
        directory,
        input_name(input_file) + "_" + method + "_" + chromosome_part + context,
    )
    mods_file = base_name + ".mods" + (".gz" if compress else "")
    exit_if_exists(
        [base_name + ".mods", base_name + ".mods.gz"],
        "Mods file with this name exists. Please rename before rerunning.",
    )
    settings = dict(
        single_end=single_end,
        ignore_orphans=ignore_orphans,
        library=library,
        add_indels=add_indels,
        start_clip=start_clip,
        end_clip=end_clip,
        method=method,
        no_information=no_information_value(no_information),
    )
    if add_underscores:
        write_reference_with_underscores(reference)
    chromosomes = (
        reference_names(reference) if per_chromosome is None else [per_chromosome]
    )
    keys = context_keys(context, user_defined_context)
    statistics, total_counts = (Counter(), Counter(), Counter()), Counter()
    with (
        open_alignments(input_file, N_threads) as alignments,
        open_text(mods_file, compress) as mods,
    ):
        mods.write(tab_line(MODS_HEADER))
        for chromosome in chromosomes:
            logs.info(
                "Starting modification calling on {} chromosome (sequence).".format(
                    chromosome
                )
            )
            sequence = read_reference(reference, [chromosome])[chromosome]
            positions, counts = find_cytosine_contexts(
                sequence, chromosome, keys, user_defined_context
            )
            total_counts.update(counts)
            true_variants = (
                read_vcf(known_snp, chromosome, sequence, N_threads, None, None)[0]
                if known_snp is not None
                else set()
            )
            # TODO no reference sequence is passed to the pileup, so pysam never applies BAQ
            # (--compute_baq, --redo_baq) or mapping quality capping (--adjust_capq_threshold).
            pileups = alignments.pileup(
                chromosome,
                ignore_overlaps=skip_clip_overlap,
                min_base_quality=minimum_base_quality,
                stepper="samtools",
                max_depth=max_depth,
                redo_baq=redo_baq,
                ignore_orphans=ignore_orphans,
                compute_baq=compute_baq,
                min_mapping_quality=minimum_mapping_quality,
                adjust_capq_threshold=adjust_acapq_threshold,
            )
            for record in call_chromosome(pileups, positions, true_variants, settings):
                mods.write(tab_line(record))
                add_to_statistics(statistics, record, user_defined_context)
            if zero_coverage:
                # TODO known bug: writing the uncovered positions always failed.
                raise TypeError(
                    "modification_calls_writer() got an unexpected keyword argument 'header'"
                )
    with open(base_name + ".stats", "w") as stats:
        stats.writelines(
            tab_line(row)
            for row in statistics_rows(
                statistics,
                total_counts,
                context,
                user_defined_context,
                settings["no_information"],
            )
        )
    logs.info("asTair modification finder finished running.")
