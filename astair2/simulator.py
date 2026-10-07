import gzip
import logging
import random
import re
from os import path

import click
import numpy
import pysam

from astair2.bam_file_parser import iterate_reads, open_alignments
from astair2.cigar_search import position_correction_cigar
from astair2.context_search import context_keys, find_cytosine_contexts
from astair2.output import input_name, output_directory
from astair2.safe_division import safe_rounder
from astair2.simple_fasta_parser import read_reference

logs = logging.getLogger(__name__)

TOP_FLAGS, BOTTOM_FLAGS = (99, 147), (83, 163)


@click.command()
@click.option(
    "reference",
    "--reference",
    "-f",
    required=True,
    help="Reference DNA sequence in FASTA format used for generation and modification of the sequencing reads at desired contexts.",
)
@click.option(
    "control_file",
    "--control_file",
    "-c",
    required=False,
    help="A VCF file with SNP status returned after WGS genotyping or a publicly available SNP list (dbSNP, 1000 genomes, etc).",
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
    "read_length",
    "--read_length",
    "-l",
    type=int,
    required=True,
    help="Desired length of pair-end sequencing reads.",
)
@click.option(
    "input_file",
    "--input_file",
    "-i",
    required=True,
    help="Sequencing reads as a BAM|CRAMfile or fasta sequence to generate reads.",
)
@click.option(
    "simulation_input",
    "--simulation_input",
    "-si",
    type=click.Choice(["bam"]),
    default="bam",
    required=False,
    help="Input file format according to the desired outcome. BAM|CRAM files can be generated with other WGS simulators allowing for sequencing errors and read distributions or can be real-life sequencing data.",
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
    "modification_level",
    "--modification_level",
    "-ml",
    type=int,
    required=False,
    help="Desired modification level; can take any value between 0 and 100.",
)
@click.option(
    "library",
    "--library",
    "-lb",
    type=click.Choice(["directional", "reverse"]),
    default="directional",
    required=False,
    help="Provide the correct library construction method. NB: Non-directional methods under development.",
)
@click.option(
    "modified_positions",
    "--modified_positions",
    "-mp",
    required=False,
    default=None,
    help="Provide a tab-delimited list of positions to be modified. By default the simulator randomly modifies certain positions. Please use seed for replication if no list is given.",
)
@click.option(
    "context",
    "--context",
    "-co",
    required=False,
    default="all",
    type=click.Choice(["all", "CpG", "CHG", "CHH"]),
    help="Explains which cytosine sequence contexts are to be modified in the output file. Default behaviour is all, which modifies positions in CpG, CHG, CHH contexts. (Default all).",
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
    "coverage",
    "--coverage",
    "-cv",
    required=False,
    type=int,
    help="Desired depth of sequencing coverage.",
)
@click.option(
    "region",
    "--region",
    "-r",
    nargs=3,
    type=click.Tuple([str, int, int]),
    default=(None, None, None),
    required=False,
    help="The one-based genomic coordinates of the specific region of interest given in the form chromosome, start position, end position, e.g. chr1 100 2000.",
)
@click.option(
    "overwrite",
    "--overwrite",
    "-ov",
    required=False,
    default=False,
    is_flag=True,
    help="Indicates whether existing output files with matching names will be overwritten. (Default False).",
)
@click.option(
    "per_chromosome",
    "--per_chromosome",
    "-chr",
    type=str,
    help="When used, it calculates the modification rates only per the chromosome given.",
)
@click.option(
    "GC_bias",
    "--GC_bias",
    "-gc",
    default=0.3,
    required=True,
    type=float,
    help="The value of total GC levels in the read above which lower coverage will be observed in Ns and fasta modes. (Default 0.5).",
)
@click.option(
    "sequence_bias",
    "--sequence_bias",
    "-sb",
    default=0.1,
    required=True,
    type=float,
    help="The proportion of lower-case letters in the read string for the Ns and fasta modes that will decrease the chance of the read being output. (Default 0.1).",
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
    "reverse_modification",
    "--rev",
    "-rv",
    default=False,
    is_flag=True,
    required=False,
    help="Returns possible or known modified position to their unmodified expected state. NB: Works only on files with MD tags (Default False).",
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
    "seed",
    "--seed",
    "-s",
    type=int,
    required=False,
    help="An integer number to be used as a seed for the random generators to ensure replication.",
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
    "--adjust_acapq_threshold",
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
def simulate(
    reference,
    control_file,
    model,
    read_length,
    input_file,
    simulation_input,
    method,
    modification_level,
    library,
    modified_positions,
    context,
    user_defined_context,
    coverage,
    region,
    overwrite,
    per_chromosome,
    GC_bias,
    sequence_bias,
    N_threads,
    reverse_modification,
    directory,
    seed,
    skip_clip_overlap,
    single_end,
    minimum_base_quality,
    minimum_mapping_quality,
    adjust_acapq_threshold,
    add_indels,
    redo_baq,
    compute_baq,
    ignore_orphans,
    max_depth,
):
    """Simulate TAPS/BS conversion on top of an existing bam/cram file."""
    modification_simulator(
        reference,
        control_file,
        model,
        read_length,
        input_file,
        simulation_input,
        method,
        modification_level,
        library,
        modified_positions,
        context,
        user_defined_context,
        coverage,
        region,
        overwrite,
        per_chromosome,
        GC_bias,
        sequence_bias,
        N_threads,
        reverse_modification,
        directory,
        seed,
        skip_clip_overlap,
        single_end,
        minimum_base_quality,
        minimum_mapping_quality,
        adjust_acapq_threshold,
        add_indels,
        redo_baq,
        compute_baq,
        ignore_orphans,
        max_depth,
    )


def listed_positions(positions_file, chromosome, input_file):
    """{position: [modified so far, target]} from a .mods-like file for one chromosome.

    The target is the listed modification level times the number of reads on the informative strand.
    Lines are read from the first one of the chromosome until the next chromosome starts.
    """
    opener = gzip.open if positions_file.endswith(".gz") else open
    targets = {}
    with (
        opener(positions_file, "rt") as lines,
        open_alignments(input_file) as alignments,
    ):
        started = False
        for line in lines:
            fields = line.split()
            if not fields or fields[0] != chromosome:
                if started:
                    break
                continue
            started = True
            level = fields[3]
            if re.search("[a-zA-Z]", level):
                continue
            try:
                level = float(level)
            except ValueError as error:
                raise ValueError(
                    "The modification level {} in {} is not a number; use --no_information 0 "
                    "when calling the positions.".format(level, positions_file)
                ) from error
            if level == 0:
                continue
            start, end, alternative = int(fields[1]), int(fields[2]), fields[7]
            informative = (
                TOP_FLAGS
                if alternative in ("C", "T")
                else BOTTOM_FLAGS if alternative in ("A", "G") else ()
            )
            reads = sum(
                1
                for read in alignments.fetch(contig=chromosome, start=start, stop=end)
                if read.flag in informative
            )
            targets[(chromosome, start, end)] = numpy.array(
                [0, level * reads], dtype=numpy.int8
            )
    return targets


def modification_level_transformation(modification_level, modified_positions):
    """The requested modification level as a fraction, 'user_provided_list' for a positions file, or None for 0."""
    if modified_positions is not None:
        return "user_provided_list"
    return modification_level / 100 if modification_level != 0 else None


def level_label(modification_level):
    return (
        modification_level
        if isinstance(modification_level, str)
        else int(modification_level * 100)
    )


def candidate_positions(positions, context):
    """The positions of the context to choose modifications from; 'all' means CpG, CHG and CHH."""
    wanted = ("CHG", "CHH", "CpG") if context == "all" else (context,)
    return {
        position
        for position, information in positions.items()
        if information[1] in wanted
    }


def random_positions(candidates, modification_level, seed):
    """A random sample of modification_level of the candidate positions."""
    required = int(safe_rounder(len(candidates) * (modification_level or 0), 1, False))
    if seed is not None:
        random.seed(seed)
    return set(random.sample(sorted(candidates), required))


def strand_change(read, reverse_modification):
    """(reference base whose positions may change, base written there) for a read of the given strand."""
    reference = "C" if read.flag in TOP_FLAGS else "G"
    return reference, (
        reference if reverse_modification else {"C": "T", "G": "A"}[reference]
    )


def reference_positions(read, sequences, reference):
    """Genomic positions of the reference base under a read, ignoring its CIGAR."""
    window = sequences[read.reference_name][
        read.reference_start : read.reference_start + read.query_alignment_length
    ].upper()
    return {
        (
            read.reference_name,
            match.start() + read.reference_start,
            match.start() + read.reference_start + 1,
        )
        for match in re.finditer(reference, window)
    }


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


def read_offsets(read, method, chosen, positions, reverse_modification):
    """Read offsets of the bases to change: the chosen positions, or in CtoT mode the others."""
    if has_indels_or_clips(read):
        return sorted(
            position_correction_cigar(
                read, method, chosen, positions, reverse_modification
            )
        )
    selected = chosen & positions
    if method == "CtoT" and not reverse_modification:
        selected = positions - selected
    return sorted(position[1] - read.reference_start for position in selected)


def convertible(base, flag):
    return (base in ("C", "T") and flag in TOP_FLAGS) or (
        base in ("G", "A") and flag in BOTTOM_FLAGS
    )


def modified_sequence(read, offsets, new_base, chosen_here, targets):
    """The read sequence with the bases at the offsets changed to new_base.

    With a positions file (targets), each position is changed in at most its target number of reads.
    """
    bases = list(read.query_sequence)
    if targets is None:
        for offset in offsets:
            if len(bases) > offset and convertible(bases[offset], read.flag):
                bases[offset] = new_base
    else:
        for index, position in enumerate(sorted(chosen_here)):
            target = targets[position]
            if target[0] <= target[1] and len(offsets) > index:
                # TODO known bug: the base checked is at index, not at offsets[index].
                if len(bases) > offsets[index] and convertible(bases[index], read.flag):
                    bases[offsets[index]] = new_base
                    target[0] += 1
    return "".join(bases)


def read_information_line(read):
    return "{}\t{}\t{}\t{}\n".format(
        read.query_name + ("/1" if read.is_read1 else "/2"),
        read.reference_name,
        read.reference_start,
        read.reference_start + read.query_length,
    )


def simulate_reads(
    reads,
    sequences,
    method,
    chosen,
    targets,
    reverse_modification,
    outbam,
    read_information,
):
    """Writes the reads with simulated modifications and returns the chosen positions they cover.

    Only properly paired reads are written."""
    covered = set()
    header_pending = True
    for read in reads:
        if header_pending:
            # TODO check: the header line repeats until the first properly paired, aligned read.
            read_information.write("#Read ID\treference\tstart\tend\n")
        read_information.write(read_information_line(read))
        if read.flag not in TOP_FLAGS + BOTTOM_FLAGS:
            continue
        if read.reference_length == 0:
            outbam.write(read)
            continue
        header_pending = False
        reference, new_base = strand_change(read, reverse_modification)
        positions = reference_positions(read, sequences, reference)
        chosen_here = chosen & positions
        covered |= chosen_here
        offsets = read_offsets(read, method, chosen, positions, reverse_modification)
        if offsets:
            qualities = read.query_qualities
            read.query_sequence = modified_sequence(
                read, offsets, new_base, chosen_here, targets
            )
            read.query_qualities = qualities
        outbam.write(read)
    return covered


def modified_positions_summary(
    covered, positions, chosen_from_list, modification_level, context
):
    """The lines of the modified positions summary of one chromosome."""
    covered = sorted(covered)
    if chosen_from_list:
        level = "Custom"
    else:
        candidates = (
            len(positions)
            if context == "all"
            else sum(
                1 for information in positions.values() if information[1] == context
            )
        )
        level = safe_rounder(len(covered) / candidates, 3, True)
    rule = "__________________________________________________________________________________________________\n"
    yield rule
    yield "Absolute modified positions: {}   |   Percentage to all positions of the desired context: {} %\n".format(
        len(covered), level
    )
    yield rule
    yield from ("{}\t{}\t{}\n".format(*position) for position in covered)


def output_names(
    directory,
    name,
    method,
    label,
    context,
    per_chromosome,
    reverse_modification,
    extension,
):
    """(simulated alignments, read information, modified positions summary) file names."""
    base = path.join(directory, "{}_{}_{}_{}".format(name, method, label, context))
    chromosome_suffix = "" if per_chromosome is None else "_" + per_chromosome
    alignments = (
        base
        + ("_reversed" if reverse_modification else "")
        + chromosome_suffix
        + extension
    )
    if per_chromosome is None:
        read_information = base + "_read_information.txt.gz"
    else:
        read_information = path.join(
            directory,
            "{}_{}_{}_{}_{}_read_information.txt.gz".format(
                name, method, label, per_chromosome, context
            ),
        )
    return (
        alignments,
        read_information,
        base + chromosome_suffix + "_modified_positions_information.txt.gz",
    )


def modification_simulator(
    reference,
    control_file,
    model,
    read_length,
    input_file,
    simulation_input,
    method,
    modification_level,
    library,
    modified_positions,
    context,
    user_defined_context,
    coverage,
    region,
    overwrite,
    per_chromosome,
    GC_bias,
    sequence_bias,
    N_threads,
    reverse_modification,
    directory,
    seed,
    skip_clip_overlap,
    single_end,
    minimum_base_quality,
    minimum_mapping_quality,
    adjust_acapq_threshold,
    add_indels,
    redo_baq,
    compute_baq,
    ignore_orphans,
    max_depth,
):
    """Simulates TAPS or bisulfite conversion of randomly chosen (or listed) cytosines in the reads of a BAM or CRAM file."""
    logs.info("asTair's cytosine modification simulator started running.")
    if library != "directional":
        raise click.BadParameter(
            "Only directional libraries can be simulated.", param_hint="--library"
        )
    directory = output_directory(directory)
    region = None if None in region else region
    name = input_name(input_file)
    with open_alignments(input_file) as alignments:
        extension = ".cram" if alignments.is_cram else ".bam"
        write_mode = "wc" if alignments.is_cram else "wb"
        chromosome_lengths = dict(zip(alignments.references, alignments.lengths))
    modification_level = modification_level_transformation(
        modification_level, modified_positions
    )
    label = level_label(modification_level)
    # TODO check: --overwrite has no effect, existing outputs are always overwritten.
    bam_name, information_name, summary_name = output_names(
        directory,
        name,
        method,
        label,
        context,
        per_chromosome,
        reverse_modification,
        extension,
    )
    sequences = read_reference(reference)
    keys = context_keys(context, user_defined_context)
    chromosomes = [region[0]] if region else list(sequences)
    with (
        open_alignments(input_file, N_threads) as template,
        pysam.AlignmentFile(
            bam_name, write_mode, reference_filename=reference, template=template
        ) as outbam,
        gzip.open(information_name, "wt", compresslevel=9) as read_information,
        gzip.open(summary_name, "wt", compresslevel=9) as summary,
    ):
        for chromosome in chromosomes:
            positions, _ = find_cytosine_contexts(
                sequences[chromosome], chromosome, keys, user_defined_context, region
            )
            if modified_positions:
                targets = listed_positions(modified_positions, chromosome, input_file)
                positions, chosen = targets, set(targets)
            else:
                targets = None
                chosen = random_positions(
                    candidate_positions(positions, context), modification_level, seed
                )
            fetch = (
                region if region else (chromosome, 0, chromosome_lengths[chromosome])
            )
            covered = simulate_reads(
                iterate_reads(input_file, N_threads, fetch),
                sequences,
                method,
                chosen,
                targets,
                reverse_modification,
                outbam,
                read_information,
            )
            summary.writelines(
                modified_positions_summary(
                    covered,
                    positions,
                    modified_positions is not None,
                    modification_level,
                    context,
                )
            )
    pysam.index(bam_name)
    logs.info("asTair's cytosine modification simulator finished running.")
