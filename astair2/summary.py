import logging
from os import path

import click

from astair2.bam_file_parser import iterate_reads
from astair2.context_search import context_keys, find_cytosine_contexts
from astair2.output import (
    exit_if_exists,
    input_name,
    open_text,
    output_directory,
    tab_line,
)
from astair2.simple_fasta_parser import (
    read_reference,
    reference_names,
    write_reference_with_underscores,
)
from astair2.vcf_reader import read_vcf

logs = logging.getLogger(__name__)

SUMMARY_HEADER = (
    "#CHROM",
    "START",
    "END",
    "READ_NAME",
    "MODIFICATION_STATUS",
    "FLAG",
    "CONTEXT",
    "SPECIFIC_CONTEXT",
    "BQ",
    "MAPQ",
    "STRAND",
    "FRAGMENT_LENGTH",
    "AS",
    "XS",
    "EDIT",
    "KNOWN_SNP",
    "INFO",
)

# read base that indicates a modification, per method and strand
MODIFIED_BASES = {
    ("mCtoT", "top"): "Tt",
    ("CtoT", "top"): "Cc",
    ("mCtoT", "bottom"): "Aa",
    ("CtoT", "bottom"): "Gg",
}


@click.command()
@click.option(
    "input_file",
    "--input_file",
    "-i",
    required=True,
    help="BAM|CRAM format file containing sequencing reads.",
)
@click.option(
    "reference",
    "--reference",
    "-f",
    required=True,
    help="Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.",
)
@click.option(
    "known_snp",
    "--known_snp",
    "-ks",
    default=None,
    required=False,
    help="VCF format file containing genotyped WGS high quality variants or known common variants in VCF format (dbSNP, 1000 genomes, etc.).",
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
    type=click.Choice(["directional"]),
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
    "single_end",
    "--single_end",
    "-se",
    default=False,
    is_flag=True,
    required=False,
    help="Indicates single-end sequencing reads (Default False).",
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
    "minimum_base_quality",
    "--minimum_base_quality",
    "-bq",
    required=False,
    type=int,
    default=0,
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
    default="0",
    type=click.Choice([".", "0", "*", "NA"]),
    required=False,
    help="What symbol should be used for a value where no enough quantative information is used. (Default *).",
)
def summarise(
    input_file,
    reference,
    known_snp,
    context,
    user_defined_context,
    library,
    method,
    region,
    minimum_base_quality,
    minimum_mapping_quality,
    per_chromosome,
    N_threads,
    directory,
    single_end,
    add_underscores,
    no_information,
):
    """Collects and outputs modification information per read."""
    read_summariser(
        input_file,
        reference,
        known_snp,
        context,
        user_defined_context,
        library,
        method,
        region,
        minimum_base_quality,
        minimum_mapping_quality,
        per_chromosome,
        N_threads,
        directory,
        single_end,
        add_underscores,
        no_information,
    )


def read_positions(read, positions):
    """(genomic position, read index) of the read bases aligned to the given positions."""
    pairs = read.get_aligned_pairs(matches_only=True)
    return [
        ((read.reference_name, reference_index, reference_index + 1), read_index)
        for read_index, reference_index in pairs
        if (read.reference_name, reference_index, reference_index + 1) in positions
    ]


def modification_status(read, read_index, method, strand_flags):
    """1 if the read base shows a modification, else 0."""
    top_flags, bottom_flags = strand_flags
    strand = (
        "top"
        if read.flag in top_flags
        else "bottom" if read.flag in bottom_flags else None
    )
    return int(
        strand is not None
        and read.query_sequence[read_index] in MODIFIED_BASES[(method, strand)]
    )


def read_rows(
    read, positions, possible_mods, true_variants, known_snp, method, strand_flags
):
    """The read summary rows of one read."""
    snp_status = "*"
    for position, read_index in read_positions(read, positions) + (
        read_positions(read, possible_mods) if possible_mods else []
    ):
        possible = known_snp is not None and position in possible_mods
        strand_base = (
            positions[position] if position in positions else possible_mods[position]
        )[2]
        if known_snp is not None and position in true_variants:
            # TODO check: once set, the WGS status also applies to the following positions of the read.
            snp_status = "WGS"
        specific_context, context = (
            possible_mods[position]
            if position in possible_mods
            else positions[position]
        )[:2]
        yield (
            read.reference_name,
            position[1],
            position[2],
            read.query_name,
            modification_status(read, read_index, method, strand_flags),
            read.flag,
            context,
            specific_context,
            read.query_qualities[read_index],
            read.mapping_quality,
            "+" if strand_base == "C" else "-",
            abs(read.template_length),
            read.get_tag("AS"),
            read.get_tag("XS"),
            read.get_tag("NM"),
            snp_status,
            "possible_modification" if possible else "true_modification",
        )


def read_summariser(
    input_file,
    reference,
    known_snp,
    context,
    user_defined_context,
    library,
    method,
    region,
    minimum_base_quality,
    minimum_mapping_quality,
    per_chromosome,
    N_threads,
    directory,
    single_end,
    add_underscores,
    no_information,
):
    """Writes, for every read base aligned to a cytosine of the requested contexts, whether it shows a modification."""
    logs.info("asTair's read information summary function started running.")
    output_directory(directory)
    chromosome_part = "" if per_chromosome is None else per_chromosome + "_"
    # TODO check: with per_chromosome there is no underscore before read_summary.
    file_name = path.join(
        directory,
        input_name(input_file)
        + "_"
        + method
        + "_"
        + chromosome_part
        + context
        + ("_" if per_chromosome is None else "")
        + "read_summary.txt.gz",
    )
    exit_if_exists(
        [file_name],
        "Read information summary file with this name exists. Please rename before rerunning.",
    )
    if add_underscores:
        write_reference_with_underscores(reference)
    if region != (None, None, None):
        chromosomes = [region[0]]
    else:
        chromosomes = (
            reference_names(reference) if per_chromosome is None else [per_chromosome]
        )
    strand_flags = ((0,), (16,)) if single_end else ((99, 147), (83, 163))
    keys = context_keys(context, user_defined_context)
    true_variants, possible_mods = set(), {}
    with open_text(file_name, compress=True) as output:
        output.write(tab_line(SUMMARY_HEADER))
        for chromosome in chromosomes:
            logs.info(
                "Starting read information extraction on {} chromosome (sequence).".format(
                    chromosome
                )
            )
            sequence = read_reference(reference, [chromosome])[chromosome]
            positions, _ = find_cytosine_contexts(
                sequence, chromosome, keys, user_defined_context
            )
            start, end = (
                (region[1], region[2])
                if region != (None, None, None)
                else (0, len(sequence))
            )
            if known_snp is not None:
                true_variants, possible_mods = read_vcf(
                    known_snp, chromosome, sequence, N_threads, start, end
                )
            for read in iterate_reads(input_file, N_threads, (chromosome, start, end)):
                output.writelines(
                    tab_line(row)
                    for row in read_rows(
                        read,
                        positions,
                        possible_mods,
                        true_variants,
                        known_snp,
                        method,
                        strand_flags,
                    )
                )
    logs.info("asTair read information summary function finished running.")
