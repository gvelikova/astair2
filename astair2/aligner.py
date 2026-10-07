import gzip
import logging
import os
import re
import shutil
import subprocess
import sys

import click
import pysam

from astair2.output import exit_if_exists, output_directory
from astair2.simple_fasta_parser import (
    _is_plain_gzip,
    reference_names,
    underscored_reference_path,
    write_reference_with_underscores,
)

logs = logging.getLogger(__name__)


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
    "reference",
    "--reference",
    "-f",
    required=True,
    help="Reference DNA sequence in FASTA format used for aligning of the sequencing reads.",
)
@click.option(
    "bwa_path",
    "--bwa_path",
    "-bp",
    required=False,
    help="The path to BWA for TAPS-like data and to bwameth.py for bisulfite sequencing.",
)
@click.option(
    "samtools_path",
    "--samtools_path",
    "-sp",
    required=False,
    help="The path to Samtools.",
)
@click.option(
    "directory",
    "--directory",
    "-d",
    required=True,
    help="Output directory to save files.",
)
@click.option(
    "method",
    "--method",
    "-m",
    required=False,
    default="mCtoT",
    type=click.Choice(["CtoT", "mCtoT"]),
    help="Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT)",
)
@click.option(
    "output_format",
    "--output_format",
    "-O",
    required=False,
    default="BAM",
    type=click.Choice(["CRAM", "BAM"]),
    help="Specify output format, possible options are BAM and CRAM. The Default BAM.",
)
@click.option(
    "minimum_mapping_quality",
    "--minimum_mapping_quality",
    "-mq",
    required=False,
    type=int,
    default=1,
    help="Set the minimum mapping quality for a read to be output to file (Default >=1).",
)
@click.option(
    "keep_unmapped",
    "--keep_unmapped",
    "-u",
    default=False,
    is_flag=True,
    help="Outputs the unmapped reads (Default false).",
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
    "N_threads",
    "--N_threads",
    "-t",
    default=1,
    required=True,
    help="The number of threads to spawn (Default 1).",
)
@click.option(
    "minimum_seed_length",
    "--minimum_seed_length",
    "-k",
    default=19,
    type=int,
    required=False,
    help="The minimum seed length used for alignment, see BWA manual. (Default 19).",
)
@click.option(
    "band_width",
    "--band_width",
    "-w",
    default=100,
    type=int,
    required=False,
    help="The band width for banded alignment, see BWA manual. (Default 100).",
)
@click.option(
    "dropoff",
    "--dropoff",
    "-D",
    default=100,
    type=int,
    required=False,
    help="The off-diagonal X-dropoff, see BWA manual. (Default 100).",
)
@click.option(
    "internal_seeds",
    "--internal_seeds",
    "-r",
    default=1.5,
    type=float,
    required=False,
    help="Looks for internal seeds inside a seed longer than minimum_seed_length * internal_seeds, see BWA manual. (Default 1.5).",
)
@click.option(
    "reseeding_occurence",
    "--reseeding_occurence",
    "-y",
    default=20,
    type=int,
    required=False,
    help="The seed occurrence for the 3rd round seeding, see BWA manual. (Default 20).",
)
@click.option(
    "N_skip_seeds",
    "--N_skip_seeds",
    "-c",
    default=500,
    type=int,
    required=False,
    help="Skips seeds with more than the given seed occurrences, see BWA manual. (Default 500).",
)
@click.option(
    "drop_chains",
    "--drop_chains",
    "-dc",
    default=0.5,
    type=float,
    required=False,
    help="Drops chains shorter than the specified fraction of the longest overlapping chain, see BWA manual. (Default 0.5).",
)
@click.option(
    "discard_chains",
    "--discard_chains",
    "-W",
    default=0,
    type=int,
    required=False,
    help="Discards a chain if seeded bases shorter than the specified value, see BWA manual. (Default 0).",
)
@click.option(
    "N_mate_rescues",
    "--N_mate_rescues",
    "-mr",
    default=50,
    type=int,
    required=False,
    help="Performs at most the specified rounds of mate rescues for each read, see BWA manual. (Default 50).",
)
@click.option(
    "skip_mate_rescue",
    "--skip_mate_rescue",
    "-s",
    is_flag=True,
    required=False,
    help="NB: Does not recommend unless necessary: skips mate rescue in mCtoT mode, see BWA manual. If set, orphan reads (paired reads that are not in a proper pair) will be generated. Ensure ignore_orphans in the caller is set to False.",
)
@click.option(
    "skip_pairing",
    "--skip_pairing",
    "-P",
    is_flag=True,
    required=False,
    help="NB: Does not recommend unless necessary: skips read pairing in mCtoT mode, but does rescue mates unless mate_skipping is also performed, see BWA manual. If set, orphan reads (paired reads that are not in a proper pair) will be generated. Ensure ignore_orphans in the caller is set to False.",
)
@click.option(
    "match_score",
    "--match_score",
    "-A",
    default=1,
    type=int,
    required=False,
    help="The score for a sequence match, which scales the remaing scoring options, see BWA manual. (Default 1).",
)
@click.option(
    "mismatch_penalty",
    "--mismatch_penalty",
    "-B",
    default=4,
    type=int,
    required=False,
    help="The penalty for a mismatch, see BWA manual. (Default 4).",
)
@click.option(
    "gap_open_penalty",
    "--gap_open_penalty",
    "-o",
    default="6,6",
    type=str,
    required=False,
    help="The gap open penalties for deletions and insertions, see BWA manual. (Default 6,6).",
)
@click.option(
    "gap_extension_penalty",
    "--gap_extension_penalty",
    "-E",
    default="1,1",
    type=str,
    required=False,
    help="The gap extension penalty with a cost size calculated as {-O} + {-E}*k, see BWA manual. (Default 1,1).",
)
@click.option(
    "end_clipping_penalty",
    "--end_clipping_penalty",
    "-L",
    default="5,5",
    type=str,
    required=False,
    help="The penalty for 5-prime- and 3-prime-end clipping, see BWA manual. (Default 5,5).",
)
@click.option(
    "unpaired_penalty",
    "--unpaired_penalty",
    "-U",
    default=17,
    type=int,
    required=False,
    help="The penalty for an unpaired read pair, see BWA manual. (Default 17).",
)
@click.option(
    "read_type",
    "--read_type",
    "-x",
    default="null",
    type=click.Choice(["null", "pacbio", "ont2d", "intractg"]),
    required=False,
    help="Changes multiple parameters unless overridden, see BWA manual. (Default null).",
)
@click.option(
    "smart_pairing",
    "--smart_pairing",
    "-smp",
    default=False,
    is_flag=True,
    required=False,
    help="Ignores read2, see BWA manual.",
)
@click.option(
    "read_group",
    "--read_group",
    "-rg",
    default=None,
    type=str,
    required=False,
    help='Adds the given read group line "@RG\\tID:ids\\tSM:name", see BWA manual. (Default None).',
)
@click.option(
    "header_string",
    "--header_string",
    "-hs",
    default=None,
    type=str,
    required=False,
    help="Adds the given string to header if it starts with @, or to file if it does not, see BWA manual.",
)
@click.option(
    "include_alt",
    "--include_alt",
    "-al",
    default=False,
    is_flag=True,
    required=False,
    help="Treats ALT contigs as part of the primary assembly, see BWA manual.",
)
@click.option(
    "split_alignment",
    "--split_alignment",
    "-spa",
    default=False,
    is_flag=True,
    required=False,
    help="Takes the alignment with the smallest coordinate as primary in case of split alignment, see BWA manual.",
)
@click.option(
    "supplementary_mapq",
    "--supplementary_mapq",
    "-smq",
    default=False,
    is_flag=True,
    required=False,
    help="Does not modify the MAPQ of supplementary alignments, see BWA manual.",
)
@click.option(
    "minimum_score",
    "--minimum_score",
    "-T",
    default=30,
    type=int,
    required=False,
    help="Minimum score to ouput, see BWA manual. (Default 30).",
)
@click.option(
    "alternative_score",
    "--alternative_score",
    "-h",
    default="5,200",
    type=str,
    required=False,
    help="Defines the number of reads that will be tagged as XA if they have 80 percent of the maximum score, see BWA manual. (Default 5,200).",
)
@click.option(
    "all_alignments",
    "--all_alignments",
    "-aa",
    default=False,
    is_flag=True,
    required=False,
    help="Outputs all SE or unpaired PE reads, see BWA manual.",
)
@click.option(
    "fasta_comment",
    "--fasta_comment",
    "-fc",
    default=None,
    required=False,
    help="Adds FASTA/FASTQ comment to output, see BWA manual.",
)
@click.option(
    "fasta_header",
    "--fasta_header",
    "-fh",
    default=None,
    required=False,
    help="Outputs the reference header in the XR tag, see BWA manual.",
)
@click.option(
    "clip_supplementary",
    "--clip_supplementary",
    "-cs",
    default=False,
    is_flag=True,
    required=False,
    help="Enables soft clipping on supplementary alignments, see BWA manual.",
)
@click.option(
    "mark_splitted",
    "--mark_splitted",
    "-ms",
    default=False,
    is_flag=True,
    required=False,
    help="Labels short split reads as secondary, see BWA manual.",
)
@click.option(
    "reads_distribution",
    "--reads_distribution",
    "-rd",
    default=None,
    required=False,
    help="Specifies the mean, standard deviation (10 percent of the mean if absent), max (4 sigma from the mean if absent) and min of the insert size distribution for FR orientation, see BWA manual. Must be provided as FLOAT,FLOAT,INT,INT (Default reads distribution metrics are inferred from the data).",
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
    "use_underscores",
    "--use_underscores",
    "-uu",
    default=False,
    is_flag=True,
    required=False,
    help="Uses as a reference the fasta file with added underscores in the sequence names that is afterwards used for calling. (Default False).",
)
@click.option(
    "temp_dir",
    "--temp_dir",
    "-td",
    required=False,
    help="Provides a custom directory to write temporary files. (Default the chosen directory for the output).",
)
@click.option(
    "compress",
    "--compress",
    "-z",
    is_flag=True,
    default=False,
    required=False,
    help="Should the reference FASTA be compressed after the run (Default False).",
)
@click.option(
    "sort_chunck_size",
    "--sort_chunck_size",
    "-cz",
    default="768M",
    required=False,
    help="WARNING: sorting large files >= 50 GB might require to increaase the maximum memory per thread parameter; recognised sufixes are  K/M/G (Default 768M).",
)
def align(
    fq1,
    fq2,
    reference,
    bwa_path,
    samtools_path,
    directory,
    method,
    output_format,
    minimum_mapping_quality,
    keep_unmapped,
    N_threads,
    minimum_seed_length,
    band_width,
    dropoff,
    internal_seeds,
    reseeding_occurence,
    N_skip_seeds,
    drop_chains,
    discard_chains,
    N_mate_rescues,
    skip_mate_rescue,
    skip_pairing,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    end_clipping_penalty,
    unpaired_penalty,
    read_type,
    single_end,
    smart_pairing,
    read_group,
    header_string,
    include_alt,
    split_alignment,
    supplementary_mapq,
    minimum_score,
    alternative_score,
    all_alignments,
    fasta_comment,
    fasta_header,
    clip_supplementary,
    mark_splitted,
    reads_distribution,
    add_underscores,
    temp_dir,
    use_underscores,
    compress,
    sort_chunck_size,
):
    """Align raw reads in fastq format to a reference genome. bwa is required to align TAPS reads, and bwa-meth fif you plan to process BS-seq data."""
    run_alignment(
        fq1,
        fq2,
        reference,
        bwa_path,
        samtools_path,
        directory,
        method,
        output_format,
        minimum_mapping_quality,
        keep_unmapped,
        N_threads,
        minimum_seed_length,
        band_width,
        dropoff,
        internal_seeds,
        reseeding_occurence,
        N_skip_seeds,
        drop_chains,
        discard_chains,
        N_mate_rescues,
        skip_mate_rescue,
        skip_pairing,
        match_score,
        mismatch_penalty,
        gap_open_penalty,
        gap_extension_penalty,
        end_clipping_penalty,
        unpaired_penalty,
        read_type,
        single_end,
        smart_pairing,
        read_group,
        header_string,
        include_alt,
        split_alignment,
        supplementary_mapq,
        minimum_score,
        alternative_score,
        all_alignments,
        fasta_comment,
        fasta_header,
        clip_supplementary,
        mark_splitted,
        reads_distribution,
        add_underscores,
        temp_dir,
        use_underscores,
        compress,
        sort_chunck_size,
    )


def find_executable(path, default_name, option):
    """The executable given by path, or found on PATH by name; stops with an error if there is none."""
    executable = shutil.which(path or default_name)
    if executable is None:
        logs.error(
            "{} was not found. Please install it or give its location with {}.".format(
                path or default_name, option
            )
        )
        sys.exit(1)
    return executable


def which_path(bwa_path, samtools_path, method):
    """The aligner (bwa for mCtoT, bwameth.py for CtoT) and samtools executables."""
    aligner = find_executable(
        bwa_path, "bwa" if method == "mCtoT" else "bwameth.py", "--bwa_path"
    )
    return aligner, find_executable(samtools_path, "samtools", "--samtools_path")


def _flag(enabled, option):
    return [option] if enabled else []


def _value(value, option):
    return [] if value is None else [option, str(value)]


def bwa_mem_command(aligner, reference, fastq_files, settings):
    """The bwa mem command line for TAPS-like (mCtoT) data."""
    s = settings
    return (
        [
            aligner,
            "mem",
            "-t",
            str(s["N_threads"]),
            "-k",
            str(s["minimum_seed_length"]),
            "-w",
            str(s["band_width"]),
            "-d",
            str(s["dropoff"]),
            "-r",
            str(s["internal_seeds"]),
            "-y",
            str(s["reseeding_occurence"]),
            "-c",
            str(s["N_skip_seeds"]),
            "-D",
            str(s["drop_chains"]),
            "-W",
            str(s["discard_chains"]),
            "-m",
            str(s["N_mate_rescues"]),
        ]
        + _flag(s["skip_mate_rescue"], "-S")
        + _flag(s["skip_pairing"], "-P")
        + [
            "-A",
            str(s["match_score"]),
            "-B",
            str(s["mismatch_penalty"]),
            "-O",
            s["gap_open_penalty"],
            "-E",
            s["gap_extension_penalty"],
            "-L",
            s["end_clipping_penalty"],
            "-U",
            str(s["unpaired_penalty"]),
        ]
        + ([] if s["read_type"] == "null" else ["-x", s["read_type"]])
        + _flag(s["smart_pairing"], "-p")
        + _value(s["read_group"], "-R")
        + _value(s["header_string"], "-H")
        + _flag(s["include_alt"], "-j")
        + _flag(s["split_alignment"], "-5")
        + _flag(s["supplementary_mapq"], "-q")
        + ["-T", str(s["minimum_score"]), "-h", s["alternative_score"]]
        + _flag(s["all_alignments"], "-a")
        # TODO check: bwa's -C and -V are switches, the values given to them become extra arguments.
        + _value(s["fasta_comment"], "-C")
        + _value(s["fasta_header"], "-V")
        + _flag(s["clip_supplementary"], "-Y")
        + _flag(s["mark_splitted"], "-M")
        + _value(s["reads_distribution"], "-I")
        + [reference]
        + fastq_files
    )


def alignment_commands(
    aligner,
    samtools,
    reference,
    fastq_files,
    method,
    output_format,
    settings,
    temp_name,
):
    """The aligner, samtools view and samtools sort commands of the alignment pipeline."""
    if method == "mCtoT":
        align = bwa_mem_command(aligner, reference, fastq_files, settings)
        sort_options = ["-m", settings["sort_chunck_size"]]
    else:
        align = [
            aligner,
            "-t",
            str(settings["N_threads"]),
            "--reference",
            reference,
        ] + fastq_files
        sort_options = ["-T", temp_name]
    keep_unmapped = settings["keep_unmapped"]
    view = (
        [
            samtools,
            "view",
            "-hC" if output_format == "CRAM" else "-hb",
            "-T",
            reference,
            "-q",
            str(0 if keep_unmapped else settings["minimum_mapping_quality"]),
        ]
        + ([] if keep_unmapped else ["-F", "4"])
        + ["-O", output_format]
    )
    sort = (
        [samtools, "sort"]
        + sort_options
        + ["-@", str(settings["N_threads"]), "-O", output_format]
    )
    return align, view, sort


def run_pipeline(commands, output_file):
    """Runs commands connected by pipes, writing the last one's output to a file; True if all succeed."""
    processes = []
    with open(output_file, "wb") as output:
        for index, command in enumerate(commands):
            last = index == len(commands) - 1
            processes.append(
                subprocess.Popen(
                    command,
                    stdin=processes[-1].stdout if processes else None,
                    stdout=output if last else subprocess.PIPE,
                )
            )
            if len(processes) > 1:
                processes[-2].stdout.close()
        exit_codes = [process.wait() for process in processes]
    for command, exit_code in zip(commands, exit_codes):
        if exit_code != 0:
            logs.error("{} failed with exit code {}.".format(command[0], exit_code))
    return all(exit_code == 0 for exit_code in exit_codes)


def bgzip_in_place(file_name):
    """Replaces a file by its BGZIP-compressed version, like bgzip does."""
    pysam.tabix_compress(file_name, file_name + ".gz", force=True)
    os.remove(file_name)


def gunzip_in_place(file_name):
    """Replaces a gzip-compressed file by its uncompressed version, like gunzip does."""
    target = file_name[: -len(".gz")]
    with gzip.open(file_name, "rb") as source, open(target, "wb") as output:
        shutil.copyfileobj(source, output)
    os.remove(file_name)
    return target


def check_index(
    aligner, reference, method, output_format, add_underscores, use_underscores
):
    """The reference to align to, after creating the requested copy with underscores in the sequence
    names and building the aligner index if it is missing."""
    reference_names(reference)
    if add_underscores:
        write_reference_with_underscores(reference)
    if use_underscores and os.path.isfile(underscored_reference_path(reference)):
        reference = underscored_reference_path(reference)
    index_suffix = ".bwt" if method == "mCtoT" else ".bwameth.c2t"
    if not os.path.isfile(reference + index_suffix):
        if (
            output_format == "CRAM"
            and reference.endswith(".gz")
            and _is_plain_gzip(reference)
        ):
            # samtools cannot use a plain gzip reference for CRAM
            reference = gunzip_in_place(reference)
        if subprocess.run([aligner, "index", reference]).returncode != 0:
            logs.error("Indexing the reference {} failed.".format(reference))
            sys.exit(1)
    return reference


def output_name(fq1, paired):
    name = os.path.splitext(os.path.basename(fq1))[0]
    if paired:
        name = re.sub("(_R1|_1)", "", name)
    if fq1.endswith(".gz"):
        name = "_".join(name.split(".")[:-1])
    return name


def run_alignment(
    fq1,
    fq2,
    reference,
    bwa_path,
    samtools_path,
    directory,
    method,
    output_format,
    minimum_mapping_quality,
    keep_unmapped,
    N_threads,
    minimum_seed_length,
    band_width,
    dropoff,
    internal_seeds,
    reseeding_occurence,
    N_skip_seeds,
    drop_chains,
    discard_chains,
    N_mate_rescues,
    skip_mate_rescue,
    skip_pairing,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    end_clipping_penalty,
    unpaired_penalty,
    read_type,
    single_end,
    smart_pairing,
    read_group,
    header_string,
    include_alt,
    split_alignment,
    supplementary_mapq,
    minimum_score,
    alternative_score,
    all_alignments,
    fasta_comment,
    fasta_header,
    clip_supplementary,
    mark_splitted,
    reads_distribution,
    add_underscores,
    temp_dir,
    use_underscores,
    compress,
    sort_chunck_size,
):
    """Aligns the reads to the reference with bwa (mCtoT) or bwa-meth (CtoT), and writes a sorted and indexed BAM or CRAM file."""
    settings = {key: value for key, value in locals().items()}
    logs.info("asTair genome aligner started running.")
    directory = output_directory(directory)
    paired = not single_end or fq2 is not None
    fastq_files = [fq1, fq2] if paired and fq2 else [fq1]
    output_file = os.path.join(
        directory, output_name(fq1, paired) + "_" + method + "." + output_format.lower()
    )
    exit_if_exists(
        [output_file],
        "The output files will not be overwritten. Please rename the input or the existing output files before rerunning if the input is different.",
    )
    initially_compressed = reference.endswith(".gz")
    aligner, samtools = which_path(bwa_path, samtools_path, method)
    reference = check_index(
        aligner,
        reference,
        method,
        output_format,
        add_underscores or use_underscores,
        use_underscores,
    )
    temp_name = (
        os.path.join(temp_dir, "temp") if temp_dir else os.path.join(directory, "temp")
    )
    try:
        succeeded = (
            run_pipeline(
                alignment_commands(
                    aligner,
                    samtools,
                    reference,
                    fastq_files,
                    method,
                    output_format,
                    settings,
                    temp_name,
                ),
                output_file,
            )
            and subprocess.run([samtools, "index", output_file]).returncode == 0
        )
    finally:
        if not reference.endswith(".gz") and (compress or initially_compressed):
            bgzip_in_place(reference)
    if not succeeded:
        # Do not leave an empty or truncated alignment that looks like a result.
        for partial in (output_file, output_file + ".bai", output_file + ".crai"):
            if os.path.isfile(partial):
                os.remove(partial)
        logs.error("asTair genome aligner failed; no output was written.")
        sys.exit(1)
    logs.info("asTair genome aligner finished running.")
