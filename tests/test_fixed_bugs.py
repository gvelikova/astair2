"""Regression tests for bugs inherited from asTair that have been fixed in astair2."""

import gzip
import inspect
from pathlib import Path
from unittest import mock

from click.testing import CliRunner

import astair2.aligner as aligner
from astair2.caller import strand_flags
from astair2.cli import cli

DATA = Path(__file__).resolve().parent / "test_data"
LAMBDA = str(DATA / "lambda_phage.fa")
TAPS = str(DATA / "small_real_taps_lambda_mCtoT.bam")


def run(*args):
    return CliRunner().invoke(cli, [str(a) for a in args])


def read_text(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return handle.read()


def data_rows(path):
    return [line for line in read_text(path).splitlines() if not line.startswith("#")]


def call_cpg(directory, *options):
    directory.mkdir()
    result = run(
        "call", "-i", TAPS, "-f", LAMBDA, "-co", "CpG", "-d", directory, *options
    )
    assert result.exit_code == 0, result.output
    return directory / "small_real_taps_lambda_mCtoT_mCtoT_CpG.mods"


def test_orphan_reads_can_be_included(tmp_path):
    """Including orphans used to drop every position. The test BAM has no orphans, so the output
    must match the default."""
    default = call_cpg(tmp_path / "default")
    with_orphans = call_cpg(tmp_path / "orphans", "-io", "False")
    assert data_rows(with_orphans)
    assert read_text(with_orphans) == read_text(default)


def test_orphan_flags_count_for_their_strand():
    assert strand_flags("C", False, "directional", ignore_orphans=False) == (
        (99, 147, 97, 145),
        (83, 163, 81, 161),
    )
    assert strand_flags("G", False, "directional", ignore_orphans=False) == (
        (83, 163, 81, 161),
        (99, 147, 97, 145),
    )
    assert strand_flags("C", False, "directional") == ((99, 147), (83, 163))


def test_zero_coverage_adds_uncovered_positions_in_order(tmp_path):
    """--zero_coverage used to crash."""
    default = call_cpg(tmp_path / "default")
    with_uncovered = call_cpg(tmp_path / "uncovered", "-zc")
    rows = data_rows(with_uncovered)
    assert run("find", "-f", LAMBDA, "-co", "CpG", "-d", tmp_path).exit_code == 0
    assert len(rows) == len(data_rows(tmp_path / "lambda_phage_CpG.bed"))
    starts = [int(row.split("\t")[1]) for row in rows]
    assert starts == sorted(starts)
    covered = [row for row in rows if row.split("\t")[11] != "0"]
    assert covered == data_rows(default)
    uncovered = [row.split("\t") for row in rows if row.split("\t")[11] == "0"]
    assert uncovered[0][3:8] == [
        "*",
        "0",
        "0",
        uncovered[0][6],
        {"C": "T", "G": "A"}[uncovered[0][6]],
    ]
    # uncovered positions do not change the statistics
    assert read_text(with_uncovered.with_suffix(".stats")) == read_text(
        default.with_suffix(".stats")
    )


def test_overlap_clipping_can_be_disabled(tmp_path):
    """--skip_clip_overlap used to be a flag that defaulted to True."""
    skipped = call_cpg(tmp_path / "skip", "-sc")
    kept = call_cpg(tmp_path / "keep", "-kc")
    assert read_text(skipped) == read_text(call_cpg(tmp_path / "default"))

    def informative_reads(path):
        return sum(
            int(row.split("\t")[4]) + int(row.split("\t")[5]) for row in data_rows(path)
        )

    assert informative_reads(kept) > informative_reads(skipped)


def test_find_compressed_output_matches_uncompressed(tmp_path):
    """The compressed output used to lose the CONTEXT column."""
    plain, compressed = tmp_path / "plain", tmp_path / "gz"
    plain.mkdir(), compressed.mkdir()
    assert run("find", "-f", LAMBDA, "-co", "CpG", "-d", plain).exit_code == 0
    assert (
        run("find", "-f", LAMBDA, "-co", "CpG", "--gz", "-d", compressed).exit_code == 0
    )
    expected = read_text(plain / "lambda_phage_CpG.bed")
    actual = read_text(compressed / "lambda_phage_CpG.bed.gz")
    # compare a bool: pytest's diff of two multi-megabyte strings takes minutes
    assert actual.split("\n", 1)[0] == expected.split("\n", 1)[0]
    assert actual == expected, "compressed and uncompressed outputs differ"


def test_find_per_chromosome_only_reports_that_chromosome(tmp_path):
    """--per_chromosome used to only change the output file name."""
    reference = tmp_path / "two.fa"
    reference.write_text(">chrA\nACGTTCAGCCGATCCA\n>chrB\nTTCGACCAGGCTCCGA\n")
    assert run("find", "-f", reference, "-chr", "chrB", "-d", tmp_path).exit_code == 0
    chromosomes = {
        row.split("\t")[0] for row in data_rows(tmp_path / "two_chrB_all.bed")
    }
    assert chromosomes == {"chrB"}


def test_simulate_per_chromosome_ignores_other_sequences(tmp_path):
    """--per_chromosome used to simulate every reference sequence, failing on those absent from the BAM."""
    reference = tmp_path / "lambda_plus.fa"
    reference.write_text(Path(LAMBDA).read_text() + ">chrExtra\nACGTTCAGCCGATCCA\n")
    result = run(
        "simulate",
        "-f",
        reference,
        "-l",
        "80",
        "-i",
        TAPS,
        "-ml",
        "50",
        "-co",
        "CpG",
        "-chr",
        "lambda",
        "-s",
        "1",
        "-d",
        tmp_path,
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "small_real_taps_lambda_mCtoT_mCtoT_50_CpG_lambda.bam").exists()


def test_align_short_options_are_unique():
    """-ms used to be declared for both --minimum_score and --mark_splitted."""
    options = [opt for param in cli.commands["align"].params for opt in param.opts]
    assert len(options) == len(set(options))


def test_align_minimum_score_and_mark_splitted_options():
    names = list(inspect.signature(aligner.run_alignment).parameters)
    with mock.patch.object(aligner, "run_alignment") as run_alignment:
        result = CliRunner().invoke(
            cli,
            [
                "align",
                "-1",
                "reads.fq.gz",
                "-f",
                "ref.fa",
                "-d",
                ".",
                "-ms",
                "-T",
                "20",
            ],
        )
    assert result.exit_code == 0, result.output
    arguments = dict(zip(names, run_alignment.call_args.args))
    assert arguments["mark_splitted"] is True
    assert arguments["minimum_score"] == 20


def test_simulate_region(tmp_path):
    """The region branch used to call the context search with a missing argument."""
    result = run(
        "simulate",
        "-f",
        LAMBDA,
        "-l",
        "80",
        "-i",
        TAPS,
        "-ml",
        "60",
        "-co",
        "CpG",
        "-r",
        "lambda",
        "100",
        "20000",
        "-s",
        "11",
        "-d",
        tmp_path,
    )
    assert result.exit_code == 0


def test_phred_quality_lines_starting_with_at_sign(tmp_path):
    """Quality lines starting with '@' (Phred 31) used to lose their first score."""
    fastq = tmp_path / "reads_1.fq.gz"
    fastq.write_bytes(gzip.compress(b"@read1\nTACG\n+\n@###\n@read2\nTACG\n+\n@###\n"))
    assert run("phred", "-1", fastq, "--se", "-d", tmp_path).exit_code == 0
    means = (tmp_path / "reads_1_total_Phred.txt").read_text().splitlines()[2]
    assert "thymines: 31.0" in means
    assert "adenines: 2.0" in means
