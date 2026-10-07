"""Bugs inherited from asTair, pinned as strict expected failures.

Each test describes the correct behaviour. While the bug exists the test fails
and is reported as xfail; once the bug is fixed the test passes, and because
the marker is strict, pytest reports it as XPASS(strict) so the marker gets
removed together with the fix.
"""

import csv
import gzip
from pathlib import Path

import pysam
import pytest
from click.testing import CliRunner

from astair2.cli import cli

DATA = Path(__file__).resolve().parent / "test_data"
LAMBDA = str(DATA / "lambda_phage.fa")
TAPS = str(DATA / "small_real_taps_lambda_mCtoT.bam")


def run(*args):
    return CliRunner().invoke(cli, [str(a) for a in args])


def data_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [line for line in handle if not line.startswith("#")]


@pytest.mark.xfail(
    strict=True,
    reason="with --ignore_orphans False every position is dropped (inherited TypeError in the flag expectations)",
)
def test_call_keeps_positions_when_orphans_are_included(tmp_path):
    result = run(
        "call", "-i", TAPS, "-f", LAMBDA, "-co", "CpG", "-io", "False", "-d", tmp_path
    )
    assert result.exit_code == 0
    assert data_rows(tmp_path / "small_real_taps_lambda_mCtoT_mCtoT_CpG.mods")


@pytest.mark.xfail(
    strict=True,
    reason="writing the uncovered positions always fails",
)
def test_call_zero_coverage_outputs_uncovered_positions(tmp_path):
    result = run("call", "-i", TAPS, "-f", LAMBDA, "-co", "CpG", "-zc", "-d", tmp_path)
    assert result.exit_code == 0
    rows = data_rows(tmp_path / "small_real_taps_lambda_mCtoT_mCtoT_CpG.mods")
    find = run("find", "-f", LAMBDA, "-co", "CpG", "-d", tmp_path)
    assert find.exit_code == 0
    assert len(rows) == len(data_rows(tmp_path / "lambda_phage_CpG.bed"))


@pytest.mark.xfail(
    strict=True,
    reason="--skip_clip_overlap is a flag defaulting to True, so it cannot be disabled",
)
def test_call_overlap_clipping_can_be_disabled():
    option = next(
        p for p in cli.commands["call"].params if p.name == "skip_clip_overlap"
    )
    assert not (option.is_flag and option.default is True and not option.secondary_opts)


@pytest.mark.xfail(
    strict=True, reason="the gzip writer formats six values into five placeholders"
)
def test_find_compressed_output_matches_uncompressed(tmp_path):
    plain, compressed = tmp_path / "plain", tmp_path / "gz"
    plain.mkdir(), compressed.mkdir()
    assert run("find", "-f", LAMBDA, "-co", "CpG", "-d", plain).exit_code == 0
    assert (
        run("find", "-f", LAMBDA, "-co", "CpG", "--gz", "-d", compressed).exit_code == 0
    )
    with open(plain / "lambda_phage_CpG.bed") as handle:
        expected = handle.read()
    with gzip.open(compressed / "lambda_phage_CpG.bed.gz", "rt") as handle:
        actual = handle.read()
    # Compare the headers first and then a plain bool: pytest's diff of two
    # multi-megabyte strings takes minutes.
    assert actual.split("\n", 1)[0] == expected.split("\n", 1)[0]
    assert actual == expected, "compressed and uncompressed outputs differ"


def test_simulate_region(tmp_path):
    """Fixed: the region branch used to call the context search with a missing argument."""
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


@pytest.mark.xfail(
    strict=True, reason="-ms is declared for both --minimum_score and --mark_splitted"
)
def test_align_short_options_are_unique():
    options = [opt for param in cli.commands["align"].params for opt in param.opts]
    assert len(options) == len(set(options))


@pytest.mark.xfail(
    strict=True,
    reason="positions where most reads differ from the reference and no opposite-strand "
    "read shows the reference base or A are dropped unless --no_information is 0",
)
def test_call_output_does_not_depend_on_no_information_symbol(tmp_path):
    reference = tmp_path / "ref.fa"
    sequence = "ACGTACGATTGCAAGCTTCGATCGGATCCGTAGCTAGCTAACGTTGCATGCAAGCTGACGTCAGTCAGCTAGCATGACGTAGCTAGCATCGATCGTACG"
    reference.write_text(">chrT\n" + sequence + "\n")
    pysam.faidx(str(reference))
    bam = tmp_path / "reads.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrT", "LN": len(sequence)}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        for number in range(5):
            # fully modified top-strand reads (every C read as T), no bottom-strand reads
            read = pysam.AlignedSegment()
            read.query_name, read.flag, read.reference_id, read.reference_start = (
                "r%d" % number,
                99,
                0,
                0,
            )
            read.query_sequence = sequence[:60].replace("C", "T")
            read.query_qualities = pysam.qualitystring_to_array("F" * 60)
            read.cigarstring, read.mapping_quality = "60M", 60
            read.next_reference_id, read.next_reference_start, read.template_length = (
                0,
                40,
                100,
            )
            output.write(read)
    pysam.index(str(bam))
    rows = {}
    for symbol in ("0", "*"):
        out = tmp_path / ("ni_zero" if symbol == "0" else "ni_star")
        out.mkdir()
        assert (
            run(
                "call",
                "-i",
                bam,
                "-f",
                reference,
                "-co",
                "CpG",
                "-ni",
                symbol,
                "-d",
                out,
            ).exit_code
            == 0
        )
        rows[symbol] = [
            row.split("\t")[:3] for row in data_rows(out / "reads_mCtoT_CpG.mods")
        ]
    assert rows["0"], "the fully modified CpGs should be called"
    assert rows["*"] == rows["0"]


@pytest.mark.xfail(
    strict=True, reason="the UNMOD_COUNT and MOD_COUNT columns are multiplied by 100"
)
def test_mbias_counts_are_read_counts(tmp_path):
    assert (
        run("mbias", "-f", LAMBDA, "-i", TAPS, "-l", "80", "-d", tmp_path).exit_code
        == 0
    )
    with pysam.AlignmentFile(TAPS) as bam:
        reads = sum(1 for _ in bam.fetch(until_eof=True))
    with open(tmp_path / "small_real_taps_lambda_mCtoT_Mbias.txt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    counts = [
        float(row[column]) for row in rows for column in row if "_COUNT_" in column
    ]
    assert max(counts) <= reads


@pytest.mark.xfail(
    strict=True, reason="--per_chromosome only changes the output file name"
)
def test_find_per_chromosome_only_reports_that_chromosome(tmp_path):
    reference = tmp_path / "two.fa"
    reference.write_text(">chrA\nACGTTCAGCCGATCCA\n>chrB\nTTCGACCAGGCTCCGA\n")
    assert run("find", "-f", reference, "-chr", "chrB", "-d", tmp_path).exit_code == 0
    chromosomes = {
        row.split("\t")[0] for row in data_rows(tmp_path / "two_chrB_all.bed")
    }
    assert chromosomes == {"chrB"}


def test_phred_quality_lines_starting_with_at_sign(tmp_path):
    """Fixed: quality lines starting with '@' (Phred 31) used to lose their first score."""
    fastq = tmp_path / "reads_1.fq.gz"
    fastq.write_bytes(gzip.compress(b"@read1\nTACG\n+\n@###\n@read2\nTACG\n+\n@###\n"))
    assert run("phred", "-1", fastq, "--se", "-d", tmp_path).exit_code == 0
    summary = (tmp_path / "reads_1_total_Phred.txt").read_text()
    assert "thymines: 31.0" in summary.splitlines()[2]
    assert "adenines: 2.0" in summary.splitlines()[2]
