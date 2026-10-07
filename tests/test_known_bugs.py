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
