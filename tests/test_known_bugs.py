"""Bugs inherited from asTair, pinned as strict expected failures.

Each test describes the correct behaviour. While the bug exists the test fails
and is reported as xfail; once the bug is fixed the test passes, and because
the marker is strict, pytest reports it as XPASS(strict) so the marker gets
removed together with the fix.
"""

import gzip
from pathlib import Path

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
    reason="flags_expectation calls list.extend with two arguments; "
    "clean_pileup swallows the TypeError and drops every position",
)
def test_call_keeps_positions_when_orphans_are_included(tmp_path):
    result = run(
        "call", "-i", TAPS, "-f", LAMBDA, "-co", "CpG", "-io", "False", "-d", tmp_path
    )
    assert result.exit_code == 0
    assert data_rows(tmp_path / "small_real_taps_lambda_mCtoT_mCtoT_CpG.mods")


@pytest.mark.xfail(
    strict=True,
    reason="modification_calls_writer is called with an unsupported header= argument",
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
