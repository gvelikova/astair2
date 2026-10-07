"""Runs `astair2 align` with stub bwa/bwameth.py/samtools executables that record their arguments."""

import os
import stat
from pathlib import Path

import pytest
from click.testing import CliRunner

from astair2.cli import cli

DATA = Path(__file__).resolve().parent / "test_data"

STUB = """#!/bin/sh
printf "%s %s\\n" "$(basename "$0")" "$*" >> "{log}"
case "$*" in
  *" index "*|"index "*) touch "$2.bwt" "$2.bwameth.c2t" 2>/dev/null; touch "$2.bai" ;;
  *) [ -p /dev/stdin ] && cat > /dev/null; echo "{name} output" ;;
esac
"""


@pytest.fixture
def tools(tmp_path):
    log = tmp_path / "calls.log"
    paths = {}
    for name in ("bwa", "bwameth.py", "samtools"):
        path = tmp_path / name
        path.write_text(STUB.format(log=log, name=name))
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        paths[name] = str(path)
    return paths, log


@pytest.fixture
def reference(tmp_path):
    copy = tmp_path / "lambda.fa"
    copy.write_text((DATA / "lambda_phage.fa").read_text())
    return str(copy)


def calls(log):
    return [line.split() for line in log.read_text().splitlines()]


def call(log, *prefix):
    """The one recorded call starting with prefix; pipeline stages start concurrently, so order is not fixed."""
    matching = [c for c in calls(log) if c[: len(prefix)] == list(prefix)]
    assert len(matching) == 1, calls(log)
    return matching[0]


def test_taps_pair_end(tools, reference, tmp_path):
    paths, log = tools
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        cli,
        [
            "align",
            "-1",
            str(DATA / "small_real_taps_lambda_1.fq.gz"),
            "-2",
            str(DATA / "small_real_taps_lambda_2.fq.gz"),
            "-f",
            reference,
            "-bp",
            paths["bwa"],
            "-sp",
            paths["samtools"],
            "-d",
            str(out),
            "-rg",
            r"@RG\tID:a\tSM:b",
        ],
    )
    assert result.exit_code == 0, result.output
    assert call(log, "bwa", "index") == ["bwa", "index", reference]
    bwa = call(log, "bwa", "mem")
    assert bwa[:2] == ["bwa", "mem"] and bwa[-3:] == [
        reference,
        str(DATA / "small_real_taps_lambda_1.fq.gz"),
        str(DATA / "small_real_taps_lambda_2.fq.gz"),
    ]
    assert bwa[bwa.index("-R") + 1] == r"@RG\tID:a\tSM:b"
    assert call(log, "samtools", "view") == [
        "samtools",
        "view",
        "-hb",
        "-T",
        reference,
        "-q",
        "1",
        "-F",
        "4",
        "-O",
        "BAM",
    ]
    assert call(log, "samtools", "sort") == [
        "samtools",
        "sort",
        "-m",
        "768M",
        "-@",
        "1",
        "-O",
        "BAM",
    ]
    output = out / "small_real_taps_lambda_mCtoT.bam"
    # indexing runs after the whole pipeline has finished
    assert calls(log)[-1] == ["samtools", "index", str(output)]
    assert output.read_text() == "samtools output\n"


def test_wgbs_single_end_keeps_unmapped(tools, reference, tmp_path):
    paths, log = tools
    out = tmp_path / "out"
    out.mkdir()
    fastq = str(DATA / "small_real_taps_lambda_mCtoT_SE.fq.gz")
    result = CliRunner().invoke(
        cli,
        [
            "align",
            "-1",
            fastq,
            "-se",
            "-m",
            "CtoT",
            "-u",
            "-O",
            "CRAM",
            "-f",
            reference,
            "-bp",
            paths["bwameth.py"],
            "-sp",
            paths["samtools"],
            "-d",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert call(log, "bwameth.py", "-t") == [
        "bwameth.py",
        "-t",
        "1",
        "--reference",
        reference,
        fastq,
    ]
    assert call(log, "samtools", "view") == [
        "samtools",
        "view",
        "-hC",
        "-T",
        reference,
        "-q",
        "0",
        "-O",
        "CRAM",
    ]
    assert call(log, "samtools", "sort") == [
        "samtools",
        "sort",
        "-T",
        os.path.join(str(out), "temp"),
        "-@",
        "1",
        "-O",
        "CRAM",
    ]
    assert (out / "small_real_taps_lambda_mCtoT_SE_CtoT.cram").exists()


def test_failed_aligner_is_not_indexed(tools, reference, tmp_path):
    paths, log = tools
    Path(paths["bwa"]).write_text(
        '#!/bin/sh\necho "bwa $*" >> "{}"\nexit 3\n'.format(log)
    )
    (Path(reference + ".bwt")).touch()
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        cli,
        [
            "align",
            "-1",
            str(DATA / "small_real_taps_lambda_1.fq.gz"),
            "-2",
            str(DATA / "small_real_taps_lambda_2.fq.gz"),
            "-f",
            reference,
            "-bp",
            paths["bwa"],
            "-sp",
            paths["samtools"],
            "-d",
            str(out),
        ],
    )
    assert result.exit_code == 1
    assert not any(call[:2] == ["samtools", "index"] for call in calls(log))
    assert list(out.iterdir()) == []


def taps_align_args(reference, out, *extra):
    return [
        "align",
        "-1",
        str(DATA / "small_real_taps_lambda_1.fq.gz"),
        "-2",
        str(DATA / "small_real_taps_lambda_2.fq.gz"),
        "-f",
        reference,
        "-d",
        str(out),
        *extra,
    ]


def test_failed_samtools_leaves_no_output(tools, reference, tmp_path):
    paths, log = tools
    Path(paths["samtools"]).write_text(
        '#!/bin/sh\necho "samtools $*" >> "{}"\nexit 1\n'.format(log)
    )
    (Path(reference + ".bwt")).touch()
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        cli,
        taps_align_args(reference, out, "-bp", paths["bwa"], "-sp", paths["samtools"]),
    )
    assert result.exit_code == 1
    assert list(out.iterdir()) == []


def test_missing_aligner_is_reported(tools, reference, tmp_path, caplog):
    paths, _ = tools
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        cli,
        taps_align_args(
            reference,
            out,
            "-bp",
            str(tmp_path / "no_such_bwa"),
            "-sp",
            paths["samtools"],
        ),
    )
    assert result.exit_code == 1
    assert "no_such_bwa was not found" in caplog.text
    assert "--bwa_path" in caplog.text


def test_missing_samtools_on_path_is_reported(
    tools, reference, tmp_path, caplog, monkeypatch
):
    paths, _ = tools
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        cli, taps_align_args(reference, out, "-bp", paths["bwa"])
    )
    assert result.exit_code == 1
    assert "samtools was not found" in caplog.text


def test_failed_index_is_reported(tools, reference, tmp_path):
    paths, log = tools
    Path(paths["bwa"]).write_text(
        '#!/bin/sh\necho "bwa $*" >> "{}"\nexit 2\n'.format(log)
    )
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        cli,
        taps_align_args(reference, out, "-bp", paths["bwa"], "-sp", paths["samtools"]),
    )
    assert result.exit_code == 1
    assert calls(log) == [["bwa", "index", reference]]


def test_existing_output_is_not_overwritten(tools, reference, tmp_path):
    paths, log = tools
    out = tmp_path / "out"
    out.mkdir()
    existing = out / "small_real_taps_lambda_mCtoT.bam"
    existing.write_text("previous result")
    result = CliRunner().invoke(
        cli,
        taps_align_args(reference, out, "-bp", paths["bwa"], "-sp", paths["samtools"]),
    )
    assert result.exit_code == 1
    assert existing.read_text() == "previous result"
    assert not log.exists()
