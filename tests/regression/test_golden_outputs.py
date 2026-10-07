"""Golden-output regression tests.

Every scenario in ``scenarios.py`` is run through the real ``astair2`` CLI and
all of its output files are normalised into one text snapshot, which must
match the stored copy in ``golden/`` byte for byte. Any refactoring that
changes results, file names or formatting fails here.

After an intentional change in output, regenerate the snapshots with
``pytest tests/regression --update-golden`` and review the diff.
"""

import difflib
import gzip
from pathlib import Path

import pysam
import pytest
from click.testing import CliRunner

from astair2.cli import cli
from tests.regression.scenarios import SCENARIOS

DATA = Path(__file__).resolve().parent.parent / "test_data"
GOLDEN = Path(__file__).resolve().parent / "golden"
SKIPPED_SUFFIXES = (".bai", ".crai", ".csi", ".tbi", ".fai")


def _alignment_dump(path):
    try:
        with pysam.AlignmentFile(str(path), check_sq=False) as bam:
            header = str(bam.header)
            reads = [read.to_string() for read in bam.fetch(until_eof=True)]
    except (ValueError, OSError):
        return "<unreadable alignment file>\n"
    return "\n".join([header, *reads]) + "\n"


def _render(path):
    if path.suffix in (".bam", ".cram"):
        return _alignment_dump(path)
    if path.suffix == ".pdf":
        return "<pdf, {} bytes>\n".format(
            "non-empty" if path.stat().st_size else "zero"
        )
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    return raw.decode("utf8")


def snapshot(out_dir, exit_codes):
    """Normalised, deterministic text form of everything a scenario produced."""
    parts = ["exit codes: {}\n".format(exit_codes)]
    for path in sorted(
        p
        for p in out_dir.rglob("*")
        if p.is_file() and not p.name.endswith(SKIPPED_SUFFIXES)
    ):
        parts.append("===== {} =====\n".format(path.relative_to(out_dir)))
        parts.append(_render(path))
    return "".join(parts).replace(str(out_dir), "<OUT>").replace(str(DATA), "<DATA>")


def run_scenario(name, out_dir):
    runner = CliRunner()
    exit_codes = []
    for args in SCENARIOS[name]:
        args = [a.format(data=DATA, out=out_dir) for a in args]
        exit_codes.append(runner.invoke(cli, args).exit_code)
    return snapshot(out_dir, exit_codes)


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_golden_output(name, tmp_path, request):
    actual = run_scenario(name, tmp_path)
    golden_file = GOLDEN / (name + ".txt.gz")
    if request.config.getoption("--update-golden"):
        GOLDEN.mkdir(exist_ok=True)
        golden_file.write_bytes(gzip.compress(actual.encode("utf8"), mtime=0))
        return
    if not golden_file.exists():
        pytest.fail("No golden snapshot for {}; run with --update-golden".format(name))
    expected = gzip.decompress(golden_file.read_bytes()).decode("utf8")
    if actual != expected:
        diff = difflib.unified_diff(
            expected.splitlines(),
            actual.splitlines(),
            "golden",
            "actual",
            lineterm="",
            n=2,
        )
        pytest.fail(
            "Output of scenario {} changed:\n{}".format(
                name, "\n".join(list(diff)[:80])
            ),
            pytrace=False,
        )
