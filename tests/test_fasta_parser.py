import gzip

import pysam
import pytest

from astair2 import simple_fasta_parser as sfp

RECORDS = {
    "chr1": "ACTGCTCCCTGGaaaTCGGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGG",
    "chr2": "AAACCTGCcctGttug",
    "chr3": "CTGATCGTTTAGCAGCA",
}


def fasta_text(records, width=10, newline="\n", description=" some description"):
    lines = []
    for name, sequence in records.items():
        lines.append(">" + name + description)
        lines.extend(sequence[i : i + width] for i in range(0, len(sequence), width))
    return newline.join(lines) + newline


@pytest.fixture(params=["plain", "gzip", "bgzip"])
def reference(request, tmp_path):
    plain = tmp_path / "reference.fa"
    plain.write_text(fasta_text(RECORDS))
    if request.param == "plain":
        return str(plain)
    if request.param == "gzip":
        compressed = tmp_path / "reference.fa.gz"
        compressed.write_bytes(gzip.compress(plain.read_bytes()))
        return str(compressed)
    compressed = str(tmp_path / "reference.fa.gz")
    pysam.tabix_compress(str(plain), compressed)
    return compressed


def test_reference_names(reference):
    assert sfp.reference_names(reference) == ["chr1", "chr2", "chr3"]


def test_read_whole_reference(reference):
    assert sfp.read_reference(reference) == RECORDS


def test_read_selected_sequence_of_multiline_reference(reference):
    assert sfp.read_reference(reference, ["chr2"]) == {"chr2": RECORDS["chr2"]}


def test_missing_sequence_is_an_error(reference):
    with pytest.raises(KeyError):
        sfp.read_reference(reference, ["chrX"])


def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        sfp.read_reference(str(tmp_path / "absent.fa"))


def test_records_with_windows_line_endings_and_blank_lines():
    text = fasta_text(RECORDS, newline="\r\n").replace(">chr2", "\r\n\r\n>chr2")
    assert dict(sfp.fasta_records(text.splitlines(keepends=True))) == RECORDS


def test_reference_with_underscores(tmp_path):
    plain = tmp_path / "with spaces.fa"
    plain.write_text(fasta_text(RECORDS))
    target = sfp.write_reference_with_underscores(str(plain))
    assert target == str(tmp_path / "with spaces_no_spaces.fa.gz")
    with pysam.FastaFile(target) as fasta:
        assert list(fasta.references) == [
            "chr1_some_description",
            "chr2_some_description",
            "chr3_some_description",
        ]
        assert fasta.fetch("chr3_some_description") == RECORDS["chr3"]
    # an existing copy is reused
    assert sfp.write_reference_with_underscores(str(plain)) == target
