from types import SimpleNamespace

import pytest

from astair2.read_context import cytosine_context, read_cytosine_calls

REFERENCE = "AAACAATACGCGgggCTGCATCCCGCGCCG"
# read positions of the top-strand Cs and bottom-strand Gs of REFERENCE, by context
TOP = {"CpG": {8, 10, 23, 25, 28}, "CHG": {15, 22, 27}, "CHH": {3, 18, 21}}
BOTTOM = {"CpG": {9, 11, 24, 26, 29}, "CHG": {12, 17}, "CHH": {13, 14}}


@pytest.mark.parametrize("top_strand, expected", [(True, TOP), (False, BOTTOM)])
def test_cytosine_contexts(top_strand, expected):
    base = "C" if top_strand else "G"
    positions = [
        index for index, letter in enumerate(REFERENCE) if letter.upper() == base
    ]
    found = {
        context: {
            p
            for p in positions
            if cytosine_context(top_strand, p, REFERENCE) == context
        }
        for context in expected
    }
    assert found == expected


def aligned_read(sequence):
    return SimpleNamespace(
        query_name="read",
        query_sequence=sequence,
        get_aligned_pairs=lambda: [(index, index) for index in range(len(sequence))],
    )


@pytest.mark.parametrize(
    "method, top_strand, read_sequence, modified",
    [
        ("mCtoT", True, REFERENCE.upper(), False),  # TAPS, nothing converted
        (
            "mCtoT",
            True,
            REFERENCE.upper().replace("C", "T"),
            True,
        ),  # TAPS, all converted
        ("CtoT", True, REFERENCE.upper(), True),  # bisulfite, nothing converted
        (
            "CtoT",
            False,
            REFERENCE.upper().replace("G", "A"),
            False,
        ),  # bisulfite, all converted
    ],
)
def test_read_cytosine_calls(method, top_strand, read_sequence, modified):
    expected = TOP if top_strand else BOTTOM
    calls = read_cytosine_calls(
        aligned_read(read_sequence), REFERENCE, top_strand, method
    )
    assert calls == {
        p: (context, modified)
        for context, positions in expected.items()
        for p in positions
    }


def test_read_beyond_reference_is_an_error():
    with pytest.raises(ValueError):
        read_cytosine_calls(aligned_read(REFERENCE + "CCC"), REFERENCE, True, "mCtoT")
