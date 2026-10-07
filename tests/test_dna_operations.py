import sys
import unittest

from astair2 import DNA_sequences_operations as so


class DNAOperationsTest(unittest.TestCase):
    """Tests whether the resulting DNA sequences match the expected reverse, complementary and reverse complementary output."""

    def test_reverse_DNA(self):
        """Tests whether the resulting DNA sequences match the expected reverse output."""
        self.assertEqual(
            so.reverse("AAAACTCTGCCGggAaaccttCG"), "GCttccaaAggGCCGTCTCAAAA"
        )

    def test_complementary_DNA(self):
        """Tests whether the resulting DNA sequences match the expected complementary output."""
        self.assertEqual(
            so.complementary("AAAACTCTGCCGggAaaccttCG"), "TTTTGAGACGGCccTttggaaGC"
        )

    def test_reverse_complementary_DNA(self):
        """Tests whether the resulting DNA sequences match the expected reverse complementary output."""
        self.assertEqual(
            so.reverse_complementary("AAAACTCTGCCGggAaaccttCG"),
            "CGaaggttTccCGGCAGAGTTTT",
        )


if __name__ == "__main__":
    unittest.main()


def test_cigar_with_sequence_match_operations():
    """Fixed: '=' operations used to be skipped, misaligning names and lengths."""
    from astair2.cigar_search import cigar_search, correct_positions_for_cigar

    assert cigar_search("5=2X3I4=") == (
        ["=", "X", "I", "="],
        [5, 7, 10, 14],
        [5, 2, 3, 4],
    )
    assert correct_positions_for_cigar("5=3I4=", [2, 6, 9]) == [2, 6, 12]
