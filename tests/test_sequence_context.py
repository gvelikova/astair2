import unittest

import pytest

from astair2 import context_search as context


class SequenceSearchOutputTest(unittest.TestCase):
    """Tests whether the expected cytosine contexts are found on the top and bottom strands.
    context.find_cytosine_contexts returns {(chromosome, start, end): (specific context, context,
    strand base)} for the motif keys it is given, and user-defined contexts are reported at their
    first cytosine."""

    def test_context_search_CHH_top(self):
        """Tests whether the expected CHH positions will be discovered on the top DNA strand."""
        data_context = {}
        for key in ["CHH"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AACTTCATCACT", "test_string", (key,), None
                )[0]
            )
        self.assertEqual(
            data_context,
            {
                ("test_string", 2, 3): ("CTT", "CHH", "C"),
                ("test_string", 5, 6): ("CAT", "CHH", "C"),
                ("test_string", 8, 9): ("CAC", "CHH", "C"),
            },
        )

    def test_context_search_CHH_bottom(self):
        """Tests whether the expected CHH positions will be discovered on the bottom DNA strand."""
        data_context = {}
        for key in ["CHHb"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AAGGCTTTGccc", "test_string", (key,), None
                )[0]
            )
        self.assertEqual(
            data_context,
            {
                ("test_string", 2, 3): ("CTT", "CHH", "G"),
                ("test_string", 3, 4): ("CCT", "CHH", "G"),
                ("test_string", 8, 9): ("CAA", "CHH", "G"),
            },
        )

    def test_context_search_CHG_top(self):
        """Tests whether the expected CHG positions will be discovered on the top DNA strand."""
        data_context = {}
        for key in ["CHG"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AACTTCAGCACT", "test_string", (key,), None
                )[0]
            )
        self.assertEqual(data_context, {("test_string", 5, 6): ("CAG", "CHG", "C")})

    def test_context_search_CHG_bottom(self):
        """Tests whether the expected CHG positions will be discovered on the bottom DNA strand."""
        data_context = {}
        for key in ["CHGb"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AACTTCAGCACT", "test_string", (key,), None
                )[0]
            )
        self.assertEqual(data_context, {("test_string", 7, 8): ("CTG", "CHG", "G")})

    def test_context_search_CpG(self):
        """Tests whether the expected CpG positions will be discovered."""
        data_context = {}
        for key in ["CG", "CGb"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AAGCGTTTGccc", "test_string", (key,), None
                )[0]
            )
        self.assertEqual(
            data_context,
            {
                ("test_string", 3, 4): ("CGT", "CpG", "C"),
                ("test_string", 4, 5): ("CGC", "CpG", "G"),
            },
        )

    def test_context_search_user_short(self):
        """Tests whether the expected short context positions will be discovered."""
        data_context = {}
        for key in ["user"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AAcTGCGTTTGcccTTGAC", "test_string", (key,), "AAcTG"
                )[0]
            )
        self.assertEqual(
            data_context,
            {
                ("test_string", 2, 3): ("AAcTG", "user defined context", "C"),
                ("test_string", 16, 17): ("AAcTG", "user defined context", "G"),
            },
        )

    def test_context_search_user_long(self):
        """Tests whether the expected long context positions will be discovered."""
        data_context = {}
        for key in ["user"]:
            data_context.update(
                context.find_cytosine_contexts(
                    "AATTTCCCGAcgtAAcTGTTaaagggcTgcACGTTTGcccTTGAC",
                    "test_string",
                    (key,),
                    "AATTTCCCGAcgt",
                )[0]
            )
        self.assertEqual(
            data_context,
            {
                ("test_string", 5, 6): ("AATTTCCCGAcgt", "user defined context", "C"),
                ("test_string", 28, 29): ("AATTTCCCGAcgt", "user defined context", "G"),
            },
        )

    def test_user_context_without_cytosine_is_rejected(self):
        """A user-provided context must contain a cytosine."""
        with pytest.raises(ValueError):
            context.context_keys("CpG", "AATTT")


if __name__ == "__main__":
    unittest.main()
