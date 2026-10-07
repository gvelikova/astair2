"""Locating cytosines of the requested sequence contexts in a reference sequence.

Every context is a fixed-length motif matched case-insensitively. Forward
motifs are searched on the given strand and report the position of their
first base (a C); reverse motifs are searched on the complementary strand and
report the position of their last base (the G opposite the C).
"""

import re
from collections import Counter

from astair2.DNA_sequences_operations import complementary

USER_CONTEXT = "user defined context"

# key: (regex on the upper-cased strand, context label, searched on the complement)
MOTIFS = {
    "CHG": ("C[ACT]G", "CHG", False),
    "CHGb": ("G[ACT]C", "CHG", True),
    "CHH": ("C[ACT][ACT]", "CHH", False),
    "CHHb": ("[ACT][ACT]C", "CHH", True),
    "CG": ("CG[ACTG]", "CpG", False),
    "CGb": ("[ACTG]GC", "CpG", True),
    "CN": ("CN[NACT]", "CN", False),
    "CNb": ("N[NACT]C", "CN", True),
}

CONTEXT_KEYS = {
    "all": ("CHG", "CHGb", "CHH", "CHHb", "CG", "CGb", "CN", "CNb"),
    "CpG": ("CG", "CGb"),
    "CHG": ("CHG", "CHGb"),
    "CHH": ("CHH", "CHHb"),
}


def context_keys(context, user_defined_context):
    """The motif keys searched for a context choice, in the order they are applied."""
    if user_defined_context and "C" not in user_defined_context.upper():
        raise ValueError(
            "The user-defined context {} does not contain a cytosine.".format(
                user_defined_context
            )
        )
    return CONTEXT_KEYS[context] + (("user",) if user_defined_context else ())


def _overlapping_matches(pattern, sequence):
    """Start index and text of every (possibly overlapping) match of pattern."""
    return (
        (match.start(), match.group(1))
        for match in re.finditer("(?=({}))".format(pattern), sequence)
    )


def _in_region(position, region):
    return region is None or (position >= region[1] and position + 1 <= region[2])


def _motif_hits(key, upper, complement):
    """(position, specific context, context label, reference base) for one motif."""
    pattern, label, on_complement = MOTIFS[key]
    if on_complement:
        return (
            (start + 2, found[::-1], label, "G")
            for start, found in _overlapping_matches(pattern, complement)
        )
    return (
        (start, found, label, "C")
        for start, found in _overlapping_matches(pattern, upper)
    )


def _user_hits(user_defined_context, upper, complement):
    """(position, counted key, user context, label, reference base) for a user-defined context."""
    motif = user_defined_context.upper()
    forward_offset = motif.find("C")
    reverse_offset = len(motif) - 1 - motif[::-1].find("C")
    forward = (
        (start + forward_offset, user_defined_context, USER_CONTEXT, "C")
        for start, _ in _overlapping_matches(re.escape(motif), upper)
    )
    reverse = (
        (start + reverse_offset, user_defined_context, USER_CONTEXT, "G")
        for start, _ in _overlapping_matches(re.escape(motif), complement)
    )
    return forward, reverse


def find_cytosine_contexts(
    sequence,
    sequence_name,
    keys,
    user_defined_context=None,
    region=None,
    stranded=False,
):
    """Finds the cytosines of the requested contexts in a sequence.

    Returns a dictionary {(name, start, end[, strand]): (specific context, context, reference base)}
    and a Counter of how many times each motif key and specific context occurs in the whole
    sequence, irrespective of the region.
    """
    upper = sequence.upper()
    complement = complementary(upper)
    positions, counts = {}, Counter()
    for key in keys:
        if key == "user":
            hit_groups, counted = (
                _user_hits(user_defined_context, upper, complement),
                None,
            )
        else:
            hit_groups, counted = (_motif_hits(key, upper, complement),), key
        for hits in hit_groups:
            for position, specific, label, base in hits:
                counts[counted or USER_CONTEXT] += 1
                if counted:
                    counts[specific] += 1
                if _in_region(position, region):
                    strand = ("+" if base == "C" else "-",) if stranded else ()
                    positions[(sequence_name, position, position + 1) + strand] = (
                        specific,
                        label,
                        base,
                    )
    return positions, counts
