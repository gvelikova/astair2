import re
from itertools import accumulate

_CIGAR_OPERATION = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_search(cigar_string):
    """Splits a CIGAR string into operation names, cumulative end positions and lengths.

    Zero-length operations count as one base for the cumulative positions."""
    operations = _CIGAR_OPERATION.findall(cigar_string)
    lengths = [int(length) for length, _ in operations]
    names = [name for _, name in operations]
    positions = list(accumulate(length or 1 for length in lengths))
    return names, positions, lengths


def _apply_operation(name, operation_index, end, length, positions):
    """Shifts read positions past a deletion or insertion, or drops positions inside a clip."""
    if name == "D":
        return [x if x < end else x - length for x in positions]
    if name == "I":
        return [x if x < end else x + length for x in positions]
    if name in ("S", "H"):
        return [x for x in positions if (x > end if operation_index == 0 else x < end)]
    return positions


def correct_positions_for_cigar(cigar_string, positions):
    """Corrects reference-derived read positions for the indels and clips of a CIGAR string."""
    for operation_index, (name, end, length) in enumerate(
        zip(*cigar_search(cigar_string))
    ):
        positions = _apply_operation(name, operation_index, end, length, positions)
    return positions


def position_correction_cigar(
    read, method, random_sample, positions, reverse_modification
):
    """Read offsets of the genomic positions (chrom, start, end) that the simulator changes in a read,
    corrected for its indels and clips.

    In CtoT mode the positions to change are those not in the random sample (the unmodified ones,
    which get converted), otherwise those in the sample. The genomic positions are converted to
    read offsets at the first indel or clip; that operation only shifts the offsets, it does not
    drop clipped ones."""
    names, ends, lengths = cigar_search(read.cigarstring)
    if not positions:
        return positions
    selected = random_sample.intersection(positions)
    if method == "CtoT" and not reverse_modification:
        selected = [x for x in positions if x not in selected]
    offset = abs(read.qstart - read.reference_start)
    converted = False
    for operation_index, (name, end, length) in enumerate(zip(names, ends, lengths)):
        if name not in ("D", "I", "S", "H"):
            continue
        if not converted:
            positions, converted = [x[1] - offset for x in selected], True
            if name in ("S", "H"):
                continue
        positions = _apply_operation(name, operation_index, end, length, positions)
    return positions
