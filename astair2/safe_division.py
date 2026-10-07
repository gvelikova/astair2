def non_zero_division(x, y, sign):
    """x / y, or 0 if only y is 0, or sign if both are 0."""
    if y == 0:
        return 0 if x != 0 else sign
    return x / y


def non_zero_division_NA(x, y):
    """x / y, or 'NA' if y is 0, or 0 if x is 0."""
    if y == 0:
        return "NA"
    if x == 0:
        return 0
    return x / y


def safe_rounder(data, precision, multiple):
    """Rounds numbers, as a percentage if multiple is true, and leaves strings untouched."""
    if isinstance(data, str):
        return data
    return round(data * 100 if multiple else data, precision)
