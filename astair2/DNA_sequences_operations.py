_COMPLEMENT = str.maketrans("TtCcGgAa", "AaGgCcTt")


def complementary(sequence):
    """Takes an input DNA string and gives its complementary."""
    return sequence.translate(_COMPLEMENT)


def reverse_complementary(sequence):
    """Takes an input DNA string and gives its reverse complementary."""
    return sequence[::-1].translate(_COMPLEMENT)


def reverse(sequence):
    """Takes an input DNA string and gives its reverse."""
    return sequence[::-1]
