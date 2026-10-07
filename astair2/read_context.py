"""Cytosine calls along a single read, shared by the M-bias and ID-bias reports."""

TOP_CHG = ("CAG", "CCG", "CTG")
TOP_CHH = ("CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT")
BOTTOM_CHG = ("CAG", "CGG", "CTG")
BOTTOM_CHH = ("AAG", "AGG", "ATG", "GAG", "GGG", "GTG", "TAG", "TGG", "TTG")


def strand_bases(top_strand):
    """(reference base, converted base) of the cytosines a read reports on."""
    return ("C", "T") if top_strand else ("G", "A")


def cytosine_context(top_strand, reference_index, sequence):
    """CpG, CHG or CHH for a cytosine of the read's strand at a reference position, else None."""
    if top_strand:
        if sequence[reference_index : reference_index + 2].upper() == "CG":
            return "CpG"
        trinucleotide = sequence[reference_index : reference_index + 3].upper()
        chg, chh = TOP_CHG, TOP_CHH
    else:
        if sequence[reference_index - 1 : reference_index + 1].upper() == "CG":
            return "CpG"
        trinucleotide = sequence[reference_index - 2 : reference_index + 1].upper()
        chg, chh = BOTTOM_CHG, BOTTOM_CHH
    return "CHG" if trinucleotide in chg else "CHH" if trinucleotide in chh else None


def reference_cytosines(read, sequence, reference_base):
    """(read index, reference index) of the read bases aligned to reference_base."""
    try:
        return [
            (read_index, reference_index)
            for read_index, reference_index in read.get_aligned_pairs()
            if reference_index is not None
            and read_index is not None
            and sequence[reference_index].upper() == reference_base
        ]
    except IndexError as error:
        raise ValueError(
            "Read {} is aligned beyond the end of its reference sequence.".format(
                read.query_name
            )
        ) from error


def read_cytosine_calls(read, sequence, top_strand, method):
    """{read index: (context, modified)} for the read bases at CpG, CHG and CHH cytosines."""
    reference_base, converted_base = strand_bases(top_strand)
    query = read.query_sequence
    calls = {}
    for read_index, reference_index in reference_cytosines(
        read, sequence, reference_base
    ):
        context = cytosine_context(top_strand, reference_index, sequence)
        if context:
            converted = query[read_index].upper() == converted_base
            calls[read_index] = (
                context,
                converted if method == "mCtoT" else not converted,
            )
    return calls
