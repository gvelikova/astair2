"""Reading reference sequences from FASTA files.

Plain and BGZIP-compressed files are read through their faidx index (which
pysam creates for plain files if it is missing). Files compressed with plain
gzip cannot be indexed and are parsed in full instead.
"""

import gzip
import os

import pysam


def _is_plain_gzip(fasta_file):
    """True for gzip files that are not BGZF (BGZF sets the FEXTRA flag with a 'BC' subfield)."""
    with open(fasta_file, "rb") as handle:
        header = handle.read(14)
    return header[:2] == b"\x1f\x8b" and not (header[3] & 4 and header[12:14] == b"BC")


def _indexed(fasta_file):
    """An open pysam.FastaFile, or None if the file cannot be indexed (plain gzip)."""
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(
            "The reference FASTA file {} does not exist.".format(fasta_file)
        )
    return None if _is_plain_gzip(fasta_file) else pysam.FastaFile(fasta_file)


def _open_text(fasta_file):
    return (
        gzip.open(fasta_file, "rt") if fasta_file.endswith(".gz") else open(fasta_file)
    )


def fasta_records(lines):
    """Yields (name, sequence) for each FASTA record; the name is the first word of the header."""
    name, chunks = None, []
    for line in lines:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = (line[1:].split() or [""])[0], []
        elif name is not None:
            chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def reference_names(fasta_file):
    """The sequence names in a FASTA file, in file order."""
    fasta = _indexed(fasta_file)
    if fasta is not None:
        with fasta:
            return list(fasta.references)
    with _open_text(fasta_file) as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]


def read_reference(fasta_file, names=None):
    """{name: sequence} for the requested sequence names, or for all sequences if names is None."""
    fasta = _indexed(fasta_file)
    if fasta is not None:
        with fasta:
            wanted = fasta.references if names is None else names
            missing = [name for name in wanted if name not in fasta.references]
            if missing:
                raise KeyError(
                    "Sequences {} are not in the reference FASTA file {}.".format(
                        missing, fasta_file
                    )
                )
            return {name: fasta.fetch(name) for name in wanted}
    with _open_text(fasta_file) as handle:
        sequences = dict(fasta_records(handle))
    if names is None:
        return sequences
    missing = [name for name in names if name not in sequences]
    if missing:
        raise KeyError(
            "Sequences {} are not in the reference FASTA file {}.".format(
                missing, fasta_file
            )
        )
    return {name: sequences[name] for name in names}


def underscored_reference_path(fasta_file):
    """Where the copy of a reference with underscores instead of spaces in its sequence names is kept."""
    absolute = os.path.abspath(fasta_file)
    stem = os.path.splitext(os.path.basename(os.path.splitext(absolute)[0]))[0]
    return os.path.join(os.path.dirname(absolute), stem + "_no_spaces.fa.gz")


def write_reference_with_underscores(fasta_file):
    """Writes a BGZIP-compressed copy of the reference in which spaces in the sequence names are
    replaced by underscores, unless it already exists, and returns its path."""
    target = underscored_reference_path(fasta_file)
    if not os.path.isfile(target):
        plain = target[: -len(".gz")]
        with _open_text(fasta_file) as source, open(plain, "w") as output:
            for line in source:
                line = line.rstrip("\r\n")
                output.write(
                    (line.replace(" ", "_") if line.startswith(">") else line) + "\n"
                )
        pysam.tabix_compress(plain, target, force=True)
        os.remove(plain)
    return target
