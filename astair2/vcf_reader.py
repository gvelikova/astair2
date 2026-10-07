import logging
import os
from itertools import product

import pysam

from astair2.DNA_sequences_operations import reverse_complementary

logs = logging.getLogger(__name__)

_CONTEXT_OF = {
    # later groups override earlier ones (CGN is CN, not CpG)
    **{"C" + x + z: "CHH" for x, z in product("ACT", "ACT")},
    **{"C" + x + "G": "CHG" for x in "ACT"},
    **{"CG" + x: "CpG" for x in "ACTGN"},
    **{"CN" + z: "CN" for z in "NACTG"},
    **{"C" + z + "N": "CN" for z in "NACTG"},
}


def _variant_context(variant, fasta):
    """(specific context, context, strand base) of a variant whose alternative alleles include C or G.

    Near the sequence ends, or beyond them, the specific context is shorter than three bases and
    the context is CN."""
    if "C" in variant.alts:
        base = "C"
        subcontext = (
            "C" + fasta[variant.start + 1 : variant.start + 3].upper()
            if variant.start + 3 < len(fasta)
            else "CNN"
        )
    else:
        base = "G"
        subcontext = (
            reverse_complementary(
                fasta[variant.start - 2 : variant.start].upper() + "G"
            )
            if variant.start - 2 >= 0
            else "CNN"
        )
    return subcontext, _CONTEXT_OF[subcontext] if len(subcontext) == 3 else "CN", base


def read_vcf(vcf_file, chromosome, fasta, threads, start, end):
    """Reads the variants of one chromosome (optionally restricted to start-end) from a VCF file.

    Returns the positions (chrom, start, end) of variants with a C or G reference allele, and
    {position: (specific context, context, strand base, reference allele)} for variants with a C or
    G alternative allele, which may create new cytosines. The FILTER column is not taken into account.
    """
    if os.path.isfile(vcf_file) and not os.path.isfile(vcf_file + ".tbi"):
        logs.info(
            "Creating a TABIX index for the provided VCF file. This may take a while..."
        )
        pysam.tabix_index(vcf_file, preset="vcf")
    true_variants, possible_mods = set(), {}
    with pysam.VariantFile(vcf_file, "r", threads=threads) as variants:
        name = next(
            (
                name
                for name in (chromosome, chromosome[3:])
                if variants.is_valid_reference_name(name)
            ),
            None,
        )
        if name is None:
            raise ValueError(
                "The chromosome {} is not in the VCF file {}; the chromosome names must match "
                "those in the FASTA and BAM files.".format(chromosome, vcf_file)
            )
        for variant in variants.fetch(name):
            variant_chrom = (
                "chr" + variant.chrom
                if variant.chrom.isnumeric() and not chromosome.isnumeric()
                else variant.chrom
            )
            if variant_chrom != chromosome:
                continue
            if not (start is None and end is None) and not (
                variant.start >= start and variant.start + 1 <= end
            ):
                continue
            position = (variant_chrom, variant.start, variant.start + 1)
            if variant.ref in ("C", "G"):
                true_variants.add(position)
            if variant.alts is not None and {"C", "G"} & set(variant.alts):
                possible_mods[position] = _variant_context(variant, fasta) + (
                    variant.ref,
                )
    return true_variants, possible_mods
