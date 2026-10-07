"""End-to-end CLI scenarios whose outputs are pinned by golden snapshots.

Each scenario is a list of CLI invocations (argument lists for the
``astair2`` command). ``{data}`` is replaced by the test data directory and
``{out}`` by a fresh, scenario-specific output directory.
"""

LAMBDA = "{data}/lambda_phage.fa"
TAPS = "{data}/small_real_taps_lambda_mCtoT.bam"
TAPS_SE = "{data}/small_real_taps_lambda_mCtoT_SE.bam"
TAPS_REV = "{data}/small_real_taps_lambda_mCtoT_reversed.bam"
WGBS = "{data}/small_real_wgbs_lambda_CtoT.bam"
SYNTH = "{data}/small_lambda.bam"
HG38 = "{data}/hg38_chr10_20000-60000.fa"
HG38_BAM = "{data}/small_real_taps_chr10:20000-60000_pos.bam"
VCF = "{data}/GRCh38p7_common_snps_sample.vcf.gz"
FQ1 = "{data}/small_real_taps_lambda_1.fq.gz"
FQ2 = "{data}/small_real_taps_lambda_2.fq.gz"
FQ_SE = "{data}/small_real_taps_lambda_mCtoT_SE.fq.gz"
OUT = ["-d", "{out}"]

SCENARIOS = {
    # call
    "call_taps_all": [["call", "-i", TAPS, "-f", LAMBDA, *OUT]],
    "call_taps_cpg_gz": [
        [
            "call",
            "-i",
            TAPS,
            "-f",
            LAMBDA,
            "-co",
            "CpG",
            "--gz",
            "-ni",
            "0",
            "-bq",
            "13",
            *OUT,
        ]
    ],
    "call_wgbs_chg": [
        [
            "call",
            "-i",
            WGBS,
            "-f",
            LAMBDA,
            "-m",
            "CtoT",
            "-co",
            "CHG",
            "-ni",
            "NA",
            *OUT,
        ]
    ],
    "call_se_chh": [
        ["call", "-i", TAPS_SE, "-f", LAMBDA, "-se", "-co", "CHH", "-ni", ".", *OUT]
    ],
    "call_user_context": [
        ["call", "-i", TAPS, "-f", LAMBDA, "-co", "CpG", "-uc", "CAG", *OUT]
    ],
    "call_vcf": [["call", "-i", HG38_BAM, "-f", HG38, "-ks", VCF, *OUT]],
    "call_chr_clip": [
        [
            "call",
            "-i",
            TAPS,
            "-f",
            LAMBDA,
            "-chr",
            "lambda",
            "-scl",
            "5",
            "-ecl",
            "5",
            "-co",
            "CpG",
            *OUT,
        ]
    ],
    "call_reverse_library": [
        ["call", "-i", TAPS_REV, "-f", LAMBDA, "-li", "reverse", "-co", "CpG", *OUT]
    ],
    "call_zero_coverage_orphans": [
        ["call", "-i", TAPS, "-f", LAMBDA, "-co", "CHG", "-zc", "-io", "False", *OUT]
    ],
    "call_keep_overlaps_gz": [
        [
            "call",
            "-i",
            WGBS,
            "-f",
            LAMBDA,
            "-m",
            "CtoT",
            "-co",
            "CpG",
            "-kc",
            "--gz",
            *OUT,
        ]
    ],
    "call_filters": [
        [
            "call",
            "-i",
            TAPS,
            "-f",
            LAMBDA,
            "-co",
            "CpG",
            "-mq",
            "10",
            "-md",
            "50",
            "-cbq",
            "False",
            *OUT,
        ]
    ],
    # find
    "find_all": [["find", "-f", LAMBDA, *OUT]],
    "find_cpg_gz": [["find", "-f", LAMBDA, "-co", "CpG", "--gz", *OUT]],
    "find_user_context": [["find", "-f", LAMBDA, "-co", "CHG", "-uc", "CAG", *OUT]],
    "find_chr": [["find", "-f", HG38, "-co", "CHH", "-chr", "chr10", *OUT]],
    # mbias
    "mbias_pe": [["mbias", "-f", LAMBDA, "-i", TAPS, "-l", "80", *OUT]],
    "mbias_se": [["mbias", "-f", LAMBDA, "-i", TAPS_SE, "-l", "75", "-se", *OUT]],
    "mbias_wgbs_na": [
        [
            "mbias",
            "-f",
            LAMBDA,
            "-i",
            WGBS,
            "-l",
            "80",
            "-m",
            "CtoT",
            "-ni",
            "NA",
            "-chr",
            "lambda",
            *OUT,
        ]
    ],
    "mbias_plot": [
        [
            "mbias",
            "-f",
            LAMBDA,
            "-i",
            TAPS,
            "-l",
            "80",
            "-p",
            "-c",
            "red,blue,green",
            *OUT,
        ]
    ],
    # idbias
    "idbias_pe": [["idbias", "-f", LAMBDA, "-i", TAPS, "-l", "80", *OUT]],
    "idbias_se": [["idbias", "-f", LAMBDA, "-i", TAPS_SE, "-l", "75", "-se", *OUT]],
    "idbias_wgbs_plot": [
        [
            "idbias",
            "-f",
            LAMBDA,
            "-i",
            WGBS,
            "-l",
            "80",
            "-m",
            "CtoT",
            "-ni",
            "*",
            "-p",
            *OUT,
        ]
    ],
    # phred
    "phred_pe": [["phred", "-1", FQ1, "-2", FQ2, *OUT]],
    "phred_pe_absolute_plot": [
        ["phred", "-1", FQ1, "-2", FQ2, "-cm", "absolute", "-p", *OUT]
    ],
    "phred_se": [["phred", "-1", FQ_SE, "--se", *OUT]],
    # filter
    "filter_taps": [["filter", "-f", LAMBDA, "-i", TAPS, *OUT]],
    "filter_wgbs": [
        ["filter", "-f", LAMBDA, "-i", WGBS, "-m", "CtoT", "--bases_noncpg", "2", *OUT]
    ],
    "filter_se": [["filter", "-f", LAMBDA, "-i", TAPS_SE, "--se", *OUT]],
    # simulate (always seeded, otherwise the output is random)
    "simulate_chg_taps": [
        [
            "simulate",
            "-f",
            LAMBDA,
            "-l",
            "75",
            "-i",
            SYNTH,
            "-ml",
            "50",
            "-co",
            "CHG",
            "-s",
            "7",
            *OUT,
        ]
    ],
    "simulate_cpg_wgbs": [
        [
            "simulate",
            "-f",
            LAMBDA,
            "-l",
            "75",
            "-i",
            SYNTH,
            "-m",
            "CtoT",
            "-ml",
            "30",
            "-co",
            "CpG",
            "-s",
            "7",
            *OUT,
        ]
    ],
    "simulate_all_reverse": [
        [
            "simulate",
            "-f",
            LAMBDA,
            "-l",
            "80",
            "-i",
            TAPS,
            "-ml",
            "100",
            "-co",
            "all",
            "-rv",
            "-s",
            "3",
            *OUT,
        ]
    ],
    "simulate_region": [
        [
            "simulate",
            "-f",
            LAMBDA,
            "-l",
            "80",
            "-i",
            TAPS,
            "-ml",
            "60",
            "-co",
            "CpG",
            "-r",
            "lambda",
            "100",
            "20000",
            "-s",
            "11",
            *OUT,
        ]
    ],
    "simulate_chr_wgbs_reverse": [
        [
            "simulate",
            "-f",
            LAMBDA,
            "-l",
            "80",
            "-i",
            WGBS,
            "-m",
            "CtoT",
            "-ml",
            "40",
            "-co",
            "CHH",
            "-chr",
            "lambda",
            "-rv",
            "-s",
            "5",
            *OUT,
        ]
    ],
    "simulate_from_list": [
        ["call", "-i", TAPS, "-f", LAMBDA, "-co", "CpG", "-bq", "13", "-ni", "0", *OUT],
        [
            "simulate",
            "-f",
            LAMBDA,
            "-l",
            "80",
            "-i",
            TAPS,
            "-co",
            "all",
            "-mp",
            "{out}/small_real_taps_lambda_mCtoT_mCtoT_CpG.mods",
            *OUT,
        ],
    ],
    # summarise
    "summarise_vcf": [
        ["summarise", "-i", HG38_BAM, "-f", HG38, "-ks", VCF, "-co", "CpG", *OUT]
    ],
    "summarise_wgbs_region": [
        [
            "summarise",
            "-i",
            HG38_BAM,
            "-f",
            HG38,
            "-m",
            "CtoT",
            "-co",
            "all",
            "-r",
            "chr10",
            "20000",
            "30000",
            *OUT,
        ]
    ],
    "summarise_lambda_chh": [
        ["summarise", "-i", TAPS, "-f", LAMBDA, "-co", "CHH", "-bq", "20", *OUT]
    ],
}
