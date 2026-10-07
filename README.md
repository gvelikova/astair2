# asTair2

_`asTair2` is a toolchain to process DNA modification sequencing data. It supports both methods that convert modified cytosines to thymines, such as [TET-Assisted Pyridine Borane Sequencing (TAPS)](https://www.nature.com/articles/s41587-019-0041-2) (`--method mCtoT`, the default), and methods that convert unmodified cytosines to thymines, such as whole-genome bisulfite sequencing (WGBS) (`--method CtoT`)._

`asTair2` is an updated and modernised continuation of [asTair](#authors-and-history). It runs on current Python (3.12+) and current versions of its dependencies, while keeping asTair's sub-commands, options and output formats. Further development, including refactoring beyond the original code base, happens here.

## Who asTair2 is for

asTair2 is intended for:

- **Existing asTair users** who need their pipelines to keep working on current Python and dependency versions, with the same commands and output formats.
- **Researchers working in Python**, who want to install with `pip`, run the tool from scripts or notebooks, and read or extend the code in a familiar language.
- **Automated and agent-based workflows** built on Python tooling, where a pip-installable package with a stable command-line interface is easiest to integrate.

## Related tools

[rastair](https://github.com/bsbludwig/rastair) ([rastair.com](https://www.rastair.com)) is a separate command-line tool, written in Rust, for the simultaneous detection of genetic variants and methylated positions from short-read data produced by TAPS, Illumina 5-Base sequencing or similar methods that convert 5mC to T. asTair2 flags likely homozygous C-to-T variants at the cytosine positions it reports, so that they are not mistaken for modifications (see [Interpret results](#3-interpret-results)); users who need general variant calling alongside methylation calling, or who prefer a compiled Rust tool, may want to consider rastair.

asTair2 and rastair are independent projects with different goals. asTair2 exists to keep the original asTair maintained, usable on current software and openly available under the GPL; it is not intended to replace or compete with rastair. Results from the two tools are not guaranteed to be identical, so please use one tool consistently within a study and cite the tool you used.

# Basic usage
## 0. Installation 

`asTair2` requires Python 3.12 or newer. We recommend installing it into a virtual environment, which keeps it separate from the system Python (many Linux distributions and Homebrew do not allow `pip` to install into the system Python):

```bash
git clone https://github.com/gvelikova/astair2.git
cd astair2
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
```

To also enable the plotting options (`--plot`), install the optional matplotlib dependency:

```bash
python -m pip install ".[plot]"
```

Activate the environment (`source .venv/bin/activate`) in each new terminal before running `astair2`.

All Python dependencies are installed by `pip`. The alignment step (`astair2 align`) additionally calls external programs, which must be on your `PATH`; all other sub-commands work without them:

| Data | Required for `astair2 align` |
| ---- | ---------------------------- |
| TAPS (`--method mCtoT`) | [`bwa`](https://github.com/lh3/bwa), [`samtools`](http://www.htslib.org/) |
| Bisulfite (`--method CtoT`) | [`bwa-meth`](https://github.com/brentp/bwa-meth) (which uses `bwa`), [`samtools`](http://www.htslib.org/) |

One way to install all of them is from [Bioconda](https://bioconda.github.io/):

```bash
conda install -c conda-forge -c bioconda bwa samtools bwameth
```

Alternatively, point `astair2 align` to existing installations with `--samtools_path` and `--bwa_path` (the path to `bwa` for TAPS data, or to `bwameth.py` for bisulfite data).

You should now be able to call `astair2`:

```bash
astair2 --help
```
```text
Usage: astair2 [OPTIONS] COMMAND [ARGS]...

  asTair2 (tools for processing cytosine modification sequencing data)

Options:
  --help  Show this message and exit.

Commands:
  align      Align raw reads in fastq format to a reference genome.
  call       Call modified cytosines from a bam or cram file.
  filter     Look for sequencing reads with more than N CpH modifications.
  find       Output positions of Cs from fasta file per context.
  idbias     Generate indel count per read length information (IDbias).
  mbias      Generate modification per read length information (Mbias).
  phred      Calculate per base (A, C, T, G) Phred scores for each strand.
  simulate   Simulate TAPS/BS conversion on top of an existing bam/cram...
  summarise  Collects and outputs modification information per read.

  __________________________________About__________________________________
  asTair2 is developed by Gergana V. Velikova. It is based on asTair, which
  was written by Gergana V. Velikova and Benjamin Schuster-Boeckler. This code
  is made available under the GNU General Public License v3, see LICENSE.txt
  for more details.

                                                           Version: 1.x.x
```

In general, you can use `--help` on all `astair2` sub-commands to get detailed instructions on the available options.

The tutorial below uses __TAPS__ paired-end sequencing reads. All analyses can also be run on single-end data with the `-se` option. For __bisulfite__ (WGBS) data, the same steps apply with `--method CtoT`; see [Analysis of WGBS data](#analysis-of-wgbs-data-or-other-unmodified-cytosine-to-thymine-conversion-methods).

## 1. Align reads

We will assume that you have generated paired-end sequencing data, which is stored in two fastq files. For this brief tutorial, we use the files `lambda.phage_test_sample_1.fq.gz` and `lambda.phage_test_sample_2.fq.gz` (about 80 MB in total), which you can download here:

```bash
# Or use curl -O if wget is not available
wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_1.fq.gz
wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_2.fq.gz
```

The raw reads need to be aligned. asTair2 contains a command to help with this. It requires [`bwa`](https://github.com/lh3/bwa) and [`samtools`](http://www.htslib.org/) (see [Installation](#0-installation)). (If you prefer to use a different aligner, [skip to step 2](#2-call-methylation).)

You will also need an indexed reference genome to align to, which can be given as a bgzip compressed file. For this example we are using the lambda phage genome, which you can download with

```bash
wget https://zenodo.org/record/2582855/files/lambda_phage.fa
wget https://zenodo.org/record/2582855/files/lambda_phage.fa.fai
```

Now, you are ready to align:
```bash
mkdir -p output_dir
astair2 align -f lambda_phage.fa -1 lambda.phage_test_sample_1.fq.gz -2 lambda.phage_test_sample_2.fq.gz -d output_dir
```
This writes the sorted, indexed alignment `output_dir/lambda_phage_test_sample_mCtoT.bam`. (Dots in the input name are replaced by underscores. Use `-O CRAM` for CRAM output.) If a required program is missing or any alignment step fails, `astair2 align` exits with a non-zero status and does not leave a partial output file. Existing output files are never overwritten.

If the reference FASTA file contains spaces in the header, alignment and calling will proceed using only the first word in the description unless the parameters `--add_underscores` and `--use_underscores` (aligner only) are used.

## 2. Call methylation

Once your fastq files are aligned and sorted (done automatically by `astair2 align`), you can run `astair2 call` to create a list of putative modified positions:

```bash
astair2 call -i output_dir/lambda_phage_test_sample_mCtoT.bam -f lambda_phage.fa --context CpG --minimum_base_quality 13 -d output_dir/
```
You can skip positions from the 5' or 3' of the reads if they seem to show Mbias by the `--start_clip` and `--end_clip` options. NB: In case positions are poorly covered or are covered only by reads' start and end positions, the usage of `--start_clip` and `--end_clip` can can alter the modified/unmodified ratio or show the genomic position as uncovered.
The `--no_information` option is also critical and we recommend to use 0, followed by a filtering of positions by the sum of modified and unmodified positions to be greater than 0.

## 3. Interpret results
After calling methylation, you will find two additional files in `output_dir`:

1. `lambda_phage_test_sample_mCtoT_mCtoT_CpG.stats`
2. `lambda_phage_test_sample_mCtoT_mCtoT_CpG.mods`

The `.stats` file contains global statistics on the modification rate in different sequence contexts. You can use this to get an idea of the overall level of modification in your sample. Here you will find information about how many cytosine positions of certain context are in the reference, how many of them were covered, and how many reads were modified/unmodified at the covered positions on the relevant strand assuming directionality. In our example here, we used a 1:1 mixture of in-vitro modified and unmodified lambda phage, so the results show a methylation rate of approximately 50%. (The exact numbers depend slightly on the `bwa` version used for alignment; the values below were obtained with `bwa` 0.7.19 and `samtools` 1.21.)

| CONTEXT | SPECIFIC_CONTEXT | MEAN_MODIFICATION_RATE_PERCENT | TOTAL_POSITIONS | COVERED_POSITIONS | MODIFIED | UNMODIFIED |
| ------- | ---------------- | ------------------------------ | --------------- | ----------------- | -------- | ---------- |
| CpG     | *                | **47.356**                     | 6225            | 6172              | 353133   | 392570     |
| *       | CGA              | 44.159                         | 1210            | 1203              | 64527    | 81596      |
| *       | CGC              | 47.74                          | 1730            | 1716              | 97016    | 106202     |
| *       | CGG              | 47.811                         | 1847            | 1826              | 106299   | 116033     |
| *       | CGT              | 49.009                         | 1438            | 1427              | 85291    | 88739      |

The `.mods` file contains per-position information on your sample:

| CHROM  | START | END | MOD_LEVEL | MOD | UNMOD | REF | ALT | SPECIFIC_CONTEXT | CONTEXT | SNV | TOTAL_DEPTH |
| ------ | ----- | --- | --------- | --- | ----- | --- | --- | ---------------- | ------- | --- | ----------- |
| lambda | 3     | 4   | 1.0       | 23  | 0     | C   | T   | CGG              | CpG     | No  | 57          |
| lambda | 4     | 5   | 0.0       | 0   | 41    | G   | A   | CGC              | CpG     | No  | 71          |
| lambda | 6     | 7   | 1.0       | 46  | 0     | C   | T   | CGA              | CpG     | No  | 104         |
| lambda | 7     | 8   | 1.0       | 72  | 0     | G   | A   | CGC              | CpG     | No  | 127         |
| lambda | 12    | 13  | 1.0       | 101 | 0     | C   | T   | CGC              | CpG     | No  | 240         |
| lambda | 13    | 14  | 0.0       | 0   | 145   | G   | A   | CGA              | CpG     | No  | 250         |

The header should be mostly self-explanatory. `MOD` and `UNMOD` refer to the number of reads covering that base that showed evidence of modification/no modification, and were of the right orientation to be meaningful for modification calling. The total coverage, including reads that were oriented in a way that no modification information can be extracted, is shown in `TOTAL_DEPTH`. `SNV` indicates whether the position may be a genetic C-to-T variant in the genome of the sample rather than a modified base:

- `No`: no evidence of a variant.
- `homozygous`: at least 80% of the informative reads from the opposite strand, which are not affected by the conversion, show the alternative base, indicating a likely homozygous C-to-T (or G-to-A) variant.
- `WGS_known`: the position is listed in the VCF file of known variants given with `--known_snp`.

Positions flagged as `homozygous` or `WGS_known` are excluded from the totals in the `.stats` file. Heterozygous variants and other variant types are not detected.

# Other useful information

## Recommendations for data pre-processing

1. Do quality control of the sequencing reads and do quality trimming before mapping and dispose of very short reads, using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or similar tools.
2. In most cases, it will be best to remove PCR duplicates before running the modification caller, unless your reads are non-randomly fragmented (e.g. enzymatically digested).
3. Check the fragment (insert) size distribution and decide on an overlap removal method for paired-end reads. The simplest option is the default removal of overlaps handled by `astair2 call`, which will randomly select one of two overlapping reads. This behaviour can be disabled by the `-sc` option, in case you are using a more sophisticated overlap-clipping tool.
4.  For speed and convenience we recommend using the `--per_chromosome` option, if possible, in order to run multiple processes in parallel. This also reduces the memory requirement when asTair2 is run on a desktop machine.

## Analysis of WGBS data (or other unmodified cytosine to thymine conversion methods)

The analysis pipeline for bisulfite sequencing data follows the same steps as TAPS data analysis, with `--method CtoT`. To avoid Bismark-style double alignments, `astair2 align --method CtoT` uses [`bwa-meth`](https://github.com/brentp/bwa-meth) (see [Installation](#0-installation)).

A small lambda phage WGBS example is included in this repository under `tests/test_data/`. From the repository root:

```bash
mkdir -p output_dir
astair2 align -f data/lambda_phage.fa -1 tests/test_data/small_real_wgbs_lambda_1.fq.gz -2 tests/test_data/small_real_wgbs_lambda_2.fq.gz --method CtoT -d output_dir/
```

You can now use `astair2 call` with `--method CtoT` for the modification calling:
```bash
astair2 call -i output_dir/small_real_wgbs_lambda_CtoT.bam -f data/lambda_phage.fa --method CtoT --context CpG --minimum_base_quality 13 -d output_dir/
```

All other sub-commands (`filter`, `mbias`, `idbias`, `summarise`, `simulate`) accept `--method CtoT` in the same way.
# Running the tests

With the virtual environment from the [installation](#0-installation) step activated:

```bash
python -m pip install -e ".[plot,test]"
python -m pytest
```

# Local quality checks

These are the checks run by CI. With the virtual environment activated:

```bash
python -m pip install -e ".[dev]"
black --check astair2 tests
mypy astair2
```

To run them automatically before each commit, together with a [gitleaks](https://github.com/gitleaks/gitleaks) scan for accidentally committed secrets:

```bash
pre-commit install
pre-commit run --all-files
```

pre-commit downloads and builds gitleaks itself the first time it runs, so no separate installation is needed. Commit from a terminal where the virtual environment is active, so that the hooks find `black` and `mypy`.

# Citation

If you use asTair2 in your work, please cite:

> Velikova, G. V. (2020). *Development and application of computational methods to study DNA modifications* [DPhil thesis]. University of Oxford. https://doi.org/10.5287/ora-x5zmqda6d

```bibtex
@phdthesis{velikova2020astair,
  author = {Velikova, Gergana V.},
  title  = {Development and application of computational methods to study {DNA} modifications},
  school = {University of Oxford},
  year   = {2020},
  type   = {{DPhil} thesis},
  doi    = {10.5287/ora-x5zmqda6d},
  url    = {https://ora.ox.ac.uk/objects/uuid:9c7e449a-f0c1-439d-82c6-1b64e52d2dd6}
}
```

# Authors and history

`asTair2` is developed and maintained by Gergana V. Velikova.

It is based on asTair, which was written by **Gergana V. Velikova** and **Benjamin Schuster-Boeckler** and released under the GNU General Public License v3. asTair2 is distributed under the same license.

asTair2 is not affiliated with the [rastair](https://github.com/bsbludwig/rastair) project (see [Related tools](#related-tools)).

# License

This software is made available under the terms of the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

