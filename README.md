# asTair2

_`asTair2` is a toolchain to process DNA modification sequencing data. It was designed primarily for handling [TET-Assisted Pyridine Borane (TAPS) sequencing](https://www.nature.com/articles/s41587-019-0041-2) output, but also contains functions that are useful for Bisulfite Sequencing (BS) data._

`asTair2` is an updated and modernised continuation of [asTair](#authors-and-history). It runs on current Python (3.12+) and current versions of its dependencies, while keeping asTair's sub-commands, options and output formats. Further development, including refactoring beyond the original code base, happens here.

# Basic usage
## 0. Installation 

`asTair2` requires Python 3.12 or newer. From a clone of this repository, install it with `pip`:

```bash
pip install .
```

To also enable the plotting options (`--plot`), install the optional matplotlib dependency:

```bash
pip install ".[plot]"
```

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

All of the examples in the main part of the current tutorial are based on the assumption that the input sequencing data are __TAPS__  pair-end sequencing reads, however, asTair2 analyses can be run in single-end mode (`--se`).  Also, asTair2 enables you to run analyses on __WGBS__ data, which requires a running installation of  [`bwa-meth`](https://github.com/brentp/bwa-meth) for the alignment step. For more information on WGBS analyses you may check the section [Analysis of WGBS data (or other unmodified cytosine to thymine conversion methods)](#analysis-of-wgbs-data-or-other-unmodified-cytosine-to-thymine-conversion-methods).

## 1. Align reads

We will assume that you have generated paired-end sequencing data, which is stored in two fastq files. For this brief tutorial, we assume the files are called `lambda.phage_test_sample_R1.fq.gz` and `lambda.phage_test_sample_R2.fq.gz`. If you want to follow this tutorial, you can download the files here:

```bash
# Or use curl -O if wget is not available
wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_1.fq.gz
wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_2.fq.gz
```

The raw reads need to be aligned. asTair2 contains a command to help with this. It assumes that [`bwa`](https://github.com/lh3/bwa) and [`samtools`](http://www.htslib.org/) are available on your system. (If you prefer to use a different aligner, [skip to step 2](#2-call-methylation).)

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
If the reference FASTA file contains spaces in the header, algnment and calling will proceed using only the first word in the description unless the parameters '--add_undescores' and '--use_underscores' (aligner only) are used.

## 2. Call methylation

Once your fastq files are aligned and sorted (done automatically by `astair2 align`), you can run `astair2 call` to create a list of putative modified positions:

```bash
astair2 call -i output_dir/lambda.phage_test_sample_mCtoT.cram -f lambda_phage.fa --context CpG --minimum_base_quality 13 -d output_dir/
```
You can skip positions from the 5' or 3' of the reads if they seem to show Mbias by the `--start_clip` and `--end_clip` options. NB: In case positions are poorly covered or are covered only by reads' start and end positions, the usage of `--start_clip` and `--end_clip` can can alter the modified/unmodified ratio or show the genomic position as uncovered.
The `--no_information` option is also critical and we recommend to use 0, followed by a filtering of positions by the sum of modified and unmodified positions to be greater than 0.

## 3. Interpret results
After calling methylation, you will find two additional files in `output_dir`:

1. `lambda.phage_test_sample_mCtoT_mCtoT_CpG.stats`
2. `lambda.phage_test_sample_mCtoT_mCtoT_CpG.mods`

The `.stats` file contains global statistics on the modification rate in different sequence contexts. You can use this to get an idea of the overall level of modification in your sample. Here you will find information about how many cytosine positions of certain context are in the reference, how many of them were covered, and how many reads were modified/unmodified at the covered positions on the relevant strand assuming directionality. In our example here, we used a 1:1 mixture of in-vitro modified and unmodified lambda phage, so the results show a methylation rate of approximately 50% :


| CONTEXT | SPECIFIC_CONTEXT | MEAN_MODIFICATION_RATE_PERCENT | TOTAL_POSITIONS | COVERED_POSITIONS | MODIFIED | UNMODIFIED |
| ------- | ---------------- | ------------------------------ | --------------- | ----------------- | ----------------- | ----------------- |
| CpG     |  *                | **48.225**                     | 6225            | 6225  | 356153           | 382377 |
|  *       | CGA              | 44.647                         | 1210            | 1210   | 64160          | 79545 |
|  *       | CGC              | 48.595                         | 1730            | 1730    | 97842           | 103499 |
|  *       | CGG              | 48.936                         | 1847            | 1847    | 108283            | 112991 |
|  *      | CGT              | 49.862                         | 1438            | 1438     | 85868            | 86342          |

The `.mods` file contains per-position information on your sample:

| CHROM | START | END   | MOD_LEVEL | MOD    | UNMOD   | REF   | ALT  | SPECIFIC_CONTEXT | CONTEXT  | SNV     | TOTAL_DEPTH |
| ----- | ----- | ----- | --------- | ------ | ------- | ----- | ---- | --------| ------- | ----------------- | ----------- |
| lambda |3 | 4 | 1.0 | 23 | 0 | C | T | CGG | CpG | No | 57 |
| lambda |4 | 5 | 0.0 | 0 | 34 | G | A | CGC | CpG | No | 71 |
| lambda  |6 | 7 | 1.0 | 38 | 0 | C | T | CGA | CpG | No | 104 |
| lambda  |7 | 8 | 1.0 | 58 | 0 | G | A | CGC | CpG | No | 127 |
| lambda | 12 | 13 | 1.0 | 88 | 0 | C | T | CGC | CpG | No | 240 |
| lambda  |13 | 14 | 0.0 | 0 | 139 | G | A | CGA | CpG | No | 250 |

The header should be mostly self-explanatory. `MOD` and `UNMOD` refer to the number of reads covering that base that showed evidence of modification/no modification, and were of the right orientation to be meaningful for modification calling. The total coverage, including reads that were oriented in a way that no modification information can be extracted, is shown in `TOTAL_DEPTH`. `SNV` gives a heuristic indication whether the position is indeed a modified base, or a genetic C to T variant in the genome of the sample.

# Other useful information

## Recommendations for data pre-processing

1. Do quality control of the sequencing reads and do quality trimming before mapping and dispose of very short reads, using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or similar tools.
2. In most cases, it will be best to remove PCR duplicates before running the modification caller, unless your reads are non-randomly fragmented (e.g. enzymatically digested).
3. Check the fragment (insert) size distribution and decide on an overlap removal method for paired-end reads. The simplest option is the default removal of overlaps handled by `astair2 call`, which will randomly select one of two overlapping reads. This behaviour can be disabled by the `-sc` option, in case you are using a more sophisticated overlap-clipping tool.
4.  For speed and convenience we recommend using the `--per_chromosome` option, if possible, in order to run multiple processes in parallel. This also reduces the memory requirement when asTair2 is run on a desktop machine.

## Analysis of WGBS data (or other unmodified cytosine to thymine conversion methods)

The analysis pipeline for bisulfite sequencing data does follows the same steps as TAPS data analysis, but requires different options. We again start from fastq files. To avoid Bismark-style double-alignments, we prefer to use `bwa meth`, which can be used directly through `astair2 align` when you choose the `--method CtoT` option.

```bash
mkdir -p output_dir
astair2 align -f lambda_phage.fa -1 lambda.phage_test_sample_BS_1.fastq.gz -2 lambda.phage_test_sample_BS_2.fastq.gz --method CtoT -d output_dir/
```

You can now use `astair2 call` with `--method CtoT` for the modifcation calling:
```bash
astair2 call -i output_dir/lambda.phage_test_sample_BS_CtoT.cram -f lambda_phage.fa --method CtoT --context CpG --minimum_base_quality 13 -d output_dir/
```
# Running the tests

```bash
pip install ".[plot,test]"
pytest
```

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

# License

This software is made available under the terms of the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

