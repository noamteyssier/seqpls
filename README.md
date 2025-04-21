# seqpls

This is a paired FASTQ grep tool for sequence analysis.

It's built using the same matching algorithm as [`bqtools grep`](https://github.com/arcinstitute/bqtools), and was developed to showcase the difference between BINSEQ and FASTQ formats.

It accepts FASTQ files (compressed or uncompressed) and lets you match on regular expressions or fixed strings on either the R1, R2 or both.

## Installation

```bash
cargo install seqpls
```

## Usage

```bash
# See full help menu
seqpls --help

# Search for a fixed string in an unpaired FASTQ file
seqpls -e "ACGT" <some_fastq>

# Search for a string in the R1 and a regex in the R2
seqpls -T3 -e "ACGT" -R "[AC][TG][AC][TG]" <some_r1> <some_r2>

# Filter sequences without string in either using 3 threads
seqpls -T3 -v -F "ACGT" <some_r1> <some_r2>
```
