# indelConverter

## Prerequisite
* samtools > 1.3
* Python  > 3.7
* corresponding referenece genome fasta

## Installation
``` shell
git clone https://github.com/shanghungshih/indelConverter.git
```

## Acceptable format
- `txt`: tabular format (ex. database file from Annovar, note: it allows contig name without `chr`, ex. `1 10144 10145 TA T` or `chr1 10144 10145 TA T`)
- `vcf`: vcf format

## Parameters
* `--in_file` - input file
* `--out_file` - output file
* `--in_reference` - input corresponding reference fasta
* `--type` - input file type
* `--to_dash` - if true, convert indel in input file to format with '-' (ex. `chr1 10144 10145 TA T`  to  `chr1 10144 10145 A -`), else on the contrary
* `cmd_samtools` local samtools binary path

## Quick start
### Convert dash format to non-dash format for `txt`
```
python3 indelConverter.py --in_file data/dash.txt --in_reference /path/to/reference/ucsc.hg19.fasta --out_file data/out_dash.txt --type txt --cmd_samtools samtools
```

### Convert non-dash format to dash format for `vcf` (example file will be added in future, please use your own vcf file)
```
python3 indelConverter.py --in_file data/dash.vcf --in_reference /path/to/reference/ucsc.hg19.fasta --out_file data/out_dash.vcf --type vcf --to_dash
```
