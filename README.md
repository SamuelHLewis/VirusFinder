# VirusFinder
## Purpose
A pipeline to identify viruses in paired-end RNA-Seq data
## Requirements
Written in python3.

Requires:

Bowtie2

Trinity v2.2.0 or later

TransDecoder v3.0.0 or later

DIAMOND

# Usage
Basic usage is:
```bash
python3 VirusFinder.py -l leftreads.fq -r rightreads.fq -g hostgenome -v virusdatabase -n nrdatabase
```
VirusFinder accepts the following additional arguments (all of which have defaults already set):

-c (number of cores)

-m (minimum length of viral contigs)
