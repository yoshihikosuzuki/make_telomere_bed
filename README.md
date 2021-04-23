# make_telomere_bed

**Input**:

* Single `.fasta` file of contigs/scaffolds
* Unit sequence of telomeres (cf. http://telomerase.asu.edu/sequences_telomere.html)

**Output**: 

* Single `.bed` file each line of which is as follows:

```txt
<contig_name> <start_position> <end_position> <copy_number>
```

* (TRF's original output file for the fasta file and filtered TRF output file containing only the telomere-related lines)

**Dependencies**:

* [Tandem Repeat Finder](https://github.com/Benson-Genomics-Lab/TRF)
* Python 3

## How to install

```bash
$ git clone https://github.com/yoshihikosuzuki/make_telomere_bed
$ cd make_telomere_bed
$ python setup.py install
```

After this, a command `make_telomere_bed` should be available.

## How to run

```txt
usage: make_telomere_bed [-h] [-t TRF_PATH] [-s] [-I]
                        contig_fasta unit_sequence

positional arguments:
  contig_fasta          A fasta file name.
  unit_sequence         Unit sqeuence like `TTAGGG`. Both capital and small
                        are OK.

optional arguments:
  -h, --help            show this help message and exit
  -t TRF_PATH, --trf_path TRF_PATH
                        Path to TRF's executable. [trf]
  -s, --split_contigs   Split contigs in `contig_fasta` into 1 Mbp
                        subsequences. TRF sometimes freezes for long contigs.
                        Use this option if TRF takes forever. [False]
  -I, --ignore_exist    Ignore existing files and generate them again (except
                        the final .bed file). [False]
```

**NOTES**:

* Results are written to files, not stdin.
* This program does not run TRF if the result file already exists. Therefore, you can try different telomere unit sequences without running TRF every time. If you wish to ignore the current TRF result and re-generate it, then use the `-I` option.
* This program can be actually used for any tandem repeats other than telomeres. What this program does is to extract outputs of TRF whose consensus unit sequence is exactly same as the given unit sequence while permitting reverse complement sequence and cyclic alignment. Therefore, this program is suitable for tandem repeats with short units (e.g. <10 bp) but not for those with long units (e.g. >50 bp) because longer units are less likely to have exact matches with other units in the same tandem repeat array.
