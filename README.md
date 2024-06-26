# filter_reads.py
Filters BAM files created by bowtie2 for better read mapping. Use for genomes from metagenomes.
```
usage: filter_reads.py [-h] [-m MISMATCH_THRESHOLD] [-q MIN_MAPQ]
                       [-l MAX_INSERT_LENGTH] [-u MIN_INSERT_LENGTH]
                       [-w WRITE] [-g] [--log LOG]
                       bam fasta

Reports read statistics for a BAM mapping file. Do you have a SAM? Not to fear, just run these three commands:
 
samtools view -S -b sample.sam > sample.bam

samtools sort sample.bam -o sample.sorted.bam

samtools index sample.sorted.bam
 in that order!

positional arguments:
  bam                   Sorted .bam file
  fasta                 Fasta file the bam is mapped to

optional arguments:
  -h, --help            show this help message and exit
  -m MISMATCH_THRESHOLD, --mismatch_threshold MISMATCH_THRESHOLD
                        Minimum percent identity of read pairs to consensus to use the reads - default is to run at 0.97.
  -q MIN_MAPQ, --min_mapq MIN_MAPQ
                        Minimum mapq score of EITHER read in a pair to use that pair. Default: 2.
  -l MAX_INSERT_LENGTH, --max_insert_length MAX_INSERT_LENGTH
                        Maximum insert size between two reads - default is to use 3x median insert size.
  -u MIN_INSERT_LENGTH, --min_insert_length MIN_INSERT_LENGTH
                        Minimum insert size between two reads - default is to use 50 bp.
  -w WRITE, --write WRITE
                        File name to write read statistics to.
  -g, --generate_sam    Include to create a new filtered SAM or BAM to write to (will be a sam if ends in .sam, bam if ends in .bam).
  --log LOG             File to log results to.

```

### New version: Filter reads where pairs can map to different scaffolds.

By default, `filter_reads.py` will remove read pairs that map to different scaffolds. However, in some cases (e.g., HiC), 
read pairs can map to different scaffolds legitimately. 

`filter_reads_all_scaffolds.py` has this filter (and insert length filters) removed. 

Usage:
```
usage: filter_reads_all_scaffolds.py [-h] [-m MISMATCH_THRESHOLD] [-q MIN_MAPQ] [-w WRITE] [-g GENERATE_BAM] [--log LOG]
                                     bam fasta
```

Example:


```
python filter_reads_all_scaffolds.py example.bam example.fasta -g example_filtered.bam
```

By default, this will remove any reads that map to the same starting position as any other (only 1 read per start position allowed). 
(If you have two identical reads, this removes one of them)

to turn off this behavior, pass `--allow_dups`.


