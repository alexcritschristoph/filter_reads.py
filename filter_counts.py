import os
import sys
import glob
import pysam
import logging
import argparse
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict

def filter_paired_reads(bam, genes, min_ani = 0.95, min_mapq = 2):
    '''
    Filter reads from a .bam file
    Returns:
        pair2info - dictionary of read pair -> (mismatches, insert distance, mapq score, combined length)
    '''
    samfile = pysam.AlignmentFile(bam)

    
    total = mapped_pairs = mapq_good = mm_good = 0

    for scaff in genes:
        read_count = 0
        for read in samfile.fetch(scaff[0], scaff[1], scaff[2]):
            total += 1
            if read.has_tag('NM'):
                readmm = float(read.get_tag('NM')) #number of mismatches in pair
                read_length = float(read.infer_query_length()) #total length of pair
                read_ani =  1 - (readmm / read_length) #pair %ANI to reference
                # Final filter
                if read.mapping_quality >= min_mapq:
                    mapq_good += 1
                    if readmm >= min_ani:
                        read_count += 1
                        mm_good += 1

        if read_count > 0:
            print(bam.split("/")[-1].split(".")[0] + "\t" + scaff[-1] + "\t" + str(read_count))

    samfile.close()
    


def get_genes(genes_file):
    genes = []
    for record in SeqIO.parse(genes_file, 'fasta'):
        scaf = "_".join(record.id.split("_")[:-1])
        start = int(record.description.split("#")[1].strip())
        stop = int(record.description.split("#")[2].strip())
        genes.append([scaf, start, stop, record.id])
    return genes

if __name__ == '__main__':
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description="""Reports read statistics for a BAM mapping file, can also produce a filtered BAM.\n
        Input should be a BAM and a FASTA file of some or all sequences in that BAM.\t
        Output includes read statistics and an optional filtered BAM file.""", formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument("genes_file", help=".hits genes file")
    # Optional arguments
    parser.add_argument("-m", "--mismatch_threshold", action="store", default=0.95, \
        help='Minimum percent identity of read pairs to consensus to use the reads - default is to run at 0.95.')
    parser.add_argument("-q", "--min_mapq", action="store", default=2, \
        help='Minimum required mapq score of EITHER read in a pair to use that pair. Default: 2.')
    parser.add_argument('--log', action='store', default='filter.log', \
        help ="File to log results to.")

    args = parser.parse_args()
#    setup_logger(args.log)
    
    import glob
    for fn in glob.glob('./*.bam'):

        genes = get_genes(args.genes_file)
        filter_paired_reads(fn, genes, float(args.mismatch_threshold), int(args.min_mapq))

