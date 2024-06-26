import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pysam
from Bio import SeqIO
from tqdm import tqdm


def get_fasta(fasta_file = None):
    positions = []
    total_length = 0
    for rec in SeqIO.parse(fasta_file, "fasta"):
        start = 0
        total_length += len(rec.seq)
        positions.append([str(rec.id),start, len(rec.seq)])

    return [positions, total_length]

def filter_reads(bam, positions, fasta_length, filter_cutoff = 0.97, min_mapq = 2, write_data = None, write_bam = False, log=None, allow_dups = False):

    # read sets
    observed_read1s = set()
    observed_read2s = set()
    mapped_pairs = set()
    final_reads = set()
    

    # counters
    total_read_count = 0
    total_read_pairs = 0
    total_mapped_pairs = 0
    mapped_read_lengths = 0
    duplicate_positions = 0

    # storing data
    read_data = {}
    pair_mapqs = {}
    pair_mismatch = {}
    pair_inserts = {}
    read_positions = defaultdict(list)


    samfile = pysam.AlignmentFile(bam)

    #for printing out a new bam file
    if write_bam:
        print("Copying header for new bam...")
        samfile_out = pysam.AlignmentFile(write_bam + "_temp", "wb", template=samfile)
        reads_all = defaultdict(list)

    print("READING BAM: " + bam.split("/")[-1])
    print("Using reads with >" + str(filter_cutoff) + "% PID to consensus reference.")
    if log:
        log_file = open(log, 'a+')
        log_file.write("READING BAM: " + bam.split("/")[-1] + "\n")
        log_file.write("Using reads with >" + str(filter_cutoff) + "% PID to consensus reference." + "\n")


    ## STEP 1: collect paired reads and their information
    for gene in tqdm(positions, desc='Getting read pairs: '):
        for read in samfile.fetch(gene[0], gene[1], gene[2]):
            total_read_count += 1
 
            #store all reads if we're going to write them back to a new bam file
            if write_bam:
                reads_all[read.query_name].append(read)
   
            ## If we've seen this read's pair before
            if (read.is_read2 and read.query_name in observed_read1s) or (read.is_read1 and read.query_name in observed_read2s):

                #But if we haven't already seen this complete pair, then we can complete the pair and store the information
                #Also check that the pair is on the same scaffold
                if read.query_name not in mapped_pairs:
                    total_read_pairs += 1

                    if read.get_reference_positions() != []:
                        total_mapped_pairs += 1
                        mapped_pairs.add(read.query_name) #add to found 

                        read_positions[read.query_name].append(read.reference_name + ":" + str(read.get_reference_positions()[0]))

                        #for calculating mean read length
                        mapped_read_lengths += float(read_data[read.query_name]['len'])

                        #set mismatch percentage
                        pair_mismatch[read.query_name] = 1- ( ( float(read_data[read.query_name]['nm']) + float(read.get_tag('NM')) ) / ( float(read_data[read.query_name]['len']) + read.infer_query_length()) )
                        #set insert size
                        if read.get_reference_positions()[-1] > read_data[read.query_name]['start']:
                            pair_inserts[read.query_name] = read.get_reference_positions()[-1] - read_data[read.query_name]['start']
                        else:
                            pair_inserts[read.query_name] = read_data[read.query_name]['stop'] - read.get_reference_positions()[0]
                        #set mapq
                        pair_mapqs[read.query_name] = read.mapping_quality
                        if read_data[read.query_name]['mapq'] > read.mapping_quality:
                            pair_mapqs[read.query_name] = read_data[read.query_name]['mapq'] 

            #this is the first time we see a read from this pair and don't double count 
            elif (read.is_read1 and read.query_name not in observed_read1s) or (read.is_read2 and read.query_name not in observed_read2s):
                if read.get_reference_positions() != []: # don't use unmapped reads
                    if read.is_read1:
                        observed_read1s.add(read.query_name)
                    else:
                        observed_read2s.add(read.query_name)
                    #record the data for this read
                    read_data[read.query_name] = {"nm": read.get_tag('NM'), "len": read.infer_query_length(), "mapq": read.mapping_quality, "start": read.get_reference_positions()[0], 'stop': read.get_reference_positions()[-1], 'scaf': gene[0]}


    ## STEP 2: INSERT SIZE CUTOFF, MAPQ CUTOFF, AND MISMATCH CUTOFF
    mapped_read_lengths = mapped_read_lengths / total_mapped_pairs

    good_length = 0.0
    mapq_good = 0.0
    filter_cutoff_good = 0.0
    
    print("Filtering reads...")

    positions_observed = set()
    
    for read_pair in mapped_pairs:
        good_length += 2
        if pair_mapqs[read_pair] > min_mapq:
            mapq_good += 2

            # Which set does this read go into?
            if pair_mismatch[read_pair] > filter_cutoff:
                filter_cutoff_good += 2
                final_reads.add(read_pair)

                #write out to new bam file if option selected
                if write_bam:

                    found_pos = False
                    for pos in read_positions[read_pair]:
                        if not allow_dups and pos in positions_observed:
                            found_pos = True
                            duplicate_positions += 1
                        else:
                            positions_observed.add(pos)
                    
                    if not found_pos:
                        for read in reads_all[read_pair]:
                            samfile_out.write(read)
    
    print("**READ STATSTICS**")
    print("total reads found: " + str(total_read_count))
    print("average mapped read length: " + str(mapped_read_lengths))
    print("total fasta length: " + str(fasta_length))
    print("expected possible coverage: " + str(float(total_read_count)*mapped_read_lengths / fasta_length))
    print("total paired reads: " + str(total_read_pairs*2) + " (" + str(int(100*total_read_pairs*2.0 / total_read_count)) + "%)")
    print("total same scaffold mapped paired reads: " + str(total_mapped_pairs*2) + " (" + str(int(100*total_mapped_pairs*2.0 / total_read_count)) + "%)")
    print("")
    print("reads which also pass both pair insert size filters: " + str(good_length) + " (" + str(int(100*float(good_length) / total_read_count)) + "%)")
    print("reads which pass minimum mapq threshold of " + str(min_mapq) + ": " + str(mapq_good) + " (" + str(int(100*float(mapq_good) / total_read_count)) +  "%)")
    print("(final) reads which also pass read pair PID >" + str(filter_cutoff - duplicate_positions) + "%: " + str(filter_cutoff_good) + " (" + str(int(100*float(filter_cutoff_good) / total_read_count)) + "%)")
    print("Duplicate reads: " +str(duplicate_positions))
    print("(final) expected coverage: " + str(float(filter_cutoff_good - duplicate_positions) * mapped_read_lengths / fasta_length))
    if log:
        log_file.write("\n**READ STATSTICS**")
        log_file.write("\ntotal reads found: " + str(total_read_count))
        log_file.write("\naverage mapped read length: " + str(mapped_read_lengths))
        log_file.write("\ntotal fasta length: " + str(fasta_length))
        log_file.write("\nexpected possible coverage: " + str(float(total_read_count)*mapped_read_lengths / fasta_length))
        log_file.write("\ntotal paired reads: " + str(total_read_pairs*2) + " (" + str(int(100*total_read_pairs*2.0 / total_read_count)) + "%)")
        log_file.write("\ntotal same scaffold mapped paired reads: " + str(total_mapped_pairs*2) + " (" + str(int(100*total_mapped_pairs*2.0 / total_read_count)) + "%)")
        log_file.write("\n")
        log_file.write("\nreads which also pass both pair insert size filters: " + str(good_length) + " (" + str(int(100*float(good_length) / total_read_count)) + "%)")
        log_file.write("\nreads which pass minimum mapq threshold of " + str(min_mapq) + ": " + str(mapq_good) + " (" + str(int(100*float(mapq_good) / total_read_count)) +  "%)")
        log_file.write("\n(final) reads which also pass read pair PID >" + str(filter_cutoff) + "%: " + str(filter_cutoff_good) + " (" + str(int(100*float(filter_cutoff_good) / total_read_count)) + "%)")
        log_file.write("\nDuplicate reads: " +str(duplicate_positions))
        log_file.write("\n(final) expected coverage: " + str(float(filter_cutoff_good - duplicate_positions) * mapped_read_lengths / fasta_length))

    ## STEP 3: WRITE DATA IF NEEDED
    if write_data:
        f = open(write_data, 'w+')
        for read_pair in mapped_pairs:
            f.write(read_pair + "\t" + "\t" + str(pair_inserts[read_pair]) + "\t" + str(pair_mapqs[read_pair]) + "\t" + str(pair_mismatch[read_pair]) + "\n")
        f.close()
    ## STEP 4: WRITE NEW BAM IF NEEDED (TODO)

    samfile.close()
    if write_bam:
        samfile_out.close()
        print("sorting new bam")
        pysam.sort("-o", write_bam, write_bam + "_temp")
        os.system('rm ' + write_bam + "_temp")

    if log:
        log_file.close()
    return final_reads

        
            


if __name__ == '__main__':

    """ This is executed when run from the command line """

    parser = argparse.ArgumentParser(description="""Reports read statistics for a BAM mapping file. Do you have a SAM? Not to fear, just run these three commands:\n 
samtools view -S -b sample.sam > sample.bam\n
samtools sort sample.bam -o sample.sorted.bam\n
samtools index sample.sorted.bam\n in that order!""", formatter_class=argparse.RawTextHelpFormatter)


    # Required positional arguments
    parser.add_argument("bam", help="Sorted .bam file")
    parser.add_argument("fasta", help="Fasta file the bam is mapped to")


    # Optional arguments
    parser.add_argument("-m", "--mismatch_threshold", action="store", default=0.97, \
        help='Minimum percent identity of read pairs to consensus to use the reads - default is to run at 0.97.')

    parser.add_argument("-q", "--min_mapq", action="store", default=2, \
        help='Minimum mapq score of EITHER read in a pair to use that pair. Default: 2.')

    parser.add_argument("-w", "--write", action="store", default=None, \
        help='File name to write read statistics to.')

    parser.add_argument("-g", "--generate_bam", default=None, \
        help='Include to create a new filtered BAM to write to. (end your file name in .bam for BAM, .sam for SAM) ')

    parser.add_argument("-d", "--allow_dups", action="store_true", default=False, \
        help='Allow duplicate positions in the new bam file. Default: False means NO duplicate positions allowed by default.')
    
    parser.add_argument('--log', action='store', default=None, \
        help ="File to log results to.")




    # Parse
    args = parser.parse_args()
    positions = get_fasta(args.fasta)

    filter_reads(args.bam, positions[0], positions[1], float(args.mismatch_threshold),  int(args.min_mapq), write_data = args.write, write_bam=args.generate_bam, log=args.log, allow_dups=args.allow_dups)
