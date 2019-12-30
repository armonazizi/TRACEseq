# Armon Azizi (aazizi@stanford.edu)
#
# This script takes as input a merged fastq file derived from amplicon sequencing data.
# It aligns each read to a master amplicon and bins each read based on a number of heuristics.
# Each read is determined to be derived from a Homologous Recombination (HR) event, a Non-Homologous End Joining Event (NHEJ),
# an unmodified event. If the read does not align to the expected master sequence, it is binned as being non-HBB or contaminating DNA.
#
# The order of operations is as follows:
# 1. Find the expected cut site on the master read after alignment with amplicon read. If cut site is not found, the read is likely junk, or non-HBB.
# 2. If cut site is found, and if there are insertions or deletions flanking the cut site, classify read as NHEJ.
# 3. If read is not NHEJ, and the anchor bases (PAM associated bases changed after successful HR) are correctly modified, classify read as HR.
# 4. If read does not have changed anchor bases, and not NHEJ, classify as unmodified.
#
# Requires skbio.
#
# Usage: python align_and_filter_pre_tuba.py <input_fastq> <output_HR_fastq> <output_NHEJ_reads> <output_WT_reads> <non-HBB>



from skbio.alignment import StripedSmithWaterman
from skbio.alignment import local_pairwise_align_ssw
import skbio.sequence as seq
import sys

# Read 4 lines from fastq file and return as array
def get_read(fastq):
    result=[]
    for i in range(4):
        result.append(fastq.readline().strip())

    return result

# Write 4 lines to fastq file from array
def write_read(fastq, read):

    for line in read:
        fastq.write(line+'\n')

# Check next line without iterating
def peek_line(f):
    pos = f.tell()
    line = f.readline()
    f.seek(pos)
    return line


# Open all files for reading and writing
input_fastq=open(sys.argv[1], 'r')
hr_fastq=open(sys.argv[2], 'w')
nhej_fastq=open(sys.argv[3], 'w')
wt_fastq=open(sys.argv[4], 'w')
non_hbb=open(sys.argv[5], 'w')


# HBB Query sequence
master = "TCGATAGCAATTCGTTCACTAGCAACCTCAAACAGACACCATGGTNCANNTNACNCCNGANGANAANNNNGCAGTCACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA"
wildtype="TCGATCATGCTTAGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA"

# Define HBB cut site position on master sequence
cut_site_position=77

# Make smith_waterman_alignment_object
master_query = StripedSmithWaterman(master)
wildtype_query = StripedSmithWaterman(wildtype)

total_reads=0
hr_reads=0
nhej_reads=0
wt_reads=0

# Iterate over all reads and align
# Bin into HR, NHEJ, or WT Reads
while 1:
    # get next read
    read=get_read(input_fastq)

    # check for EOF
    if read[0]=="":
        break

    total_reads+=1
    seq=read[1]

    # Align read to master sequence
    master_alignment = master_query(seq)


    #### Check for NHEJ ####

    # Get aligned master and query sequences
    read_alignment = master_alignment.aligned_target_sequence
    ref_alignment = master_alignment.aligned_query_sequence

    # find cut site on query sequence
    cut_site_position=ref_alignment.find("CTGCCCTGTGGGG")

    # if cut site is not found, bin as non-HBB
    if cut_site_position==-1:
        write_read(non_hbb, read)
        continue

    # Check for insertions/deletions flanking cut site.
    read_region_of_interest=read_alignment[cut_site_position-20:cut_site_position+20]
    ref_region_of_interest=read_alignment[cut_site_position-20:cut_site_position+20]

    deletions=read_region_of_interest.count('-')
    insertions=ref_region_of_interest.count('-')

    # If insertions/deletions found, bin read as NHEJ
    if deletions > 0 or insertions > 0:
        write_read(nhej_fastq, read)
        nhej_reads+=1
        continue

    #### Check for HR and unmodified ####

    master_anchor_region="CAGTCACTGCCCTGTGGGG"

    # Find anchor base position on master read
    anchor_position=ref_alignment.find(master_anchor_region)

    # If anchor bases on query read are correctly modified, bin as hr event, otherwise, bin as WT/unmodified.
    if read_alignment[anchor_position+1]=='A' and read_alignment[anchor_position+4]=='C':
        write_read(hr_fastq, read)
        hr_reads+=1
    else:
        write_read(wt_fastq, read)
        wt_reads+=1

#print "file\ttotal_reads\twt_reads\tnhej_reads\thr_reads\twt_fraction\tnhej_fraction\thr_fraction"
#print sys.argv[1]+'\t'+str(total_reads)+'\t'+str(wt_reads)+'\t'+str(nhej_reads)+'\t'+str(hr_reads)+'\t'+str(float(wt_reads)/float(total_reads))+'\t'+str(float(nhej_reads)/float(total_reads))+'\t'+str(float(hr_reads)/float(total_reads))

# Close all output fastq files.
input_fastq.close()
wt_fastq.close()
nhej_fastq.close()
hr_fastq.close()
