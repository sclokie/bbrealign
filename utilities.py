import sys
import re
from collections import deque
from collections import defaultdict
import subprocess

def passes_conditions(column_5, column_6, nm, nm_threshold):
    h_or_s_count = column_6.count('H') + column_6.count('S')
    m_count = column_6.count('M')
    numbers = [int(num) for num in re.findall(r'\d+', column_6)]

    return column_5 > 50 and h_or_s_count == 1 and m_count == 1 and all(num > 5 for num in numbers) and nm <= nm_threshold

def filter_sam(reads, nm_threshold=1):
    stored_data = {}
    keys_queue = deque(maxlen=1000)
    filtered_reads = []

    for read in reads:
        sam_line = read.to_string()
        columns = sam_line.split('\t')

        if len(columns) >= 6:
            key = columns[0]
            column_5 = int(columns[4])
            column_6 = columns[5]

            # Check if the NM tag is present before attempting to access it
            if read.has_tag('NM'):
                nm = read.get_tag('NM')
            else:
                continue

            if key not in stored_data:
                stored_data[key] = []
                keys_queue.append(key)

                if len(keys_queue) == 1000:
                    oldest_key = keys_queue.popleft()
                    if oldest_key in stored_data:
                        del stored_data[oldest_key]

            stored_data[key].append(read)

            if passes_conditions(column_5, column_6, nm, nm_threshold):
                matching_reads = stored_data[key]
                for matching_read in matching_reads:
                    # Add a check for one M only
                    if matching_read.cigarstring.count('M') == 1:
                        filtered_reads.append(matching_read)

    return filtered_reads

import os
import pysam
from collections import deque



def filter_bam_file(input_bam, output_bam, window_size=100000):
    input_sam = input_bam.replace('.bam', '.sam')
    output_sam = output_bam.replace('.bam', '.sam')

    # Convert BAM to SAM
    pysam.view('-h', '-o', input_sam, input_bam, catch_stdout=False)

    # Read input SAM file and create output SAM file
    with pysam.AlignmentFile(input_sam, 'r') as infile, pysam.AlignmentFile(output_sam, 'wh', template=infile) as outfile:
        stored_reads = {}
        keys_queue = deque(maxlen=window_size)

        for read in infile:
            read_name = read.query_name

            if read_name not in stored_reads:
                stored_reads[read_name] = []
                keys_queue.append(read_name)

                if len(keys_queue) == window_size:
                    oldest_key = keys_queue.popleft()
                    if oldest_key in stored_reads:
                        del stored_reads[oldest_key]

            stored_reads[read_name].append(read)

            if read.cigarstring and 'D' in read.cigarstring:
                matching_reads = stored_reads[read_name]
                for matching_read in matching_reads:
                    outfile.write(matching_read)

    # Convert filtered SAM to BAM
    pysam.view('-bS', '-o', output_bam, output_sam, catch_stdout=False)

    # Remove temporary SAM files
    os.remove(input_sam)
    os.remove(output_sam)



def split_bam_by_deletion_length(input_bam, output_bam_gt1000, output_bam_lt1000):
    # Read input BAM file and create output BAM files
    # need to add in coverage filter too:

    # @transform(["*.bbmap.cigar.filtered.bam"], suffix(".bbmap.cigar.filtered.bam"), ".bbmap.roi.bed")
    # def calculate_depth(infile, outfile):
    #     # create genome file using a custom function that summarises a given fasta on the fly -
    #     # that is the same one used by bbmap in this case
    #     summarise_fasta(config_dict['bbmap_genome'], 'genome_sizes.txt')
    #     cmd = f"bedtools genomecov -bga -ibam {infile} \
    #     | awk '$4 > 5' \
    #     | bedtools merge -i - \
    #     | bedtools slop -i - -g genome_sizes.txt -b 300 > {outfile}"
    #     os.system(cmd)
    #
    with pysam.AlignmentFile(input_bam, 'rb') as infile, \
         pysam.AlignmentFile(output_bam_gt1000, 'wb', template=infile) as outfile_gt1000, \
         pysam.AlignmentFile(output_bam_lt1000, 'wb', template=infile) as outfile_lt1000:

        for read in infile:
            if read.cigarstring and 'D' in read.cigarstring:
                deletions = [int(x) for x in re.findall(r'(\d+)D', read.cigarstring)]
                max_deletion = max(deletions)

                if max_deletion >= 1000:
                    outfile_gt1000.write(read)
                else:
                    outfile_lt1000.write(read)
    # Filter by depth
    for output_bam in [output_bam_gt1000, output_bam_lt1000]:
        print("PROCESSING:",output_bam)
        outfile = output_bam.replace('.bam', '.roi.bed')
        print(outfile)
        cmd = f"bedtools genomecov -bga -ibam {output_bam} \
                | awk '$4 > 5' \
                | bedtools merge -i - \
                | bedtools slop -i - -g genome_sizes.txt -b 300 \
                | bedtools intersect -wa -a {output_bam} -b - > {outfile}"

        os.system(cmd)
        os.system(f"mv {outfile} {output_bam}")


def merge_and_count_deletions(in_bam,out_bed):
    import pandas as pd
    from pybedtools import BedTool
    #convert bam to bed
    command = f"""bedtools bamtobed -cigar -tag NM -i {in_bam} > {in_bam}.bamtobed.tmp"""
    subprocess.run(command, shell=True, check=True)

    # Load the BED file into a pandas DataFrame
    df = pd.read_table(f'{in_bam}.bamtobed.tmp', header=None)


    # Define column names based on the format
    df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'custom']
    print(df)
    # Group by the 'score' column and apply operation to each group
    grouped = df.groupby('score')
    print ('grouped:',grouped)

    # Open the output file
    with open('output.bed', 'w') as f:

        # For each group (i.e., rows with the same 'score')
        for score, group in grouped:
            print(score,'-->',group)

            # Skip if group is empty
            if group.empty:
                continue

            # Create a temporary BED file
            group[['chrom', 'start', 'end', 'name', 'score', 'strand', 'custom']].to_csv('temp.bed', sep='\t',
                                                                                         header=False, index=False)

            # Create a BedTool object from the temporary BED file
            bt = BedTool('temp.bed')

            # Merge the features in the BedTool object
            merged = bt.merge(c='4,5,6,7', o='collapse,collapse,collapse,collapse')

            # Write the merged features to the output file
            for interval in merged:
                fields = str(interval).strip().split('\t')
                names = fields[3].split(',')
                counts = len(names)
                if len(fields) > 6:
                    num_operations = len(re.findall('[0-9]+[MIDNSHP=X]', fields[6]))
                    f.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{','.join(names)}\t{counts}\t{num_operations}\n")
                else:
                    f.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{','.join(names)}\t{counts}\n")

            # Delete temp.bed after processing each group
            os.remove('temp.bed')

from Bio import SeqIO
def summarise_fasta(fasta_file, output_file):
    with open(fasta_file, 'r') as f, open(output_file, 'w') as out:
        for record in SeqIO.parse(f, 'fasta'):
            sequence_name = record.id.split()[0]  # Extract the first part of the header
            sequence_length = len(record.seq)
            out.write(f"{sequence_name}\t{sequence_length}\n")




def update_cigar_distribution(cigar_distribution, cigar_string):
    if cigar_string in cigar_distribution:
        cigar_distribution[cigar_string] += 1
    else:
        cigar_distribution[cigar_string] = 1

    return cigar_distribution

#A way to anotate a bam!
# Could be handy for OncoVir work

# import pysam
#     from pybedtools import BedTool
#
#     # Load the BAM file
#     bamfile = pysam.AlignmentFile(f'{infile}', 'rb')
#
#     # Load the BED file
#     bedfile = BedTool('hg19.genes.bed')
#
#     # Create a new BAM file with the gene name tag ('GN')
#     with pysam.AlignmentFile(f'{outfile}', 'wb', header=bamfile.header) as outputfile:
#         for read in bamfile:
#             # Get the read's region in BED format
#             read_bed = BedTool('{}\t{}\t{}'.format(read.reference_name, read.reference_start, read.reference_end),
#                                from_string=True)
#
#             # Intersect the read's region with the BED file
#             intersection = bedfile.intersect(read_bed)
#
#             # If there's an intersection, add a tag with the gene name
#             for gene in intersection:
#                 read.set_tag('GN', gene.name)
#
#             # Write the read to the new BAM file
#             outputfile.write(read)
