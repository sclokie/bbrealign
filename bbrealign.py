
import argparse
from glob import glob
import sys
import os
import os.path
import re
import shutil
from argparse import RawTextHelpFormatter
from ruffus import *
import time

from utilities import filter_sam, update_cigar_distribution, filter_bam_file
from utilities import summarise_fasta, split_bam_by_deletion_length, merge_and_count_deletions

parser = argparse.ArgumentParser(description='\
Process bwa2 aligned data to extract interesting reads that could be realigned using bbmap.',
                                 epilog="Notes:\n\
The script assumes bams are already split by contig.", \
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-c', '--config', action="store", dest="config_file", \
                    help=' REQUIRED. A file containing paths of all programs to use,\n CPUs to use etc. This is system specific. [Required]', \
                    required=True)

parser.add_argument('-n', '--named-directory', action="store", dest="named_directory", \
                    help=' Required. A name to use for run summary. [Required]', \
                    required=True)

parser.add_argument('-chemistry', '--chemistry', action="store", dest="chemistry", \
                    help=' REQUIRED. The capture chemistry used. [Required]', \
                    required=True)

parser.add_argument('-debug', '--debug', action="store", dest="debug", \
                    help=' NOT REQUIRED. If set to Y, then large intermediate files will not be zeroed, to save disk space', \
                    required=False)

args = parser.parse_args()

# Set global dictionaries:
samplesheetDict = {}
sampleDict = {}
config_dict = {}
#config_dict['cpu_count'] = 12

# Will want to read pipeline version in from command?
pipeline_version = '0.1.0'

# Need a way to get the bam files in.
# This script is intended to be used as part of a typical bioinformatics pipeline - accepting a bam file
# and extracting and realigning interesting reads from it. The bam file has to be
# split by chromosome to function.


# Requirements:
# config files of genomic resources - need to get from command line now.
# bam file input: get from regexsimply state: "just docker run <bbrealign> -args <config_file> in this dir and it'll work"


def readConfigFile(config_file):
    """
    Read in the config file that lists the locations of genomic resources such as genome.fasta,
    and indexes for bwa and computational resources available on the host machine
    """
    with open(config_file, 'r') as infile:
        for line in infile:
            if line.startswith('#') or line.startswith('\n'):
                pass
            else:
                line = line.split('=', 1)
                config_dict[line[0]] = line[1][:-1]


readConfigFile(args.config_file)
cur_dir = os.getcwd()

# @jobs_limit(12)
# @transform(["/data/*.split.bam"],suffix(".split.bam"),".splitreads.bam")
# def extract_split_reads(infile,outfile):
#    print(infile,'-->',outfile)
#    name = re.sub(".split.bam","",infile)
#    os.system(f"/app/sambamba-0.8.2-linux-amd64-static index -t 4 {infile}")
#    # Run python script to extract split reads with little filtering:
#    # Edit distance of 1 or more
#    os.system(f"python realignbam.py \
#                -i {infile} \
#                -o {outfile} \
#                -t {os.getcwd()}")
#
# @follows(extract_split_reads)
# @subdivide(["/data/*splitreads.bam"],formatter(r".splitreads.bam$"),"/data/{basename[0]}.fastq.gz")
# def splitreads_to_fastq(infile,outfile):
#     print(infile,'-->',outfile)
#     cmd = f"samtools fastq -o {outfile} -0 /dev/null {infile}"
#     os.system(cmd)
#
# splitreadsregex = "(?P<FILESTRING>.+).splitreads.R[12].fastq.gz$"
# @jobs_limit(5)
# @follows(splitreads_to_fastq)
# @transform("/data/*.splitreads.fastq.gz", suffix(".splitreads.fastq.gz"),".bbmap.bam")
# def bbmap_split_reads(infile,outfile):
#     name = re.sub(".bbmap.bam","",outfile)
#     os.system(f"bbmap.sh \
#            -Xmx24g \
#            ref={config_dict['bbmap_genome']} \
#            nodisk=t \
#            rgid={name} \
#            rgpl=illumina \
#            rgsm={name} \
#            in={infile} \
#            out={outfile} \
#            outu={name}.bbmap.unmapped.fastq \
#            unpigz=t \
#            threads=6 \
#            maxindel=1000k \
#            rescuedist=1000000 \
#            usejni=t \
#            rgid={name} \
#            mdtag=t \
#            nhtag=t \
#            xmtag=t \
#            nmtag=t \
#            stoptag=t \
#            lengthtag=t \
#            idtag=t \
#            inserttag=t \
#            scoretag=t \
#            usemodulo=t \
#            statsfile={name}.bbmap.align.stats.log \
#            machineout=t")

#@jobs_limit(4)
#@follows(bbmap_split_reads)
# @transform(["/data/*.bbmap.bam"], suffix(".bbmap.bam"), ".bbmap.sorted.bam")
# def sort_bbmap_bam(infile, outfile):
#      print(infile, '-->', outfile)
#      """
#      sort the bam using sambamba, to result in a coordinate sorted bam, suitable for GATK
#      This is CONSIDERABLY faster than letting picard FixMateInformation do it and is more stable in reagards to system resource.
#      """
#      #temp_dir=os.getcwd()
#      temp_dir=os.getcwd()
#      os.system(f"./sambamba-0.8.2-linux-amd64-static sort -m 4G -t {config_dict['sambamba_sort_threads']} --tmpdir {temp_dir} -o {outfile} {infile} &> {outfile}.sambamba-coordinate-sort.log")
#      time.sleep(2)
#      if config_dict['debugmode'] == 'T':
#          pass
#      else:
#          os.system(f"> {infile}")

# @follows(sort_bbmap_bam)
@collate("/data/*.bbmap.sorted.bam", formatter("([^/]+).chr([0-9]|[0-9][0-9]|X|Y|MT).bbmap.sorted.bam$"),"{path[0]}/{1[0]}.merged.bbmap.bam")
def merge_bbmap(infiles,outfile):
    """
    Sort the chromosomes in proper order so increase speed.
    i.e. the MergeSamFiles will then not have to sort the file.
    Note sambamba failed here with a segmentation fault or a 'Read reference ID is out of range' error. Picard works.
    """
    print(infiles,'-->',outfile,'\n')
    log_name   = outfile[:-10]
    m_x_y_bams = list(infiles[-3:])           # slice the M, X and Y
    cwd = os.getcwd()
    m_x_y_bams_sorted = [m_x_y_bams[1],m_x_y_bams[2],m_x_y_bams[0]] # get into order: X,Y,M
    split_bams = list(infiles)                # convert ruffus tuple to list
    split_bams = (split_bams[:-3])            # remove X and Y for sorting
    split_bams.sort(key=lambda x: int(x.split('chr')[1].split('.')[0])) # sort on the chr value within the file name
    sorted_bams = split_bams + m_x_y_bams_sorted  # Create the file order for vcf concatenation
    bam_files = ' INPUT='.join(sorted_bams)
    cmd = f"picard MergeSamFiles \
            SO=coordinate \
            INPUT={bam_files} \
            OUTPUT={outfile} \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR={cwd} \
            CREATE_INDEX=true 2>{log_name}.merge.log"

    os.system(cmd)


@follows(merge_bbmap)
@transform(["/data/*.merged.bbmap.bam"],suffix(".merged.bbmap.bam"),".merged.bbmap.bam.bai")
def index_merged_bbmap(infile,outfile):
    print(infile,'-->',outfile)
    os.system(f"./sambamba-0.8.2-linux-amd64-static index -t 12 {infile}")
@follows(index_merged_bbmap)
@transform(["/data/*.merged.bbmap.bam"],suffix(".merged.bbmap.bam"),".merged.bbmap.cigar.filtered.bam")
def filter_bbmap_on_cigar(infile,outfile):
    outfile_temp = f"{outfile}.temp"
    filter_bam_file(infile,outfile_temp)
    os.system(f"samtools sort -o {outfile} {outfile}.temp")
    os.system(f"samtools index {outfile}")
    os.system(f"rm {outfile_temp}")
#
# @follows(filter_bbmap_on_cigar)
# @transform(["*.bbmap.cigar.filtered.bam"],suffix(".bbmap.cigar.filtered.bam"),".bbmap.roi.bed")
# def calculate_depth(infile,outfile):
#     #create genome file using a custom function that summarises a given fasta on the fly -
#     # that is the same one used by bbmap in this case
#     summarise_fasta(config_dict['bbmap_genome'],'genome_sizes.txt')
#     cmd = f"bedtools genomecov -bga -ibam {infile} \
#     | awk '$4 > 5' \
#     | bedtools merge -i - \
#     | bedtools slop -i - -g genome_sizes.txt -b 300 > {outfile}"
#     os.system(cmd)
#
# @follows(calculate_depth)
# @transform(["*.bbmap.roi.bed"], suffix(".bbmap.roi.bed"), ".bbmap.roi.bam")
# def filter_bbmap_bam_for_roi(infile, outfile):
#     # filter bam
#     bam_name = re.sub(".bbmap.roi.bed",".bbmap.cigar.filtered.bam",infile)
#     filter_cmd = f"intersectBed -wa -a {bam_name} -b {infile} > {outfile}"
#     os.system(filter_cmd)
#
# @follows(filter_bbmap_bam_for_roi)
# @transform(["*.bbmap.roi.bam"],suffix(".bbmap.roi.bam"),".bbmap.roi.bam.bai")
# def index_filtered_bbmap(infile,outfile):
#     print(infile,'-->',outfile)
#     os.system(cmd.sambamba_index(infile,6))
#
# @follows(index_filtered_bbmap)
# @transform(["*.bbmap.roi.bam"],suffix(".bbmap.roi.bam"),".bbmap_realigned.roi.bed")
# def calculate_depth_for_bbmap_realigned(infile,outfile):
#     #create genome file using a custom function that summarises a given fasta on the fly -
#     # that is the same one used by bbmap in this case
#     summarise_fasta(config_dict['bbmap_genome'],'genome_sizes.txt')
#     cmd = f"bedtools genomecov -bga -ibam {infile} \
#     | awk '$4 > 5' \
#     | bedtools merge -i - \
#     | bedtools slop -i - -g genome_sizes.txt -b 300 > {outfile}"
#     os.system(cmd)
#
# @follows(calculate_depth_for_bbmap_realigned)
# @transform(["*.bbmap.roi.bam"],suffix(".bbmap.roi.bam"),".annotated.bed")
# def create_annotated_bed(infile,outfile):
#     print(infile,'-->',outfile)
#     # download a list of genes, remove header and merge
#     # os.system(f"mysql \
#     #             --user=genome \
#     #             --host=genome-mysql.soe.ucsc.edu \
#     #             -A -e \"select chrom, txStart, txEnd, name2 from hg19.refGene\" \
#     #             | sed \'1d\' \
#     #             | bedtools sort -i - \
#     #             | bedtools merge -c 4 -o distinct | sed 's/chr//' > hg19.genes.bed")
#     #merge_and_count_deletions(infile,outfile)
#     # command = f"""bedtools bamtobed -cigar -tag NM -i {infile} \
#     #                | bedtools merge -c 5,5,7 -o first,count,first \
#     #                | intersectBed -loj -a stdin -b hg19.genes.bed \
#     #                | awk '{{OFS="\\t"}} {{printf "%s\\t%s\\t%s\\t%.0f\\t%s\\t%s\\t%s\\n", $1, $2, $3, $4, $5, $10, $6}}' \
#     #                > {outfile}"""
#
#     command = f"""bedtools bamtobed -cigar -tag NM -i {infile} \
#                   | intersectBed -loj -a stdin -b hg19.genes.bed \
#                   | awk '{{OFS="\\t"}} {{printf "%s\\t%s\\t%s\\t%.0f\\t%s\\t%s\\t%s\\n", $1, $2, $3, $5, $7, $10, $11}}' \
#                   > {outfile}"""
#
#     subprocess.run(command, shell=True, check=True)
#
# @follows(create_annotated_bed)
# @transform(["*.annotated.bed"],suffix(".annotated.bed"),".annotated.summary.bed")
# def create_summary_bed(infile,outfile):
#     print(infile,'-->',outfile)
#     import re
#     from collections import defaultdict
#
#     bed_dict = defaultdict(list)
#
#     with open(f'{infile}', 'r') as f:
#         for line in f:
#             fields = line.strip().split("\t")  # strip() removes newline characters at the end of lines
#             if len(fields) < 7:  # Ignore lines with fewer than 7 columns
#                 continue
#
#             # If column 7 is '.', replace it with a short name
#             if fields[6] == '.':
#                 fields[6] = f"{fields[0]}_{fields[1][-3:]}_{fields[2][-3:]}"
#
#             # Count the number of operations in column 5
#             num_operations = len(re.findall('[0-9]+[MIDNSHP=X]', fields[4]))
#
#             if num_operations > 4:  # Ignore rows with more than 4 operations
#                 continue
#
#             key = (fields[3], fields[6])  # Creates a tuple from the 4th and 7th column
#             bed_dict[key].append(fields)
#
#     # Now bed_dict contains all the rows grouped by the 4th and 7th column
#     # Let's create a summary
#     summary = {}
#
#     for key, values in bed_dict.items():
#         score = int(key[0])
#         if score <= 40:  # Ignore rows with score <= 20
#             continue
#         chromosome = values[0][0]
#         start_positions = [int(v[1]) for v in values]
#         end_positions = [int(v[2]) for v in values]
#         min_start = min(start_positions)
#         max_end = max(end_positions)
#         count = len(values)
#         if count >= 5:  # Only include rows with count >= 3
#             summary[key] = (chromosome, min_start, max_end, count)  # Store chromosome, min start, max end, count
#
#     # Now, 'summary' is a dictionary with keys as (column 4, column 7) tuples and values as summarized events
#     # Let's write it to a new BED file
#     with open(f'{outfile}', 'w') as f:
#         for key, value in summary.items():
#             f.write(f"{value[0]}\t{value[1]}\t{value[2]}\t{key[0]}|{key[1]}|{value[3]}\n")
#
#
# @follows(create_summary_bed)
# @transform(["*.annotated.summary.bed"],suffix(".annotated.summary.bed"),".annotated.summary.txt")
# def summarise_deletions(infile, outfile):
#     import pandas as pd
#
#     # Load the bed file into a pandas DataFrame
#     df = pd.read_csv(f'{infile}', sep='\t', header=None)
#
#     # Split the 4th column into separate columns
#     df[['score', 'name', 'count']] = df[3].str.split('|', expand=True)
#
#     # Convert 'score' and 'count' columns to int for proper sorting
#     df[['score', 'count']] = df[['score', 'count']].apply(pd.to_numeric)
#
#     # Drop the original 4th column
#     df = df.drop(3, axis=1)
#
#     # Rename the columns for clarity
#     df.columns = ['chromosome', 'start', 'end', 'score', 'name', 'count']
#
#     # Sort by 'score' then 'count'
#     df = df.sort_values(['score', 'count'], ascending=False)
#
#     # Save the sorted dataframe to a new bed file
#     df.to_csv(f'{outfile}', sep='\t', index=False, header=False)
#
#
# @follows(summarise_deletions)
# @transform(["*annotated.bed"], suffix("annotated.bed"), "distribution.png")
# def plot_distributions(infile,outfile):
#
#     # plot the data
#     import matplotlib.pyplot as plt
#
#     # Read the BED file and extract deletion sizes
#     deletion_sizes = []
#     with open(f'{infile}', 'r') as bedfile:
#         for line in bedfile:
#             fields = line.strip().split('\t')
#             deletion_size = float(fields[3].split('|')[0])
#             deletion_sizes.append(deletion_size)
#
#     # Create the histogram
#     plt.hist(deletion_sizes, bins=30, edgecolor='black')
#     plt.title(f'{infile}')
#     plt.xlabel('Deletion Size')
#     plt.ylabel('Frequency')
#
#     # Set y-axis to log scale
#     plt.yscale('log')
#
#     # Save the figure as PNG
#     plt.savefig(f'{outfile}', format='png')

pipeline_run(multiprocess = int(config_dict['cpu_count']),verbose=1)