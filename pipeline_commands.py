
import re
import os
import pandas as pd
import io
import pickle
import subprocess
#import allel

#update
# A collection of functions useful for NGS pipelines
# Implemented on VM1-CentOS7-binfz-wdc-1 on 01/09/2020
# Major changes are bwa-mem to bwa-mem2
# They expect already cleaned file names and assume as little as possible, for maximum reuse.

# Program paths:
Trimmomatic='/opt/programs/Trimmomatic-0.38/trimmomatic-0.38.jar'
Trimmomatic_adapters='/opt/programs/Trimmomatic-0.38/adapters/NexteraPE-PE.fa'

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
                line = line.split('=')
                config_dict[line[0]] = line[1][:-1]

def create_poly_list(poly_excel):
    """
    Create a poly list from an excel file on the S drive
    should remove spaces and replace with unserscore, also removing any trailing spaces.
    """
    from openpyxl import load_workbook
    wb = load_workbook(filename="%s" % poly_excel)
    sheet_ranges=wb['variants_for_pipeline']
    var_sheet=wb.get_sheet_by_name('variants_for_pipeline')
    rows= var_sheet.get_highest_row()
    variants_used_for_filtering = 0
    variants_not_checked = 0
    with open("poly_list_file.txt", 'w+') as poly_list_file:
        for i in range(2, rows+1):
            CHROM_cell='D'+str(i)
            POS_cell='E'+str(i)
            REF_cell='F'+str(i)
            ALT_cell='G'+str(i)
            filter_list_cell='BE'+str(i)
            check1_cell = 'BF'+str(i)
            check2_cell = 'BG'+str(i)
            
            #Only record variants that have been second checked - i.e. if either colums report None, skip them.
            check1_result = sheet_ranges[check1_cell].value
            check2_result = sheet_ranges[check2_cell].value
            if check1_result and check2_result != None:
                variants_used_for_filtering+=1
                variant="%s\t%s\t%s\t%s\t%s\t%s\n" % \
                ((i-1),sheet_ranges[CHROM_cell].value, sheet_ranges[POS_cell].value,sheet_ranges[REF_cell].value,sheet_ranges[ALT_cell].value,sheet_ranges[filter_list_cell].value)
                poly_list_file.write(variant)
            else:
                variants_not_checked+=1
                
    print("Number of variants that will be used for filtering: %s\nNumber of variants still to be checked by a clinical scientist: %s" \
          %(variants_used_for_filtering,variants_not_checked))    

def bcl2fastq(config_dict,illumina_folder):
    #Configure bacl2fastq qith many threads
    cmd = "{bcl2fastq} --no-lane-splitting --fastq-compression-level 1 -r 12 -p 64 -w 12 -R {illumina_folder}".format(bcl2fastq=config_dict['bcl2fastq'], illumina_folder=illumina_folder)
    return cmd
    
def trim_adapters_from_paried_end_reads(in_fwd_FASTQ,in_rev_FASTQ,out_fwd_FASTQ,out_fwd_unpaired,out_rev_FASTQ,out_rev_unpaired,log_name):
    """
    Remove bases with low quality using Trimmomatic
    Add in illumina adapter removal too
    ILLUMINACLIP:/opt/programs/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10
    """
    cmd = "java -Xmx12g -XX:ParallelGCThreads=4 -jar {Trimmomatic} PE -threads 4 -phred33 \
          {in_fwd_FASTQ} {in_rev_FASTQ} \
          {out_fwd_FASTQ} {out_fwd_unpaired} \
          {out_rev_FASTQ} {out_rev_unpaired} \
          ILLUMINACLIP:{adapters}:2:30:10 \
          SLIDINGWINDOW:4:25 \
          LEADING:3 \
          TRAILING:3 \
          MINLEN:25 2>{log}".format(
                            Trimmomatic      = Trimmomatic,
                            in_fwd_FASTQ     = in_fwd_FASTQ,
                            in_rev_FASTQ     = in_rev_FASTQ,
                            out_fwd_FASTQ    = out_fwd_FASTQ,
                            out_rev_FASTQ    = out_rev_FASTQ,
                            out_fwd_unpaired = out_fwd_unpaired,
                            out_rev_unpaired = out_rev_unpaired,
                            adapters         = Trimmomatic_adapters,
                            log              = log_name)
          
    return cmd
    
    
    
def fastp_trim_paired_end_reads(config_dict,in_fwd_FASTQ,in_rev_FASTQ,out_fwd_FASTQ,out_fwd_unpaired,out_rev_FASTQ,out_rev_unpaired,log_name):
    """
    Remove bases with low quality using fastp
    https://github.com/OpenGene/fastp/releases
    Removes illumina adapter by default
    
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 25
    
    replicates Trimmomatic sliding window function
    
    Consider Q20! We lost ~12 million reads using 25
    """
    cmd = "{fastp} \
          --thread {fastp_threads} \
          -i {in_fwd_FASTQ}  -I {in_rev_FASTQ} \
          -o {out_fwd_FASTQ} -O {out_rev_FASTQ} \
          --unpaired1 {out_fwd_unpaired} --unpaired2 {out_rev_unpaired} \
          --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
          --json {log}.json \
          --html {log}.html \
          2>{log}".format(
                            fastp            = config_dict['fastp'],
                            fastp_threads    = config_dict['fastp_threads'],
                            in_fwd_FASTQ     = in_fwd_FASTQ,
                            in_rev_FASTQ     = in_rev_FASTQ,
                            out_fwd_FASTQ    = out_fwd_FASTQ,
                            out_rev_FASTQ    = out_rev_FASTQ,
                            out_fwd_unpaired = out_fwd_unpaired,
                            out_rev_unpaired = out_rev_unpaired,
                            log              = log_name)
          
    return cmd

def fastp_trim_paired_end_reads_split(config_dict,in_fwd_FASTQ,in_rev_FASTQ,out_fwd_FASTQ,out_fwd_unpaired,out_rev_FASTQ,out_rev_unpaired,log_name):
    """
    Remove bases with low quality using fastp
    https://github.com/OpenGene/fastp/releases
    Removes illumina adapter by default
    
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 25
    
    replicates Trimmomatic sliding window function
    
    Consider Q20! We lost ~12 million reads using 25
    """
    cmd = "{fastp} \
          --thread 12 \
          -i {in_fwd_FASTQ}  -I {in_rev_FASTQ} \
          -o {out_fwd_FASTQ} -O {out_rev_FASTQ} \
          --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
          --split=48 \
          --split_prefix_digits=4 \
          --json {log}.json \
          --html {log}.html \
          2>{log}".format(
                            fastp            = config_dict['fastp'],
                            fastp_threads    = config_dict['fastp_threads'],
                            in_fwd_FASTQ     = in_fwd_FASTQ,
                            in_rev_FASTQ     = in_rev_FASTQ,
                            out_fwd_FASTQ    = out_fwd_FASTQ,
                            out_rev_FASTQ    = out_rev_FASTQ,
                            log              = log_name)
          
    return cmd


def align_FASTQ_mem2(config_dict,worksheet,illumina_folder,chemistry,trimmed_fwd_FASTQ,trimmed_rev_FASTQ,out_primary_bam,root_name):
    """
    Align quality filtered FASTQ files using bwa mem. Use samblaster to directly pipe to duplicate marker.
    SAMBLASTER marks duplicates in a single pass over a SAM file in which all alignments for the same read-id are grouped together...and ...outputs read-id grouped SAM
    To do: pipe into sambamba for sorting (multicore).

    Read group configuration: RG -
                                    ID = worksheet.flowcellID.Lane#             
                                    CN = centre (WMRGL)
                                    DS = description (TSO)
                                    DT = date and time (ISO)
                                    SM = Sample name (sliced infile name)
                                    LB = Library (sliced infile name)
                                    PL = Platform key word (Illumina)
                                    PN =
    ./bam polishBam  --in testFiles/sortedSam.sam --out results/updatedSam.sam \
    --log results/updated.log \
    --checkSQ \
    --fasta testFiles/testFasta.fa \
    --AS my37 \
    --UR testFasta.fa \
    --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none"

    """
    #readConfigFile(config_file)
    
    import datetime
    i = datetime.datetime.now()
    date_time_iso = "%s" % i.isoformat()
    
    #lane = re.search(r"(L00[1-8])",trimmed_fwd_FASTQ)
    lane = re.search(r"_([1-8])_R[12]_001",trimmed_fwd_FASTQ)
    if lane is not None:
        lane = list(filter(None,lane.groups()[0]))[0]
    else:
        lane = 'premerged'
    original_FASTQ = re.sub(r".qfilter.fastq.gz",".fastq.gz",trimmed_fwd_FASTQ)
    path           = os.path.realpath("%s" % original_FASTQ)
    sample_name     = out_primary_bam#.split("/")[-1][:-12]
    # example: 200622_NB551645_0039_AHNYJJBGXF/
    #try/except added 27-03-2021
    try:
        split_folder = illumina_folder.split("_")[3]
    except AttributeError:
        split_folder = None#'FC_UNKNOWN'
    if split_folder is not None:
        flowcell = split_folder
    else:
        flowcell = 'FC_UNKNOWN'
    
    #sed \'s/-R @RG.*//\' - | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    # Capture a cleaned up version of the command line to put in the bam header
    tidy_cmd_line = "{bwa2} mem -t {threads} -M -k 18 \
                     {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                     {samblaster} --removeDups -u {sample_name}.{lane}.unmapped.fastq 2> \
                     {sample_name}.{lane}.duplicate.Metrics.txt | \
                     {samtools} view -Sb -o {out_primary_bam} 2>{sample_name}.{lane}.samtoolsBamConvertErrors.log".format(\
                     bwa2              = config_dict['bwa2'],
                     threads           = config_dict['bwa2_threads'],
                     samtools          = config_dict['samtools'],
                     samblaster        = config_dict['samblaster'],
                     worksheet         = worksheet,
                     flowcell          = flowcell,
                     lane              = lane,
                     trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                     trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                     date              = date_time_iso,
                     sample_name       = sample_name,
                     chemistry         = chemistry,
                     genome            = config_dict['bwa2_genome_fasta'],
                     out_primary_bam   = out_primary_bam)
                          
     
    align_cmd = "{bwa2} mem -t {threads} -M -k 18 -R '@RG\\tID:{worksheet}.{flowcell}.{lane}\\tCN:WMRGL\\tDS:{chemistry}\\tDT:{date}\\tSM:{sample_name}\\tLB:{worksheet}\\tPL:ILLUMINA\\tCL:{cmd_line}' \
                {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                sed \'s/-R @RG.*//\' - | \
                {samblaster} --removeDups -u {sample_name}.{lane}.unmapped.fastq 2> \
                {sample_name}.{lane}.duplicate.Metrics.txt | \
                {samtools} view -Sb -o {out_primary_bam} 2>{sample_name}.{lane}.samtoolsBamConvertErrors.log".format(\
                    bwa2              = config_dict['bwa2'],
                    threads           = config_dict['bwa2_threads'],
                    samtools          = config_dict['samtools'],
                    samblaster        = config_dict['samblaster'],
                    worksheet         = worksheet,
                    flowcell          = flowcell,
                    lane              = lane,
                    trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                    trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                    date              = date_time_iso,
                    sample_name       = sample_name,
                    chemistry         = chemistry,
                    genome            = config_dict['bwa2_genome_fasta'],
                    out_primary_bam   = out_primary_bam,
                    cmd_line          = tidy_cmd_line)
                    
    return align_cmd
    
    
def exome_align_mem2_pipe_to_sort_single(config_dict,worksheet,illumina_folder,chemistry,trimmed_fwd_FASTQ,trimmed_rev_FASTQ,out_primary_bam,root_name):
    """
    Align quality filtered FASTQ files using bwa mem. Use samblaster to directly pipe to duplicate marker.
    SAMBLASTER marks duplicates in a single pass over a SAM file in which all alignments for the same read-id are grouped together...and ...outputs read-id grouped SAM
    To do: pipe into sambamba for sorting (multicore).

    Read group configuration: RG -
                                    ID = worksheet.flowcellID.Lane#             
                                    CN = centre (WMRGL)
                                    DS = description (TSO)
                                    DT = date and time (ISO)
                                    SM = Sample name (sliced infile name)
                                    LB = Library (sliced infile name)
                                    PL = Platform key word (Illumina)
                                    PN =
    ./bam polishBam  --in testFiles/sortedSam.sam --out results/updatedSam.sam \
    --log results/updated.log \
    --checkSQ \
    --fasta testFiles/testFasta.fa \
    --AS my37 \
    --UR testFasta.fa \
    --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none"

    """
    #readConfigFile(config_file)
    
    import datetime
    i = datetime.datetime.now()
    date_time_iso = "%s" % i.isoformat()
    
    lane = re.search(r"(L00[1-8])",trimmed_fwd_FASTQ)
    #lane = re.search(r"_([1-8])_R[12]_001",trimmed_fwd_FASTQ)
    
    if lane is not None:
        #lane = list(filter(None,lane.groups()[0]))
        lane = lane.groups()[0]
        #print ('*'*50,lane,'*'*50)
    else:
        lane = 'premerged'
    original_FASTQ = re.sub(r".qfilter.fastq.gz",".fastq.gz",trimmed_fwd_FASTQ)
    path           = os.path.realpath("%s" % original_FASTQ)
    sample_name     = out_primary_bam#.split("/")[-1][:-12]
    # example: 200622_NB551645_0039_AHNYJJBGXF/
    #try/except added 27-03-2021
    try:
        split_folder = illumina_folder.split("_")[3]
    except AttributeError:
        split_folder = None#'FC_UNKNOWN'
    if split_folder is not None:
        flowcell = split_folder
    else:
        flowcell = 'FC_UNKNOWN'
    
    #sed \'s/-R @RG.*//\' - | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    # Capture a cleaned up version of the command line to put in the bam header
    
    tidy_cmd_line = "{bwa2} mem -t {threads} -M -k 18 \
                     {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                     {samblaster} --removeDups -u {sample_name}.{lane}.unmapped.fastq 2> \
                     {sample_name}.{lane}.duplicate.Metrics.txt | \
                     {samtools} view -Sb - 2>{sample_name}.{lane}.samtoolsBamConvertErrors.log | \
                     {samtools} sort -m 1024M -@ 2 -T {sample_name}.temp -O bam - > \
                     {out_primary_bam}".format(\
                     bwa2              = config_dict['bwa2'],
                     threads           = config_dict['bwa2_threads'],
                     samtools          = config_dict['samtools'],
                     samblaster        = config_dict['samblaster'],
                     worksheet         = worksheet,
                     flowcell          = flowcell,
                     lane              = lane,
                     trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                     trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                     date              = date_time_iso,
                     sample_name       = sample_name,
                     chemistry         = chemistry,
                     genome            = config_dict['bwa2_genome_fasta'],
                     out_primary_bam   = out_primary_bam)
                          
     
    align_cmd = "{bwa2} mem -t {threads} -M -k 18 -R '@RG\\tID:{worksheet}.{flowcell}.{lane}\\tCN:WMRGL\\tDS:{chemistry}\\tDT:{date}\\tSM:{sample_name}\\tLB:{worksheet}\\tPL:ILLUMINA\\tCL:{cmd_line}' \
                {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                sed \'s/-R @RG.*//\' - | \
                {samblaster} --removeDups -u {sample_name}.{lane}.unmapped.fastq 2> \
                {sample_name}.{lane}.duplicate.Metrics.txt | \
                {samtools} view -Sb - 2>{sample_name}.{lane}.samtoolsBamConvertErrors.log | \
                {samtools} sort -m 1024M -@ 2 -T {sample_name}.temp -O bam - > \
                {out_primary_bam}".format(\
                    bwa2              = config_dict['bwa2'],
                    threads           = config_dict['bwa2_threads'],
                    samtools          = config_dict['samtools'],
                    samblaster        = config_dict['samblaster'],
                    worksheet         = worksheet,
                    flowcell          = flowcell,
                    lane              = lane,
                    trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                    trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                    date              = date_time_iso,
                    sample_name       = sample_name,
                    chemistry         = chemistry,
                    genome            = config_dict['bwa2_genome_fasta'],
                    out_primary_bam   = out_primary_bam,
                    cmd_line          = tidy_cmd_line)
                    
    return align_cmd
    
#####
def wgs_align_mem2_pipe_to_sambamba_sort_single(config_dict,worksheet,illumina_folder,chemistry,trimmed_fwd_FASTQ,trimmed_rev_FASTQ,out_primary_bam,root_name):
    """
    Align quality filtered FASTQ files using bwa mem. Use samblaster to directly pipe to duplicate marker.
    SAMBLASTER marks duplicates in a single pass over a SAM file in which all alignments for the same read-id are grouped together...and ...outputs read-id grouped SAM
    To do: pipe into sambamba for sorting (multicore).

    Read group configuration: RG -
                                    ID = worksheet.flowcellID.Lane#             
                                    CN = centre (WMRGL)
                                    DS = description (TSO)
                                    DT = date and time (ISO)
                                    SM = Sample name (sliced infile name)
                                    LB = Library (sliced infile name)
                                    PL = Platform key word (Illumina)
                                    PN =
    ./bam polishBam  --in testFiles/sortedSam.sam --out results/updatedSam.sam \
    --log results/updated.log \
    --checkSQ \
    --fasta testFiles/testFasta.fa \
    --AS my37 \
    --UR testFasta.fa \
    --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none"

    """
    #readConfigFile(config_file)
    
    import datetime
    i = datetime.datetime.now()
    date_time_iso = "%s" % i.isoformat()
    
    lane = re.search(r"(L00[1-8])",trimmed_fwd_FASTQ)
    #lane = re.search(r"_([1-8])_R[12]_001",trimmed_fwd_FASTQ)
    
    if lane is not None:
        #lane = list(filter(None,lane.groups()[0]))
        lane = lane.groups()[0]
        #print ('*'*50,lane,'*'*50)
    else:
        lane = 'premerged'
    original_FASTQ = re.sub(r".qfilter.fastq.gz",".fastq.gz",trimmed_fwd_FASTQ)
    path           = os.path.realpath("%s" % original_FASTQ)
    sample_name     = root_name#.split("/")[-1][:-12]
    # example: 200622_NB551645_0039_AHNYJJBGXF/
    #try/except added 27-03-2021
    try:
        split_folder = illumina_folder.split("_")[3]
    except AttributeError:
        split_folder = None#'FC_UNKNOWN'
    if split_folder is not None:
        flowcell = split_folder
    else:
        flowcell = 'FC_UNKNOWN'
    
    #sed \'s/-R @RG.*//\' - | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    # Capture a cleaned up version of the command line to put in the bam header
     #time /opt/programs/bwa-mem2-2.2.1_x64-linux-hot-fix/bwa-mem2 mem -t 48 -M -k 18 /data/genomic-resources/human/gatk_bundle/hs37d5/bwa-mem2/hs37d5.fa 109JH_R1.fastq.gz 109JH_R2    .fastq.gz | /opt/programs/sambamba-0.7.1-linux-static view -t 2 -S -f bam -l 0 /dev/stdin | /opt/programs/sambamba-0.7.1-linux-static sort -t 12 -m 12GB /dev/stdin --tmpdir=/data/wgs/temp -o 109JH-test.bam
    tidy_cmd_line = "{bwa2} mem -t {threads} -M -k 18 \
                     {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                     {samblaster} --removeDups -u {sample_name}.{lane}.unmapped.fastq 2> \
                     {sample_name}.{lane}.duplicate.Metrics.txt | \
                     {sambamba} view -l 0 -t 4 -S -f bam -l 0 /dev/stdin 2>{sample_name}.{lane}.sambambaConvertErrors.log | \
                     {sambamba} sort -t 4 -m 12GB /dev/stdin --tmpdir=/mnt/fire/scratch -o {out_primary_bam} ".format(\
                     bwa2              = config_dict['bwa2'],
                     threads           = config_dict['bwa2_threads'],
                     sambamba          = config_dict['sambamba'],
                     samblaster        = config_dict['samblaster'],
                     worksheet         = worksheet,
                     flowcell          = flowcell,
                     lane              = lane,
                     trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                     trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                     date              = date_time_iso,
                     sample_name       = sample_name,
                     chemistry         = chemistry,
                     genome            = config_dict['bwa2_genome_fasta'],
                     out_primary_bam   = out_primary_bam)
                          
     
    align_cmd = "{bwa2} mem -t {threads} -M -k 18 -R '@RG\\tID:{worksheet}.{flowcell}.{lane}\\tCN:WMRGL\\tDS:{chemistry}\\tDT:{date}\\tSM:{sample_name}\\tLB:{worksheet}\\tPL:ILLUMINA\\tCL:{cmd_line}' \
                {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                sed \'s/-R @RG.*//\' - | \
                {samblaster} --removeDups -u {sample_name}.{lane}.unmapped.fastq 2> \
                {sample_name}.{lane}.duplicate.Metrics.txt | \
                {sambamba} view -l 0 -t 4 -S -f bam -l 0 /dev/stdin 2>{sample_name}.{lane}.sambambaConvertErrors.log | \
                {sambamba} sort -t 4 -m 12GB /dev/stdin --tmpdir=/mnt/fire/scratch -o {out_primary_bam} ".format( \
                    bwa2              = config_dict['bwa2'],
                    threads           = config_dict['bwa2_threads'],
                    sambamba          = config_dict['sambamba'],
                    samblaster        = config_dict['samblaster'],
                    worksheet         = worksheet,
                    flowcell          = flowcell,
                    lane              = lane,
                    trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                    trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                    date              = date_time_iso,
                    sample_name       = sample_name,
                    chemistry         = chemistry,
                    genome            = config_dict['bwa2_genome_fasta'],
                    out_primary_bam   = out_primary_bam,
                    cmd_line          = tidy_cmd_line)
                    
    return align_cmd

#####


def extract_split_reads(config_dict,bam,bam_out,cwd,log_name):
    """
    extract useful reads for bbmap realignment
    """

    cmd = f"python {config_dict['extract_splits']} \
            -i {bam} \
            -o {bam_out} \
            -t {cwd} \
             2>{log_name}.extract_splits.log"

    return cmd

#####
    
def merge_bams(config_dict,bams,merged_bam,cwd,log_name):
    """
    Merge sorted bam files from each lane
    """
    bam_files = ' --INPUT '.join(bams)
    cmd = "java -Xmx4g -XX:ParallelGCThreads=1 -Djava.io.tmpdir=/tmp -jar /home/sclokie/programs/picard.jar MergeSamFiles \
    --SO coordinate \
    --INPUT {bams} \
    --OUTPUT {merged_bam} \
    --VALIDATION_STRINGENCY LENIENT \
    --TMP_DIR {cur_dir} \
    --CREATE_INDEX true 2>{name}.merge.log".format(
                              picard_mergesam = config_dict['picard_mergesam'],
                              bams       = ' --INPUT '.join(bams),
                              merged_bam      = merged_bam,
                              cur_dir         = cwd,
                              name            = log_name)
    return cmd 

def samtools_merge_bams(bams,merged_bam):
    """
    Merge sorted bam files from each lane
    """
    bam_files = ' '.join(bams)
    cmd = "samtools merge  {merged_bam} {bam_files}".format(merged_bam=merged_bam,bam_files=bam_files)
    
    return cmd 

def sambamba_merge(config_dict,threads,merged_bam,in_bams,log_name):
    """
    Index bam file using sambamba (fast!)
    """
    
    cmd = f"{config_dict['sambamba']} merge -l 9 -t {threads} {merged_bam} {' '.join(in_bams)} 2> {log_name}.sambamba_merge_split_bams.log"
    
    return cmd     


def samtools_index(bam):
    """
    Index bam file using samtools
    """
    cmd = "samtools index {bam}".format(bam=bam)
    
    return cmd 

def sambamba_index(bam,threads):
    """
    Index bam file using sambamba (fast!)
    default to 6
    """
    cmd = "sambamba-0.8.0-linux-amd64-static index -t {threads} {bam}".format(bam=bam,threads=threads)
    
    return cmd     

def sambamba_sort(config_dict,unsortedbam,sortedbam):
    """
    sort bam file using sambamba (fast!)
    """
    temp_dir=os.getcwd()
    cmd = f"/opt/programs/sambamba/sambamba-0.8.0-linux-amd64-static sort -m 4G -t {config_dict['sambamba_sort_threads']} --tmpdir {temp_dir} -o {sortedbam} {unsortedbam} &> {sortedbam}.sambamba-coordinate-sort.log"
    return cmd    
    
def sambamba_sort_picard(unsortedbam,sortedbam):
    """
    sort bam file using sambamba (fast!)
    I found it mfaster to give a very large number of threads to sambamba and then jobs_limit for Ruffus.
    ... large # Ruffus jobs and low threads here used a lot of system resource and much less user.
    """
    #temp_dir=os.getcwd()
    temp_dir='/wgs-scratch'
    cmd = f"/opt/programs/sambamba-0.7.1-linux-static sort --sort-picard -m 8G -t 4 --tmpdir {temp_dir} -o {sortedbam} {unsortedbam} &> {sortedbam}.sambamba-picard-sort.log"
    
    return cmd  
    
def mosdepth_coverage(bam,name):
    """
    Run mosdepth to calculate some quality coverage metrics
    /opt/programs/miniconda3/bin/conda config --add channels defaults
    /opt/programs/miniconda3/bin/conda config --add channels bioconda
    /opt/programs/miniconda3/bin/conda config --add channels conda-forge
    /opt/programs/miniconda3/bin/conda install mosdepth
    """
    cmd = "/opt/programs/miniconda3/bin/mosdepth \
           -t 3 \
           --by 100000 \
           --quantize 0:1:2 \
           {name} \
           {bam}".format(name=name,bam=bam)
    return cmd


def abra_realign_bam(config_dict,input_bam,realigned_bam,bwa_genome_fasta,abra_min_edge_pruning,abra_threads,cur_dir,name):
    """
    Realign sorted bam using ABRA - https://github.com/mozack/abra and PMID:24907369
    """
    cmd = ("java -Xmx4G -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -jar {abra} \
    --in {input_bam}\
    --out {realigned_bam}\
    --ref {bwa_genome_fasta} \
    --mer {abra_min_edge_pruning} \
    --target-kmers {kmers} \
    --threads 2 \
    --no-debug \
    --working {cur_dir} 2>{name}.stderr.abra.log 1>{name}.stout.abra.log".format(
                                                               abra                  = config_dict['abraRealignbam'],
                                                               input_bam             = input_bam,
                                                               realigned_bam         = realigned_bam,
                                                               bwa_genome_fasta      = bwa_genome_fasta,
                                                               abra_min_edge_pruning = abra_min_edge_pruning,
                                                               abra_threads          = int(abra_threads),
                                                               kmers                 = config_dict['abra_kmers'],
                                                               cur_dir               = cur_dir,
                                                               name                  = name))
    return cmd

def abra2_realign_bam(config_dict,input_bam,realigned_bam,bwa_genome_fasta,abra_min_edge_pruning,abra_threads,cur_dir,name):
    """
    Realign sorted bam using ABRA - https://github.com/mozack/abra and PMID:24907369
    """
    cmd = ("/opt/programs/amazon-corretto-8.265.01.1-linux-x64/bin/java -Xmx8G -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -jar {abra2} \
    --in {input_bam}\
    --out {realigned_bam}\
    --ref {bwa_genome_fasta} \
    --mer {abra_min_edge_pruning} \
    --target-kmers {kmers} \
    --threads {abra_threads} \
    --skip \
    --mcr 400000 \
    --tmpdir {cur_dir} 2>{name}.stderr.abra.log 1>{name}.stout.abra.log".format(
                                                               abra2                 = config_dict['abra2Realignbam'],
                                                               input_bam             = input_bam,
                                                               realigned_bam         = realigned_bam,
                                                               bwa_genome_fasta      = bwa_genome_fasta,
                                                               abra_min_edge_pruning = abra_min_edge_pruning,
                                                               kmers                 = config_dict['abra_kmers'],
                                                               abra_threads          = abra_threads,
                                                               cur_dir               = cur_dir,
                                                               name                  = name))
    return cmd


    
def wgs_abra_realign_split_bam(config_dict,input_bam,realigned_bam,bwa_genome_fasta,abra_min_edge_pruning,abra_threads,cur_dir,name):
    """
    Realign sorted bam using ABRA - https://github.com/mozack/abra and PMID:24907369
    Add the --nosort option, as it uses samtools == very slow!
    """
    cmd = ("java -Xmx8G -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2 -jar {abra} \
    --in {input_bam}\
    --out {realigned_bam}\
    --ref {bwa_genome_fasta} \
    --mer {abra_min_edge_pruning} \
    --nosort \
    --threads {abra_threads} \
    --tmpdir {cur_dir} 2>{name}.stderr.abra.log 1>{name}.stout.abra.log".format(
                                                               abra                  = config_dict['abra2Realignbam'],
                                                               input_bam             = input_bam,
                                                               realigned_bam         = realigned_bam,
                                                               bwa_genome_fasta      = bwa_genome_fasta,
                                                               abra_min_edge_pruning = abra_min_edge_pruning,
                                                               abra_threads          = abra_threads,
                                                               cur_dir               = cur_dir,
                                                               name                  = name))
    return cmd

def wgs_abra2_realign_bam(config_dict,input_bam,realigned_bam,bwa_genome_fasta,abra_min_edge_pruning,abra_threads,cur_dir,name):
    """
    Realign sorted bam using ABRA - https://github.com/mozack/abra and PMID:24907369
    --targets  /data/genomic-resources/human/selected_regions/NCBI_RefSeq+-50000.merged.bed \
    --target-kmers /data/genomic-resources/human/selected_regions/hg19-ensemble-31-07-2018+-5.merged.primaryAssemblyOnly.nochrs.kmers.bed \
    --target-kmers /data/genomic-resources/human/selected_regions/wgs_RefSeq_with_50k_padding_kmers.bed \
    --target-kmers /data/genomic-resources/human/selected_regions/wgs_RefSeq_with_1k_padding_kmers.bed \
    -Dsamjdk.intel_deflator_so_path=/opt/programs/abra2/abra2-2.23/lib/libIntelDeflater.so
    -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4
    """
    cmd = ("/home/sclokie/programs/amazon-corretto-8.292.10.1-linux-x64/bin/java -Xmx12G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -jar {abra} \
    --in {input_bam}\
    --out {realigned_bam}\
    --ref {bwa_genome_fasta} \
    --mer {abra_min_edge_pruning} \
    --target-kmers /mnt/fire/genomes/human/selected_regions/wgs_RefSeq_with_1k_padding_kmers.bed \
    --threads {abra_threads} \
    --nosort \
    -cl 1 \
    --gkl \
    --skip \
    --mcr 400000 \
    --tmpdir {cur_dir} 2>{name}.stderr.abra.log 1>{name}.stout.abra.log".format(
                                                               abra                  = config_dict['abra2Realignbam'],
                                                               input_bam             = input_bam,
                                                               realigned_bam         = realigned_bam,
                                                               bwa_genome_fasta      = bwa_genome_fasta,
                                                               abra_min_edge_pruning = abra_min_edge_pruning,
                                                               abra_threads          = abra_threads,
                                                               cur_dir               = cur_dir,
                                                               name                  = name))
    return cmd

    

def step1recalibrate_gatk4_BQSR_big_bam(config_dict,realigned_bam,BQSR_recal_data_table,log_name):
    """
    Recalibrate using GATK v4
    These files have the wrong chromosome
    --known-sites /data/process/genomes/human/gatk_bundle/dbsnp_137.hg19.vcf \
    --known-sites /data/process/genomes/human/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf \
    --known-sites /data/process/genomes/human/gatk_bundle/1000G_phase1.indels.hg19.vcf \
           
           #remove this part on 14-11-2019:
           --intervals {intervals} \
           intervals             = config_dict['ROI_ensembl_pad10'],
    Observation on 05/11/2020:
    Giving more than 8GB RAM used a lot more of the 'system' resources, rather than 'user'. This made it slower
    """
    #cmd = "{gatk4} --java-options \"-XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -Xmx48g\" BaseRecalibrator \
    cmd = "{gatk4} --java-options \"-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx20g\" BaseRecalibrator \
           -I {realigned_bam} \
           -R {reference_fasta} \
           --known-sites {dbsnp} \
           -L {targets} \
           -O {BQSR_recal_data_table} \
           --use-jdk-inflater TRUE \
           2>{log_name} ".format(
                                       gatk4                 = config_dict['gatk4'],
                                       dbsnp                 = config_dict['dbsnp'],
                                       targets               = config_dict['ROI_ensembl_pad5_nochrs'],
                                       realigned_bam         = realigned_bam,
                                       BQSR_recal_data_table = BQSR_recal_data_table,
                                       reference_fasta       = config_dict['bwa_genome_fasta'],
                                       log_name              = log_name)    
    return cmd


def step1recalibrate_gatk4_BQSR_small_bam(config_dict,realigned_bam,BQSR_recal_data_table,log_name):
    """
    Recalibrate using GATK v4
    These files have the wrong chromosome
    --known-sites /data/process/genomes/human/gatk_bundle/dbsnp_137.hg19.vcf \
    --known-sites /data/process/genomes/human/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf \
    --known-sites /data/process/genomes/human/gatk_bundle/1000G_phase1.indels.hg19.vcf \
           
           #remove this part on 14-11-2019:
           --intervals {intervals} \
           intervals             = config_dict['ROI_ensembl_pad10'],
    Observation on 05/11/2020:
    Giving more than 8GB RAM used a lot more of the 'system' resources, rather than 'user'. This made it slower
    
    -L {targets} \
    targets               = config_dict['ROI_ensembl_pad5_nochrs'],
    """
    #cmd = "{gatk4} --java-options \"-XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -Xmx48g\" BaseRecalibrator \
    cmd = "{gatk4} --java-options \"-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx6g\" BaseRecalibrator \
           -I {realigned_bam} \
           -R {reference_fasta} \
           --known-sites {dbsnp} \
           -O {BQSR_recal_data_table} \
           --use-jdk-inflater TRUE \
           2>{log_name} ".format(
                                       gatk4                 = config_dict['gatk4'],
                                       dbsnp                 = config_dict['dbsnp'],
                                       realigned_bam         = realigned_bam,
                                       BQSR_recal_data_table = BQSR_recal_data_table,
                                       reference_fasta       = config_dict['bwa_genome_fasta'],
                                       log_name              = log_name)    
    return cmd



    
def wgs_step1recalibrate_gatk4_BQSR(config_dict,realigned_bam,BQSR_recal_data_table,log_name):
    """
    Recalibrate using GATK v4
    These files have the wrong chromosome
    --known-sites /data/process/genomes/human/gatk_bundle/dbsnp_137.hg19.vcf \
    --known-sites /data/process/genomes/human/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf \
    --known-sites /data/process/genomes/human/gatk_bundle/1000G_phase1.indels.hg19.vcf \
           
    """

    cmd = "{gatk4} --java-options \"-XX:ParallelGCThreads=1 -XX:ConcGCThreads=2 -Xmx4g\" BaseRecalibrator \
           -I {realigned_bam} \
           -R {reference_fasta} \
           --known-sites /network/processed/HiSeq_data/Gnomad/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz \
           -O {BQSR_recal_data_table} \
           --use-jdk-inflater TRUE \
           2>{log_name} ".format(
                                       gatk4                 = config_dict['gatk4'],
                                       realigned_bam         = realigned_bam,
                                       BQSR_recal_data_table = BQSR_recal_data_table,
                                       reference_fasta       = config_dict['bwa_genome_fasta'],
                                       log_name              = log_name)    
    return cmd


def step2recalibrate_gatk4_BQSR(config_dict,recalibrated_bam,bam_file_ext,name,temp_dir):    
    """
    Apply the covariates determined during the previous step.
    
    --intervals {intervals} \
    intervals        = config_dict['ROI_ensembl_pad10'],
    """
    cmd = "{gatk4} --java-options \"-XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx6G\" ApplyBQSR \
            -R {reference_fasta} \
            -I {realigned_bam} \
            --emit-original-quals \
            --bqsr-recal-file {name}.split.BQSR_recal_data_table \
            -O {recalibrated_bam} \
            --create-output-bam-index \
            --use-jdk-inflater TRUE \
            --use-jdk-deflater TRUE \
            --tmp-dir {temp_dir} \
            2>{name}.BQSR_step1_ApplyBQSR.log".format(
                                  gatk4            = config_dict['gatk4'],
                                  realigned_bam    = name + bam_file_ext,
                                  recalibrated_bam = recalibrated_bam,
                                  reference_fasta  = config_dict['bwa_genome_fasta'],
                                  name             = name,
                                  temp_dir         = temp_dir)
    return cmd

def wgs_step2recalibrate_gatk4_BQSR(config_dict,recalibrated_bam,bbmap_bam,bam_file_ext,name,temp_dir):    
    """
    Apply the covariates determined during the previous step.
    """
    cmd = "{gatk4} --java-options \"-XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -Xmx4g\" ApplyBQSR \
            -R {reference_fasta} \
            -I {bbmap_bam} \
            --emit-original-quals \
            --bqsr-recal-file {name}.split.BQSR_recal_data_table \
            -O {recalibrated_bam} \
            --create-output-bam-index \
            --use-jdk-inflater TRUE \
            --use-jdk-deflater TRUE \
            --tmp-dir {temp_dir} \
            2>{name}.BQSR_step1_ApplyBQSR.log".format(
                                  gatk4            = config_dict['gatk4'],
                                  bbmap_bam        = bbmap_bam,
                                  recalibrated_bam = recalibrated_bam,
                                  reference_fasta  = config_dict['bwa_genome_fasta'],
                                  name             = name,
                                  temp_dir         = temp_dir)
    return cmd


     
def plot_covariates(config_dict,BQSR_recal_data_table,out_pdf,log_name):     
    """
    Plot both covariate tables in the same graph
    https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
    """
    cmd = " {gatk4} --java-options \"-XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -Xmx4g\" AnalyzeCovariates \
            -before {BQSR_recal_data_table} \
            -after  {BQSR_recal_data_table}_post\
            -plots {out_pdf} \
            --intermediate-csv-file {name}.csv \
            2>{log_name}".format(
                                     gatk4                 = config_dict['gatk4'],
                                     BQSR_recal_data_table = BQSR_recal_data_table,
                                     out_pdf               = out_pdf,
                                     log_name              = log_name,
                                     name                  = re.sub(".log","",log_name))

    return cmd 
     
     
def create_samtools_stats(config_dict,in_bam,out_stats):
    """
    Run samtools stats that generates a text file for processing with plot-bamstats
    :return: A text file
    """
    cmd = "{samtools} stats {in_bam} > {out_stats}".format(samtools=config_dict['samtools'],in_bam=in_bam,out_stats=out_stats)
    
    return cmd
    
def plot_samtools_stats(config_dict,stats,html):
    """
    Run samtools stats that generates a text file for processing with plot-bamstats
    :return: A text file
    """
    cmd = "{bamstats} {stats} -p {html}".format(bamstats=config_dict['plotbamstats'],stats = stats,html=html)
     
    return cmd    
    
def HaploTypeCaller_v4_single_vcf(config_dict,ROI,bam,output_file_name,output_file_name_root,ploidy,log_name):
    """
    GATKv4 now does read backed phasing by default ()
    Pertinant flags example:
    GT:AD:DP:GQ:PGT:PID:PL  0/1:1,2:3:64:0|1:108185301_C_G:81,0,64
    ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
    ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
    16-05-2019:
    remove the chr names from here:
    --dbsnp /data/process/genomes/human/gatk_bundle/dbsnp_137.hg19.excluding_sites_after_129_primary_assembly.nochr.vcf \
    """
    
    cmd = " {gatk} --java-options \"-XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx8G\" HaplotypeCaller \
             -R {bwa_genome_fasta} \
             --dbsnp {dbsnp} \
             --intervals {ROI} \
             --interval-padding 50 \
             --native-pair-hmm-threads 4 \
             -I {infile} \
             --output {outfile} \
             -bamout {output_file_name_root}.bamout.bam \
             --create-output-bam-md5 \
             --create-output-variant-md5 \
             --create-output-bam-index \
             --create-output-variant-index \
             -G StandardAnnotation \
             -G StandardHCAnnotation \
             --use-jdk-inflater TRUE \
             --use-jdk-deflater TRUE \
             --smith-waterman FASTEST_AVAILABLE \
             --sample-ploidy {ploidy} \
             2> {log_name}.HaploTypeCaller4.log\n \
             ".format( gatk                  = config_dict['gatk4'],
                       ROI                   = ROI,
                       dbsnp                 = config_dict['dbsnp'],
                       bwa_genome_fasta      = config_dict['bwa2_genome_fasta'],
                       infile                = bam,
                       outfile               = output_file_name,
                       output_file_name_root = output_file_name_root,
                       ploidy                = ploidy,
                       log_name              = log_name)
    
    return cmd


def wgs_HaploTypeCaller_v4_single_vcf(config_dict,ROI,bam,output_file_name,output_file_name_root,ploidy,log_name,temp_dir):
    """
    GATKv4 now does read backed phasing by default ()
    Pertinant flags example:
    GT:AD:DP:GQ:PGT:PID:PL  0/1:1,2:3:64:0|1:108185301_C_G:81,0,64
    ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
    ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
    16-05-2019:
    remove the chr names from here:
    --dbsnp /data/process/genomes/human/gatk_bundle/dbsnp_137.hg19.excluding_sites_after_129_primary_assembly.nochr.vcf \
    
    MAKE SURE TO SET the vm threads
    
    ulimit -u 8192
     
    --pair-hmm-implementation LOGLESS_CACHING \
     --pair-hmm-implementation LOGLESS_CACHING \
    """
    
    cmd = " {gatk} --java-options \"-XX:+UseParallelGC -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -Xmx10G\" HaplotypeCaller \
             -R {bwa_genome_fasta} \
             --dbsnp {dbsnp} \
             --native-pair-hmm-threads 4 \
             -I {infile} \
             --intervals {ROI} \
             --output {outfile} \
             -bamout {output_file_name_root}.bamout.bam \
             --create-output-bam-md5 \
             --create-output-variant-md5 \
             --create-output-bam-index \
             --create-output-variant-index \
             -G StandardAnnotation \
             -G StandardHCAnnotation \
             --use-jdk-inflater TRUE \
             --use-jdk-deflater TRUE \
             --smith-waterman FASTEST_AVAILABLE \
             --sample-ploidy {ploidy} \
             --tmp-dir {temp_dir} \
             2> {log_name}.HaploTypeCaller4.log\n \
             ".format( gatk                  = config_dict['gatk4'],
                       ROI                   = ROI,
                       dbsnp                 = config_dict['dbsnp'],
                       bwa_genome_fasta      = config_dict['bwa2_genome_fasta'],
                       infile                = bam,
                       outfile               = output_file_name,
                       output_file_name_root = output_file_name_root,
                       ploidy                = ploidy,
                       log_name              = log_name,
                       temp_dir              = temp_dir)
    
    return cmd


def wgs_HaploTypeCaller_v4_GVCF(config_dict,ROI,bam,output_file_name,output_file_name_root,ploidy,log_name,temp_dir):
    """
    GATKv4 now does read backed phasing by default ()
    Pertinant flags example:
    GT:AD:DP:GQ:PGT:PID:PL  0/1:1,2:3:64:0|1:108185301_C_G:81,0,64
    ##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
    ##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
    16-05-2019:
    remove the chr names from here:
    --dbsnp /data/process/genomes/human/gatk_bundle/dbsnp_137.hg19.excluding_sites_after_129_primary_assembly.nochr.vcf \
    
    MAKE SURE TO SET the vm threads
    
    ulimit -u 8192
     
    --pair-hmm-implementation LOGLESS_CACHING \
     --pair-hmm-implementation LOGLESS_CACHING \
    """
    
    cmd = " {gatk} --java-options \"-XX:+UseParallelGC -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -Xmx10G\" HaplotypeCaller \
             -R {bwa_genome_fasta} \
             --dbsnp {dbsnp} \
             --native-pair-hmm-threads 2 \
             -I {infile} \
             --intervals {ROI} \
             --output {outfile} \
             -bamout {output_file_name_root}.bamout.bam \
             --create-output-bam-md5 \
             --create-output-variant-md5 \
             --create-output-bam-index \
             --create-output-variant-index \
             -G StandardAnnotation \
             -G StandardHCAnnotation \
             -ERC GVCF \
             --use-jdk-inflater TRUE \
             --use-jdk-deflater TRUE \
             --smith-waterman FASTEST_AVAILABLE \
             --sample-ploidy {ploidy} \
             --tmp-dir {temp_dir} \
             2> {log_name}gvcf.HaploTypeCaller4.log\n \
             ".format( gatk                  = config_dict['gatk4'],
                       ROI                   = ROI,
                       dbsnp                 = config_dict['dbsnp'],
                       bwa_genome_fasta      = config_dict['bwa2_genome_fasta'],
                       infile                = bam,
                       outfile               = output_file_name,
                       output_file_name_root = output_file_name_root,
                       ploidy                = ploidy,
                       log_name              = log_name,
                       temp_dir              = temp_dir)
    
    return cmd

def wgs_add_gvcf_to_db(config_dict,ROI,gvcf,output_file_name_root,log_file,temp_dir,chr_specific_db):
    """
    For initial set up need to use the parameter:

    --genomicsdb-workspace-path
    --genomicsdb-update-workspace-path
    Once created, use:

    --genomicsdb-update-workspace-path

    """
    cmd = "{gatk} --java-options \"-Xms2G -Xmx2G -XX:ParallelGCThreads=2\" GenomicsDBImport  \
            --genomicsdb-update-workspace-path \
            {chr_specific_db} \
           -R /mnt/fire/genomes/human/bwa2/human_g1k_v37_decoy.fasta \
           --batch-size 1  \
           --tmp-dir /mnt/fire/scratch/ \
           --variant {gvcf} \
           --intervals {ROI} 2> {log_file} ".format(
                                                               gatk                  = config_dict['gatk4'],
                                                               ROI                   = ROI,
                                                               dbsnp                 = config_dict['dbsnp'],
                                                               bwa_genome_fasta      = config_dict['bwa2_genome_fasta'],
                                                               gvcf                  = gvcf,
                                                               outfile               = output_file_name_root,
                                                               output_file_name_root = output_file_name_root,
                                                               log_file              = log_file,
                                                               temp_dir              = temp_dir,
                                                               chr_specific_db       = chr_specific_db)

    return cmd

def genotype_wgs_gvcf(config_dict,ROI,chrom,out_vcf,output_file_name_root,log_file,temp_dir,chr_specific_db):
    """
    For initial set up need to use the parameter:

    --genomicsdb-workspace-path


   -G AlleleSpecificAnnotation \
    Once created, use:

    --genomicsdb-update-workspace-path
       -V gendb:///zpool/ngs/projects/nicky/wgs/batch2/wgs-db-chr{chrom}_gdb \
    """
    cmd = "{gatk} --java-options \"-Xms2G -Xmx2G -XX:ParallelGCThreads=2\" GenotypeGVCFs  \
           -V gendb://{chr_specific_db} \
           -R /mnt/fire/genomes/human/bwa2/human_g1k_v37_decoy.fasta \
           --tmp-dir /mnt/fire/scratch/ \
           --dbsnp {dbsnp} \
           --call-genotypes true \
           --intervals {ROI} \
           --annotate-with-num-discovered-alleles \
           -G AS_StandardAnnotation \
           -G GenotypeAnnotation \
           -G StandardAnnotation \
           -O {outfile} 2> {log_file} ".format(
                                                               gatk                  = config_dict['gatk4'],
                                                               ROI                   = ROI,
                                                               chrom                 = chrom,
                                                               dbsnp                 = config_dict['dbsnp'],
                                                               bwa_genome_fasta      = config_dict['bwa2_genome_fasta'],
                                                               outfile               = out_vcf,
                                                               output_file_name_root = output_file_name_root,
                                                               log_file              = log_file,
                                                               temp_dir              = temp_dir,
                                                               chr_specific_db       = chr_specific_db)
    return cmd


def VQSR_stage1_SNPs(config_dict,merged_vcf_in,root_name,out_file,output_file_name_root,log_file,temp_dir):
    """
    Run this on the genotyped merged vcf -
    """

    cmd = "{gatk} --java-options \"-Xms4G -Xmx4G -XX:ParallelGCThreads=2\" VariantRecalibrator \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 \
           -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           -R /mnt/fire/genomes/human/bwa2/human_g1k_v37_decoy.fasta \
           -V {infile} \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /mnt/fire/genomes/human/gatk/build37/hapmap_3.3.b37.vcf.gz  \
           --resource:omni,known=false,training=true,truth=false,prior=12.0  /mnt/fire/genomes/human/gatk/build37/1000G_omni2.5.b37.vcf.gz \
           --resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/fire/genomes/human/gatk/build37/1000G_phase1.snps.high_confidence.b37.vcf.gz \
           -an QD \
           -an MQ \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR  \
           -mode SNP \
           -O {outfile} \
           --tranches-file {output_file_name_root}.vqsr.stage1.SNP.tranches \
           --rscript-file {output_file_name_root}.vqsr.stage1.SNP.plots.R \
           2> {log_file}".format(gatk=config_dict['gatk4'],
                                 dbsnp=config_dict['dbsnp'],
                                 bwa_genome_fasta=config_dict['bwa2_genome_fasta'],
                                 infile=merged_vcf_in,
                                 outfile=out_file,
                                 output_file_name_root=root_name,
                                 log_file=log_file,
                                 temp_dir=temp_dir)

    return cmd




def VQSR_stage1_INDELs(config_dict,merged_vcf_in,root_name,out_file,output_file_name_root,log_file,temp_dir):
    """
    Run this on the genotyped merged vcf -
    """

    cmd = "{gatk} --java-options \"-Xms4G -Xmx4G -XX:ParallelGCThreads=2\" VariantRecalibrator \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
           -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           -R /mnt/fire/genomes/human/bwa2/human_g1k_v37_decoy.fasta \
           -V {infile} \
           --resource:mills,known=false,training=true,truth=true,prior=12.0 /mnt/fire/genomes/human/gatk/build37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/fire/genomes/human/gatk/build37/dbsnp_138.b37.vcf.gz \
           -an QD \
           -an MQ \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR  \
           -an DP \
           -mode INDEL \
           -O {outfile} \
           --tranches-file {output_file_name_root}.vqsr.stage1.INDEL.tranches \
           --rscript-file {output_file_name_root}.vqsr.stage1.INDEL.plots.R \
           2> {log_file}".format( gatk                  = config_dict['gatk4'],
                                  dbsnp                 = config_dict['dbsnp'],
                                  bwa_genome_fasta      = config_dict['bwa2_genome_fasta'],
                                  infile                = merged_vcf_in,
                                  outfile               = out_file,
                                  output_file_name_root = root_name,
                                  log_file              = log_file,
                                  temp_dir              = temp_dir)

    return cmd


def VQSR_stage2_SNPs(config_dict,merged_vcf_in,root_name,out_file,output_file_name_root,log_file,temp_dir):
    """
    """

    cmd = "{gatk} --java-options \"-Xms2G -Xmx2G -XX:ParallelGCThreads=2\" ApplyVQSR \
           -V {root_name}.merged.vcf.gz \
           --recal-file {output_file_name_root}.vqsr_SNP.recal.txt \
           -mode SNP \
           --truth-sensitivity-filter-level 99.9 \
           --create-output-variant-index true \
           --tranches-file {root_name}.vqsr.stage1.SNP.tranches \
           -O {outfile} \
           2> {log_file}".format(gatk=config_dict['gatk4'],
                                 output_file_name_root = output_file_name_root,
                                 dbsnp=config_dict['dbsnp'],
                                 bwa_genome_fasta=config_dict['bwa2_genome_fasta'],
                                 infile=merged_vcf_in,
                                 outfile=out_file,
                                 log_file=log_file,
                                 temp_dir=temp_dir,
                                 root_name=root_name)

    return cmd


def VQSR_stage2_INDELs(config_dict,merged_vcf_in,root_name,out_file,output_file_name_root,log_file,temp_dir):
    """  -V SNP.recalibrated_99.9.vcf.gz \
  -mode INDEL \
  --recal-file merged_indel1.recal \
  --tranches-file output_indel1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O indel.SNP.recalibrated_99.9.vcf.gz
    """

    cmd = "{gatk} --java-options \"-Xms2G -Xmx2G -XX:ParallelGCThreads=2\" ApplyVQSR \
           -V {infile} \
           --recal-file {root_name}.vqsr_INDEL.recal.txt \
           -mode INDEL \
           --truth-sensitivity-filter-level 99.9 \
           --create-output-variant-index true \
           --tranches-file {root_name}.vqsr.stage1.INDEL.tranches \
           -O {outfile} \
           2> {log_file}".format(gatk=config_dict['gatk4'],
                                 output_file_name_root = output_file_name_root,
                                 dbsnp=config_dict['dbsnp'],
                                 bwa_genome_fasta=config_dict['bwa2_genome_fasta'],
                                 infile=merged_vcf_in,
                                 outfile=out_file,
                                 log_file=log_file,
                                 temp_dir=temp_dir,
                                 root_name=root_name)

    return cmd
def funcotate_merged_vcf(config_dict,merged_vcf_in,root_name,out_file,output_file_name_root,log_file,temp_dir):
    """
    """

    cmd = "{gatk} --java-options \"-Xms2G -Xmx2G -XX:ParallelGCThreads=2\" Funcotator \
           -V {infile} \
           -R /mnt/fire/genomes/human/bwa2/human_g1k_v37_decoy.fasta \
           -O {outfile} \
           --data-sources-path {funcotator_resources} \
           --ref-version hg19 \
           --output-file-format VCF \
           2> {log_file}".format(gatk=config_dict['gatk4'],
                                 output_file_name_root = output_file_name_root,
                                 dbsnp=config_dict['dbsnp'],
                                 funcotator_resources=config_dict['funcotator_resources'],
                                 bwa_genome_fasta=config_dict['bwa2_genome_fasta'],
                                 infile=merged_vcf_in,
                                 outfile=out_file,
                                 log_file=log_file,
                                 temp_dir=temp_dir,
                                 root_name=root_name)

    return cmd

def vep_anotate_merged_vcf(config_dict,merged_vcf_in,root_name,vep_out,log_file,temp_dir):
    """
    docker run -t -i -v /mnt/fire/genomes/human/:/genomics:Z -v /mnt/fire/projects/nicky/wgs/wgs-db:/input:Z -v /zpool/genomic-resources/human/vep_cache:/data ensemblorg/ensembl-vep vep --cache --offline --format vcf --vcf --input_file /input/batch2.SNP_INDEL.recalibrated_99.9.vcf.gz --output_file /data/batch2.SNP_INDEL.recalibrated_99.9.vep.vcf.gz --everything --refseq --assembly GRCh37 --fasta /genomics/human_g1k_v37_decoy.fasta --compress_output bgzip
    """
    # Return just the .fasta file, to ensure correct mouting within the docker container
    fasta = config_dict['bwa2_genome_fasta'].split('/')[-1]

    cmd = "docker run -t -i \
           -v /mnt/fire/genomes/human/:/genomics:Z \
           -v /mnt/fire/projects/nicky/wgs/wgs-db:/input:Z \
           -v /zpool/genomic-resources/human/vep_cache:/data ensemblorg/ensembl-vep \
           vep \
           --fork 30 \
           --cache \
           --offline \
           --format vcf \
           --vcf --force_overwrite \
           --compress_output bgzip \
           --input_file /input/{merged_vcf_in} \
           --output_file /input/{vep_out} \
           --everything \
           --refseq \
           --assembly GRCh37 \
           --fasta /genomics/{bwa_genome_fasta} \
           2> {log_file}".format(gatk=config_dict['gatk4'],
                                 dbsnp=config_dict['dbsnp'],
                                 funcotator_resources=config_dict['funcotator_resources'],
                                 bwa_genome_fasta=fasta,
                                 merged_vcf_in=merged_vcf_in,
                                 vep_out=vep_out,
                                 log_file=log_file,
                                 temp_dir=temp_dir,
                                 root_name=root_name)

    return cmd

def wgs_preprocess_intervals_for_cnv(gatk4,bwa2_genome_fasta,output_file_name):
    """
    Procedure taken from:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531092
    Just run this once and save a copy with other reference data
    /opt/programs/gatk-4.2.2.0/gatk PreprocessIntervals \
     -R /data/genomic-resources/human/gatk_bundle/hs37d5/bwa-mem2/hs37d5.fa \
     --bin-length 1000 \
     --interval-merging-rule OVERLAPPING_ONLY \
     --padding 0 \
     -O /data/genomic-resources/human/gatk_bundle/hs37d5/bwa-mem2/hs37d5.targets_C.preprocessed.interval_list
    """
    cmd = f"{gatk4} PreprocessIntervals \
           -R {bwa2_genome_fasta} \
           --bin-length 1000 \
           --interval-merging-rule OVERLAPPING_ONLY \
           --padding 0 \
           -O {output_file_name}.preprocessed.interval_list" 
    
    return cmd
 
    
def wgs_collect_counts(gatk4,in_bam,intervals,exclude_intervals,outfile,log_name):
    """
    Create intervals as above. Use all for now.
    """
    cmd = f"{gatk4} CollectReadCounts \
         -I {in_bam} \
         -L {intervals} \
         --exclude-intervals {exclude_intervals} \
         --format TSV \
         --use-jdk-inflater \
         --use-jdk-deflater \
         --interval-merging-rule OVERLAPPING_ONLY \
         -O {outfile} \
         2> {log_name}.gatk-cnv.CollectReadCounts.log \n "
    
    return cmd    
    
def wgs_annotate_cnv_counts(gatk,in_tsv,segmental_dups,genome_fa,output,name):
    cmd = f"{gatk} \
    -L {in_tsv} \
    --segmental-duplication-track {segmental_dups} \
    -R {genome_fa} \
    -imr OVERLAPPING_ONLY \
    -O {output} "    
    
def wgs_filter_intervals(gatk4,interval_list,annotated_intervals,in_tsvs,exclude_intervals,outfile,log_name):
    """
    
    
    """
    cmd = f"{gatk4} FilterIntervals \
          -L {interval_list} \
          --annotated-intervals {annotated_intervals} \
          {in_tsvs} \
          -XL {exclude_intervals} \
          -imr OVERLAPPING_ONLY \
          -O {outfile} \
          2> {log_name}.gatk-cnv.filter_intervals.log \n"
    
    return cmd
    

def wgs_determine_contig_ploidy(gatk4,filtered_interval_list,in_tsvs,ploidy_priors,outfile,log_name):
    cmd = f"{gatk4} DetermineGermlineContigPloidy \
           -L {filtered_interval_list} \
           --interval-merging-rule OVERLAPPING_ONLY \
           {in_tsvs} \
           --contig-ploidy-priors {ploidy_priors} \
           --output . \
           --output-prefix {outfile} \
           --verbosity DEBUG \
           2> {log_name}.gatk-cnv.determine.ploidy.log > LOG\n"
    
    return cmd

def wgs_germline_cohort_cnv_caller(gatk4,filtered_interval_list,in_tsvs,ploidy_priors,outfile,log_name):
    cmd = f"{gatk4} GermlineCNVCaller \
           --run-mode COHORT \
           -L {filtered_interval_list} \
           {in_tsvs} \
           --contig-ploidy-calls ploidy-calls \
            \
           --interval-merging-rule OVERLAPPING_ONLY \
           {in_tsvs} \
           --contig-ploidy-priors {ploidy_priors} \
           --output . \
           --output-prefix {outfile} \
           --verbosity DEBUG \
           2> {log_name}.gatk-cnv.determine.ploidy.log > LOG\n"

def wgs_split_interval_list(gatk4):
    cmd = f"{gatk4} IntervalListTools \
           --INPUT {gc.filtered.interval} \
           --SUBDIVISION_MODE INTERVAL_COUNT \
           --scatter_content 5000 \
           --OUTPUT scatter"
    
def wgs_germline_case_cnv_caller(gatk4,in_tsv,ploidy_priors,outfile,log_name):
    cmd = f"{gatk4} GermlineCNVCaller \
           --run-mode CASE \
           {in_tsv} \
           --contig-ploidy-calls ploidy-calls \
           --model {cohortX_model} \
           --output {case_vers_cohortX} \
           --output-prefix {outfile} \
           --verbosity DEBUG \
           2> {log_name}.gatk-cnv.CNV-case.calling.log > LOG\n"
    
    
    
    
    
    
    
#cancer functions    
def wgs_denoise_counts(gatk4,in_counts,panel_of_normals,out_standardised,out_denoised,log_name):
    """
    Run gatk4 on each count file and create a denoised file
    """
    cmd = f"{gatk4} --java-options \"-Xmx12g\" DenoiseReadCounts \
             -I {in_counts} \
             --count-panel-of-normals {panel_of_normals} \
             --standardized-copy-ratios {out_standardised} \
             --denoised-copy-ratios {out_denoised} \
             2> {log_name}.gatk-cnv.denoise_counts.log \n"
    
    return cmd

def wgs_plot_denoised_data(gatk4,standardized_counts,denoised_counts,genome_dict,plot_dir,out_prefix,log_name):
    """
    plot the data
    set min contig length to ignore alt contigs
    """
    cmd = f"{gatk4} PlotDenoisedCopyRatios \
            --standardized-copy-ratios {standardized_counts} \
            --denoised-copy-ratios {denoised_counts} \
            --sequence-dictionary {genome_dict} \
            --minimum-contig-length 46709983 \
            --output {plot_dir} \
            --output-prefix {out_prefix} \
            2> {log_name}.gatk-cnv.plot.denoise_counts.log \n"
    
    return cmd

def reheader_vcf(infile,outfile):
    
    with open("%s" % infile,'r+') as inputfile:
        for line in inputfile:
            line = line.strip().split
    return cmd

def call_mutect2(gatk,bam,mito_out_vcf,fasta,regions,output_file_name):
    cmd = f"{gatk} --java-options \"-XX:+UseParallelGC -XX:ParallelGCThreads=2 -Xmx4G\" Mutect2 \
            -I {bam} \
            -O {mito_out_vcf} \
            -R {fasta} \
            -L {regions} \
            --interval-padding 50 \
            -bamout {output_file_name}.bamout.bam \
            --mitochondria-mode true \
            --disable-bam-index-caching false \
            -G AlleleSpecificAnnotation \
            -G AS_StandardAnnotation \
            -G ReducibleAnnotation \
            -G StandardAnnotation \
            --use-jdk-inflater TRUE \
            --use-jdk-deflater TRUE \
            --create-output-bam-md5 \
            --create-output-variant-md5 \
            --create-output-bam-index \
            --create-output-variant-index \
            2> {output_file_name}.Mutect2onMito.log\n "
    return cmd 



def merge_vcfs(config_dict,in_vcfs,merged_vcf,log_name):
    cmd = "{java18} -jar {picard2} MergeVcfs \
    I={in_vcfs}\
    O={merged_vcf} \
    2>{log_name}.merged_vcfs.log".format(
                                  java18     = config_dict['java18'],
                                  picard2    = config_dict['picard2'],
                                  in_vcfs    = ' I='.join(in_vcfs),
                                  merged_vcf = merged_vcf,
                                  log_name   = log_name)
    return cmd


def merge_cohort_vcfs(config_dict,in_vcfs,merged_vcf,log_name):
    cmd = "{java18} -jar {picard2} MergeVcfs \
    I={in_vcfs}\
    O={merged_vcf} \
    2>{log_name}.merged_vcfs.log".format(
                                  java18     = config_dict['java18'],
                                  picard2    = config_dict['picard2'],
                                  in_vcfs    = ' I='.join(in_vcfs),
                                  merged_vcf = merged_vcf,
                                  log_name   = log_name)
    return cmd


def merge_vcfs(config_dict,in_vcfs,merged_vcf,log_name,comment):
    cmd = "{java18} -jar {picard2} MergeVcfs \
    I={in_vcfs}\
    O={merged_vcf} \
    CO={comment} \
    2>{log_name}.merged_vcfs.log".format(
                                  java18     = config_dict['java18'],
                                  picard2    = config_dict['picard2'],
                                  in_vcfs    = ' I='.join(in_vcfs),
                                  merged_vcf = merged_vcf,
                                  log_name   = log_name,
                                  comment    = comment)
    return cmd
       
def merge_bamout(config_dict,in_bams,merged_bam,log_name):    
    cmd =  "{gatk4} --java-options \" -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx32g \" MergeSamFiles \
            --INPUT {bams} \
            --USE_THREADING TRUE \
            --USE_JDK_DEFLATER \
            --USE_JDK_INFLATER \
            --COMPRESSION_LEVEL 5 \
            --VALIDATION_STRINGENCY LENIENT \
            --CREATE_INDEX TRUE \
            --CREATE_MD5_FILE TRUE \
            --OUTPUT {merged_bam} \
            2> {log_name}.mergeBamout.log".format(
                                                  gatk4      = config_dict['gatk4'],
                                                  bams       = ' --INPUT '.join(in_bams),
                                                  merged_bam = merged_bam,
                                                  log_name   = log_name)    
    return cmd
    
    
    
def SingleSampleGenomicsDBImport(config_dict,ROI,infiles,chrom,log_name,outfile):    
    cmd = "/opt/programs/gatk-4.1.2.0/gatk --java-options \"-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" \
              GenomicsDBImport \
              --genomicsdb-workspace-path {output} \
              --intervals {chromosome} \
              --batch-size 3 \
              -V {vcfs} \
              --reader-threads 1 \
              2> {log_name}GenomicsDBImport.log".format(
                                                        vcfs       = ' -V '.join(infiles),
                                                        chromosome = chrom,
                                                        output     = outfile,
                                                        log_name   = log_name,
                                                        cwd        = os.getcwd())
    return cmd
    

def vep_annotate_gatk_vcf(config_dict,infile,outfile,log_name,outtype):    
    #print(infile,'-->',outfile)
    """
    is --phased needed?
    --plugin GeneSplicer,/opt/programs/GeneSplicer/bin/linux/genesplicer,/opt/programs/GeneSplicer/human/ \
    """
    
    #log_name = outfile[:-4]
    cmd = "{vep} \
           --{outtype} \
           --buffer_size 8000 \
           --flag_pick \
           --exclude_predicted \
           --var_synonyms \
           --hgvsg \
           --hgvs \
           --everything \
           --dont_skip \
           --fork 8 \
           --force_overwrite \
           --offline \
           --assembly GRCh37 \
           --fasta {vep_genome_fasta} \
           --db_version {vepdb_version} \
           --cache_version {vepcache_version} \
           --refseq  \
           --variant_class \
           --transcript_version \
           --protein \
           --symbol \
           --phased \
           --pubmed \
           --canonical \
           --allele_number \
           --dir {vep_dir} \
           --dir_cache {vep_cache_dir} \
           --dir_plugins {vep_plugins} \
           -custom /data/genomic-resources/human/vep-related/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
           -custom /data/genomic-resources/human/vep-related/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
           --plugin SingleLetterAA \
           -plugin CADD,/data/genomic-resources/human/CADD/whole_genome_SNVs.tsv.gz,/data/genomic-resources/human/CADD/InDels_inclAnno.tsv.gz \
           --plugin Blosum62 \
           --plugin ExACpLI,/data/genomic-resources/human/vep-related/ExACpLI_values.txt \
           --plugin LoFtool,/data/genomic-resources/human/vep-related/LoFtool_scores.txt \
           --plugin SpliceRegion \
           --plugin MaxEntScan,/data/genomic-resources/human/vep-related/maxEntscan/fordownload \
           --plugin MPC,/data/genomic-resources/human/vep-related/fordist_constraint_official_mpc_values_v2.txt.gz \
           --plugin FATHMM_MKL,/data/genomic-resources/human/vep-related/fathmm-MKL_Current.tab.gz  \
           --plugin dbNSFP,/data/genomic-resources/human/dbNSFPv4.0a/dbNSFP4.0a.gz,aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc,Uniprot_entry,HGVSc_ANNOVAR,HGVSp_ANNOVAR,HGVSc_snpEff,HGVSp_snpEff,HGVSc_VEP,HGVSp_VEP,APPRIS,GENCODE_basic,TSL,VEP_canonical,cds_strand,refcodon,codonpos,codon_degeneracy,Ancestral_allele,AltaiNeandertal,Denisova,VindijiaNeandertal,SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,MVP_score,MVP_rankscore,MPC_score,MPC_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,DEOGEN2_score,DEOGEN2_rankscore,DEOGEN2_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-raw_coding_rankscore,Eigen-pred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-PC-phred_coding,GenoCanyon_score,GenoCanyon_rankscore,integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_rankscore,HUVEC_confidence_value,LINSIGHT,LINSIGHT_rankscore,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian,phyloP30way_mammalian_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons30way_mammalian,phastCons30way_mammalian_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,bStatistic,bStatistic_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,UK10K_AC,UK10K_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,ExAC_nonTCGA_AC,ExAC_nonTCGA_AF,ExAC_nonTCGA_Adj_AC,ExAC_nonTCGA_Adj_AF,ExAC_nonTCGA_AFR_AC,ExAC_nonTCGA_AFR_AF,ExAC_nonTCGA_AMR_AC,ExAC_nonTCGA_AMR_AF,ExAC_nonTCGA_EAS_AC,ExAC_nonTCGA_EAS_AF,ExAC_nonTCGA_FIN_AC,ExAC_nonTCGA_FIN_AF,ExAC_nonTCGA_NFE_AC,ExAC_nonTCGA_NFE_AF,ExAC_nonTCGA_SAS_AC,ExAC_nonTCGA_SAS_AF,ExAC_nonpsych_AC,ExAC_nonpsych_AF,ExAC_nonpsych_Adj_AC,ExAC_nonpsych_Adj_AF,ExAC_nonpsych_AFR_AC,ExAC_nonpsych_AFR_AF,ExAC_nonpsych_AMR_AC,ExAC_nonpsych_AMR_AF,ExAC_nonpsych_EAS_AC,ExAC_nonpsych_EAS_AF,ExAC_nonpsych_FIN_AC,ExAC_nonpsych_FIN_AF,ExAC_nonpsych_NFE_AC,ExAC_nonpsych_NFE_AF,ExAC_nonpsych_SAS_AC,ExAC_nonpsych_SAS_AF,gnomAD_exomes_flag,gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_exomes_nhomalt,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AN,gnomAD_exomes_AFR_AF,gnomAD_exomes_AFR_nhomalt,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AN,gnomAD_exomes_AMR_AF,gnomAD_exomes_AMR_nhomalt,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AN,gnomAD_exomes_ASJ_AF,gnomAD_exomes_ASJ_nhomalt,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AF,gnomAD_exomes_EAS_nhomalt,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AN,gnomAD_exomes_FIN_AF,gnomAD_exomes_FIN_nhomalt,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AN,gnomAD_exomes_NFE_AF,gnomAD_exomes_NFE_nhomalt,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AN,gnomAD_exomes_SAS_AF,gnomAD_exomes_SAS_nhomalt,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AN,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_controls_AC,gnomAD_exomes_controls_AN,gnomAD_exomes_controls_AF,gnomAD_exomes_controls_nhomalt,gnomAD_exomes_controls_AFR_AC,gnomAD_exomes_controls_AFR_AN,gnomAD_exomes_controls_AFR_AF,gnomAD_exomes_controls_AFR_nhomalt,gnomAD_exomes_controls_AMR_AC,gnomAD_exomes_controls_AMR_AN,gnomAD_exomes_controls_AMR_AF,gnomAD_exomes_controls_AMR_nhomalt,gnomAD_exomes_controls_ASJ_AC,gnomAD_exomes_controls_ASJ_AN,gnomAD_exomes_controls_ASJ_AF,gnomAD_exomes_controls_ASJ_nhomalt,gnomAD_exomes_controls_EAS_AC,gnomAD_exomes_controls_EAS_AN,gnomAD_exomes_controls_EAS_AF,gnomAD_exomes_controls_EAS_nhomalt,gnomAD_exomes_controls_FIN_AC,gnomAD_exomes_controls_FIN_AN,gnomAD_exomes_controls_FIN_AF,gnomAD_exomes_controls_FIN_nhomalt,gnomAD_exomes_controls_NFE_AC,gnomAD_exomes_controls_NFE_AN,gnomAD_exomes_controls_NFE_AF,gnomAD_exomes_controls_NFE_nhomalt,gnomAD_exomes_controls_SAS_AC,gnomAD_exomes_controls_SAS_AN,gnomAD_exomes_controls_SAS_AF,gnomAD_exomes_controls_SAS_nhomalt,gnomAD_exomes_controls_POPMAX_AC,gnomAD_exomes_controls_POPMAX_AN,gnomAD_exomes_controls_POPMAX_AF,gnomAD_exomes_controls_POPMAX_nhomalt,gnomAD_genomes_flag,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_nhomalt,gnomAD_genomes_AFR_AC,gnomAD_genomes_AFR_AN,gnomAD_genomes_AFR_AF,gnomAD_genomes_AFR_nhomalt,gnomAD_genomes_AMR_AC,gnomAD_genomes_AMR_AN,gnomAD_genomes_AMR_AF,gnomAD_genomes_AMR_nhomalt,gnomAD_genomes_ASJ_AC,gnomAD_genomes_ASJ_AN,gnomAD_genomes_ASJ_AF,gnomAD_genomes_ASJ_nhomalt,gnomAD_genomes_EAS_AC,gnomAD_genomes_EAS_AN,gnomAD_genomes_EAS_AF,gnomAD_genomes_EAS_nhomalt,gnomAD_genomes_FIN_AC,gnomAD_genomes_FIN_AN,gnomAD_genomes_FIN_AF,gnomAD_genomes_FIN_nhomalt,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AN,gnomAD_genomes_NFE_AF,gnomAD_genomes_NFE_nhomalt,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,gnomAD_genomes_controls_AC,gnomAD_genomes_controls_AN,gnomAD_genomes_controls_AF,gnomAD_genomes_controls_nhomalt,gnomAD_genomes_controls_AFR_AC,gnomAD_genomes_controls_AFR_AN,gnomAD_genomes_controls_AFR_AF,gnomAD_genomes_controls_AFR_nhomalt,gnomAD_genomes_controls_AMR_AC,gnomAD_genomes_controls_AMR_AN,gnomAD_genomes_controls_AMR_AF,gnomAD_genomes_controls_AMR_nhomalt,gnomAD_genomes_controls_ASJ_AC,gnomAD_genomes_controls_ASJ_AN,gnomAD_genomes_controls_ASJ_AF,gnomAD_genomes_controls_ASJ_nhomalt,gnomAD_genomes_controls_EAS_AC,gnomAD_genomes_controls_EAS_AN,gnomAD_genomes_controls_EAS_AF,gnomAD_genomes_controls_EAS_nhomalt,gnomAD_genomes_controls_FIN_AC,gnomAD_genomes_controls_FIN_AN,gnomAD_genomes_controls_FIN_AF,gnomAD_genomes_controls_FIN_nhomalt,gnomAD_genomes_controls_NFE_AC,gnomAD_genomes_controls_NFE_AN,gnomAD_genomes_controls_NFE_AF,gnomAD_genomes_controls_NFE_nhomalt,gnomAD_genomes_controls_POPMAX_AC,gnomAD_genomes_controls_POPMAX_AN,gnomAD_genomes_controls_POPMAX_AF,gnomAD_genomes_controls_POPMAX_nhomalt,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_hgvs,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain,GTEx_V7_gene,GTEx_V7_tissue,Geuvadis_eQTL_target_gene \
           -i {input} \
           -o {output} \
           --compress_output bgzip \
           2> {log_name}".format(
                                     vep=config_dict['vep'],
                                     vep_genome_fasta=config_dict['vep_genome_fasta'],
                                     vep_dir=config_dict['vep_dir'],
                                     vep_cache_dir=config_dict['vep_cache_dir'],
                                     vepdb_version=config_dict['vepdb_version'],
                                     vepcache_version=config_dict['vepcache_version'],
                                     vep_plugins=config_dict['vep_plugins'],
                                     input=infile,
                                     output=outfile,
                                     log_name=log_name,
                                     outtype=outtype)
    return cmd
    
def vep_minimal_annotate_gatk_vcf(config_dict,infile,outfile,log_name,outtype):    
    #print(infile,'-->',outfile)
    """
    is --phased needed?
    """
    
    #log_name = outfile[:-4]
    cmd = "{vep} \
           --{outtype} \
           --buffer_size 100000 \
           --flag_pick \
           --exclude_predicted \
           --hgvsg \
           --hgvs \
           --dont_skip \
           --fork 10 \
           --force_overwrite \
           --offline \
           --assembly GRCh37 \
           --fasta {vep_genome_fasta} \
           --db_version {vepdb_version} \
           --cache_version {vepcache_version} \
           --refseq  \
           --variant_class \
           --transcript_version \
           --protein \
           --symbol \
           --phased \
           --pubmed \
           --canonical \
           --allele_number \
           --dir {vep_dir} \
           --dir_cache {vep_cache_dir} \
           --dir_plugins {vep_plugins} \
           -custom /data/genomes/human/annotations/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
           -custom /data/genomes/human/annotations/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
           --plugin SingleLetterAA \
           --plugin ExACpLI,/data/genomes/human/annotations/ExACpLI_values.txt \
           -i {input} \
           -o {output} \
           --compress_output bgzip \
           2> {log_name}.vep.log".format(
                                     vep=config_dict['vep'],
                                     vep_genome_fasta=config_dict['vep_genome_fasta'],
                                     vep_dir=config_dict['vep_dir'],
                                     vep_cache_dir=config_dict['vep_cache_dir'],
                                     vepdb_version=config_dict['vepdb_version'],
                                     vepcache_version=config_dict['vepcache_version'],
                                     vep_plugins=config_dict['vep_plugins'],
                                     input=infile,
                                     output=outfile,
                                     log_name=log_name,
                                     outtype=outtype)
    return cmd    
    


#Read the table into a df and write to the sqlite db
def write_df_to_sql(df,vcf_file_name, table_name, engine):
    print("Writing:", vcf_file_name, "to:", table_name)
    try:
        df.to_sql(table_name, con=engine)
    except ValueError:
        print("Tables already exist, made today; therefore passing")
        pass
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df.head()) # A method to print df nicely
    
def vcf_to_df(vcf):
    #
    #rename(columns={'#CHROM':'CHROM'})
    with open(vcf,'r') as f:
        lines = [l for l in f if not l.startswith('##')]
        #print(lines,"*"*100,"\n")
        #callset = allel.read_vcf(infile,fields="*")
    #print(callset.keys())
    # Save vcf as zarr
    #allel.vcf_to_zarr(infile,outfile,group='22', fields='*',alt_number=8,log=sys.stdout,
    #compressor=numcodecs.Blosc(cname='zstd',clevel=1,shuffle=False))
    
        data = io.StringIO(''.join(lines))
        print(data)
        df = pd.read_csv(data,
                      dtype={'#CHROM':str,'POS':int,'ID':str,'REF':str,'ALT':str,'QUAL':str,'FILTER':str,'INFO':str},sep="\t")
    
import pandas as pd
import vcf

def read_vcf_to_df(f):
    reader = vcf.Reader(open(f))
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    return out    
#def add_tags_to_vcf(vcf_in,vcf_out,tag1_name,tag2_name,tag3_name,tag4_name,tag5_name,poly_link,freq_link,sanger_file,worksheet):
def add_tags_to_vcf(vcf_in,vcf_out,tag1_name,tag2_name,tag4_name,tag5_name,poly_link,freq_link,worksheet):
    """
    write arbitary new fields to a VCF file, to the INFO field
    """
    ######################################################################
    # A better method, but the raw ASCII introdfuced by Annovar means PyVCF misses information:
    #vcf_reader = vcf.Reader(filename='%s' % vcf_in)
    #vcf_writer = vcf.Writer(open('%s' %(vcf_out), 'w'), vcf_reader)
    #for record in vcf_reader:
    #    vcf_writer.write_record(record)        
    ######################################################################
    
    # Parse the poly list
    poly_db = {}
    def parse_WMRGL_polymorphism_data(poly_link):
        poly_set = set()
        with open('%s' % poly_link, 'r') as db_file:
            #next(db_file) # removed on 04-01-2016 - no need to skip first line.
            for line in db_file:
                line=line.strip()
                line = line.split('\t')
                (key, variant_key) = (line[0],"{0}|{1}|{2}|{3}".format(line[1],line[2],line[3],line[4]))  # tuple = (bin,'chr+pos+ref+alt')
                poly_db[variant_key] = line[5]
                poly_set.add(variant_key)
    
        return poly_db
        for k,v in poly_db.items():
            print ("%"*200,k,v)
    
    
    parse_WMRGL_polymorphism_data(poly_link)
    #def data_formatter(poly_anno,frequency,AF,sanger_result,worksheet):
    def data_formatter(poly_anno,frequency,AF):
        """ print the variants in VCF format, adding optional fields as required
        
        deleted:
        
        {tag3_name}={tag3_value}
        
        tag3_name  = tag3_name,\
        tag3_value = sanger_result,\
        
        
        ... on 13/10/20
        """
        print(v)
        data = "{Chrom}\t{Pos}\t{ID}\t{REF}\t{ALT}\t{Qual}\t{Filter}\t{INFO};{tag1_name}={tag1_value};{tag2_name}={tag2_value};{tag4_name}={tag4_value};{tag5_name}={tag5_value}\t{FORMAT}\t{Genotype}\n".format(\
               Chrom      = v[0],\
               Pos        = v[1],\
               ID         = v[2],\
               REF        = v[3],\
               ALT        = v[4],\
               Qual       = v[5],\
               Filter     = v[6],\
               INFO       = v[7],\
               FORMAT     = v[8],\
               Genotype   = v[9],\
               tag1_name  = tag1_name,\
               tag1_value = poly_anno,\
               tag2_name  = tag2_name,\
               tag2_value = frequency,\
               tag4_name  = tag4_name,\
               tag4_value = worksheet,\
               tag5_name  = tag5_name,\
               tag5_value = AF)
        
        return data

    with open(freq_link,'rb') as pickled_dict:
        master_database = pickle.load(pickled_dict)
    #master_database = shelve.open(freq_link)
    all_vars = master_database.keys()
    total_samples_set = set()
    for var in all_vars:
	      for run_sample in master_database[var]:
	          total_samples_set.add(run_sample[1])
    total_samples = len(total_samples_set)

    #print ('ADD TAGS TO VCF',vcf_in)
    with open("%s" % (vcf_in), 'r') as input_vcf:
        with open("%s" % (vcf_out), 'w+') as output_vcf:
            for line in input_vcf:
                
                if line.startswith("##"):
                    output_vcf.write(line)
                    continue
                if line.startswith("#"):
                    output_vcf.write(line)
                    continue
                if line.startswith("#CHROM"):
                    output_vcf.write(line)
                    continue
                if len(line.strip()) == 0:                 # deal with empty vcfs
                    print("this line starts with nothing")
                    output_vcf.write(line)
                    break
                else:
                    v = line.strip().split('\t')
                    variant_key = "{0}|{1}|{2}|{3}".format(v[0],v[1],v[3],v[4])
                    if master_database[variant_key]:
                        #Calculate: # runs, # samples ever, # variants ever, # percentage variant over all runs ever.
                        # use list comprehension to interogate the list of tuples identified by the variant_key
                        total_runs = len(set([x[0] for x in master_database[variant_key]])) # set added on 02-09-15
                        total_vars = len(set([x[1] for x in master_database[variant_key]]))
                        percent = (float(total_vars)/float(total_samples))*100
                        pct = round(percent,1)
                        data = '%s|%s|%s|%s' % (total_vars,total_samples,pct,total_runs)
                        AF = (float(total_vars)/float(total_samples))
                    freq_result_tuple = [data,AF]
                    if v[0].startswith('chr'): # ignores variants found on unplaced scaffold sequences as these break the bedtools intersect
                        sanger_result = add_sanger_primers_to_vcf(sanger_file,variant_key)
                    if variant_key in poly_db.keys():
                        output_vcf.write(data_formatter(poly_db[variant_key],freq_result_tuple[0],freq_result_tuple[1],sanger_result,worksheet))
                    else:
                        #output_vcf.write(data_formatter(".",freq_result_tuple[0],freq_result_tuple[1],sanger_result,worksheet))
                        output_vcf.write(data_formatter(".",freq_result_tuple[0],freq_result_tuple[1]))

    #master_database.close()

    
    
#NOT USED/NEEDED
#@transform(["*.split.g.vcf.gz"],suffix(".split.g.vcf.gz"),".split.g.fixed.vcf.gz")
def rename_vcf_header(infile,outfile):
    old_name = re.sub(".bwa.drm.realn.recal.chr([0-9]|[0-9][0-9]|M|X|Y).split.g.vcf.gz","",infile)
    fixed_name = re.sub(".split.g.vcf.gz","",infile)
    os.system("bgzip -d {vcf}".format(vcf=infile))
    os.system("sed \'s/FORMAT\t{old_name}/FORMAT\t{fixed_name}/\' {infile} > {outfile}".format(old_name=old_name,fixed_name=fixed_name,infile=infile[:-3],outfile=outfile[:-3]))
    os.system("bgzip {vcf}".format(vcf=fixed_name + '.split.g.fixed.vcf'))
    os.system("tabix {vcf}".format(vcf=fixed_name + '.split.g.fixed.vcf.gz'))
    
####################################################################################################################################################
################################################   NOT USED , as it breaks the genomeVariantDB import   ############################################
####################################################################################################################################################
####################################################################################################################################################
#@jobs_limit(8)
#@subdivide(["*bwa.drm.realn.recal.bam"],formatter(r".bwa.drm.realn.recal.bam$"),
#("{basename[0]}.chr1.split.g.vcf.gz","{basename[0]}.chr2.split.g.vcf.gz","{basename[0]}.chr3.split.g.vcf.gz",
#"{basename[0]}.chr4.split.g.vcf.gz","{basename[0]}.chr5.split.g.vcf.gz","{basename[0]}.chr6.split.g.vcf.gz",
#"{basename[0]}.chr7.split.g.vcf.gz","{basename[0]}.chr8.split.g.vcf.gz","{basename[0]}.chr9.split.g.vcf.gz",
#"{basename[0]}.chr10.split.g.vcf.gz","{basename[0]}.chr11.split.g.vcf.gz","{basename[0]}.chr12.split.g.vcf.gz",
#"{basename[0]}.chr13.split.g.vcf.gz","{basename[0]}.chr14.split.g.vcf.gz","{basename[0]}.chr15.split.g.vcf.gz",
#"{basename[0]}.chr16.split.g.vcf.gz","{basename[0]}.chr17.split.g.vcf.gz","{basename[0]}.chr18.split.g.vcf.gz",
#"{basename[0]}.chr19.split.g.vcf.gz","{basename[0]}.chr20.split.g.vcf.gz","{basename[0]}.chr21.split.g.vcf.gz",
#"{basename[0]}.chr22.split.g.vcf.gz","{basename[0]}.chrX.split.g.vcf.gz","{basename[0]}.chrY.split.g.vcf.gz",
#"{basename[0]}.chrM.split.g.vcf.gz"))
def run_HaploTypeCaller4(infile,outfiles):
    """
    Use parallel, as a quick work around. Other possibilities:
    Create a seperate @transform function for each subdivide
    Create a seperate function somewhere else and exploit Ruffus oo to run it.
    """     
    # create a text file of all the commands, to give to parallel
    with open("GATKv4_HC_jobs",'w+') as outputfile:
        for i,outfilename in enumerate(outfiles):
            output_file_name = outfiles[i]
            output_file_name_root = outfiles[i][:-9]
            chrom = re.search("(.chr)([0-9a-zA-Z]+)(.)(.+)",outfilename).group(2)
            command = cmd.HaploTypeCaller_v4(config_dict           = config_dict,
                                             ROI                   = config_dict['ROI_chr%s' % chrom],
                                             bam                   = infile,
                                             output_file_name      = output_file_name,
                                             output_file_name_root = output_file_name_root  )
            outputfile.write(command)
    os.system("parallel -j 8 :::: GATKv4_HC_jobs")
    
#@follows(run_HaploTypeCaller4)
#@collate("*.split.g.vcf.gz", formatter("([^/]+).chr([0-9]|[0-9][0-9]|X|Y|M).split.g.vcf.gz$"),"{path[0]}/{1[0]}.HC.g.vcf.gz")
def merge_split_gvcfs(infiles,outfile):
    #print (infiles,outfile)
    log_name   = outfile[:-7]
    m_x_y_vcfs = list(infiles[-3:])           # slice the M, X and Y
    m_x_y_vcfs_sorted = [m_x_y_vcfs[1],m_x_y_vcfs[2],m_x_y_vcfs[0]] # get into order: X,Y,M
    split_vcfs = list(infiles)                # convert ruffus tuple to list
    split_vcfs = (split_vcfs[:-3])            # remove X and Y for sorting
    split_vcfs.sort(key=lambda x: int(x.split('chr')[1].split('.')[0])) # sort on the chr value within the file name
    newlist    = split_vcfs + m_x_y_vcfs_sorted  # Create the file order for vcf concatenation
    
    print(cmd.merge_gvcfs(config_dict = config_dict,
                               in_vcfs    = newlist,
                               merged_gvcf = outfile,
                               log_name   = log_name))
    

#@follows(run_HaploTypeCaller4)
#@collate("*.bamout.bam", formatter("([^/]+).chr([0-9]|[0-9][0-9]|X|Y|M).split.bamout.bam$"),"{path[0]}/{1[0]}.HC.bamout.bam")
#def merge_bamout(infiles,outfile):
#    """
#    Sort the chromosomes in proper order so increase speed.
#    i.e. the MergeSamFiles will then not have to sort the file.
#    """
#    
#    log_name   = outfile[:-4]
#    m_x_y_bams = list(infiles[-3:])           # slice the M, X and Y
#    m_x_y_bams_sorted = [m_x_y_bams[1],m_x_y_bams[2],m_x_y_bams[0]] # get into order: X,Y,M
#    split_bams = list(infiles)                # convert ruffus tuple to list
#    split_bams = (split_bams[:-3])            # remove X and Y for sorting
#    split_bams.sort(key=lambda x: int(x.split('chr')[1].split('.')[0])) # sort on the chr value within the file name
#    newlist    = split_bams + m_x_y_bams_sorted  # Create the file order for vcf concatenation
#     
#    os.system(cmd.merge_bamout(config_dict = config_dict,
#                               in_bams    = newlist,
#                               merged_bam = outfile,
#                               log_name   = log_name))
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################    
#@transform(["*.bwa.drm.realn.recal.bam"],suffix(".bwa.drm.realn.recal.bam"),".HC.g.vcf.gz")
def runGATKv4(infile,outfile):
    os.system(cmd.HaploTypeCaller_v4(config_dict          = config_dict,
                                    ROI                   = config_dict['ROI_ensembl_pad5'],
                                    bam                   = infile,
                                    output_file_name      = outfile,
                                    output_file_name_root = outfile[:-7]))    
    
    
    
    
    



#bbmap - on FASTQs straight from bcl2fastq
def bbmap_default_wgs_single_read(config_dict,worksheet,chemistry,FASTQ,out_primary_bam,root_name):
    """
    Run bbmap with settings that caught the very long INDEL in FGF14
    """
    cmd = f"{config_dict['bbmap']} \
           -Xmx24g \
           ref={config_dict['bbmap_genome']} \
           nodisk=t \
           rgid={root_name} \
           rgpl=illumina \
           rgsm={root_name} \
           in={FASTQ} \
           out={out_primary_bam} \
           outu={root_name}.bbmap.unmapped.grch37.fastq \
           unpigz=t \
           threads=6 \
           maxindel=1000k \
           rescuedist=1000000 \
           usejni=t \
           rgid={root_name} \
           mdtag=t \
           nhtag=t \
           xmtag=t \
           nmtag=t \
           stoptag=t \
           lengthtag=t \
           idtag=t \
           inserttag=t \
           scoretag=t \
           usemodulo=t \
           statsfile={root_name}.bbmap.align.stats.log \
           machineout=t"
    
    return cmd
    
    
    
    
#bbmap
def bbmap_align_FASTQ(config_dict,worksheet,chemistry,trimmed_fwd_FASTQ,trimmed_rev_FASTQ,out_primary_bam,root_name):
    """
    Align quality filtered FASTQ files using bwa mem. Use samblaster to directly pipe to duplicate marker.
    SAMBLASTER marks duplicates in a single pass over a SAM file in which all alignments for the same read-id are grouped together...and ...outputs read-id grouped SAM
    To do: pipe into sambamba for sorting (multicore).

    Read group configuration: RG -
                                    ID = worksheet.flowcellID.Lane#             
                                    CN = centre (WMRGL)
                                    DS = description (TSO)
                                    DT = date and time (ISO)
                                    SM = Sample name (sliced infile name)
                                    LB = Library (sliced infile name)
                                    PL = Platform key word (Illumina)
                                    PN =
    ./bam polishBam  --in testFiles/sortedSam.sam --out results/updatedSam.sam \
    --log results/updated.log \
    --checkSQ \
    --fasta testFiles/testFasta.fa \
    --AS my37 \
    --UR testFasta.fa \
    --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none"

    """
    #readConfigFile(config_file)
    
    import datetime
    i = datetime.datetime.now()
    date_time_iso = "%s" % i.isoformat()
    
    
    
    lane = re.search(r"(L00[1-8])",trimmed_fwd_FASTQ)
    #lane = '%s' % re.findall("(L00[1-8])",trimmed_fwd_FASTQ)
    if lane is not None:
        lane = list(filter(None,lane.groups()))
    else:
        lane = 'merged'
    original_FASTQ = re.sub(r".qfilter.fastq.gz",".fastq.gz",trimmed_fwd_FASTQ)
    path           = os.path.realpath("%s" % original_FASTQ)
    sample_name    = re.sub(r".qfilter.fastq.gz","",trimmed_fwd_FASTQ[:-24])
    sample_name    = re.sub(r"_L00[1-8]","",sample_name)
    #sample_name_for_bam = re.sub(r"%s/" % args.named_directory,"",sample_name[:-29])
    
    flowcell_match = re.search(r"(_|.)([A-Z0-9]{9})(_|.)",path)
    if flowcell_match is not None:
       flowcell = filter(None, flowcell_match.groups()[1])
    else:
       flowcell = 'FC_UNKNOWN'
    #sed \'s/-R @RG.*//\' - | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    #grep -f /network/sequenced/genomic_resources/gatk_bundle/chrs.txt | \
    
    
    #maxindel=200k
    #--removeDups 22-11-2019
    align_cmd = "{{bbmap} in={trimmed_fwd_FASTQ} in2={trimmed_frev_FASTQ} ref={genome} maxindel=2k".format(
                    bbmap             = config_dict['bbmap'],
                    trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                    trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                    genome            = config_dict['bwa_genome_fasta'])
    
    
    align_cmd = "{bwa} mem -t {threads} -M -k 18 -R '@RG\\tID:{worksheet}.{flowcell}.{lane}\\tCN:WMRGL\\tDS:{chemistry}\\tDT:{date}\\tSM:{sample_name}\\tLB:{worksheet}\\tPL:ILLUMINA' \
                {genome} {trimmed_fwd_FASTQ} {trimmed_rev_FASTQ} 2> {sample_name}.{lane}.bwa.log | \
                sed \'s/-R @RG.*//\' - | \
                samblaster -u {sample_name}.{lane}.unmapped.fastq 2> \
                {sample_name}.{lane}.duplicate.Metrics.txt | \
                {samtools} view -Sb - 2>{sample_name}.{lane}.samtoolsBamConvertErrors.log | \
                {samtools} sort -m 512M -@ 2 -T {sample_name}.{lane}.temp -O bam - > \
                {out_primary_bam}".format(
                    bwa               = config_dict['bwa'],
                    threads           = config_dict['bwa_threads'],
                    samtools          = config_dict['samtools'],
                    worksheet         = worksheet,
                    flowcell          = flowcell,
                    lane              = lane[0],
                    trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                    trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                    date              = date_time_iso,
                    sample_name       = sample_name,
                    chemistry         = chemistry,
                    genome            = config_dict['bwa_genome_fasta'],
                    out_primary_bam   = out_primary_bam)
                    
    return align_cmd


def bbmap_align_FASTQ_for_recalibration_matrix(config_dict,trimmed_fwd_FASTQ,trimmed_rev_FASTQ,out_primary_bam,root_name):
    """
    Run a quick bbmap to get matrix file for the run 
    """
    #readConfigFile(config_file)
    
    import datetime
    i = datetime.datetime.now()
    date_time_iso = "%s" % i.isoformat()
    
    lane = re.search(r"(L00[1-8])",trimmed_fwd_FASTQ)
    #lane = '%s' % re.findall("(L00[1-8])",trimmed_fwd_FASTQ)
    if lane is not None:
        lane = list(filter(None,lane.groups()))
    else:
        lane = 'merged'
    original_FASTQ = re.sub(r".qfilter.fastq.gz",".fastq.gz",trimmed_fwd_FASTQ)
    path           = os.path.realpath("%s" % original_FASTQ)
    sample_name    = re.sub(r".qfilter.fastq.gz","",trimmed_fwd_FASTQ[:-24])
    sample_name    = re.sub(r"_L00[1-8]","",sample_name)
    #sample_name_for_bam = re.sub(r"%s/" % args.named_directory,"",sample_name[:-29])
    
    
    #maxindel=200k
    align_cmd = "{bbmap} in={trimmed_fwd_FASTQ} in2={trimmed_rev_FASTQ} ref={genome} maxindel=2k out={name}.bam t=4".format(
                    bbmap             = config_dict['bbmap'],
                    trimmed_fwd_FASTQ = trimmed_fwd_FASTQ,
                    trimmed_rev_FASTQ = trimmed_rev_FASTQ,
                    genome            = config_dict['bwa_genome_fasta'],
                    name              = root_name)
    
 
                    
    return align_cmd
    

##################
# this was code for merging bams, removed after I decided to merge FASTQs instead
#def regex_for_merging():
    """
    A method to merge if needed, but don't if only one available
    """
    with open('mergeinfo.temp', 'r+') as inputfile:
        for line in inputfile:
            line = line.strip().split('\t')
            if line[0] == "merge=yes":
                regex = "([^/]+)_L00[1-2].bwa.drm.bam$"
                #regex = "([^/]+)_[1-8].bwa.drm.bam$" # possible not needed for non-illumina data
                print (regex)
            if line[0] == "merge=no":
                regex = "([^/]+).bwa.drm.bam$"
                
                
    return regex
  
#@follows(align_FASTQ)
#@collate("*.bwa.drm.bam", formatter(regex_for_merging()),"{path[0]}/.bwa.drm.merged.bam")
#@collate(["*bwa.drm.bam"], formatter("([^/]+)_L00[12].bwa.drm.bam$"),"{1[0]}.bwa.drm.merged.bam")
#@collate("*.bwa.drm.bam", formatter(regex_for_merging()),"{path[0]}/{1[0]}.bwa.drm.merged.bam") # needed for non-merge operation
def merge_lane_bams(infiles,outfile):
    """
    3 - Use samtools to merge the bam from any number of lanes
    """
    print(infiles,'-->',outfile)
    name = infiles[0]
    
    with open('mergeinfo.temp') as f:
        merge = first_line = f.readline().strip()
    
    if merge == "merge=yes":
        print ("Actually run a merge operation")
        os.system(cmd.merge_bams(config_dict = config_dict,
                             bams        = infiles,
                             merged_bam  = outfile,
                             name        = name,
                             cur_dir     = os.getcwd()))
    if merge == "merge=no":
        print ("Just rename files to include merged in the outfile name")
        os.system("mv {infile} {outfile}".format(infile=infiles[0],outfile=outfile))
    
    #name = outfile[:-19]        


####NOTNEEDED@transform(["*.bwa.drm.realn.recal.bam"],suffix(".bwa.drm.realn.recal.bam"),".mean_depth_data")
def create_mean_depth_data(infile,outfile):
    """
    Uses the depth calculator from Bedtools and calculates the depth for each base
    """
    print (infile,'-->',outfile)
    coverageBed="/opt/programs/bedtools-2.29.0/bedtools2/bin/coverageBed"
    bed = 'FullTargetFileUsedforCoverage_2.sorted.bed'
    # Generate mean coverage per exon
    # Use the nonamecheck to make it not break with GL0000 etc in the bam file.
    os.system(f"{coverageBed} -sorted -nonamecheck -mean -a {bed} -b {infile} > {outfile}")
    
    
#@follows(create_mean_depth_data)















################
def COPYannotate_gatk_vcf(infile,outfile,log_name,outtype):    
    #print(infile,'-->',outfile)
    """
    is --phased needed?
    """
    log_name = outfile[:-4]
    cmd = "/opt/programs/ensembl-vep-release-93.2/vep \
           --{outtype} \
           --buffer_size 10000 \
           --flag_pick \
           --exclude_predicted \
           --hgvsg \
           --hgvs \
           --dont_skip \
           --fork 32 \
           --force_overwrite \
           --offline \
           --assembly GRCh37 \
           --fasta /data/process/genomes/human/gatk_bundle/hs37d5_chr/hs37d5.fa \
           --db_version 97 \
           --cache_version 97 \
           --refseq  \
           --variant_class \
           --transcript_version \
           --protein \
           --symbol \
           --phased \
           --pubmed \
           --canonical \
           --allele_number \
           --dir /opt/programs/ensembl-vep-release-93.2/ \
           --dir_cache /opt/programs/ensembl-vep-97/ensembl-vep/ \
           --dir_plugins /opt/programs/ensembl-vep-97/ensembl-vep/Plugins/ \
           -custom /network/processed/vep_data/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
           -custom /network/processed/vep_data/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
           --plugin SingleLetterAA \
           -plugin CADD,/network/sequenced/genomic_resources/CADD/whole_genome_SNVs.tsv.gz,/network/sequenced/genomic_resources/CADD/InDels_inclAnno.tsv.gz \
           --plugin Blosum62 \
           --plugin ExACpLI,/network/sequenced/genomic_resources/ExAC/ExACpLI_values.txt \
           --plugin GeneSplicer,/opt/programs/GeneSplicer/bin/linux/genesplicer,/opt/programs/GeneSplicer/human/ \
           --plugin LoFtool,/network/sequenced/genomic_resources/ExAC/ExAC_LoFtool_scores.txt \
           --plugin SpliceRegion \
           --plugin MaxEntScan,/network/sequenced/genomic_resources/fordownload/ \
           --plugin MPC,/network/processed/vep_data/fordist_constraint_official_mpc_values_v2.txt.gz \
           --plugin FATHMM_MKL,/data/process/genomes/human/FATHM/fathmm-MKL_Current.tab.gz  \
           --plugin dbNSFP,/data/process/genomes/human/dbNSFPv3.5a/dbNSFP_hg19.gz,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Uniprot_acc_Polyphen2,\
           Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,\
           Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,\
           MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_UniprotID,MutationAssessor_variant,MutationAssessor_score,\
           MutationAssessor_score_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,\
           PROVEAN_pred,Transcript_id_VEST3,Transcript_var_VEST3,VEST3_score,VEST3_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,\
           MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,\
           MutPred_protID,MutPred_AAchange,MutPred_Top5features,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,\
           fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,Eigen_coding_or_noncoding,Eigen-raw,Eigen-phred,Eigen-PC-raw,\
           Eigen-PC-phred,Eigen-PC-raw_rankscore,GenoCanyon_score,GenoCanyon_score_rankscore,integrated_fitCons_score,integrated_fitCons_score_rankscore,\
           integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_score_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,\
           H1-hESC_fitCons_score_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_score_rankscore,HUVEC_confidence_value,GERP++_NR,GERP++_RS,\
           GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,\
           phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,SiPhy_29way_pi,\
           SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,\
           1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,\
           clinvar_rs,clinvar_clnsig,clinvar_trait,clinvar_golden_stars,Interpro_domain,GTEx_V6p_gene,GTEx_V6p_tissue \
           -i {input} \
           -o {output} \
           2> {log_name}.log".format(input=infile,output=outfile,log_name=log_name,outtype=outtype)
    return cmd
    




###############
#
# This function was removed from exome pipeline on 05-11-2020 and replaced with a seperate function for each chromsome,
# to allow Ruffus to handle the multi jobs instead.
#
#@subdivide(["*bwa.drm.realn.recal.bam"],formatter(r".bwa.drm.realn.recal.bam$"),
#("{basename[0]}.chr1.split.vcf.gz","{basename[0]}.chr2.split.vcf.gz","{basename[0]}.chr3.split.vcf.gz",
# "{basename[0]}.chr4.split.vcf.gz","{basename[0]}.chr5.split.vcf.gz","{basename[0]}.chr6.split.vcf.gz",
# "{basename[0]}.chr7.split.vcf.gz","{basename[0]}.chr8.split.vcf.gz","{basename[0]}.chr9.split.vcf.gz",
# "{basename[0]}.chr10.split.vcf.gz","{basename[0]}.chr11.split.vcf.gz","{basename[0]}.chr12.split.vcf.gz",
# "{basename[0]}.chr13.split.vcf.gz","{basename[0]}.chr14.split.vcf.gz","{basename[0]}.chr15.split.vcf.gz",
# "{basename[0]}.chr16.split.vcf.gz","{basename[0]}.chr17.split.vcf.gz","{basename[0]}.chr18.split.vcf.gz",
# "{basename[0]}.chr19.split.vcf.gz","{basename[0]}.chr20.split.vcf.gz","{basename[0]}.chr21.split.vcf.gz",
# "{basename[0]}.chr22.split.vcf.gz","{basename[0]}.chrX.split.vcf.gz","{basename[0]}.chrY.split.vcf.gz"))
def run_HaploTypeCaller4_vcf(infile,outfiles):
    """
    Use parallel, as a quick work around. Other possibilities:
    Create a seperate @transform function for each subdivide
    Create a seperate function somewhere else and exploit Ruffus oo to run it.
    changed this:
    ROI                   = config_dict['ROI_chr%s' % chrom],
    to:
    ROI                   = config_dict['ROI_%s' % chrom],
    """     
    # create a text file (GATKv4_HC_jobs) of all the commands, to give to parallel
    #print (infile,'-->',outfiles,'\n')
    ploidy = 2
    for sample in infile:
        with open("%s.GATKv4_HC_jobs" % infile,'w+') as outputfile:
            for i,outfilename in enumerate(outfiles):
                log_name   = outfiles[i][:-7]
                output_file_name = outfiles[i]
                output_file_name_root = outfiles[i][:-7]
                chrom = re.search("(.chr)([0-9a-zA-Z]+)(.)(.+)",outfilename).group(2)
                command = cmd.HaploTypeCaller_v4_single_vcf(config_dict = config_dict,
                                                            ROI                   = config_dict['ROI_%s' % chrom],
                                                            bam                   = infile,
                                                            output_file_name      = output_file_name,
                                                            output_file_name_root = output_file_name_root,
                                                            ploidy                = ploidy,
                                                            log_name              = log_name  )
                outputfile.write(command)
                
    
    sample_count  = len(glob("*GATK*"))                                                                #calculate number of jobs to run
    parallel_jobs =  int(int(100) / int(sample_count))
    os.system("parallel -j {jobs} :::: {name}.GATKv4_HC_jobs".format(jobs=4,name=infile))  #run HC in parallel, given text file of commands: 1 sample per set of contigs
    
    
#################    
#################    
# TSO500
#################
#################    




#@collate("*.bwa.drm.sorted.bam", formatter(regex_bamfilter),"{FILESTRING[0]}.bwa.drm.primary.bam")
def merge_or_rename_bam(infiles,outfile):
    """
    If sequencing Exomes, we might put 32 samples across two lanes. We therefore need to merge the bam.
    If using Xp workflow, we need to keep them seperate. Will need a different pipeline for cancer+ exome (get data from SampleSheet)
    By naming the output file the same, means I can use all the same following functions.
    """
    #print(infiles[0],'+',infiles[1],'-->',outfile)
    
    if args.merge_bams == 'T':
        #print("merging...")
        os.system(cmd.sambamba_merge(config_dict = config_dict,
                                 threads     = 6,  
                                 merged_bam  = outfile,
                                 in_bams     = infiles,
                                 log_name    = outfile))
    else:
        #THis won't work!!! If the two seperate lanes have different names then they won;t be picked by this regex?
        print ("renaming...")
        #print("mv )    