import argparse
import sys
import pysam
from utilities import filter_sam

def main(args):
    input_bam_file = args.input
    output_bam_file = args.output
    nm_threshold = args.nm

    # Open input BAM
    if input_bam_file == '-':
        input_bam = pysam.AlignmentFile(sys.stdin.buffer, "rb")
    else:
        input_bam = pysam.AlignmentFile(input_bam_file, "rb")

    filtered_reads = filter_sam(input_bam, nm_threshold)

    # Write output BAM
    if output_bam_file == '-':
        output_bam = pysam.AlignmentFile(sys.stdout.buffer, "wb", header=input_bam.header)
    else:
        output_bam = pysam.AlignmentFile(output_bam_file, "wb", header=input_bam.header)

    for read in filtered_reads:
        output_bam.write(read)

    input_bam.close()
    output_bam.close()

    if output_bam_file != '-':
        sorted_bam_file = output_bam_file.replace(".bam", "_sorted.bam")
        pysam.sort("-o", sorted_bam_file, output_bam_file)
        pysam.index(sorted_bam_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, help="Input BAM file, use '-' for stdin")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file, use '-' for stdout")
    parser.add_argument("-n", "--nm", type=int, default=1, help="NM threshold (default: 1)")

    args = parser.parse_args()
    main(args)
