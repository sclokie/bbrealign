import argparse
import sys
import os
import tempfile
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

    # Create a temporary file for storing unsorted output BAM
    with tempfile.NamedTemporaryFile(mode='wb', delete=False) as tmp_file:
        tmp_output_bam_file = tmp_file.name

    # Write output BAM
    with pysam.AlignmentFile(tmp_output_bam_file, "wb", header=input_bam.header) as output_bam:
        for read in filtered_reads:
            output_bam.write(read)

    input_bam.close()

    # Sort and index the BAM file
    pysam.sort("-o", output_bam_file, tmp_output_bam_file)
    pysam.index(output_bam_file)

    # Clean up the temporary file
    os.remove(tmp_output_bam_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, help="Input BAM file, use '-' for stdin")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file, use '-' for stdout")
    parser.add_argument("-n", "--nm", type=int, default=1, help="NM threshold (default: 1)")

    args = parser.parse_args()
    main(args)
