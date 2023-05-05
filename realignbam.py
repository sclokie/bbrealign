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
    temp_dir = args.tempdir
    max_reads = args.max_reads

    # Open input BAM
    if input_bam_file == '-':
        input_bam = pysam.AlignmentFile(sys.stdin.buffer, "rb")
    else:
        input_bam = pysam.AlignmentFile(input_bam_file, "rb")

    # Create a temporary file for storing unsorted output BAM
    with tempfile.NamedTemporaryFile(mode='wb', delete=False, dir=temp_dir) as tmp_file:
        tmp_output_bam_file = tmp_file.name

    # Write output BAM
    with pysam.AlignmentFile(tmp_output_bam_file, "wb", header=input_bam.header) as output_bam:
        read_counter = 0
        read_buffer = []
        for read in input_bam.fetch():
            read_buffer.append(read)
            read_counter += 1

            if read_counter >= max_reads:
                filtered_reads = filter_sam(read_buffer, nm_threshold)
                for filtered_read in filtered_reads:
                    output_bam.write(filtered_read)
                read_buffer = []
                read_counter = 0

        # Process the remaining reads in the buffer
        if read_buffer:
            filtered_reads = filter_sam(read_buffer, nm_threshold)
            for filtered_read in filtered_reads:
                output_bam.write(filtered_read)

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
    parser.add_argument("-t", "--tempdir", default=None, help="Temporary directory for storing intermediate files (default: system temporary directory)")
    parser.add_argument("-m", "--max_reads", type=int, default=100000, help="Maximum number of reads to process per chunk (default: 100000)")

    args = parser.parse_args()
    main(args)
