import os
import pysam
import argparse
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

            if 'D' in read.cigarstring:
                matching_reads = stored_reads[read_name]
                for matching_read in matching_reads:
                    outfile.write(matching_read)

    # Convert filtered SAM to BAM
    pysam.view('-bS', '-o', output_bam, output_sam, catch_stdout=False)

    # Remove temporary SAM files
    os.remove(input_sam)
    os.remove(output_sam)







def main():
    parser = argparse.ArgumentParser(description="Filter BAM file based on CIGAR string")
    parser.add_argument("input_bam", help="Path to input BAM file")
    parser.add_argument("output_bam", help="Path to output BAM file")
    parser.add_argument("-w", "--window_size", type=int, default=100000, help="Sliding window size (default: 100000)")

    args = parser.parse_args()

    filter_bam_file(args.input_bam, args.output_bam, window_size=args.window_size)

if __name__ == "__main__":
    main()