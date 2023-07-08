import sys
import pysam
from pybedtools import BedTool


def filter_sam_by_neighboring_reads(input_sam, output_sam, min_distance=50):
    # Read input SAM file
    with pysam.AlignmentFile(input_sam, "r", check_sq=False) as samfile:
        # Create a list of intervals for each read
        intervals = []
        for read in samfile:
            if not read.is_unmapped:
                intervals.append(
                    [samfile.get_reference_name(read.reference_id), read.reference_start, read.reference_end])

        # Create a BedTool object and merge intervals that are closer than min_distance
        bedtool = BedTool(intervals).sort().merge(d=min_distance - 1)
        merged_intervals = {f"{iv[0]}:{iv[1]}-{iv[2]}" for iv in bedtool}

        # Filter input SAM file by checking if the read is within any of the merged intervals
        with pysam.AlignmentFile(output_sam, "wh", header=samfile.header) as outfile:
            for read in samfile:
                if not read.is_unmapped:
                    read_interval = f"{samfile.get_reference_name(read.reference_id)}:{read.reference_start}-{read.reference_end}"
                    if any(read_interval in merged_interval for merged_interval in merged_intervals):
                        outfile.write(read)


if __name__ == "__main__":
    input_sam_file = sys.argv[1]
    output_sam_file = sys.argv[2]
    min_distance = int(sys.argv[3])

    filter_sam_by_neighboring_reads(input_sam_file, output_sam_file, min_distance)
