import sys
import re
from collections import deque
import pysam
from typing import List

def passes_conditions(column_5, column_6, nm, nm_threshold):
    if column_6 is None:
        return False

    h_or_s_count = column_6.count('H') + column_6.count('S')
    m_count = column_6.count('M')

    if nm is None:
        return column_5 > 50 and h_or_s_count == 1 and m_count == 1
    else:
        return column_5 > 50 and h_or_s_count == 1 and m_count == 1 and nm <= nm_threshold



def filter_sam(bam: pysam.AlignmentFile, nm_threshold: int = 2) -> List[pysam.AlignedSegment]:
    stored_data = {}
    keys_queue = deque(maxlen=1000)
    filtered_reads = []

    for read in bam.fetch():
        line = read.to_string()
        columns = line.split('\t')

        if len(columns) >= 6:
            key = columns[0]
            column_5 = read.mapping_quality
            column_6 = read.cigarstring
            nm = read.get_tag('NM') if read.has_tag('NM') else None

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
                    if matching_read.cigarstring and matching_read.cigarstring.count('M') == 1:
                        filtered_reads.append(matching_read)

    return filtered_reads


