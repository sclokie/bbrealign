import sys
import re
from collections import deque

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
