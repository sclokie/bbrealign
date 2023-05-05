import sys
import re
from collections import deque

def passes_conditions(column_5, column_6):
    h_or_s_count = column_6.count('H') + column_6.count('S')
    m_count = column_6.count('M')
    numbers = [int(num) for num in re.findall(r'\d+', column_6)]

    return column_5 > 50 and h_or_s_count == 1 and m_count == 1 and all(num > 5 for num in numbers)

def main():
    stored_data = {}
    keys_queue = deque(maxlen=1000)

    for line in sys.stdin:
        line = line.strip()

        # Print header lines to stdout and continue to the next line
        if line.startswith("@SQ") or line.startswith("@RG") or line.startswith("@PG"):
            print(line)
            continue

        columns = line.split('\t')

        if len(columns) >= 6:
            key = columns[0]
            column_5 = int(columns[4])
            column_6 = columns[5]

            if key not in stored_data:
                stored_data[key] = []
                keys_queue.append(key)

                if len(keys_queue) == 1000:
                    oldest_key = keys_queue.popleft()
                    if oldest_key in stored_data:
                        del stored_data[oldest_key]

            stored_data[key].append(columns)

            if passes_conditions(column_5, column_6):
                matching_lines = stored_data[key]
                for matching_line in matching_lines:
                    print('\t'.join(matching_line))

if __name__ == '__main__':
    main()
