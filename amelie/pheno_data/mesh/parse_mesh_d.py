#! /usr/bin/env python3

import sys, os
from collections import defaultdict


def print_record(rec):
    assert len(rec['MH']) == 1, rec['MH']
    assert len(rec['UI']) == 1, rec['UI']
    print('\t'.join([rec['UI'][0][0], rec['MH'][0][0]]))
    for sy in rec['ENTRY']:
        print('\t'.join([rec['UI'][0][0], sy[0]]))

def main():
    last_record = defaultdict(lambda: [])
    with open(sys.argv[1]) as f:
        for line in f:
            if len(line.strip()) == 0:
                print_record(last_record)
            elif line.strip() == '*NEWRECORD':
                last_record = defaultdict(lambda: [])
            else:
                # print(line)
                header = line.split(' = ')[0].strip()
                value = line.split(' = ')[1].strip().split('|')
                last_record[header].append(value)
    print_record(last_record)


if __name__ == "__main__":
    main()
