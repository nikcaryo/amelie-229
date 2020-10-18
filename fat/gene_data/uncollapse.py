#! /usr/bin/env python3

import sys
import codecs
import io

if __name__ == "__main__":
    delimiter = sys.argv[1]
    column = int(sys.argv[2]) - 1
    input_stream = io.TextIOWrapper(sys.stdin.buffer, encoding='ISO-8859-1')
    for line in input_stream:
        line = line.strip().split('\t')
        splitCol = line[column]
        split = splitCol.strip().split(delimiter)
        for s in split:
            newline = [w for w in line]
            newline[column] = s
            print('\t'.join(newline))
