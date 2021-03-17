#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/17/21

@author: moinSebi

"""

import sys

def readPaf(paf_file, header_pos):
    paf_file_new = []

    with open(paf_file) as f:
        for line in f.readlines():
            lsplit = line.split()
            if lsplit[0] in header_pos and lsplit[5] in header_pos:
                paf_file_new.append(line)
    for x in paf_file_new:
        print(x.replace("\n", ""))




def readFasta(fasta_file):
    entries = []
    for x in fasta_file:
        with open(x) as f:
            for line in f.readlines():
                if line.startswith(">"):
                    entries.append(line[1:].replace("\n", ""))
    return  set(entries)

if __name__ == "__main__":
    j = readFasta(sys.argv[2:])
    readPaf(sys.argv[1], j)
