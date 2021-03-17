#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/4/21

@author: moinSebi

"""
import paftv
import sys
import itertools

"""
TODO: 
- start on zero (or the latest connection)
- iterate one ahead - > check if the same
- if not -> iterate one futher
- check sequence difference to that one 
- optimize minium 

"""


def readGFA(file_gfa: str):

    paths = dict()
    nodes = dict()
    with open(file_gfa) as f:
        for line in f.readlines():
            if line.startswith("P"):
                lsplit = line.split()
                if paths.get(lsplit[1]) == None:
                    paths[lsplit[1]] = lsplit[2].split(",")
            if line.startswith("S"):
                lsplit = line.split()
                nodes[lsplit[1]]  = len(lsplit[2])
    return paths


def find1(dict1):
    o = []
    o2 = []
    for k,v in dict1:
        v2 = []
        o.append(v)
        for x in v:
            if v.count(x) < 0:
               v2.append(x)

        o2.append(v2)

    j2 = dict()
    for i, j in itertools.combinations(o2,2):
        j2[(i,j)] = set(set(i) - set(j))

    for i,j in itertools.combinations(o,2)
        i1 = 0
        i2 = 0



if __name__ == "__main__":
    readGFA(sys.argv[1])
