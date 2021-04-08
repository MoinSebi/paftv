#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/4/21

@author: moinSebi

"""
import sys

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
    return paths, nodes

def plus(s):
    return "+" == s[-1]

def nodePos(paths: dict(), nodes = dict()):
    pos = 0
    n2 = dict()
    for k, v in paths.items():
        pos = 0
        n  = dict()
        n2[k] = n

        for x in v:
            if n.get(x[:-1]) == None:
                n[x[:-1]] = [[plus(x), pos, pos+int(nodes[x[:-1]])]]
            else:
                n[x[:-1]].append([plus(x), pos, pos+int(nodes[x[:-1]])])
            pos += nodes[x[:-1]]
    return n2

def writePaf()




if __name__ == "__main__":
    print(sys.argv[1])
    p, n = readGFA(sys.argv[1])
    o = nodePos(p, n)
    print(len(o))

