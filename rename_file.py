#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 1/15/21

@author: moinSebi

"""
import sys

def readandrename(filename):

    with open(filename) as f:

        for x in f.readlines():
            if x.startswith(">"):
                l = x.split()
                print(l[0] + "_1")
            else:
                print(x.replace("\n", ""))



if __name__ == "__main__":
    readandrename(sys.argv[1])