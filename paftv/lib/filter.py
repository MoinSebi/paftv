#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2/2/21

@author: moinSebi

"""


def filter_region(alg_groups: dict(), region: list, maxit: int = 1000) -> dict():
    """
    :param alg_groups:
    :param region:
    :return:
    """
    r = region
    new_alg = dict()
    count = 0

    delete = []
    while (True):
        nothing_new = True
        p = []
        for key, value in alg_groups.items():
            if new_alg.get(key) == None:
                new_alg[key] = []
            for alignment in value:
                # print(type(x))
                for entry in region:
                    if entry[0] == alignment.query_name:
                        # print("hit1")
                        if getOverlap([alignment.query_start, alignment.query_end], entry):
                            j = True
                            for d in delete:
                                if d == alignment:
                                    j = False
                            if j:
                                p.append([alignment.target_name, alignment.target_start, alignment.target_end])
                                new_alg[key].append(alignment)
                                nothing_new = False
                                delete.append(alignment)


                    elif entry[0] == alignment.target_name:
                        # print("hit2")
                        if getOverlap([alignment.target_start, alignment.target_end], entry):
                            #print(alignment)
                            j = True
                            for d in delete:
                                if d == alignment:
                                    j = False
                            if j:
                                nothing_new = False
                                p.append([alignment.query_name, alignment.query_start, alignment.query_end])
                                new_alg[key].append(alignment)
                                delete.append(alignment)
        print("Iteration: {}".format(count))
        print("New Alignments: {}".format(len(p)))
        print("Total Alignments: {}".format(len(delete)))
        print()
        region = []
        region.extend(p)
        count += 1
        if nothing_new or count == maxit+1:
            break


    return new_alg

def reduceAlg(alg_groups: dict) -> dict:

    new_alg_temp = dict()

    for k, v in alg_groups.items():
        new_alg_temp[k] = sorted(v, key = lambda x: (x.query_name, x.target_name, x.target_start, x.target_end, x.query_start, x.query_end))






def getOverlap(start_end: list, region: list ) -> bool:
    return max(0, min(start_end[1], region[2]) - max(start_end[0], region[1]))

def filter_transpon(t_list: list, gen_dict: dict, alg_groups: dict) -> dict:

    new_alg = dict()

    for k, v in alg_groups.items():
        new_alg[k] = []
        for alignment in v:
            if (alignment.query_name == "TAIR10_Chr1_1"):
                #print("hjdahskjda")
                pass
            if alignment.target_name in list(gen_dict.keys()) and alignment.query_name in list(gen_dict.keys()):
                #print(alignment.target_name)
                if alignment.target_name not in t_list[gen_dict[alignment.query_name]]:
                    new_alg[k].append(alignment)
            else:
                new_alg[k].append(alignment)

    return new_alg


def getOverlap1(start_end: list, region: list ) -> bool:
    return max(0, min(start_end[1], region[1]) - max(start_end[0], region[0]))

def overlaps(alg_dict: dict()):
    """

    :param alg_dict:
    :return:
    """

    
    """"
    TODO:
     - This will only work with 2 genomes
     - Care of k in alg_dict
     - Maybe also to able to give list


     """

    # if overlap
    check = dict()

    for k,v in alg_dict.items():
        for x in v:
            if check.get((x.query_name, x.target_name)) == None:
                check[(x.query_name, x.target_name)]  = [[x.target_start, x.target_end]]
            else:
                check[(x.query_name, x.target_name)].append([x.target_start, x.target_end])

    for k, v in check.items():
        #print(k)
        #print(v)
        v2 = sorted(v, key = lambda x: x[0])
        check[k] = v2

    odict = dict()
    for k,v in check.items():
            print(len(v))
            print(v)
            overlapList = []
            for index in range(len(v)):
                if len(overlapList) > 0:
                    if getOverlap1(overlapList[-1], v[index]):
                        overlapList[-1] = [min(overlapList[-1][0], v[index][0]), max(overlapList[-1][1], v[index][1])]
                    else:
                        overlapList.append(v[index])
                else:
                    overlapList.append(v[index])
            odict[k] = overlapList

    sdict = dict()
    for k, v in odict.items():
        p = 0
        for x in v:
            p += x[1] -x[0]
        sdict[k] = p
        print(p)

    allfiles = set()
    for k, v in sdict.items():
        allfiles.add(k[0])
    for k, v in sdict.items():
        if k[0].startswith("10015"):
            print(k)
            print(v)

    maxes = dict()
    for x in allfiles:
        op = []
        for k,v in sdict.items():
            if k[0] == x:
                op.append((k, v))
        op = sorted(op, key = lambda d: d[1])
        maxes[x] = (op[-1][0][1], op[-1][1])
    #print(allfiles)











