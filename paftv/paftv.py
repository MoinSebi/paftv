#!/usr/bin/python3


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import sys
from paftv.lib import alitv, filter
#from collections import defaultdict
#import itertools
#import json
#import logging
from enum import Enum
from pathlib import Path
from Bio import SeqIO
import argparse


"""
TODO: 
- Option to automatic translocations detection? 


"""



class Alignment:
    """
    A single alignment region between two genomes.
    This class implements values required to be converted to AliTV Link objects.
    """

    def __init__(self, query_name, query_start, query_end,
                 target_name, target_start, target_end,
                 identity, strand, length):
        """
        :param query_name: Query sequence name
        :type query_name: str
        :param query_start: Query start (0-based)
        :type query_start: int
        :param query_end: Query end (0-based)
        :type query_end: int
        :param target_name: Target sequence name
        :type target_name: str
        :param target_start: Target start (0-based)
        :type target_start: int
        :param target_end: Target end (0-based)
        :type target_end: int
        :param identity: Percent identity of a sequence, represented as a float between 0 and 100.
        :type identity: float
        :param strand: Relative strand
        :type strand: Strand
        :param length: Alignment length
        :type length: int
        """

        self.query_name = query_name
        self.query_start = query_start
        self.query_end = query_end

        self.target_name = target_name
        self.target_start = target_start
        self.target_end = target_end


        self.identity = identity
        if isinstance(strand, Strand):
            self.strand = strand
        else:
            raise ValueError(
                "{} is not a valid Alignment.strand value. Alignment objects need to be instantiated with a member of"
                " minitv.align.Strand.".format(type(strand)))
        self.length = length

    def __len__(self):
        return self.length

class PAFAlignment(Alignment):
    """
    A single alignment between two genomes in PAF format.
    """

    def __init__(self, query_name, query_len, query_start, query_end, strand, target_name, target_len, target_start,
                 target_end, matches, length, mapping_quality, id, tags):
        """
        :param query_name: Query sequence name
        :type query_name: str
        :param query_len: Query sequence length
        :type query_len: int
        :param query_start: Query start (0-based)
        :type query_start: int
        :param query_end: Query end (0-based)
        :type query_end: int
        :param strand: Relative strand: "+" or "-"
        :type strand: str
        :param target_name: Target sequence name
        :type target_name: str
        :param target_len: Target sequence length
        :type target_len: int
        :param target_start: Target start on original strand (0-based)
        :type target_start: int
        :param target_end: Target end on original strand (0-based)
        :type target_end: int
        :param matches: Number of residue matches
        :type matches: int
        :param length: Alignment block length
        :type length: int
        :param mapping_quality: Mapping quality (0-255; 255 for missing)
        :type mapping_quality: int
        :param tags: A dictionary of SAM-like typed key-value pairs
        :type tags: dict
        """

        if strand == '+':
            strand = Strand.FORWARD
        elif strand == '-':
            strand = Strand.REVERSE
        else:
            raise ValueError("{} is an invalid PAF alignment strand value.".format(strand))


        ## This is the identity
        if id != 0:
            identity = id
        else:
            identity = (matches / length) * 100
        super().__init__(query_name, query_start, query_end, target_name, target_start, target_end, identity,
                         strand, length)

        self.query_len = query_len
        self.target_len = target_len
        self.matches = matches

        if mapping_quality < 0 or mapping_quality > 255:
            raise ValueError("'{}' is an invalid mapping quality for the PAF format. Mapping quality must be between 0 "
                             "and 255.".format(mapping_quality))
        else:
            self.mapping_quality = mapping_quality

        self._tags = tags

    def get_tag(self, key):
        return self._tags[key]

    def set_tag(self, key, value):
        self._tags[key] = value


    # Compare
    def __eq__(self, other):
        if not isinstance(other, PAFAlignment):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.query_end == other.query_end and \
               self.query_start == other.query_start and \
                self.target_start == other.target_start and \
                self.target_end == other.target_end and \
                self.target_name == other.target_name and \
                self.query_name == other.query_name and \
                self.target_name == other.target_name


def parse_output_minimap2(output, id):
    """
    :param output: minimap2 stdout
    :type output: io.TextIOBase or io.StringIO
    :rtype: list[Alignment]
    """

    alignments = []
    for line in output:

        split_line = line.rstrip().split("\t")
        query_name = split_line[0]
        query_len = int(split_line[1])
        query_start = int(split_line[2])
        query_end = int(split_line[3])
        strand = split_line[4]
        target_name = split_line[5]
        target_len = int(split_line[6])
        target_start = int(split_line[7])
        target_end = int(split_line[8])
        matches = int(split_line[9])
        length = int(split_line[10])
        mapping_quality = int(split_line[11])
        id_numb = 0
        if id:
            if len(split_line) > 12 and len(set([x for x in split_line[12:] if x.startswith("id")])):
                id_numb = float(split_line[12].split(":")[-1])
            else:
                print("Ignoring id tag - not existing")
                id = False
        tags = {}

        if len(split_line) > 12:
            # parse SAM tags
            key_value_pairs = split_line[12:]
            for key_value_pair in key_value_pairs:
                key, type_char, value = key_value_pair.split(':')
                if type_char == "i":
                    tags[key] = int(value)
                elif type_char == 'f':
                    tags[key] = float(value)
                else:
                    tags[key] = value

        new_alignment = PAFAlignment(query_name, query_len, query_start, query_end, strand, target_name,
                                     target_len, target_start, target_end, matches, length, mapping_quality, id_numb, tags)
        alignments.append(new_alignment)

    return alignments

def parseGroups(group_file, xornot):
    """

    :param group_file:

    name, fastaid1, fastfa2
    :return:
    """

    data = dict()
    count = 0
    with open(group_file) as f:
        for lines in f.readlines():
            lsplit = [x.replace("\n","").replace(" ", "") for x in lines.split(",")]
            if xornot:
                p = [x + "_1" for x in lsplit]
                data[count] = lsplit + p
            else:
                data[count] = lsplit
            count += 1
    o = dict()
    for k,v in data.items():
        for x in v:
            o[x] = k
    return o, data



class FASTAFile:
    """
    Represents a FASTA file of one or more sequences.
    """

    def __init__(self, path):
        """
        :param path: path to the FASTA file
        :type path: pathlib.Path
        """

        self.path = path
        self.name = path.stem

        self.records = list(SeqIO.parse(str(path), 'fasta'))

    def __len__(self):
        """
        :return: number of records in the FASTA file
        """

        return len(self.records)


class Strand(Enum):
    """ File schema independent method to encode alignment strand. """
    FORWARD = 1
    UNKNOWN = 0
    REVERSE = -1

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    main()






def main():

    ## ARGUMENT PARSING

    parser = argparse.ArgumentParser();
    parser.add_argument(
        "-p",
        "--paf",
        metavar="paf",
        help="Alignment file in paf format",
        required=True
    )

    parser.add_argument(
        "-q",
        "--query",
        metavar="query.fa",
        nargs='+',
        help="FASTA file of sequence(s). This argument can instead take a single file of file names containing "
             "paths relative to the working directory.",
        required=True
    )

    # parser.add_argument(
    #     "-r",
    #     "--region",
    #     help="Bounds of region to view, in standard samtools region format (refSeqID:start-end)"
    # )
    #
    # parser.add_argument(
    #     "--min_link_identity",
    #     help="Minimum link identity for AliTV to initially display",
    #     type=float,
    #     default=0.0
    # )
    #
    # parser.add_argument(
    #     "--min_link_length",
    #     help="Minimum link length for AliTV to initially display",
    #     type=int,
    #     default=0
    # )
    #
    # parser.add_argument(
    #     "--min_aln_cov",
    #     help="Minimum percent coverage, from 0 to 100, by alignments that a sequence needs to be initially visible in "
    #          "AliTV.",
    #     type=float,
    #     default=0.0
    # )
    #
    # parser.add_argument(
    #     "--min_ref_cov",
    #     help="Minimum percent coverage, from 0 to 100, of a query sequence by alignments to the reference needed for"
    #          "the query sequence to be initially visible in AliTV.",
    #     type=float,
    #     default=0.0
    # )

    parser.add_argument(
        "-o",
        "--output",
        help="Output name",
        required=True,
    )

    parser.add_argument(
        "-X",
        "--allvsall",
        help="Query and the reference are the same (all-vs-all alignment) - only one genome [default: off]",
        # required=True,
        action="store_true"

    )
    parser.add_argument(
        "-c",
        "--color",
        help="Min, mid and max identity colors in hex format (without #) example 'D21414,FFEE05,1DAD0A (also possible in AliTV)",
        type=str
    )

    # parser.add_argument(
    #     "-t",
    #     "--transposon",
    #     help = "Only show 'jumping' alignments",
    #     action="store_true"
    # )

    # parser.add_argument(
    #     "--maxiteration",
    #     help = "Number of maximal iteration if region linking is used (only works with -r, ignored otherwise)",
    #     type = int
    # )

    # parser.add_argument(
    #     "-g",
    #     "--group_list",
    #     help = "CSV format - each 'group' (e.g. same chromosome) in one line. fasta_entry1, fasta_entry2, fasta_entry3 --> Check examples"
    # )

    parser.add_argument(
        "--mmid",
        help="Minimum and maximum identity range (e.g. 0,100) [default: dynamic min/max of all alignments]"
    )

    args = parser.parse_args()

    # Make a fasta file list
    paths = tuple([Path(x) for x in args.query])
    fasta_files = [FASTAFile(x) for x in paths]

    alignment_groups = {}  # key = Tuple(target, query)
    alignment_lists = {}  # This one saves the alignments in the same style as alignment groups, but holds a list of strings (lines) as values
    paf = args.paf

    if args.allvsall:
        if len(fasta_files) != 1:
            print("Only provide one fasta file with the -X flag")
            sys.exit()
        else:
            p = FASTAFile(Path(args.query[0]))
            p.name = p.name + "_1"
            for x in p.records:
                x.name = x.name + "_1"
                x.id = x.id + "_1"
            fasta_files.append(p)

            with open(paf) as f:

                for line1 in f.readlines():
                    split_line = line1.rstrip().split("\t")
                    split_line[0] = split_line[0] + "_1"
                    line2 = "\t".join(split_line)
                    # print(alignment_lists.get((split_line[5], split_line[0])))
                    if alignment_lists.get((fasta_files[0].name, fasta_files[1].name)) != None:
                        # print(line1.rstrip(),
                        alignment_lists[fasta_files[0].name, fasta_files[1].name].append(line2.rstrip())
                    else:
                        # print(len(line1))
                        alignment_lists[fasta_files[0].name, fasta_files[1].name] = [line2.rstrip()]
            for k, v in alignment_lists.items():
                alignment_groups[k] = parse_output_minimap2(v, False)



    # If you have two distict input files
    else:
        '''
        # This is to find the right FASTA file afterwards
        Dict[entry] --> Fasta file

        '''

        entry2fasta = dict()
        entrycheck = []
        for x in fasta_files:
            for y in x.records:
                entrycheck.append(x)
                entry2fasta[y.id] = x.name

        # Check if we have double entry -> I
        if len(entrycheck) != len(set(entrycheck)):
            print("Duplicated fasta entries")

        all_combinations = set()
        with open(paf) as f:

            for line1 in f.readlines():
                split_line = line1.rstrip().split("\t")
                # print(alignment_lists.get((split_line[5], split_line[0])))
                if entry2fasta.get(split_line[0]) != None:
                    if entry2fasta.get(split_line[5]) != None:
                        if entry2fasta[split_line[0]] != entry2fasta[split_line[5]]:

                            query = entry2fasta[split_line[0]]
                            target = entry2fasta[split_line[5]]

                            all_combinations.add((target, query))

                            # Sort out self alignments
                            if query == target:
                                continue;
                            else:
                                if alignment_lists.get((target, query)) != None:
                                    alignment_lists[(target, query)].append(line1.rstrip())
                                else:
                                    alignment_lists[(target, query)] = [line1.rstrip()]

        for k, v in alignment_lists.items():
            alignment_groups[k] = parse_output_minimap2(v, False)

    filter.reduceAlg(alignment_groups)
    ## If you define a region, make a new alignment
    new_alg = {}
    count = 0

    # if args.region != None:
    #     region = [[args.region.split(":")[-2], int(args.region.split(":")[-1].split("-")[0]),
    #                int(args.region.split(":")[-1].split("-")[1])]]
    #     if args.maxiteration == None:
    #         new_alg = filter.filter_region(alignment_groups, region)
    #     else:
    #         new_alg = filter.filter_region(alignment_groups, region, args.maxiteration)
    #     alignment_groups = new_alg

    # if args.group_list != None:
    #     new_alg = {}
    #     data, gen = parseGroups(args.group_list, args.allvsall)
    #     new_alg = filter.filter_transpon(gen, data, alignment_groups)
    #     alignment_groups = new_alg

    # We calculate the max and min identity
    minId = 100
    maxId = 0
    if args.mmid != None:
        minId = int(args.mmid.split(",")[0])
        maxId = int(args.mmid.split(",")[1])
    else:
        for k, v in alignment_groups.items():
            for x in v:
                if minId > x.identity:
                    minId = x.identity
                if maxId < x.identity:
                    maxId = x.identity

    # -----------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------

    # Try new functions here:
    # filter.overlaps(alignment_groups)

    """ 
    This is Alitv - try not to touch this
    Changes:
        - min, max ID in AliTv (maybe change to function)
        - changeColors function 
        - If regions given - use new alignment group  #
        - Removed Order and orient sequence 
    """
    ali_tv = alitv.AliTV(fasta_files, minId, maxId)

    if args.color != None:
        ali_tv.changeColors(args.color.split(","))

    ali_tv.load_links(alignment_groups)

    # ali_tv.set_soft_filters(args.min_link_identity, args.min_link_length)

    # if args.min_aln_cov:
    #     print('Hiding chromosomes by alignment coverage')
    #     ali_tv.hide_chromosomes_by_alignment_coverage(args.min_aln_cov)

    # if args.min_ref_cov:
    #     print('Hiding chromosomes by reference alignment coverage')
    #     ali_tv.hide_chromosomes_by_reference_coverage(args.min_ref_cov)

    # ali_tv.order_and_orient_sequences(alignment_groups)

    ali_tv.optimize_configuration()

    # logging.info('Generating JSON output')

    # Leave it like it is
    a = ali_tv.get_json(indent=1)

    with open(args.output, "w") as l:
        l.write(a)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
