# paftv


This is a fork form https://github.com/weigelworld/minitv.  
It is now now possible to use your own alignment in the PAF format. 

**Requirements:** 
- python3.5 or higher
- biopython

```bash
git clone https://github.com/MoinSebi/paftv
pip install biopython
```

**Script only - no installation needed**


### usage 

#### Basic + Help: 
```bash
minimap2 seq1.fasta seq2.fasta > aln.paf
paftv.py -p aln.paf -q seq1.fasta seq2.fasta -o out.json
```

```bash
paftv.py -h 
usage: paftv.py [-h] -p paf -q query.fa [query.fa ...] [-r REGION]
                [--min_link_identity MIN_LINK_IDENTITY] [--min_link_length MIN_LINK_LENGTH]
                [--min_aln_cov MIN_ALN_COV] [--min_ref_cov MIN_REF_COV] -o OUTPUT [-X] [-c COLOR] [-i]
                [--maxiteration MAXITERATION] [-g GROUP_LIST]

optional arguments:
  -h, --help            show this help message and exit
  -p paf, --paf paf     Alignment file in paf format
  -q query.fa [query.fa ...], --query query.fa [query.fa ...]
                        FASTA file of sequence(s). This argument can instead take a single file of file
                        names containing paths relative to the working directory.
  -r REGION, --region REGION
                        Bounds of region to view, in standard samtools region format (refSeqID:start-
                        end)
  --min_link_identity MIN_LINK_IDENTITY
                        Minimum link identity for AliTV to initially display
  --min_link_length MIN_LINK_LENGTH
                        Minimum link length for AliTV to initially display
  --min_aln_cov MIN_ALN_COV
                        Minimum percent coverage, from 0 to 100, by alignments that a sequence needs to
                        be initially visible in AliTV.
  --min_ref_cov MIN_REF_COV
                        Minimum percent coverage, from 0 to 100, of a query sequence by alignments to
                        the reference needed forthe query sequence to be initially visible in AliTV.
  -o OUTPUT, --output OUTPUT
                        Output name
  -X, --allvsall        Query and the reference are the same (all-vs-all alignment) - only one genome
  -c COLOR, --color COLOR
                        Name min, mid and max identity colors in hex format (without #) example
                        'D21414,FFEE05,1DAD0A
  -i, --identity        Use the identity flag in the paf format [default: off] (normally computed by
                        alignment length)
  --maxiteration MAXITERATION
                        Number of maximal iteration if region linking is used (only works with -r,
                        ignored otherwise)
  -g GROUP_LIST, --group_list GROUP_LIST
                        CSV format - each 'group' (e.g. same chromosome) in one line. fasta_entry1,
                        fasta_entry2, fasta_entry3 --> Check examples

```


#### Example 1 - All-vs-all 
```bash
cat AAA.fasta AAB.fasta YCT.fasta CAA.fasta > seq.fasta
minimap2 seq.fasta seq.fasta -X > aln.paf
paftv.py -p aln.paf -X -q seq.fasta -o out.json
```

#### Example 2 - All-vs-all but with a subsetted view
Comment: You are only interested in AAA.fasta and AAB.fasta (or you need a smaller file)
```bash
cat AAA.fasta AAB.fasta YCT.fasta CAA.fasta > seq.fasta
minimap2 seq.fasta seq.fasta -X > aln.paf
paftv.py -p aln.paf -q AAA.fasta AAB.fasta -o out.json
```

#### Example 3 - Check a specific region
```bash
minimap2 AAA_CAA.cat.fasta YCT_AAB.cat.fasta  > aln.paf
paftv.py -p aln.paf -q AAA_CAA.cat.fasta YCT_AAB.cat.fasta  -o out.json -r AAA_6:10000-20000 --maxiteration 5 
```


#### Exmaple 4 - Only show moving alignments (to other chromosomes)
```bash
minimap2 AAA_CAA.cat.fasta YCT_AAB.cat.fasta  > aln.paf
paftv.py -p aln.paf -q AAA_CAA.cat.fasta YCT_AAB.cat.fasta  -o out.json -g group_example.txt
```


    

### Vizualization: 
    
Upload your JSON here: [https://alitvteam.github.io/AliTV/d3/AliTV.html](https://alitvteam.github.io/AliTV/d3/AliTV.html)
 

