# paftv


This is a fork form https://github.com/weigelworld/minitv.  

**Changes**:  
> It is now possible to use your own alignment in the PAF format. 

**Requirements:** 
- python3.5 or higher
- biopython

**Installation**  
Source: 
```
git clone https://github.com/MoinSebi/paftv
pip install biopython
python3 setup.py install 
paftv -h 
```



### Usage 

#### Basic + Help: 
```
minimap2 seq1.fasta seq2.fasta > aln.paf
paftv.py -p aln.paf -q seq1.fasta seq2.fasta -o out.json
```

```
usage: paftv.py [-h] -p paf -q query.fa [query.fa ...] -o OUTPUT [-X] [-c COLOR] [--mmid MMID]

optional arguments:
  -h, --help            show this help message and exit
  -p paf, --paf paf     Alignment file in paf format
  -q query.fa [query.fa ...], --query query.fa [query.fa ...]
                        FASTA file of sequence(s). This argument can instead take a single file of file names containing paths relative to the working directory.
  -o OUTPUT, --output OUTPUT
                        Output name
  -X, --allvsall        Query and the reference are the same (all-vs-all alignment) - only one genome [default: off]
  -c COLOR, --color COLOR
                        Min, mid and max identity colors in hex format (without #) example 'D21414,FFEE05,1DAD0A (also possible in AliTV)
  --mmid MMID           Minimum and maximum identity range (e.g. 0,100) [default: dynamic min/max of all alignments]

```
#### Example 1 - Two genomes
```
minimap2 AAA.fasta AAB.fasta > aln.paf
paftv -p aln.paf -q AAA.fasta AAB.fasta -o out.json
```

#### Example 2 - All-vs-all 
```
cat AAA.fasta AAB.fasta YCT.fasta CAA.fasta > seq.fasta
minimap2 seq.fasta seq.fasta -X > aln2.paf
paftv -p aln.paf -X -q seq.fasta -o out.json
```

#### Example 3 - All-vs-all but only look at specific genomes
Comment: You are only interested in AAA.fasta and AAB.fasta (or you need a smaller file)
```
cat AAA.fasta AAB.fasta YCT.fasta CAA.fasta > seq.fasta
minimap2 seq.fasta seq.fasta -X > aln2.paf
paftv -p aln.paf -q AAA.fasta AAB.fasta -o out.json
```
  
 

#### Alternative mapping approach 
- Change parameters in pan-minimap2 
```
ls > examples/*fasta > genomes.txt
chmod +x scripts/pan-minimap2
./scripts/pan-minimap2 $(cat genomes.txt)  > output.paf
```

#### Example sequences from the [1011 genomes project](https://www.nature.com/articles/s41586-018-0030-5)
    

### Visualization: 
     
Upload your JSON here: [https://alitvteam.github.io/AliTV/d3/AliTV.html](https://alitvteam.github.io/AliTV/d3/AliTV.html)  
If the file is too big, this might crash you browser. 

**Example data in the git repository**

### Filter PAF file with [**fpa**](https://github.com/natir/fpa)  

Example:
 >  cat aln.paf | fpa drop -l 200 - > aln.bigger200.paf 
