# paftv

This is a fork form https://github.com/weigelworld/minitv.  
It is now now possible to use your own alignment in the PAF format. 

Requirements: 
- biopython: pip install biopython

Usage: 
- git clone https://github.com/MoinSebi/paftv.git
- paftv.py -h for help message
- fixed flags
    - -r/--region --> Now displays all connected alignments (and their children)
- added flags
    - -p paf file
    - -q/--query query files
    - -X for allvsall alignment  
    - --maxiterations: Maximal number iterations if --region is selected
    - -c/--color: Colors for min,mid and max identity. Hex format without "#"
    - -g/--group_list: Add this flag with file to a group.csv file. Check 
    - -o output file
    
    
    
Upload the JSON here: https://alitvteam.github.io/AliTV/d3/AliTV.html
 

