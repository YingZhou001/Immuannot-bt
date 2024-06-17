import sys,re

overlap_cut = float(sys.argv[1])
match_cut = float(sys.argv[2])

def map_overlap(cigar_string):
    regex = re.compile(r'(\d+?)([A-Z])')
    ret = regex.findall(cigar_string)
    tot = 0
    match = 0
    if ret :
        for l,tag in ret :
            if tag in ['N', 'I'] :
                tot += float(l)
            if tag == 'M' :
                tot += float(l)
                match += float(l)
        return(match/tot)
    else :
        return(0)

for line in sys.stdin:
    if not line: break
    llst = line.split('\t')
    ctg = llst[0]
    v = llst[3]
    d = llst[4]
    j = llst[5]
    v_cigar = llst[7]
    j_cigar = llst[9]
    v_overlap = 0
    j_overlap = 0
    if v_cigar : v_overlap = map_overlap(v_cigar)
    if j_cigar : j_overlap = map_overlap(j_cigar)
    c1 = v_overlap > overlap_cut and j_overlap > overlap_cut

    v_match_rate = 0
    j_match_rate = 0
    if llst[18] :
        v_match_rate = float(llst[18])
    if llst[19] :
        j_match_rate = float(llst[19])
    c2 = v_match_rate > match_cut and j_match_rate > match_cut

    if c1 and c2 : 
            sys.stdout.write(line)
