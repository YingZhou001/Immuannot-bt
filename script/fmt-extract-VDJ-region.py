import sys, re
import pprint as pp

def rev_complement(seq) :
    rc_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
              'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    ret = []
    for x in seq[::-1] :
        if x in rc_map :
            ret.append(rc_map[x])
        else :
            ret.append(x)
    return(''.join(ret))

def read_id(line: str) :
    line[0:2] == 'AC'
    regex = re.compile(r'AC\s+(.*);')
    ret = regex.search(line)
    sid = ''
    if ret :
        sid = ret.group(1)
    return(sid)

def read_seq(buf: str) :
    llst = buf.split('\n')
    tag = False
    seq = []
    for line in llst :
        if tag :
            seq.append("".join(filter(lambda x: x.isalpha(), line)))
        if line[0:2] == 'SQ' : tag = True
    return("".join(seq))

def find_allele(annot: str):
    regex = re.compile(r'IMGT_allele="(.*?)"')
    ret = regex.search(annot)
    if ret :
        return(ret.group(1))
    else :
        return('')

def find_func(annot: str) :
    if 'pseudo' in annot : return('P')
    elif 'ORF' in annot : return('ORF')
    elif 'functional' in annot : return('F')
    else :
        return('U')

def find_codon_start(annot: str) :
    regex = re.compile(r'codon_start=(\d+)')
    ret = regex.search(annot)
    if ret :
        return(ret.group(1))
    else :
        return('')

def check_complete_data(arr : list) :
    # can also shrink the interval size base on the elements it includes
    fro = sys.maxsize
    to = 0
    for x in arr[5] :
        if fro > x[1] : fro = x[1]
        if to < x[2] : to = x[2]
    tag = False
    allele_tag = True if arr[2] != '.' else False
    seg_length = to - fro + 1

    if arr[0] == "D-GENE" and seg_length > 5 : tag = True
    elif arr[0] == "V-GENE" and seg_length > 100 : tag = True
    elif arr[0] == "J-GENE" and seg_length > 15 : tag = True


    if not allele_tag :
        tag = False
    if tag :
        return([True, [fro, to]])
    else :
        return([False, []])


def read_FT(buf: str) :
    lines = buf.split('\n')
    key1s = ["V-GENE", "D-GENE", "J-GENE"]
    key2s = ["L-PART1", "L-PART2",
             "V-HEPTAMER", "V-SPACER", "V-NONAMER",
             "1st-CYS", "CONSERVED-TRP", "2nd-CYS",
             "J-NONAMER", "J-SPACER", "J-HEPTAMER",
             "5'D-NONAMER", "5'D-SPACER", "5'D-HEPTAMER",
             "3'D-NONAMER", "3'D-SPACER", "3'D-HEPTAMER"] #,
    key3s = ["V-REGION","J-REGION", "D-REGION"]
    L = len(lines)
    i = 0
    out = []
    regex = re.compile(r'<*(\d+)..>*(\d+)')
    while i < L :
        if 'FT' in lines[i] :
            ret = regex.search(lines[i])
            key = lines[i][3:21].split()
            if key : key = key[0]
            else : 
                i += 1
                continue
            rg = []
            annot = []
            rc_tag = '+'

            if ret :
                rg = [int(ret.group(1)), int(ret.group(2))]
                if 'complement' in lines[i] : rc_tag = '-'
            else :
                llst = lines[i].split()

            if key in key1s :
                i += 1
                #ret = regex.search(lines[i])
                while not lines[i][3:21].split() :
                    annot.append(lines[i])
                    i += 1
                i -= 1
                annot_str = ','.join(annot)
                attr = []
                allele = find_allele(annot_str)
                func = find_func(annot_str)
                out.append([key, rg, rc_tag, allele, func, attr])

            if out and key in key3s :
                #print(out[-1], file=sys.stderr)
                last_rg = out[-1][1]
                last_allele = out[-1][3]
                if rg[0] >= last_rg[0] and rg[1] <= last_rg[1] :
                    out[-1][5].append([key, rg[0], rg[1]])
                    if not last_allele :
                        i += 1
                        while not lines[i][3:21].split():
                            last_allele = find_allele(lines[i])
                            if last_allele : break
                            i += 1
                        i -= 1
                        if last_allele :
                            out[-1][3] = last_allele

        i += 1
    #out = sorted(out, key = lambda x : int(x[1][0]))
    out2 = []
    for x in out :
        is_complete, rg = check_complete_data(x)
        if is_complete :
            x[1] = rg
            out2.append(x)
    return(out2)

def split_fasta(seq: str, cut: int) :
    ret = ''
    i = 0
    while i < len(seq) :
        ret += seq[i:(i+cut)] + '\n'
        i += cut
    return(ret)

buf = ''
for line in sys.stdin :
    if not line : break
    if line[0:2] == 'AC' :
        sid = read_id(line).replace(' ', '')
    buf += line
    if '//' in line :
        # processing buf
        seq = read_seq(buf)
        meta = read_FT(buf)
        for x in meta :
            total_length = x[1][1] - x[1][0] + 1
            rc_tag = x[2]
            ostr = '>' + x[3].replace(' ', '_')
            ostr += '|' + sid + ':' + str(x[1][0]) + '-' + str(x[1][1])
            ostr += ' ' + x[4]
            fro = x[1][0]-1
            to = x[1][1]+1
            if rc_tag == '+' :
                x[5] = sorted(x[5], key = lambda y : y[1])
                for y in x[5] :
                    ostr += ' ' + y[0] + '=' + str(y[1]-fro) 
                    ostr += '..' + str(y[2]-fro)
            else :
                x[5] = sorted(x[5], key = lambda y : y[1], reverse=True)
                for y in x[5] :
                    ostr += ' ' + y[0] + '=' + str(to - y[2])
                    ostr += '..' + str(to - y[1])
            ostr += ' tot_length=' + str(total_length)
            ostr += '\n'
            if not x[2] or not sid or not x[3] or not x[4] or not x[5]:
                sys.stderr.write(ostr)
                print(x, file=sys.stderr)
                #sys.stdout.write(buf)
                #exit(-1)
            else :
                oseq = seq[x[1][0]-1:x[1][1]]
                if rc_tag == '-' :
                    ostr += split_fasta(rev_complement(oseq), 80)
                else :
                    ostr += split_fasta(oseq, 80)
                sys.stdout.write(ostr)
        buf = ''
