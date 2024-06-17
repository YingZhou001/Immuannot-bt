import sys


for line in sys.stdin:
    if not line : break
    llst = line.split('\t')
    strand = llst[3][0]
    ctg = llst[0]
    fro = int(llst[1])
    to = int(llst[2])
    tot_len = int(llst[7])
    ostr = ctg
    if to >= tot_len - 50 or fro < 50 : continue
    if strand == '+' :
        ostr += '\t' + str(fro)
        to += 200
        if to >= tot_len : to = tot_len - 1
        ostr += '\t' + str(to)
    elif strand == '-' :
        if fro < 200 :
            ostr += '\t' + '0'
        else :
            fro -= 200
            ostr += '\t' + str(fro)
        ostr += '\t' + str(to)
    print(ostr + '\t' + llst[3].split('|')[0] + llst[5])

