import sys, gzip
import pprint as pp


def overlap(a,b, c,d) :
    # cmpare [a,b], and [c,d]
    if b < c or a > d :
        return(False)
    else :
        return(True)

def parsing_paf(line) :
    if not line : return({})
    llst = line.split('\t')
    ret = {}
    ret['ctg'] = llst[0]
    ret['ctg_len'] = int(llst[1])
    ret['ctg_fro'] = int(llst[2])
    ret['ctg_to'] = int(llst[3])
    ret['strand'] = llst[4]
    ret['allele'] = llst[5]
    ret['allele_len'] = int(llst[6])
    ret['allele_mlen'] = int(llst[8]) - int(llst[7])
    ret['nm'] = llst[12].replace('NM:i:', '')
    ret['pafline'] = line.replace('\n', '')
    # add filter here
    f1 = float(ret['allele_mlen'])/ret['allele_len'] > 0.97
    f2 = int(ret['nm']) <= 10
    f3 = float(ret['nm'])/ret['allele_mlen'] < 0.1
    if f1 and (f2 or f3) :
        return(ret)
    else :
        return({})

pafdct = {} # define as pafdct[ctg]=[[], [], []]


for line in sys.stdin:
    if not line: break
    retpaf = parsing_paf(line)
    if not retpaf : continue
    ctg = retpaf['ctg']
    if ctg not in pafdct :
        pafdct[ctg] = []
    if retpaf: 
        pafdct[ctg].append(retpaf)


# define clusters

clustdct = {}
for ctg in pafdct :
    if ctg not in clustdct :
        clustdct[ctg] = []
    sorted_pafdct = sorted(pafdct[ctg], key=lambda x: int(x['ctg_fro']))
    for paf in sorted_pafdct :
        fro = int(paf['ctg_fro'])
        to = int(paf['ctg_to'])
        if not clustdct[ctg] :
            clustdct[ctg].append([fro, to])
        else :
            L = len(clustdct[ctg])
            overlap_tag = 0
            for i in range(L) :
                old_fro, old_to = clustdct[ctg][i]
                if overlap(fro, to, old_fro, old_to) :
                    clustdct[ctg][i] = [min(fro, old_fro), max(to, old_to)]
                    overlap_tag = 1
                    break
            if overlap_tag == 0:
                clustdct[ctg].append([fro, to])

# search gene pattern


for ctg in clustdct :
    output = []
    for clust in clustdct[ctg] :
        newpaf = []
        for paf in pafdct[ctg] :
            fro = int(paf['ctg_fro'])
            to = int(paf['ctg_to'])
            if overlap(fro, to, clust[0], clust[1]) :
                newpaf.append(paf)
        # search the best allele for each cluster
        ## sorted paf
        newpaf = sorted(newpaf, key=lambda x: float(x['allele_mlen']),reverse=True)
        newpaf = sorted(newpaf, key=lambda x: float(x['nm']))
        #print(clust)
        #pp.pprint(newpaf[0])
        selpaf=newpaf[0]
        output.append([selpaf['strand']+selpaf['allele'], selpaf['ctg_fro'],
                       selpaf['ctg_to'], str(selpaf['ctg_len']),
                       str(selpaf['allele_mlen']*100/selpaf['allele_len']), 
                       str(selpaf['nm']), selpaf['pafline']])
    output = sorted(output, key=lambda x: int(x[2]))
    for x in output :
        ostr=ctg + '\t' + str(x[1]) + '\t' +  str(x[2])
        ostr+= '\t' + x[0] + '\t' + x[3]
        ostr+= '\tmapping_rate=' + x[4] + "%"
        ostr+= ',nm=' + x[5]
        ostr = ostr.replace('"', '')
        ostr+= '\t@' + x[6]
        print(ostr)


