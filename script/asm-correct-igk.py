import sys, gzip
import pprint as pp

igklist = []

with open(sys.argv[1], 'rt') as fp:
    while True:
        line = fp.readline()
        if not line : break
        gene = line.replace('\n', '')
        if gene not in igklist :
            igklist.append(gene)

igkdct = {}
with gzip.open(sys.argv[2], 'rt') as fp:
    while True:
        line = fp.readline()
        if not line : break
        llst = line.split('\t')
        gene = llst[3][1:].split('*')[0]
        if gene in igklist :
            ctg = llst[0]
            strand = llst[3][0]
            tag = '.'
            if 'D-' in gene : tag = 'D'
            else : tag = 'N'
            if ctg not in igkdct : igkdct[ctg] = {}
            key = gene.replace('D-', '-')
            if key not in igkdct[ctg]:
                igkdct[ctg][key] = {}
            igkdct[ctg][key][strand] = tag


out = {}
for ctg in igkdct :
    x = igkdct[ctg]
    out[ctg] = {}
    out[ctg]['+'] = [0, 0] # count for D, N
    out[ctg]['-'] = [0, 0] # count for D, N
    for key in x :
        for strand in x[key] :
            if x[key][strand] == 'D' : out[ctg][strand][0] += 1
            else : out[ctg][strand][1] += 1

for ctg in out :
    x = out[ctg]
    tag = ''
    for strand in x :
        if x[strand][0] > x[strand][1] : tag = 'D'
        else : tag = 'N'
        out[ctg][strand] = tag

with gzip.open(sys.argv[2], 'rt') as fp:
    while True:
        line = fp.readline()
        if not line : break
        llst = line.split('\t')
        gene = llst[3][1:].split('*')[0]
        allele = llst[3].split('|')[0]
        allele = allele.split('_or_')[0]

        if gene in igklist :
            ctg = llst[0]
            strand = llst[3][0]
            tag = out[ctg][strand]
            if 'D-' in gene and tag == 'N' :
                allele = allele.replace('D-', '-')
            elif 'D-' not in gene and tag == 'D' :
                allele = allele[1:]
                allele = allele.replace('-', 'D-')
                allele = strand + allele

        llst[3] = allele
        line = '\t'.join(llst)
        sys.stdout.write(line)
