import sys


out = {}
pairs = []

for line in sys.stdin:
    if not line: break
    if line[0] == '>' :
        alleles = line[1:].split('|')[0].split('_or_')
        pair = []
        for allele in alleles :
            gene = allele.split('*')[0]
            tag = gene[0:4]
            if tag not in out :
                out[tag] = []
            if gene not in out[tag] :
                 out[tag].append(gene)
            if len(alleles) >=2 :
                pair.append(gene)
        if pair :
            pairs.append(pair)

for tag in out :
    out[tag] = sorted(out[tag])

#print(pairs)

out_order = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", 'IGLV', "IGLJ",
             "TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", 
             "TRDV", "TRDD", "TRDJ", "TRGV", "TRGJ"]

index = 0
for tag in out_order :
    for gene in out[tag] :
        index +=1
        print(gene + '\tIGG' + str(index).zfill(6) + '\tIGT' + str(index).zfill(6))
