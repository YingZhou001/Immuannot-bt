import sys,gzip
import pprint as pp

pref = sys.argv[1]

dat = {'V': [], 'D' : [], 'J': []}

hdr = ''
seq = ''

for line in sys.stdin :
    if not line : break
    if line[0] == '>' :
        if hdr :
            dat[hdr[4]].append([hdr, seq])
            hdr = ''
            seq = ''
        hdr = line
    else:
        seq += line



for region in dat:
    ofile = pref + '.'+ region + '.fa.gz'
    with gzip.open(ofile, 'wt') as fp : 
        for x in dat[region] :
            fp.write(''.join(x))
