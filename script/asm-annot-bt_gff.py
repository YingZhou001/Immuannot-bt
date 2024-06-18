import sys, gzip, re

import pprint as pp
import paftools as pt
import algntools as at
from datetime import date


genefile = sys.argv[1]
metafile = sys.argv[2]


# load gene map
gene_indexmap = {}
meta_info = ''
with open(genefile, 'rt') as fp:
    while True:
        line = fp.readline()
        if not line: break
        if line[0] == '#' :
            meta_info += line
        else :
            gene,gene_id,trans_id = line.split()
            gene_indexmap[gene] = [gene_id,0,trans_id,0]


# load meta dictionary
def split_line(line) :
    allele, meta = line[1:].split('|')
    llst = meta.replace('\n', '').split( )
    imgt_ori = llst[0]
    imgt_type = llst[1]
    ret = {}
    for x in llst :
        if '=' in x : 
            tag,value = x.split('=')
            # value into 0-based
            if '..' in value :
                a,b = value.split('..')
                value = a + '..' + b
            ret[tag] = value
    ret['imgt_allele'] = allele
    ret['imgt_type'] = imgt_type
    return(ret)

meta_dct = {}
seq_buf = ''
with gzip.open(metafile, 'rt') as fp:
    while True:
        line = fp.readline()
        if not line : break
        if line[0] == '>' :
            if seq_buf :
                meta_dct[tag]['seq'] = seq_buf
                seq_buf = ''
            tag = line[1:].split()[0]
            meta_dct[tag] = split_line(line)
        else :
            seq_buf += line.replace('\n', '')

if not seq_buf :
    meta_dct[tag]['seq'] = seq_buf



def extract_seq(seq, sel_rg, full_rg, strand) :
    #1-based coordinate
    frg0, frg1 = full_rg.split('..')
    srg0, srg1 = sel_rg.split('..')
    frg0 = int(frg0)
    frg1 = int(frg1)
    srg0 = int(srg0)
    srg1 = int(srg1)
    if strand == '-' :
        srg0, srg1 = [frg1 - srg1, frg1 - srg0 + 1]
    else :
        srg0, srg1 = [srg0 - frg0, srg1 - frg0 + 1]
    #print(srg0, srg1, file=sys.stderr)
    return(seq[srg0:srg1])

def check_v_region(ret_rec) :
    regex = re.compile(r'seq "(.*?)"')
    l1 = ''
    l2 = ''
    cys1 = ''
    cys2 = ''
    v = ''
    for x in ret_rec :
        tag = x[2]
        attr = x[8]
        if tag == 'L-PART1' :
            ret = regex.search(attr)
            if ret :
                l1 = ret.group(1).upper()
            else:
                l1 = ''
        elif tag == 'L-PART2' :
             ret = regex.search(attr)
             if ret :
                 l2 = ret.group(1).upper()
             else:
                l2 = ''
        elif tag == 'V-REGION' :
            ret = regex.search(attr)
            if ret :
                v = ret.group(1).upper()
            else:
                v = ''
        elif tag == '1st-CYS' :
            ret = regex.search(attr)
            if ret :
                cys1 = ret.group(1).upper()
            else :
                cys1 = ''
        elif tag == '2st-CYS' :
            ret = regex.search(attr)
            if ret :
                cys2 = ret.group(1).upper()
            else :
                cys2 = ''

    seq = l1 + l2 + v
    is_atg = True if seq[0:3] == 'ATG' else False
    is_cys = True if at.translate(cys1) == 'C' and at.translate(cys2) == 'C' else False
    is_early_stop = False
    prot=''
    frame_l2 = '.'
    frame_v = '.'
    if is_atg :
        prot = at.translate(seq[0:int(len(seq)/3)*3], 1)
        if prot.find('_')  >= 0 :
            is_early_stop = True
        frame_l2 = len(l1) % 3 + 1
        frame_v = len(l1 + l2) % 3 + 1
    
    #print(len(seq), prot.find('_'), len(prot), is_atg, seq[0:10], prot, file=sys.stderr)
    new_ret_rec = []
    for x in ret_rec :
        tag = x[2]
        new_ret_rec.append(x)
        if tag == 'L-PART2' :
            new_ret_rec[-1][7] = str(frame_l2)
        if tag == 'V-REGION' :
            new_ret_rec[-1][7] = str(frame_v)
        if tag == 'gene' :
            if is_atg and is_early_stop :
                new_ret_rec[-1][5] = 'P'
        if False and tag in ["L-PART2"] :
            print(frame_l2, frame_v, is_atg and is_early_stop, file=sys.stderr)
            print(new_ret_rec[-1], file=sys.stderr)
    #pp.pprint(new_ret_rec, stream=sys.stderr)
    return(new_ret_rec)



# gtf headline, TBD
headline = '##gff-version 2\n'
headline += "##data_generation_date " +  str(date.today()) + "\n"
headline += meta_info


source = "Immuannot-ig"
out = []
for line in sys.stdin :
    pafline = line.split('@')[1]
    pafrec = pt.parsing_paf(pafline) # output 1-based coordinate
    #print(pafrec)
    ctg = pafrec['tstr']
    strand = pafrec['strand']
    tfro = pafrec['tfro']
    tto = pafrec['tto']
    template_allele = pafrec['qstr']
    qfro = pafrec['qfro']
    qto = pafrec['qto']
    qlen = pafrec['qlen']
    cg = pafrec['cg']
    cs = pafrec['cs']
    mapping_rate = (float(qto)+1 - float(qfro))/float(qlen)
    adjNM = pt.adj_nm(cs)

    template_rec = meta_dct[template_allele]
    tseq = at.recoverTargetSeqFromQuery(template_rec['seq'], cs,
                                        qfro+'..'+qto, strand)

    template_alelle = template_rec['imgt_allele']
    template_type = template_rec['imgt_type']
    gene_name = template_alelle.split('*')[0]
    tmp_map = gene_indexmap[gene_name]
    gene_id = tmp_map[0]
    if tmp_map[1] > 0 : gene_id += '.'+str(tmp_map[1])
    tmp_map[1] += 1
    trans_id = tmp_map[2]
    if tmp_map[3] > 0 : trans_id += '.'+str(tmp_map[3])
    tmp_map[3] += 1

    out_type = '.'
    frame = '.'

    attr0 = ['gene_id "'+ gene_id +'"', 'gene_name "' + gene_name]
    attr = attr0.copy()
    attr.append('template_source ' + template_allele.split('|')[1])
    attr.append('mapping_rate ' + str(round(mapping_rate, 2)))
    attr.append('adjust_nm ' + str(adjNM))

    gene_rec = [ctg, source, 'gene', tfro, tto, 
                '.', strand, frame, ' ; '.join(attr)]

    attr1 = ['gene_id "'+ gene_id +'"', 'transcript_id "' + trans_id + '"', 
             'gene_name "' + gene_name + '"']
    attr = attr1.copy()
    attr.append('imgt_allele "' + template_alelle + '"')
    attr.append('imgt_type "' + template_type + '"')
    trans_rec = [ctg, source, 'transcript', tfro, tto, 
                 '.', strand, frame, ' ; '.join(attr)]

    ret_rec = [gene_rec, trans_rec]
    for tag in template_rec :
        value = template_rec[tag]
        if '..' in value :
            msg, rg = at.intervalMapQry2Ctg(value, cg, strand, qfro, qto, tfro, tto)
            sel_seq = extract_seq(tseq, rg, tfro + '..' + tto, strand)
            fro, to = rg.split('..')
            attr = attr1.copy()
            if msg :
                attr.append('warning "'+msg+'"')
            attr.append('seq "'+sel_seq+'"')
            ret_rec.append([ctg, source, tag, fro, to, 
                            '.', strand, frame, '; '.join(attr)])
    if gene_name[3] == 'V' : 
        ret_rec = check_v_region(ret_rec)
    ostr = ''
    for x in ret_rec :
        ostr += '\t'.join(x) + '\n'
    out.append([ctg, int(tfro), int(tto), ostr])

out = sorted(out, key=lambda x: x[1])
out = sorted(out, key=lambda x: x[0])

sys.stdout.write(headline)
for x in out:
    sys.stdout.write(x[3])
