

ASM=$1
REF=$2
OUTPREF=$3
THREAD=$4

scriptpath=$(dirname $0)

#dev-options
run_mm2=true
run_annot=true

#running alignment

#ref_ctg=${REF}/human.template.alleles.fa.gz
ref_ctg=${REF}/human.all.alleles.fa.gz
out_paf=${OUTPREF}/asm.mm2.paf.gz

if ${run_mm2}; 
then
  N=1
  f=$((N * 100000))
  U=$((N * 10000000))
  N=$((N * 20))
  options="-t${THREAD} -c -k11 -n2 -w5 -g100 -A1 -B4 -O6,26 -E2,1 -s20 -f${f}
  -U1,${U} -N${N} --cs --end-bonus=10"
  ${MINIMAP2} ${options} ${ASM} ${ref_ctg}  | gzip -c > ${out_paf}

fi

#initial annotation

annot0=${OUTPREF}/tmp.asm.annot.gz
annot1=${OUTPREF}/asm.annot.gz
genelist=${OUTPREF}/IG.gene.list
gff=${OUTPREF}.gff.gz

if ${run_annot}
then
  zcat ${out_paf} | python3 ${scriptpath}/asm-annot-bt_paf.py | gzip -c > ${annot0}
  python3 ${scriptpath}/asm-correct-igk.py \
    ${scriptpath}/sup-data/igk.gene.list ${annot0} \
    | gzip -c > ${annot1}

  #generate gff

  zcat ${ref_ctg} \
    | python3 ${scriptpath}/add_gene_id.py > ${genelist}
  zcat ${annot0} \
    | python3 ${scriptpath}/asm-annot-bt_gff.py ${genelist} ${ref_ctg} \
    | gzip -c > ${gff}

fi
