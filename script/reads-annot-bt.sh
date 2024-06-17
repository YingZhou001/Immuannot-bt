

READS=$1
REF=$2
OUTPREF=$3
THREAD=$4

scriptpath=$(dirname $0)

# dev opt
run_minimap2=true
run_seqtk=true
run_trust4=true

options="-c -k11 -n2 -w5 -U1,1000000 -g100 -A1 -B4 -O6,26 -E2,1 -s20 -f100000"
v_ref=${REF}/human.V.fa.gz
d_ref=${REF}/human.D.fa.gz
j_ref=${REF}/human.J.fa.gz
vdj_ref=${REF}/human.imgt+c.fa.gz

paf=${OUTPREF}/reads.paf.gz
tmp_v_annot=${OUTPREF}/tmp.reads.v_annot.gz
tmp_bed=${OUTPREF}/tmp.reads.vdj.bed
tmp_reads=${OUTPREF}/tmp.reads.vdj.fa.gz
tmp_vdj_annot=${OUTPREF}/tmp.reads.vdj_annot.tsv.gz
out_annot=${OUTPREF}.annot.tsv.gz

if ${run_minimap2}
then
  echo "map V segments to the reads"
  ${MINIMAP2} -t ${THREAD} ${options} ${v_ref} ${READS} | gzip -c >  ${paf}
  zcat ${paf} | python3 ${scriptpath}/reads-annot-bt_paf.py | gzip -c > ${tmp_v_annot}
  zcat ${tmp_v_annot} | python3 ${scriptpath}/reads-gen-vdj_bed.py > ${tmp_bed}
fi

if ${run_seqtk}
then
  echo "extract informative reads"
  ${SEQTK} subseq ${READS} ${tmp_bed} | ${SEQTK} seq -a | gzip -c > ${tmp_reads}
fi

if ${run_trust4}
then
  echo "annotate reads"
  opts="--fasta -t ${THREAD} --needReverseComplement --noImpute --outputFormat 1"

  ${TRUST4_annotator} -f ${vdj_ref} -a ${tmp_reads} ${opts} \
    | gzip -c > ${tmp_vdj_annot}
  zcat ${tmp_vdj_annot} | python3 ${scriptpath}/reads-vdj-filter.py \
    | gzip -c > ${out_annot}
fi
