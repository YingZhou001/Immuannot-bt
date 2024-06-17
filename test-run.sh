
scriptdir=Immuannot-bt/script
refdir=ref-imgt.2023-06-23
outpref=test-out

asm=Immuannot-bt/test-data/sel-asm.fa.gz
reads=Immuannot-bt/test-data/sel-reads.fa.gz

if [ ${asm} == unset ]
then
  echo "target assembly input is missing"
else 
  mm2thread=3
  echo "[msg:] testing assembly-based annotation"
  time bash ${scriptdir}/immuannot-bt.sh -a ${asm} -d ${refdir} -o ${outpref} -t ${mm2thread}
fi

if [ ${reads} == unset ]
then
  echo "target reads input is missing"
else
  mm2thread=10
  echo "[msg:] testing annotating HiFi reads for V(D)J events"
  time bash ${scriptdir}/immuannot-bt.sh -r ${reads} -d ${refdir} -o ${outpref} -t ${mm2thread}
fi
