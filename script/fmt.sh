
DAT=$1
PREF=$2

scriptpath=$(dirname $0)

#download_link=https://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgt.dat.Z
#dat=dat.Z
#wget ${download_link} -O ${dat}

date_tag=$(stat ${DAT} | grep Modify | cut -f2 -d' ')
refdir=${PREF}.${date_tag}

mkdir -p ${refdir}
# prepare ref for assembly annotation
if true; then
  zcat ${DAT} | python3 ${scriptpath}/fmt-subdat.py \
    | python3 ${scriptpath}/fmt-extract-full-seq.py \
    | gzip -c > ${refdir}/human.template.alleles.fa.gz
  zcat ${DAT} | python3 ${scriptpath}/fmt-subdat.py \
    | python3 ${scriptpath}/fmt-extract-all-seq.py \
    | gzip -c > ${refdir}/human.all.alleles.fa.gz
  zcat ${DAT} | python3 ${scriptpath}/fmt-extract-all-seq.py \
    | gzip -c > ${refdir}/all.alleles.fa.gz
fi

#prepare for the reads based annotation
if true; then
  zcat ${DAT} | python3 ${scriptpath}/fmt-subdat.py \
    | python3 ${scriptpath}/fmt-extract-VDJ-region.py \
    | python3 ${scriptpath}/fmt-split-vdj.py ${refdir}/human
  perl ${scriptpath}/dep-tools/bin/BuildImgtAnnot.pl Homo_sapien \
    | gzip -c > ${refdir}/human.imgt+c.fa.gz
fi

