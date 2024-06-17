#downloading

scriptdir=Immuannot-bt/script
dat=imgt.dat.Z
refpref=ref-imgt


if true
then
  echo "[msg:] git clone immuannot-bt"
  git clone https://github.com/YingZhou001/Immuannot-bt.git

  echo "[msg: ~10mins] download imgt reference database"
  wget https://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgt.dat.Z

  echo "[msg:] prepare for the reference data set"
  bash ${scriptdir}/fmt.sh ${dat} ${refpref}

fi
