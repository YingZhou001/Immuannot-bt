#!/bin/bash

SCRIPTPATH=$(dirname $0)
# Environment variable for tool location

export TRUST4_annotator=${SCRIPTPATH}/dep-tools/bin/annotator
export MINIMAP2=${SCRIPTPATH}/dep-tools/bin/minimap2
export SEQTK=${SCRIPTPATH}/dep-tools/bin/seqtk
export TRUST4_ref_formatter=${SCRIPTPATH}/dep-tools/bin/BuildImgtAnnot.pl

# Set some default values:
ASM=unset
READS=unset
REF=unset

# optional parameters
OUTPREF=immuannot-bt-out
THREAD=3


usage()
{
  echo "
  Usage: bash ${SCRIPTPATH}/immuannot-bt.sh  [OPTIONS] value
                           [ -a | --assembly  target assembly       ] 
                           [ -r | --reads  target reads       ] 
                           [ -d | --data-ref  references                          ] 
                           [ -o | --outpref output prefix (optional)            ] 
                           [ -t | --thread  num of thread (optional, default 3) ] 
                           "
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n immuannot-bt \
    -o a:r:d:o:t: \
    --long assembly:,reads:,data-ref:,outpref:,thread:\
    -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
  exit 1
fi

#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -a | --assembly)  ASM="$2"    ; shift 2 ;;
    -r | --reads)  READS="$2"    ; shift 2 ;;
    -d | --data-ref)  REF="$2"    ; shift 2 ;;
    -o | --outpref) OUTPREF="$2"   ; shift 2 ;;
    -t | --thread)  THREAD="$2"    ; shift 2 ;;
    --overlap)  OVERLAP="$2"    ; shift 2 ;;
    --diff)  DIFF="$2"    ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1."
       usage
       ;;
  esac
done

if [ $ASM == unset  &&  $READS == unset ]; then
  echo "Error: no aassembly or reads input"
  usage
fi

if [ $REF == unset ]; then
  echo "Error: reference data set is required."
  usage
fi

if [ $OUTPREF == unset ]; then
  echo "Error: output prefix is required."
  usage
else
  mkdir -p $OUTPREF
  if [ $? != "0" ]; then
    echo "${OUTPREF} (-o) is not valid."
    usage
  fi
fi

start_second=`date +%s`
start=`date +%D-%H:%M:%S`

echo "########################################"
echo "##Welcome###############################"
echo "########################################"
echo ""
echo "Starting time: ${start}"
echo "#####parameters:########################"
if [ $ASM != unset ]; then
  echo "ASSEMBLY(-a)         : $ASM"
fi
if [ $READS != unset ]; then
  echo "READS(-r)            : $READS"
fi
echo "REF(-d)              : $REF"
echo "OUTPREF(-o)          : $OUTPREF"
echo "THREAD(-t)           : $THREAD"
#echo "Parameters remaining are: $@"




if true; then
  if [ $ASM != unset ]; then
    echo ""
    echo "#### annotate B/T cell receptor genes from assemblies#############"
    bash ${SCRIPTPATH}/asm-annot-bt.sh ${ASM} ${REF} ${OUTPREF} ${THREAD}
    if [ $? != "0" ]; then
      echo "ERROR: Failed to annotate B/T cell receptor genes (assemblies)"
      exit 1
    else 
      echo "Successfully annotate B/T cell receptor genes (assemblies)"
    fi
  fi
fi


if true; then
  if [ $READS != unset ]; then
    echo ""
    echo "#### annotate germline B/T cell receptor genes from reads#############"
    bash ${SCRIPTPATH}/reads-annot-bt.sh ${READS} ${REF} ${OUTPREF} ${THREAD}
    if [ $? != "0" ]; then
      echo "ERROR: Failed to annotate B/T cell receptor genes (reads)"
      exit 1
    else 
      echo "Successfully annotate B/T cell receptor genes (reads)"
    fi
  fi
fi




end_second=`date +%s`
end=`date +%D-%H:%M:%S`
runtime=$((end_second - start_second))
echo "Ending time: ${end}, Wallclock time :${runtime} seconds"
echo "####done####"

