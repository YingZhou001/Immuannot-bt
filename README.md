# Immuannot-bt pipeline introduction

Another important immunological regions in human are B cell receptor and T cell
receptor gene regions. With immuannot-bt, we are able to annotate the **V**ariable, 
**D**iversity and **J**oining segments and the flanking RSS regions for human assembly. 
Additionally, immuannot-bt can also be used to pick up the V(D)J recombination 
events from HiFi reads.


A manuscript is in preparing.

Last update date : 6/11/2024

Content
--------

- [Detection Strategy](#detection-strategy)
- [Installation](#installation)
    - [Download](#download)
    - [Dependency](#dependency)
- [Input and Output](#inputs-and-outputs)
- [An Example](#an-example)
- [Limitation](#limitation)
- [Todo list](#todo-list)

[\[top\]](#content)

# Detection Strategy

The major difficulty for annotating V, D, and J coding regions from whole genome
assembly is the short length of these genes especially for the D regions (~10-20 bp) 
and the J regions (~50bp).

| Segment | Length | +flanking area |
| --- | --- | --- |
| V | ~300bp | ~500bp |
| D | ~10-20bp | ~80bp |
| J | ~50bp | ~100bp |

One solution is to perform alignment to the restricted region that extracted
from whole genome seqeunce, together with careful evaluation to get the most
reliable annotation on those coding segments. (cite Digger)
The other way is to extend the coding seqeunces with flanking RSS regions so
that extra ~50-200bp will be added which directly increases the alignment
specitivity.
Immuannot-bt employs the second way, by generating the V,D,J reference
sequences from IMGT/LIGM-DB, it can achieve high efficiency and accuracy at the
same time when annotating human assembly.

[\[top\]](#content)

# Installation

## Download

1) Download Immuannot-bt:

```bash
git clone https://github.com/YingZhou001/Immuannot-bt.git
```

2) Download IMGT reference data set

```bash
wget https://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgt.dat.Z
```

3) Prepare the reference data set for Immuannot-bt

```bash
scriptdir=Immuannot-bt/script
refpref=ref-imgt
bash ${scriptdir}/fmt.sh ${dat} ${refpref}
## reference will be outputing into a folder named as ${refpref}.${date}
```

[\[top\]](#content)

## Dependency

Immuannot-bt depends on:

* TRUST4 [REF1]
* minimap2 (2.27-r1193) 
* seqtk (1.4-r130-dirty)

All source codes and pre-build binaries can be found in folder 'dep-tools', 
if any pre-build binary does work, you may need to re-build it
and copy it to the folder 'dep-tools/bin'.

[\[top\]](#content)

# Input and output

Immuannot-bt takes fasta or fastq format as input, and generate a gff file for
assembly-based annotation and a csv file in [AIRR format](https://docs.airr-community.org/en/v1.3.1/datarep/rearrangements.html#) for picking up V(D)J recombination from reads.

[\[top\]](#content)

# An example

```bash
scriptdir=Immuannot-bt/script
refdir=ref-imgt.2023-06-23 # need to modify this to the lastet IMGT database
outpref=test-out

asm=Immuannot-bt/test-data/sel-asm.fa.gz
reads=Immuannot-bt/test-data/sel-reads.fa.gz

mm2thread=3
echo "[msg:] testing assembly-based annotation"
time bash ${scriptdir}/immuannot-bt.sh -a ${asm} -d ${refdir} -o ${outpref} -t ${mm2thread}

echo "[msg:] testing annotating HiFi reads for V(D)J events"
time bash ${scriptdir}/immuannot-bt.sh -r ${reads} -d ${refdir} -o ${outpref} -t ${mm2thread}
```

[\[top\]](#content)

# Limitation

[\[top\]](#content)

# Todo list

[\[top\]](#content)
