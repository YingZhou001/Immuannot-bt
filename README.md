# Immuannot-bt pipeline introduction

Another important immunological regions in human are B cell receptor and T cell
receptor gene regions. With immuannot-bt, we are able to annotate the **V**ariable, 
**D**iversity and **J**oining segments and the flanking RSS regions for human assembly. 
Additionally, immuannot-bt can also be used to pick up the V(D)J recombination 
events from HiFi reads.


A manuscript is in preparing.

Last update date : 6/17/2024

Content
--------

- [Detection Strategy](#detection-strategy)
- [Installation](#installation)
    - [Download](#download)
    - [Dependency](#dependency)
- [Input and Output](#input-and-output)
- [An Example](#an-example)
- [Limitation](#limitation)
- [Todo list](#todo-list)

# Detection Strategy

The major difficulty for annotating V, D, and J coding regions from whole genome
assembly is the short sizes of those genes, especially for the D coding region which is about ~10-20bp. 
One solution is to perform basepair-to-basepair alignment to the restricted region (
extracted from whole genome seqeunce), together with careful evaluation including the flanking sequences
to get the most reliable annotation on those coding segments.[1]
The other way is to extend the coding seqeunces with flanking RSS regions so
that extra ~50-200bp will be added which directly increases the alignment
specitivity.
Immuannot-bt employs the second way, by generating the V,D,J reference
sequences from IMGT/LIGM-DB, it can achieve high efficiency and accuracy at the
same time when annotating human assembly.

| Segment | Length | +flanking area |
| --- | --- | --- |
| V | ~300bp | ~500bp |
| D | ~10-20bp | ~80bp |
| J | ~50bp | ~100bp |

For long reads data, we use TRUST4's annotator to resolve the CDR3 sequences. [2]

[1] Lees, W.D., Saha, S., Yaari, G. and Watson, C.T., 2024. Digger: directed annotation of immunoglobulin and T cell receptor V, D, and J gene sequences and assemblies. Bioinformatics, 40(3), p.btae144. \
[2] Song, L., Cohen, D., Ouyang, Z., Cao, Y., Hu, X. and Liu, X.S., 2021. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nature methods, 18(6), pp.627-630.

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
dat=imgt.dat.Z
bash ${scriptdir}/fmt.sh ${dat} ${refpref}
## reference will be outputing into a folder named as ${refpref}.${date}
```

[\[top\]](#content)

## Dependency

Immuannot-bt depends on:

* TRUST4
* minimap2 (2.27-r1193) 
* seqtk (1.4-r130-dirty)

All source codes and pre-build binaries can be found in folder 'dep-tools'. 
If any pre-build binary does work, you may need to re-build it from the source
and copy it to the folder 'dep-tools/bin'.

[\[top\]](#content)

# Input and output

Immuannot-bt takes fasta or fastq format as input, and generates a gff file for
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

