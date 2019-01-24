#!/bin/sh

infile=$1
prefix=$2

echo "input format:"
echo "     CHR POS REF ALT RSID N Z"
echo ""
echo "annotating $infile..."

path_epacts=/net/snowwhite/home/corbinq/statgen/EPACTS

anno_bin=$path_epacts/src/anno
gfile=$path_epacts/data/hg19_refFlat.txt.gz
cfile=$path_epacts/data/codon.txt
pfile=$path_epacts/data/priority.txt
rfile=$path_funcwas/data/human_g1k_v37.fasta

$anno_bin -i $infile -o $prefix.tmp -g $gfile --inputFormat plain -c $cfile -p $pfile -r $rfile

cut -f1-8 $prefix.tmp | bgzip -c > $prefix.az.gz

echo "annotated gwas output written to $prefix.az.gz"

rm -f $prefix.tmp*

tabix -p vcf $prefix.az.gz
