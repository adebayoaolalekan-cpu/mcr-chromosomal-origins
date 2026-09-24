#!/bin/bash
IQ=/tmp/iqtree-2.3.6-Linux-intel/bin/iqtree2
MSET="LG,WAG,JTT,VT,Blosum62,Dayhoff,cpREV,rtREV,PMB"
$IQ -s aln_automated1.faa -m MFP -mset $MSET -B 1000 -alrt 1000 -T 2 --prefix ml --seed 12345 --quiet
echo "STEP1_DONE"
$IQ -s og_aln_automated1.faa -m MFP -mset $MSET -B 1000 -alrt 1000 -T 2 \
   -o "OpgB_ECOLI_P39401,AslA_ECOLI_P25549,YidJ_ECOLI_P31447,YdeN_ECOLI_P77318" \
   --prefix og --seed 12345 --quiet
echo "STEP2_DONE"
