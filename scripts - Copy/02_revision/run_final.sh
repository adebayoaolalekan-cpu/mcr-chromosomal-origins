#!/bin/bash
IQ=/tmp/iqtree-2.3.6-Linux-intel/bin/iqtree2
cd /home/claude/pet/analysis/rev
$IQ -s aln_full.faa -te ml.treefile -m LG+F+I+R6 -asr -T 2 --prefix asr1 --seed 12345 --quiet
echo ASR1_DONE
$IQ -s v2_aln_automated1.faa -m LG+F+R4 -B 1000 -alrt 1000 -T 2 \
  -o "OpgB_ECOLI_P39401,AslA_ECOLI_P25549,YidJ_ECOLI_P31447,YdeN_ECOLI_P77318" \
  --prefix v2 --seed 12345 --quiet
echo V2_DONE
$IQ -s v2_aln_full.faa -te v2.treefile -m LG+F+R4 -asr -T 2 --prefix asr2 --seed 12345 --quiet
echo ASR2_DONE
