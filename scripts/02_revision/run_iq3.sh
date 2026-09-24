#!/bin/bash
IQ=/tmp/iqtree-2.3.6-Linux-intel/bin/iqtree2
cd /home/claude/pet/analysis/rev
$IQ -s aln_full.faa -te ml.treefile -m LG+F+I+R6 -asr -T 2 --prefix asrmid --seed 12345 --quiet
echo ASR_MID_DONE
