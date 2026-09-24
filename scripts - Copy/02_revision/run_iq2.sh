#!/bin/bash
IQ=/tmp/iqtree-2.3.6-Linux-intel/bin/iqtree2
cd /home/claude/pet/analysis/rev
# constrained ML tree enforcing mcr monophyly
$IQ -s aln_automated1.faa -m LG+F+I+R6 -g constraint_mcr.nwk -T 2 --prefix con --seed 12345 --quiet
echo CON_DONE
cat ml.treefile con.treefile > trees_ml_con.nwk
$IQ -s aln_automated1.faa -m LG+F+I+R6 -z trees_ml_con.nwk -n 0 -zb 10000 -zw -au -T 2 --prefix au --seed 12345 --quiet
echo AU_DONE
