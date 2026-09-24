#!/bin/bash
CODEML=/tmp/paml/src/codeml
cd /home/claude/pet/analysis/rev/fam
for f in mcr-1 mcr-3 mcr-2 mcr-4 mcr-5 mcr-8 mcr-10; do
for M in 0 1 2 7 8; do
  d="${f}_M${M}"; mkdir -p $d
  cat > $d/codeml.ctl <<CTL
seqfile = ../${f}.phy
treefile = ../${f}.nwk
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
model = 0
NSsites = $M
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
CTL
  (cd $d && nice -n 15 $CODEML codeml.ctl > codeml.stdout 2>&1)
  echo "done $f M$M"
done
done
echo CODEML_ALL_DONE
