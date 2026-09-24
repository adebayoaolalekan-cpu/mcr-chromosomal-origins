#!/bin/bash
CODEML=/tmp/paml/src/codeml
cd /home/claude/pet/analysis/rev/fam
run(){ f=$1; M=$2; d="${f}_M${M}"; mkdir -p $d
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
(cd $d && nice -n 19 $CODEML codeml.ctl > codeml.stdout 2>&1); echo "done $f M$M"; }
for f in mcr-2 mcr-4 mcr-5 mcr-8 mcr-10 mcr-1 mcr-3; do run $f 0; done
for f in mcr-1 mcr-3; do for M in 1 2 7 8; do run $f $M; done; done
echo CODEML2_DONE
