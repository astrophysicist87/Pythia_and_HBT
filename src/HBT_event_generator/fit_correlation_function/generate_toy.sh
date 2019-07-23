#! /usr/bin/env bash

echo '#                          
#' > results/HBT_pipiCF.dat && awk 'BEGIN {for (KT=0.1; KT<1.0; KT+=0.2)
for (qo=-0.125; qo<=0.126; qo+=0.025)
for (qs=-0.125; qs<=0.126; qs+=0.025)
for (ql=-0.25; ql<=0.26; ql+=0.05)
print KT, 0.0, 0.0, qo, qs, ql, 0.0, 0.0, 0.0, 1+exp(-(3**2*(qo**2+qs**2+0.5*ql**2)/0.19733**2)-((qo**4+qs**4+0.5*ql**4)/0.19733**4)), 0.001 }' | column -t >> results/HBT_pipiCF.dat
