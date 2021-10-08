#!/bin/bash -ue
mkdir results
mv 109707.pop1popout.log.weir.fst 123563.pop1pop2.log.weir.fst 109707.pop1pop2.log.weir.fst 123563.pop1popout.log.weir.fst 109707.pop2popout.log.weir.fst 123563.pop2popout.log.weir.fst mart_export.txt pop1.frq pop3.frq pop2.frq pop1vspop2.csv pop2vsout.csv pop1vsout.csv results
Rscript --vanilla pbs_calculator.R ./results pbs.png
