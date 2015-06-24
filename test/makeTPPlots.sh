#!/bin/bash

echo "make eta/pt-reweighted files"
#-1 option means all three scenarios, but you can spcify one particular scneraio with 1,2,3 (PH1_A0, PH1_A1k, and HGC respectively).
#root -b -q -l makeTrainTrees.cc\(-1\)
root -b -q -l makeTrainTrees.cc\(3\)

echo "retrain MVAs"
echo " --- PH1_A0_PU50"
#root -b -q -l PhotonIDMVA_Training_PH1.cc+\(\"PH1_A0_PU50\"\)

echo " --- PH1_A0_PU50"
#root -b -q -l PhotonIDMVA_Training_PH1.cc+\(\"PH1_A1k_PU140\"\)

echo " --- HGC"
root -b -q -l PhotonIDMVA_Training.cc+

echo "add new MVA score to trees"
#-1 option means all three scenarios, but you can spcify one particular scneraio with 1,2,3 (PH1_A0, PH1_A1k, and HGC respectively).
#root -b -q -l makeMVATrees.cc+\(-1\)
root -b -q -l makeMVATrees.cc+\(3\)

echo "make plots"
root -b -q -l plots.C

echo "transfer to web page"
#scp *pdf lc1113@lx05.hep.ph.ic.ac.uk:/home/hep/lc1113/public_html/HGC_Pt-eta_reweight/.


echo "done"




