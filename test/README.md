# Getting set up

You will need to get CMSSW_6_2_0_SLHC25_patch2 (or higher?).
(instruction untested, contact me if any problems!)

```
cmsrel CMSSW_6_2_0_SLHC25_patch2
cd CMSSW_6_2_0_SLHC25_patch2/src
git cms-init
git cms-addpkg RecoEcal/EgammaClusterAlgos
git cms-addpkg RecoEgamma/EgammaIsolationAlgos
git cms-merge-topic -u ldcorpe:topic-HGC-phoID-TP
git clone git@github.com:ldcorpe/HGCanalysis UserCode/HGCanalysis
cd UserCode/HGCanalysis
git checkout topic-HGC-phoID-TP
cd ../../
scram b -j9
```

# Introduction

This documentation is intended to help re-run the Photon ID study for the HGC TP.
The idea is to compare three scenarios:
* PH1_A0_50PU - Phase 1 detector without ageing and with 50PU
* PH1_A1k_140PU - Phase 1 detector with 1000/fb ageing and 140PU
* HGC_140PU - HGC with 140PU

The general strategy is to take the RECO simulation of gamma+jets in each case and run it through a simple analyzer to get a tree with a certain number of ID variables included.
The output is therefore a tree with a bunch of ID variables in each case.
Signal is the true gamma (ie a reco::Photon (PH1) or reco::SuperCluster (HGC)) which is gen-matched using a dR<0.1 to the genPhoton).
Background is any other photon sufficiently far away from the true photon (say dR>1), to avoid truth in the fake sample.

The next step is to train an MVA for each scenario to separate signal and background. One subtlety here is that for the training we want the signal and background distributions to be flat in eta/pt, since we want to avoid low stats bins being overlooked by the MVA. So you need to do a 2D eta/pt reweighing first.

Once you have trained the MVA and set up a new set of weights, the idea is to calculate the new MVA value for each event on the fly.

The figure of merit for comparison between scenarios will be: what is the overall fake rate (#fake photons passing some MVA cut/ #considered events) for a similar overall efficiency (#true photons found/#gen photons within the pt/eta regions we consider). So, the next step is to pick a Working Point (WP) to compare the scenarios. For example, I pick two WPs - one where all three scenarios have ~85% efficiency within the |eta| 1.6 to 2.5 region, and one with all three have ~90% efficiency. The WP is just a cut on the MVA value in each scenario.

Then we make the fake rate/ efficient plots and we are done.

# Running the ID variable trees from the RECO gamma+jet samples.

The first step is to get some basic variables to work with from each tree.
The code for the PH1 scenarios is here:
`../plugins/HoverEAnalyzer_Phase1.cc`

The code for the HGC scenarios is here:
`../plugins/HoverEAnalyzer_MVA.cc`

They differ slightly because for the HGC you do not have ready-made photons and need to use reco::SuperClusters instead.
And you need to therefore recalculate certain id variables.

You can run them each locally using the config files:
```
hoverEHggConfig_PH1_A0_PU50_MVA.py
hoverEHggConfig_PH1_A1k_PU140_MVA.py
hoverEHggConfig_MVA.py 
``` 
(the last one is for the HGC).

You can submit the jobs using crab 3 using the config files (+NB you need to add in some output directory LFN for this to work).

```
crab3_PH1_A0_PU50.py
crab3_PH1_A1k_PU140.py
crab3_HGC.py
```

To submit the jobs, you can do:

```
source initCmsEnv.sh #set up crab3 environment
cmsenv
crab submit <config file> 
```

And then you can check the status of your jobs by doing:

``crab status crab_projects/<relevant subdir> ``

Once your jobs are done you can should hadd the output files into files of the corresponding names in your local working dir:
```
HoverE_PH1_A0_PU50.root
HoverE_PH1_A1k_PU140.root
HoverE_HGC_015.root
```
(015 represents HOE cone size).

These file names are hard-coded in the next step (sorry)...
You could add some automation here if needed.

# Training the ID trees.

The first thing you need to do once you have your hadded files is re-weight the distributions in eta/pt as described above.
There is a compiled root macro to do this:

``makeTrainTrees.cc``

Which takes as input a number as an argument :
*-1 for all three samples to be reweighted
* 1 for just PH1 A0 PU50
* 2 for just PH1 A1k PU140
* 3 for just HGC

so run it like this  `` root -b -q -l makeTrainTrees.cc\(-1\) ``

It basically fills a TH2D of pt:eta with all events (separately for sig/bkg) and applies to each event in any bin a weight of 1/nEventsInBin
It will create separate output trees with the reweighted distributions:


```
HoverE_PH1_A0_PU50_weight.root
HoverE_PH1_A1k_PU140_weight.root
HoverE_HGC_weight.root
```


Next you want to do the actual training using TMVA.
I have compiled root macros that does this for PH1 and HGC scenarios.
The PH1 macro `PhotonIDMVA_Training_PH1.cc` takes as argument the name of the sample.

```
root -b -q -l PhotonIDMVA_Training_PH1.cc+\(\"PH1_A0_PU50\"\)
root -b -q -l PhotonIDMVA_Training_PH1.cc+\(\"PH1_A1k_PU140\"\)
```

And here is the HGC one (`PhotonIDMVA_Training.cc`) which needs no argument since it just does the HGC:
``root -b -q -l PhotonIDMVA_Training.cc+``

This will dump new weight files in the `weights` folder.

# Evaluating the MVAs

Once you have done training, you need to re-evaluate teh MVA on the fly and save it to a new tree.

Thankfully you can do this easy by using this root macro `makeMVATrees.cc`.

Usage is identical to `makeTrainTrees.cc`:

*-1 for all three samples to be reweighted
* 1 for just PH1 A0 PU50
* 2 for just PH1 A1k PU140
* 3 for just HGC

Use it like this : `` root -b -q -l makeMVATrees.cc+\(-1\)``


It will dump  the required output files:

```
HoverE_PH1_A0_PU50_newMVA.root
HoverE_PH1_A1k_PU140_newMVA.root
HoverE_HGC_newMVA.root
```

# Generating the plot

Finally, you can generate all the TP plots by running:

``root -b -q -l plots.C``

This root macro will do all the final bits: pick the WPs at 85% and 95% efficiency, determine the corresponding MVA cuts to apply and make all the plots you need.


NB once you have go the hang on things, you can do all the steps from training to plots by running:
``./makeTPPlots.sh``.

Good luck!

Louie Corpe  19.06.15



