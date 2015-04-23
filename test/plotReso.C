#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TMath.h"

int plotReso() {

  TFile *fin = TFile::Open("TestAMStyleHgg_140PU_testv5Geom_LC3x3_cleanedSC.root");
  fin->cd("hgg");

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  TCanvas *myc = new TCanvas("myc","myc",1);
  myc->Print("InfoPlots.pdf[");

  std::string cut1 = "isValid1==1 && dRTrue1<0.05 && TMath::fabs(etaTrue1)>1.5 && TMath::fabs(etaTrue1)<2.7";
  std::string cut2 = "isValid2==1 && dRTrue2<0.05 && TMath::fabs(etaTrue2)>1.5 && TMath::fabs(etaTrue2)<2.7";

  const unsigned nV = 14;
  const std::string var1[nV] = {
    "xvtxTrue1","yvtxTrue1","zvtxTrue1",
    "eSR3Reco1/eTrue1","eSR5Reco1/eTrue1",
    "etaSC1-etaTrue1","etaReco1-etaTrue1",
    "phiSC1-phiTrue1","phiReco1-phiTrue1",
    "etaWidth1","phiWidth1","clustersSize1","nClusters091",
    "eSeed1/eSC1"
  };
  const std::string var2[nV] = {
    "xvtxTrue2","yvtxTrue2","zvtxTrue2",
    "eSR3Reco2/eTrue2","eSR5Reco2/eTrue2",
    "etaSC2-etaTrue2","etaReco2-etaTrue2",
    "phiSC2-phiTrue2","phiReco2-phiTrue2",
    "etaWidth2","phiWidth2","clustersSize2","nClusters092",
    "eSeed2/eSC2"
  };

  TH1F *hist[nV];

  gStyle->SetOptStat("eMRuo");
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    t->Draw(var1[iV].c_str(),cut1.c_str());
    hist[iV] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(var1[iV].c_str());
    t->Draw(var2[iV].c_str(),cut2.c_str());
    hist[iV]->Add((TH1F*)(gPad->GetPrimitive("htemp")));
    std::cout << " -- Hist " << var1[iV] << " " << hist[iV]->GetEntries() << " " << hist[iV]->GetMean() << "" << hist[iV]->GetRMS() << std::endl;
    myc->Clear();
    myc->cd();
    hist[iV]->Draw();
    myc->Update();
    myc->Print("InfoPlots.pdf");

  }//loop on variables

  myc->Print("InfoPlots.pdf]");


  return 0;
}
