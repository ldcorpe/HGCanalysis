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
#include "TProfile.h"

int plotReso() {

  //std::string outfile = "InfoPlots_140";
  std::string outfile = "InfoPlots_singleGamma";

  std::string cut1 = "isValid1==1 && converted1==0 && noPhiCrack1==1 && dRTrue1<0.05 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && showerMax1>5 && showerMax1<25 && eSR3Reco1>0 && eSR3Reco1/eTrue1>70";// && TMath::Abs(etaReco1-etaTrue1)<0.1 && eSeed1/eSC1>0.95 && eSeed1/eSC1<1.01";
  std::string cut2 = "isValid2==1 && converted2==0 && noPhiCrack2==1 && dRTrue2<0.05 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && showerMax2>5 && showerMax2<25 && eSR3Reco2>0 && eSR3Reco2/eTrue2>70";// && TMath::Abs(etaReco2-etaTrue2)<0.1 && eSeed2/eSC2>0.95 && eSeed2/eSC2<1.01";

  //TFile *fin = TFile::Open("Truth_Hgg_140pu.root");
  TFile *fin = TFile::Open("Calib_singleGamma_0pu.root");
  if (!fin) return 1;
  fin->cd("hgg");

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  const unsigned nV = 17;
  std::string var[nV] = {
    "xvtxTrue","yvtxTrue","zvtxTrue",
    "dRTrue",
    "eSR3RecoOvereTrue","eSR5RecoOvereTrue",
    "etaSCMinusetaTrue","etaRecoMinusetaTrue",
    "phiSCMinusphiTrue","phiRecoMinusphiTrue",
    "etaWidth","phiWidth","clustersSize","nClusters09",
    "eSeedOvereSC","showerMax",
    "eSR3RecoOvereTrue"
  };
  std::string var1[nV];
  std::string var2[nV];
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    size_t over=var[iV].find("Over");
    size_t minus=var[iV].find("Minus");
    size_t end=var[iV].npos;
    size_t size=var[iV].size();
    //std::cout << var[iV] << " pos(over) " << over << " pos(minus) " << minus << " end " << size << std::endl;
    if (over==end && 
	minus==end){
      var1[iV] = var[iV]+"1";
      var2[iV] = var[iV]+"2";
    }
    else {
      if (over!=end){
	var1[iV]=var[iV].substr(0,over)+"1/"+var[iV].substr(over+4,end)+"1";
	var2[iV]=var[iV].substr(0,over)+"2/"+var[iV].substr(over+4,end)+"2";
      }
      else if (minus!=end){
	var1[iV]=var[iV].substr(0,minus)+"1-"+var[iV].substr(minus+5,end)+"1";
	var2[iV]=var[iV].substr(0,minus)+"2-"+var[iV].substr(minus+5,end)+"2";
      }
    }
    std::cout << " Vars " << var[iV] << " " << var1[iV] << " " << var2[iV] << std::endl;
  }

  const unsigned nMore = 6;
  const unsigned nC = nV+nMore;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1);
  }
  myc[0]->cd();
  myc[0]->Print((outfile+".pdf[").c_str());

  gStyle->SetOptStat("e");
  gStyle->SetOptFit(1111);

  TH2F *h2 = new TH2F("h2",";egen (GeV);ereco 3x3 (mips); photons",30,0,300,1000,0,25000);
  t->Draw("eSR3Reco1:eTrue1>>h2",cut1.c_str());
  t->Draw("eSR3Reco2:eTrue2>>+h2",cut2.c_str());
  myc[0]->Clear();
  myc[0]->Divide(1,2);



  h2->Draw("colz");

  h2->ProfileX();
  myc[0]->Update();
  TProfile *h2_pfx = (TProfile*)gDirectory->Get("h2_pfx");
  h2_pfx->SetMarkerStyle(22);
  h2_pfx->SetMarkerColor(1);
  h2_pfx->Draw("PEsame");
  h2_pfx->Fit("pol1","RL","same",90,290);

  TF1 *fit3 = h2_pfx->GetFunction("pol1");
  if (!fit3) return 1;
  double slope3 = fit3->GetParameter(1);
  double offset3 = fit3->GetParameter(0);

  myc[0]->Update();
  myc[0]->Print((outfile+".pdf").c_str());

  myc[1]->cd();
  TH2F *h2eta = new TH2F("h2eta",";#eta;ereco/egen; photons",30,1.5,3,100,0,2);

  std::ostringstream cor;
  cor.str("");
  cor << "(eSR3Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1):TMath::Abs(etaTrue1)>>h2eta";
  t->Draw(cor.str().c_str(),cut1.c_str());
  cor.str("");
  cor << "(eSR3Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2):TMath::Abs(etaTrue2)>>+h2eta";
  t->Draw(cor.str().c_str(),cut2.c_str());
  h2eta->Draw("colz");

  h2eta->ProfileX();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  myc[1]->Update();
  TProfile *h2eta_pfx = (TProfile*)gDirectory->Get("h2eta_pfx");
  h2eta_pfx->SetMarkerStyle(22);
  h2eta_pfx->SetMarkerColor(1);
  h2eta_pfx->Draw("PEsame");
  h2eta_pfx->Fit("pol1","","same");

  TF1 *fiteta3 = h2eta_pfx->GetFunction("pol1");
  if (!fiteta3) return 1;

  double p03 = fit3->GetParameter(0);
  double p13 = fit3->GetParameter(1);
  //double p23 = fit3->GetParameter(2);

  myc[1]->Update();
  myc[1]->Print((outfile+".pdf").c_str());

  /*
  myc[2]->cd();

  TH2F *h25 = new TH2F("h25",";egen (GeV);ereco 5x5 (mips); photons",100,0,1000,100,0,1000);
  t->Draw("eSR5Reco1:eTrue1>>h25",cut1.c_str());
  t->Draw("eSR5Reco2:eTrue2>>+h25",cut2.c_str());
  h25->Draw("colz");

  h25->ProfileX();
  myc[2]->Update();
  TProfile *h25_pfx = (TProfile*)gDirectory->Get("h25_pfx");
  h25_pfx->SetMarkerStyle(22);
  h25_pfx->SetMarkerColor(1);
  h25_pfx->Draw("PEsame");
  h25_pfx->Fit("pol1","RL","same",80,200);

  TF1 *fit5 = h25_pfx->GetFunction("pol1");
  if (!fit5) return 1;
  double slope5 = fit5->GetParameter(1);
  double offset5 = fit5->GetParameter(0);

  myc[2]->Update();
  myc[2]->Print((outfile+".pdf").c_str());

  myc[3]->cd();

  TH2F *h25eta = new TH2F("h25eta",";#eta;ereco 5x5/egen; photons",30,1.5,3,100,0,2);
  cor.str("");
  cor << "(eSR5Reco1-" << offset5 << ")/(" << slope5 << "*eTrue1):TMath::Abs(etaTrue1)>>h25eta";
  t->Draw(cor.str().c_str(),cut1.c_str());
  cor.str("");
  cor << "(eSR5Reco2-" << offset5 << ")/(" << slope5 << "*eTrue2):TMath::Abs(etaTrue2)>>+h25eta";
  t->Draw(cor.str().c_str(),cut2.c_str());
  //  t->Draw("eSR5Reco1/eTrue1:TMath::Abs(etaTrue1)>>h25eta",cut1.c_str());
  //t->Draw("eSR5Reco2/eTrue2:TMath::Abs(etaTrue2)>>htmp5eta",cut2.c_str());
  h25eta->Draw("colz");

  h25eta->ProfileX();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  myc[3]->Update();
  TProfile *h25eta_pfx = (TProfile*)gDirectory->Get("h25eta_pfx");
  h25eta_pfx->SetMarkerStyle(22);
  h25eta_pfx->SetMarkerColor(1);
  h25eta_pfx->Draw("PEsame");
  h25eta_pfx->Fit("pol2","","same");

  TF1 *fiteta5 = h25eta_pfx->GetFunction("pol2");
  if (!fiteta5) return 1;
  double p05 = fit5->GetParameter(0);
  double p15 = fit5->GetParameter(1);
  double p25 = fit5->GetParameter(2);
  myc[3]->Update();
  myc[3]->Print((outfile+".pdf").c_str());
  */

  myc[4]->cd();
  TH2F *h2phi = new TH2F("h2phi",";#phi;ereco 3x3/egen; photons",360,-3.1416,3.1416,100,0,2);
  //t->Draw("eSR3Reco1/eTrue1:phiTrue1>>h2phi",cut1.c_str());
  //t->Draw("eSR3Reco2/eTrue2:phiTrue2>>+h2phi",cut2.c_str());
  cor.str("");
  cor << "(eSR3Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1):phiTrue1>>h2phi";
  t->Draw(cor.str().c_str(),cut1.c_str());
  cor.str("");
  cor << "(eSR3Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2):phiTrue2>>+h2phi";
  t->Draw(cor.str().c_str(),cut2.c_str());

  h2phi->Draw("colz");
  h2phi->ProfileX();
  TProfile *h2phi_pfx = (TProfile*)gDirectory->Get("h2phi_pfx");
  h2phi_pfx->SetMarkerStyle(22);
  h2phi_pfx->SetMarkerColor(1);
  h2phi_pfx->Draw("PEsame");
  myc[4]->Update();
  myc[4]->Print((outfile+".pdf").c_str());

  myc[5]->cd();
  TH2F *hetaphi = new TH2F("hetaphi",";#phi;|#eta|; photons",36,-3.1416,3.1416,15,1.5,3);
  t->Draw("TMath::Abs(etaTrue1):phiTrue1>>hetaphi",cut1.c_str());//,"isValid1==1 && converted1==0 && dRTrue1<0.05 && TMath::Abs(etaTrue1)>1.5 && TMath::Abs(etaTrue1)<2.7 && TMath::Abs(etaReco1-etaTrue1)<0.1 && eSeed1/eSC1>0.95 && eSeed1/eSC1<1.01 && eSR3Reco1/eTrue1<0.9");
  t->Draw("TMath::Abs(etaTrue2):phiTrue2>>+hetaphi",cut2.c_str());//,"isValid2==1 && converted2==0 && dRTrue2<0.05 && TMath::Abs(etaTrue2)>1.5 && TMath::Abs(etaTrue2)<2.7 && TMath::Abs(etaReco2-etaTrue2)<0.1 && eSeed2/eSC2>0.95 && eSeed2/eSC2<1.01 && eSR3Reco2/eTrue2<0.9");
  hetaphi->Draw("colz");
  myc[5]->Update();
  myc[5]->Print((outfile+".pdf").c_str());

  cor.str("");
  cor << "(eSR3Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1)";
  var1[4] = cor.str();
  //cor.str("");
  //cor << "(eSR5Reco1-" << offset5 << ")/(" << slope5 << "*eTrue1)";
  //var1[5] = cor.str();

  cor.str("");
  cor << "(eSR3Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2)";
  var2[4] = cor.str();

  var[4]+="Calib";

  //cor.str("");
  //cor << "(eSR5Reco2-" << offset5 << ")/(" << slope5 << "*eTrue2)";
  //var2[5] = cor.str();

  TH1F *hist[nV];

  gStyle->SetOptStat("eMRuo");
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    myc[nMore+iV]->cd();
    t->Draw(var1[iV].c_str(),cut1.c_str());
    hist[iV] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(var[iV].c_str());
    std::cout << " -- first photon " << hist[iV]->GetEntries() << std::endl;
    cor.str("");
    cor << var2[iV] << ">>+" << var[iV];
    std::cout << "debug " << cor.str() << std::endl;
    t->Draw(cor.str().c_str(),cut2.c_str());
    std::cout << " -- second photon " << hist[iV]->GetEntries() << std::endl;
    std::cout << " -- Hist " << var[iV] << " " << hist[iV]->GetEntries() << " " << hist[iV]->GetMean() << " " << hist[iV]->GetRMS() << std::endl;
    hist[iV]->SetTitle("");
    hist[iV]->GetXaxis()->SetTitle(var[iV].c_str());
    hist[iV]->GetYaxis()->SetTitle("Photons");
    myc[nMore+iV]->Clear();
    myc[nMore+iV]->cd();
    hist[iV]->Draw();
    gStyle->SetOptFit(1111);
    if (iV==4) hist[iV]->Fit("gaus","LR","same",0.8,1.2);
    myc[nMore+iV]->Update();
    myc[nMore+iV]->Print((outfile+".pdf").c_str());
    
    //if (iV==3) return 1;
  }//loop on variables
  
  myc[0]->Print((outfile+".pdf]").c_str());


  return 0;
}
