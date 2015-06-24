{

int n1=0;
int n2=0;
int n3=0;
float nHGCFiles=880.;
gStyle->SetOptStat(1) ;
TFile *_file1 = TFile::Open("HoverE_PH1_A0_PU50_newMVA.root"); 
TFile *_file2 = TFile::Open("HoverE_PH1_A1k_PU140_newMVA.root"); 
TFile *_file3 = TFile::Open("HoverE_HGC_newMVA.root"); 

std::cout << " file 1 open " << _file1.IsOpen() << std::endl;
std::cout << " file 2 open " << _file2.IsOpen() << std::endl;
std::cout << " file 3 open " << _file3.IsOpen() << std::endl;

TTree *ts1 = (TTree*) _file1->Get("treePhoton");
TTree *tb1 = (TTree*) _file1->Get("treeBackground");
TTree *ts2 = (TTree*) _file2->Get("treePhoton");
TTree *tb2 = (TTree*) _file2->Get("treeBackground");
TTree *ts3 = (TTree*) _file3->Get("treePhoton");
TTree *tb3 = (TTree*) _file3->Get("treeBackground");

std::cout << " ts1 : " << ts1->GetEntries() << ", tb1 : " << tb1->GetEntries() << std::endl;
std::cout << " ts2 : " << ts2->GetEntries() << ", tb2 : " << tb2->GetEntries() << std::endl;
std::cout << " ts3 : " << ts3->GetEntries() << ", tb3 : " << tb3->GetEntries() << std::endl;


TFile *tf = new TFile("HGC_photonID_TPplots.root","RECREATE");

TCanvas *c = new TCanvas("c","c",500,500);
//--------------> Input Variables <--------------//

//////---------1
ts1->Draw("sigmaIetaIeta>>s1(100,0,0.08)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && matchIndex>-1");
tb1->Draw("sigmaIetaIeta>>b1(100,0,0.08)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100");
s1->SetLineColor(kRed); s1->SetFillColor(kRed); s1->SetFillStyle(3001);
s1->DrawNormalized(); b1->DrawNormalized("same");
c->SaveAs("PH1_A0_PU50_sigmaIetaIeta.pdf");
c->Clear();
n1=s1->GetEntries();

ts1->Draw("hoe>>s1(100,0,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && matchIndex>-1");
tb1->Draw("hoe>>b1(100,0,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100");
s1->SetLineColor(kRed); s1->SetFillColor(kRed); s1->SetFillStyle(3001);
c->SetLogy(1);
s1->DrawNormalized(); b1->DrawNormalized("same");
c->SaveAs("PH1_A0_PU50_hoe.pdf");
c->SetLogy(0);
//////---------2
ts2->Draw("sigmaIetaIeta>>s2(100,0,0.08)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && matchIndex>-1");
tb2->Draw("sigmaIetaIeta>>b2(100,0,0.08)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100");
s2->SetLineColor(kRed); s2->SetFillColor(kRed); s2->SetFillStyle(3001);
s2->DrawNormalized(); b2->DrawNormalized("same");
c->SaveAs("PH1_A1k_PU140_sigmaIetaIeta.pdf");
c->Clear();
n2=s2->GetEntries();

ts2->Draw("hoe>>s2(100,0,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && matchIndex>-1");
tb2->Draw("hoe>>b2(100,0,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100");
s2->SetLineColor(kRed); s2->SetFillColor(kRed); s2->SetFillStyle(3001);
c->SetLogy(1);
s2->DrawNormalized(); b2->DrawNormalized("same");
c->SaveAs("PH1_A1k_PU140_hoe.pdf");
c->SetLogy(0);
c->Clear();
c->Clear();
//////---------3
ts3->Draw("sigmaEtaEta>>s3(100,0,0.025)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 &&  hoe<0.5 && matchIndex>-1");
tb3->Draw("sigmaEtaEta>>b3(100,0,0.025)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 &&  hoe<0.5 && pt<100");
s3->SetLineColor(kRed); s3->SetFillColor(kRed); s3->SetFillStyle(3001);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("HGC_sigmaEtaEta.pdf");
c->Clear();

ts3->Draw("hoe>>s3(100,0,1)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && hoe<0.5 && matchIndex>-1");
tb3->Draw("hoe>>b3(100,0,1)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && hoe<0.5 && pt<100");
s3->SetLineColor(kRed); s3->SetFillColor(kRed); s3->SetFillStyle(3001);
c->SetLogy(1);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("HGC_hoe.pdf");
c->SetLogy(0);
c->Clear();
ts3->Draw("lengthCompatibility>>s3(100,-10,10)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && matchIndex>-1");
tb3->Draw("lengthCompatibility>>b3(100,-10,10)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100");
s3->SetLineColor(kRed); s3->SetFillColor(kRed); s3->SetFillStyle(3001);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("HGC_lengthCompatibility.pdf");
c->Clear();
/*ts3->Draw("hcalIso>>s3(100,0,20)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && matchIndex>-1");
tb3->Draw("hcalIso>>b3(100,0,20)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100");
s3->SetLineColor(kRed); s3->SetFillColor(kRed); s3->SetFillStyle(3001);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("HGC_hcalIso.pdf");
c->Clear();
ts3->Draw("ecalIso>>s3(100,0,20000)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && matchIndex>-1");
tb3->Draw("ecalIso>>b3(100,0,20000)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100");
s3->SetLineColor(kRed); s3->SetFillColor(kRed); s3->SetFillStyle(3001);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("HGC_ecalIso.pdf");
c->Clear();
ts3->Draw("trkIso>>s3(100,0,20)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && matchIndex>-1");
tb3->Draw("trkIso>>b3(100,0,20)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100");
s3->SetLineColor(kRed); s3->SetFillColor(kRed); s3->SetFillStyle(3001);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("HGC_trkIso.pdf");*/
c->Clear();
n3=s3->GetEntries();


std::cout << " n1 " << n1 << ", n2 " << n2 << ", n3 " << n3 << std::endl;

//--------> ROC Curves <-----------//

ts1->Draw("MVA>>s1(200,-1,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && matchIndex>-1");
tb1->Draw("MVA>>b1(200,-1,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100");

ts2->Draw("MVA>>s2(200,-1,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && matchIndex>-1");
tb2->Draw("MVA>>b2(200,-1,1)","abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100");

ts3->Draw("MVA>>s3(200,-1,1)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && matchIndex>-1 && hoe<0.5");
tb3->Draw("MVA>>b3(200,-1,1)","abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && hoe<0.5");


TLegend *dummy = new TLegend();

std::cout << "s1 " << s1->GetEntries() << ", b1 : " << b1->GetEntries() << std::endl;
std::cout << "s2 " << s2->GetEntries() << ", b2 : " << b2->GetEntries() << std::endl;
std::cout << "s3 " << s3->GetEntries() << ", b3 : " << b3->GetEntries() << std::endl;

gROOT->ProcessLine(".L makeROCs.cc");

std::pair < TGraph*, std::pair < float, float> > roc1 = makeROC2(s1 , b1, dummy, 1);
std::pair < TGraph*, std::pair < float, float> > roc2 = makeROC2(s2 , b2, dummy, 2);
std::pair < TGraph*, std::pair < float, float> > roc3 = makeROC2(s3 , b3, dummy, 3);


s1->SetLineColor(kRed); s1->SetFillColor(kRed) ; s1->SetFillStyle(3001);
s1->DrawNormalized(); b1->DrawNormalized("same");
c->SaveAs("MVAOutput_PH1_A0_PU50.pdf");
		c->Clear();
s2->SetLineColor(kRed); s2->SetFillColor(kRed) ; s2->SetFillStyle(3001);
s2->DrawNormalized(); b2->DrawNormalized("same");
c->SaveAs("MVAOutput_PH1_A1k_PU140.pdf");
		c->Clear();
s3->SetLineColor(kRed); s3->SetFillColor(kRed) ; s3->SetFillStyle(3001);
s3->DrawNormalized(); b3->DrawNormalized("same");
c->SaveAs("MVAOutput_HGC_PU140.pdf");
		c->Clear();

TMultiGraph *mul  = new TMultiGraph();

mul->Add(roc1.first);
mul->Add(roc2.first);
mul->Add(roc3.first);


    mul->SetTitle( ";signal efficiency;background rejection ");
		c->Clear();
mul->Draw( "APL" );
TLegend *tl = new TLegend(0.7,0.8,0.9,0.9);
tl->AddEntry(roc1.first,"PH1, Age0, PU50","lp");
tl->AddEntry(roc2.first,"PH1, Age1k, PU140","lp");
tl->AddEntry(roc3.first,"HGC PU140","lp");

tl->Draw();

c->SaveAs("ROCs.pdf");


float ph1_a0_pu50_wp1=-999.;
float ph1_a0_pu50_wp2=-999.;
float ph1_a1k_pu140_wp1=-999.;
float ph1_a1k_pu140_wp2=-999.;
float hgc_pu140_wp1=-999.;
float hgc_pu140_wp2=-999.;
bool flag1s1=1;
bool flag2s1=1;
bool flag1s2=1;
bool flag2s2=1;
bool flag1s3=1;
bool flag2s3=1;

//---------->Determine working points<--------------//

for (float mva =0.8 ; mva<1.0 ;mva=mva+0.01){


	ostringstream cut0;
	cut0 <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 &&hoe<0.5&& MVA> " <<mva;
	ts1->Draw("abs(etaEcal)>>s1numtmp(9,1.6,2.5)",cut0.str().c_str());
	ts1->Draw("abs(etaEcal)>>s1dentmp(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
	ts2->Draw("abs(etaEcal)>>s2numtmp(9,1.6,2.5)",cut0.str().c_str());
	ts2->Draw("abs(etaEcal)>>s2dentmp(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
	ts3->Draw("abs(etaEcal)>>s3numtmp(9,1.6,2.5)",cut0.str().c_str());
	ts3->Draw("abs(etaEcal)>>s3dentmp(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
	TH1F  *s1numtmp = (TH1F*)gDirectory->Get("s1numtmp");
	TH1F  *s1dentmp = (TH1F*)gDirectory->Get("s1dentmp");
	TH1F  *s2numtmp = (TH1F*)gDirectory->Get("s2numtmp");
	TH1F  *s2dentmp = (TH1F*)gDirectory->Get("s2dentmp");
	TH1F  *s3numtmp = (TH1F*)gDirectory->Get("s3numtmp");
	TH1F  *s3dentmp = (TH1F*)gDirectory->Get("s3dentmp");

	float as1 =s1numtmp->GetEntries();
	float bs1 =s1dentmp->GetEntries();
	float as2 =s2numtmp->GetEntries();
	float bs2 =s2dentmp->GetEntries();
	float as3 =s3numtmp->GetEntries();
	float bs3 =s3dentmp->GetEntries();
	//std::cout << " DEBUG mva - "<< mva <<", "<< a/b   <<std::endl;

	if (as1/bs1 <0.9 && flag1s1 ) {
		std::cout<<" PH1_A0_PU50 90% WP -- mva "<< mva << " a/b " << as1/bs1 << std::endl;
		ph1_a0_pu50_wp1 = mva;
		flag1s1 =0;
	}

	if (as1/bs1 <0.85 && flag2s1 ) {
		std::cout<<" PH1_A0_PU50 85% WP -- mva "<< mva << " a/b " << as1/bs1 << std::endl;
		ph1_a0_pu50_wp2 = mva;
		flag2s1 =0;
	}
	if (as2/bs2 <0.9 && flag1s2 ) {
		std::cout<<" PH1_A1k_PU140 90% WP -- mva "<< mva << " a/b " << as2/bs2 << std::endl;
		ph1_a1k_pu140_wp1 = mva;
		flag1s2 =0;
	}

	if (as2/bs2 <0.85 && flag2s2 ) {
		std::cout<<" PH1_A1k_PU140 85% WP -- mva "<< mva << " a/b " << as2/bs2 << std::endl;
		ph1_a1k_pu140_wp2 = mva;
		flag2s2 =0;
	}
	if (as3/bs3 <0.9 && flag1s3 ) {
		std::cout<<" HGC_PU140 90% WP -- mva "<< mva << " a/b " << as3/bs3 << std::endl;
		hgc_pu140_wp1 = mva;
		flag1s3 =0;
	}

	if (as3/bs3 <0.85 && flag2s3 ) {
		std::cout<<" HGC_PU140 85% WP -- mva "<< mva << " a/b " << as3/bs3 << std::endl;
		hgc_pu140_wp2 = mva;
		flag2s3 =0;
	}


}
c->Clear();

TLegend *tl = new TLegend(0.54,0.79,0.88,0.89);
TH1F *dummy2 = new TH1F("d","d",100,0,1);
//dummy2->SetMarkerColor(kWhite);
//tl->AddEntry(dummy2,"WP1 - <Eff>=90%","p");
 TPad *pad1 = new TPad("pad1","",0,0,1,1);
 TPad *pad2 = new TPad("pad2","",0,0,1,1);
pad2->SetFillStyle(4000); //will be transparent
pad2->SetFrameFillStyle(0);

pad1->Draw();
pad1->cd();


///-------------> Result wp1<--------------//
ostringstream cut1;
ostringstream cut2;
ostringstream cut3;
cut1 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && hoe<0.5  && MVA>" <<ph1_a0_pu50_wp1  ;
cut2 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && hoe<0.5  && MVA>" <<ph1_a1k_pu140_wp1;
cut3 <<"abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100 && hoe <0.5  && MVA>" <<hgc_pu140_wp1    ;
tb1->Draw("abs(etaSC)>>b1num(9,1.6,2.5)",cut1.str().c_str());
tb2->Draw("abs(etaSC)>>b2num(9,1.6,2.5)",cut2.str().c_str());
tb3->Draw("abs(etaSC)>>b3num(13,1.6,2.9)",cut3.str().c_str());
c->Clear();

//TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( b1num, b1den );
//TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( b2num, b2den );
//TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( b3num, b3den );

//b1num->Scale(1./(891473.));
//b2num->Scale(1./(994522.*(791./801.)));
//b3num->Scale(1./( 585260.*( nHGCFiles/880.)));	

b1num->Scale(1./(n1));
b2num->Scale(1./(n2));
b3num->Scale(1./(n3));	

b1num->SetLineColor(1);
b2num->SetLineColor(2);
b3num->SetLineColor(3);
b1num->SetMarkerStyle(24);
b2num->SetMarkerStyle(25);
b3num->SetMarkerStyle(26);
b1num->SetMarkerColor(1);
b2num->SetMarkerColor(2);
b3num->SetMarkerColor(3);
b3num->SetTitle(";;FakeRate");
b3num->Draw("pl[ Y+ X+");
b3num->GetYaxis()->SetRangeUser(0,0.00153*3);
b3num->GetXaxis()->SetRangeUser(1.6,2.9);
b3num->GetXaxis()->SetAxisColor(kWhite);
b3num->GetXaxis()->SetLabelColor(kWhite);
b2num->Draw("Apl same");
b1num->Draw("Apl same");


pad2->Draw();
pad2->cd();

ostringstream cut1e;
ostringstream cut2e;
ostringstream cut3e;
cut1e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && hoe<0.5 && MVA>" << ph1_a0_pu50_wp1  ;
cut2e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && hoe<0.5 && MVA>" << ph1_a1k_pu140_wp1;
cut3e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.9 && pt>30 && pt<100 && matchIndex>-1 && hoe<0.5 && MVA>" << hgc_pu140_wp1    ;
ts1->Draw("abs(etaEcal)>>s1num(9,1.6,2.5)",cut1e.str().c_str());
ts2->Draw("abs(etaEcal)>>s2num(9,1.6,2.5)",cut2e.str().c_str());
ts3->Draw("abs(etaEcal)>>s3num(13,1.6,2.9)",cut3e.str().c_str());
ts1->Draw("abs(etaEcal)>>s1den(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts2->Draw("abs(etaEcal)>>s2den(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts3->Draw("abs(etaEcal)>>s3den(13,1.6,2.9)","abs(etaEcal) >1.6 && abs(etaEcal)<2.9 && pt>30 && pt<100");
pad2->Clear();

TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( s1num, s1den );
TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( s2num, s2den );
TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( s3num, s3den );

eff1->SetLineColor(1);
eff2->SetLineColor(2);
eff3->SetLineColor(3);
eff1->SetMarkerStyle(20);
eff2->SetMarkerStyle(21);
eff3->SetMarkerStyle(22);
eff1->SetMarkerColor(1);
eff2->SetMarkerColor(2);
eff3->SetMarkerColor(3);
pad2->Clear();
eff3->GetYaxis()->SetRangeUser(0,1);
eff3->GetXaxis()->SetRangeUser(1.61,2.89);
eff3->Draw("APL");
eff3->SetTitle("Working Point 1 - <Signal Efficiency> ~ 90%;|#eta|;Signal Efficiency");
eff2->SetTitle(";;");
eff1->SetTitle(";;");
eff1->Draw("same PL");
eff2->Draw("same PL");

//c->SaveAs("Eff_wp1.pdf");/i
tl->SetNColumns(2);
tl->SetFillColor(kWhite);
tl->SetLineColor(kWhite);
tl->AddEntry(eff1,"efficiency","lp");
tl->AddEntry(b1num,"fake rate (PH1, Age0, PU50)","lp");
tl->AddEntry(eff2,"efficiency","lp");
tl->AddEntry(b2num,"fake rate (PH1, Age1k, PU140)","lp");
tl->AddEntry(eff3,"efficiency","lp");
tl->AddEntry(b3num,"fake rate (HGC PU140)","lp");
/*
TGaxis *axis = new TGaxis(gPad->GetUxmax(),
		gPad->GetUymin(),
		gPad->GetUxmax(),
		gPad->GetUymax(),
		0.01,10,510,"+L");
axis->Draw();*/
TGaxis::SetMaxDigits(2);
tl->Draw();
c->SaveAs("Result_vs_eta_wp1.pdf");
c->Write("Result_vs_eta_wp1");

c->Clear();

TLegend *tl = new TLegend(0.54,0.79,0.88,0.89);
TH1F *dummy2 = new TH1F("d","d",100,0,1);
//dummy2->SetMarkerColor(kWhite);
//tl->AddEntry(dummy2,"WP1 - <Eff>=90%","p");
 TPad *pad1 = new TPad("pad1","",0,0,1,1);
 TPad *pad2 = new TPad("pad2","",0,0,1,1);
pad2->SetFillStyle(4000); //will be transparent
pad2->SetFrameFillStyle(0);

pad1->Draw();
pad1->cd();

///-------------> Result wp2<--------------//
ostringstream cut1;
ostringstream cut2;
ostringstream cut3;
cut1 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && hoe<0.5 && MVA>" <<ph1_a0_pu50_wp2  ;
cut2 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100 && hoe<0.5  && MVA>" <<ph1_a1k_pu140_wp2;
cut3 <<"abs(etaSC) >1.6 && abs(etaSC)<2.9 && pt>30 && pt<100  && hoe<0.5 && MVA>" <<hgc_pu140_wp2    ;
tb1->Draw("abs(etaSC)>>b1num(9,1.6,2.5)",cut1.str().c_str());
tb2->Draw("abs(etaSC)>>b2num(9,1.6,2.5)",cut2.str().c_str());
tb3->Draw("abs(etaSC)>>b3num(13,1.6,2.9)",cut3.str().c_str());
c->Clear();

//TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( b1num, b1den );
//TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( b2num, b2den );
//TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( b3num, b3den );

//b1num->Scale(1./(891473.));
//b2num->Scale(1./(994522.*(791./801.)));
//b3num->Scale(1./( 585260.*( nHGCFiles/880.)));	

b1num->Scale(1./(n1));
b2num->Scale(1./(n2));
b3num->Scale(1./(n3));	

b1num->SetLineColor(1);
b2num->SetLineColor(2);
b3num->SetLineColor(3);
b1num->SetMarkerStyle(24);
b2num->SetMarkerStyle(25);
b3num->SetMarkerStyle(26);
b1num->SetMarkerColor(1);
b2num->SetMarkerColor(2);
b3num->SetMarkerColor(3);
b3num->SetTitle(";;FakeRate");
b3num->Draw("pl[ Y+ X+");
b3num->GetYaxis()->SetRangeUser(0,0.00153*3);
b3num->GetXaxis()->SetRangeUser(1.6,2.9);
b3num->GetXaxis()->SetAxisColor(kWhite);
b3num->GetXaxis()->SetLabelColor(kWhite);
b2num->Draw("Apl same");
b1num->Draw("Apl same");


pad2->Draw();
pad2->cd();

ostringstream cut1e;
ostringstream cut2e;
ostringstream cut3e;
cut1e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && hoe<0.5 && matchIndex>-1 && MVA>" << ph1_a0_pu50_wp2  ;
cut2e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && hoe<0.5 && matchIndex>-1 && MVA>" << ph1_a1k_pu140_wp2;
cut3e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.9 && pt>30 && pt<100 && hoe<0.5 && matchIndex>-1 && MVA>" << hgc_pu140_wp2    ;
ts1->Draw("abs(etaEcal)>>s1num(9,1.6,2.5)",cut1e.str().c_str());
ts2->Draw("abs(etaEcal)>>s2num(9,1.6,2.5)",cut2e.str().c_str());
ts3->Draw("abs(etaEcal)>>s3num(13,1.6,2.9)",cut3e.str().c_str());
ts1->Draw("abs(etaEcal)>>s1den(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts2->Draw("abs(etaEcal)>>s2den(9,1.6,2.5)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts3->Draw("abs(etaEcal)>>s3den(13,1.6,2.9)","abs(etaEcal) >1.6 && abs(etaEcal)<2.9 && pt>30 && pt<100");
pad2->Clear();

TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( s1num, s1den );
TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( s2num, s2den );
TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( s3num, s3den );

eff1->SetLineColor(1);
eff2->SetLineColor(2);
eff3->SetLineColor(3);
eff1->SetMarkerStyle(20);
eff2->SetMarkerStyle(21);
eff3->SetMarkerStyle(22);
eff1->SetMarkerColor(1);
eff2->SetMarkerColor(2);
eff3->SetMarkerColor(3);
pad2->Clear();
eff3->GetYaxis()->SetRangeUser(0,1);
eff3->GetXaxis()->SetRangeUser(1.61,2.89);
eff3->Draw("APL");
eff3->SetTitle("Working Point 2 - <Signal Efficiency> ~ 85%;|#eta|;Signal Efficiency");
eff2->SetTitle(";;");
eff1->SetTitle(";;");
eff1->Draw("same PL");
eff2->Draw("same PL");

//c->SaveAs("Eff_wp1.pdf");/i
tl->SetNColumns(2);
tl->SetFillColor(kWhite);
tl->SetLineColor(kWhite);
tl->AddEntry(eff1,"efficiency","lp");
tl->AddEntry(b1num,"fake rate (PH1, Age0, PU50)","lp");
tl->AddEntry(eff2,"efficiency","lp");
tl->AddEntry(b2num,"fake rate (PH1, Age1k, PU140)","lp");
tl->AddEntry(eff3,"efficiency","lp");
tl->AddEntry(b3num,"fake rate (HGC PU140)","lp");
/*
TGaxis *axis = new TGaxis(gPad->GetUxmax(),
		gPad->GetUymin(),
		gPad->GetUxmax(),
		gPad->GetUymax(),
		0.01,10,510,"+L");
axis->Draw();*/
TGaxis::SetMaxDigits(2);
tl->Draw();
c->SaveAs("Result_vs_eta_wp2.pdf");
c->Write("Result_vs_eta_wp1");





///////////////////////////////////////////////////////////
//PT SPLIT
/////////////////////////////////////////////////////////





c->Clear();




TLegend *tl = new TLegend(0.54,0.4,0.88,0.6);
dummy2->SetMarkerColor(kWhite);
//tl->AddEntry(dummy2,"WP1 - <Eff>=90%","p");
 TPad *pad1 = new TPad("pad1","",0,0,1,1);
 TPad *pad2 = new TPad("pad2","",0,0,1,1);
pad2->SetFillStyle(4000); //will be transparent
pad2->SetFrameFillStyle(0);

pad1->Draw();
pad1->cd();


///-------------> Result wp1<--------------//
ostringstream cut1;
ostringstream cut2;
ostringstream cut3;
cut1 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && MVA>" <<ph1_a0_pu50_wp1  ;
cut2 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && MVA>" <<ph1_a1k_pu140_wp1;
cut3 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && MVA>" <<hgc_pu140_wp1    ;
tb1->Draw("pt>>b1num(7,30,100)",cut1.str().c_str());
tb2->Draw("pt>>b2num(7,30,100)",cut2.str().c_str());
tb3->Draw("pt>>b3num(7,30,100)",cut3.str().c_str());
c->Clear();

//TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( b1num, b1den );
//TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( b2num, b2den );
//TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( b3num, b3den );

//b1num->Scale(1./(891473.));
//b2num->Scale(1./(994522.*(791./801.)));
//b3num->Scale(1./( 585260.*( nHGCFiles/880.)));	

b1num->Scale(1./(n1));
b2num->Scale(1./(n2));
b3num->Scale(1./(n3));	

b1num->SetLineColor(1);
b2num->SetLineColor(2);
b3num->SetLineColor(3);
b1num->SetMarkerStyle(24);
b2num->SetMarkerStyle(25);
b3num->SetMarkerStyle(26);
b1num->SetMarkerColor(1);
b2num->SetMarkerColor(2);
b3num->SetMarkerColor(3);
b3num->SetTitle(";;FakeRate");
b3num->Draw("pl[ Y+ X+");
b3num->GetYaxis()->SetRangeUser(0,0.00153*3);
b3num->GetXaxis()->SetRangeUser(30,100);
b3num->GetXaxis()->SetAxisColor(kWhite);
b3num->GetXaxis()->SetLabelColor(kWhite);
b2num->Draw("Apl same");
b1num->Draw("Apl same");


pad2->Draw();
pad2->cd();

ostringstream cut1e;
ostringstream cut2e;
ostringstream cut3e;
cut1e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && MVA>" << ph1_a0_pu50_wp1  ;
cut2e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && MVA>" << ph1_a1k_pu140_wp1;
cut3e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && MVA>" << hgc_pu140_wp1    ;
ts1->Draw("pt>>s1num(7,30,100)",cut1e.str().c_str());
ts2->Draw("pt>>s2num(7,30,100)",cut2e.str().c_str());
ts3->Draw("pt>>s3num(7,30,100)",cut3e.str().c_str());
ts1->Draw("pt>>s1den(7,30,100)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts2->Draw("pt>>s2den(7,30,100)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts3->Draw("pt>>s3den(7,30,100)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
pad2->Clear();                 

TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( s1num, s1den );
TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( s2num, s2den );
TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( s3num, s3den );

eff1->SetLineColor(1);
eff2->SetLineColor(2);
eff3->SetLineColor(3);
eff1->SetMarkerStyle(20);
eff2->SetMarkerStyle(21);
eff3->SetMarkerStyle(22);
eff1->SetMarkerColor(1);
eff2->SetMarkerColor(2);
eff3->SetMarkerColor(3);
pad2->Clear();
eff3->GetYaxis()->SetRangeUser(0,1);
eff3->GetXaxis()->SetRangeUser(30,100);
eff3->Draw("APL");
eff3->SetTitle("Working Point 1 - <Signal Efficiency> ~ 90%;p_{T};Signal Efficiency");
eff2->SetTitle(";;");
eff1->SetTitle(";;");
eff1->Draw("same PL");
eff2->Draw("same PL");

//c->SaveAs("Eff_wp1.pdf");/i
tl->SetNColumns(2);
tl->SetFillColor(kWhite);
tl->SetLineColor(kWhite);
tl->AddEntry(eff1,"efficiency","lp");
tl->AddEntry(b1num,"fake rate (PH1, Age0, PU50)","lp");
tl->AddEntry(eff2,"efficiency","lp");
tl->AddEntry(b2num,"fake rate (PH1, Age1k, PU140)","lp");
tl->AddEntry(eff3,"efficiency","lp");
tl->AddEntry(b3num,"fake rate (HGC PU140)","lp");
/*
TGaxis *axis = new TGaxis(gPad->GetUxmax(),
		gPad->GetUymin(),
		gPad->GetUxmax(),
		gPad->GetUymax(),
		0.01,10,510,"+L");
axis->Draw();*/
TGaxis::SetMaxDigits(3);
tl->Draw();
c->SaveAs("Result_vs_pt_wp1.pdf");
c->Write("Result_vs_pt_wp1");


c->Clear();




TLegend *tl = new TLegend(0.54,0.4,0.88,0.6);
//dummy2->SetMarkerColor(kWhite);
//tl->AddEntry(dummy2,"WP1 - <Eff>=90%","p");
 TPad *pad1 = new TPad("pad1","",0,0,1,1);
 TPad *pad2 = new TPad("pad2","",0,0,1,1);
pad2->SetFillStyle(4000); //will be transparent
pad2->SetFrameFillStyle(0);

pad1->Draw();
pad1->cd();


///-------------> Result wp1<--------------//
ostringstream cut1;
ostringstream cut2;
ostringstream cut3;
cut1 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && MVA>" <<ph1_a0_pu50_wp2  ;
cut2 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && MVA>" <<ph1_a1k_pu140_wp2;
cut3 <<"abs(etaSC) >1.6 && abs(etaSC)<2.5 && pt>30 && pt<100  && MVA>" <<hgc_pu140_wp2    ;
tb1->Draw("pt>>b1num(7,30,100)",cut1.str().c_str());
tb2->Draw("pt>>b2num(7,30,100)",cut2.str().c_str());
tb3->Draw("pt>>b3num(7,30,100)",cut3.str().c_str());
c->Clear();

//TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( b1num, b1den );
//TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( b2num, b2den );
//TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( b3num, b3den );


//b1num->Scale(1./(891473.));
//b2num->Scale(1./(994522.*(791./801.)));
//b3num->Scale(1./( 585260.*( nHGCFiles/880.)));	

b1num->Scale(1./(n1));
b2num->Scale(1./(n2));
b3num->Scale(1./(n3));	

b1num->SetLineColor(1);
b2num->SetLineColor(2);
b3num->SetLineColor(3);
b1num->SetMarkerStyle(24);
b2num->SetMarkerStyle(25);
b3num->SetMarkerStyle(26);
b1num->SetMarkerColor(1);
b2num->SetMarkerColor(2);
b3num->SetMarkerColor(3);
b3num->SetTitle(";;FakeRate");
b3num->Draw("pl[ Y+ X+");
b3num->GetYaxis()->SetRangeUser(0,0.00153*3);
b3num->GetXaxis()->SetRangeUser(30,100);
b3num->GetXaxis()->SetAxisColor(kWhite);
b3num->GetXaxis()->SetLabelColor(kWhite);
b2num->Draw("Apl same");
b1num->Draw("Apl same");


pad2->Draw();
pad2->cd();

ostringstream cut1e;
ostringstream cut2e;
ostringstream cut3e;
cut1e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && MVA>" << ph1_a0_pu50_wp2  ;
cut2e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && MVA>" << ph1_a1k_pu140_wp2;
cut3e <<"abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100 && matchIndex>-1 && MVA>" << hgc_pu140_wp2    ;
ts1->Draw("pt>>s1num(7,30,100)",cut1e.str().c_str());
ts2->Draw("pt>>s2num(7,30,100)",cut2e.str().c_str());
ts3->Draw("pt>>s3num(7,30,100)",cut3e.str().c_str());
ts1->Draw("pt>>s1den(7,30,100)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts2->Draw("pt>>s2den(7,30,100)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
ts3->Draw("pt>>s3den(7,30,100)","abs(etaEcal) >1.6 && abs(etaEcal)<2.5 && pt>30 && pt<100");
pad2->Clear();                 

TGraphAsymmErrors *eff1 = new TGraphAsymmErrors( s1num, s1den );
TGraphAsymmErrors *eff2 = new TGraphAsymmErrors( s2num, s2den );
TGraphAsymmErrors *eff3 = new TGraphAsymmErrors( s3num, s3den );

eff1->SetLineColor(1);
eff2->SetLineColor(2);
eff3->SetLineColor(3);
eff1->SetMarkerStyle(20);
eff2->SetMarkerStyle(21);
eff3->SetMarkerStyle(22);
eff1->SetMarkerColor(1);
eff2->SetMarkerColor(2);
eff3->SetMarkerColor(3);
pad2->Clear();
eff3->GetYaxis()->SetRangeUser(0,1);
eff3->GetXaxis()->SetRangeUser(30,100);
eff3->Draw("APL");
eff3->SetTitle("Working Point 2 - <Signal Efficiency> ~ 85%;p_{T};Signal Efficiency");
eff2->SetTitle(";;");
eff1->SetTitle(";;");
eff1->Draw("same PL");
eff2->Draw("same PL");

//c->SaveAs("Eff_wp1.pdf");/i
tl->SetNColumns(2);
tl->SetFillColor(kWhite);
tl->SetLineColor(kWhite);
tl->AddEntry(eff1,"efficiency","lp");
tl->AddEntry(b1num,"fake rate (PH1, Age0, PU50)","lp");
tl->AddEntry(eff2,"efficiency","lp");
tl->AddEntry(b2num,"fake rate (PH1, Age1k, PU140)","lp");
tl->AddEntry(eff3,"efficiency","lp");
tl->AddEntry(b3num,"fake rate (HGC PU140)","lp");
/*
TGaxis *axis = new TGaxis(gPad->GetUxmax(),
		gPad->GetUymin(),
		gPad->GetUxmax(),
		gPad->GetUymax(),
		0.01,10,510,"+L");
axis->Draw();*/
TGaxis::SetMaxDigits(3);
tl->Draw();
c->SaveAs("Result_vs_pt_wp2.pdf");
c->Write("Result_vs_pt_wp2");





}
