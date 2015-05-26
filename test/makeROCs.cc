#include "TH1.h"
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TString.h"
#include <sstream>
#include <utility>

//TGraph* makeROC2( TH1 *sig, TH1 *bkg, TLegend *leg, int tc )
std::pair<TGraph*,std::pair<float,float> > makeROC2( TH1 *sig, TH1 *bkg, TLegend *leg, int tc )
{  
   float cut95=-999.;
   float cut90=-999.;

    gStyle->SetOptStat( 0 );
    sig->Sumw2( 0 );
    bkg->Sumw2( 0 );
    TGraph *g_e = new TGraph();

    g_e->SetMarkerColor( tc );
    g_e->SetLineColor( tc );
    g_e->SetMarkerSize( 1 );
    //g_e->SetMarkerStyle(21);

    bool printFlag = 0;
    bool printFlag2 = 0;
    int N0 = bkg->GetNbinsX();
    std::cout << "debug " << N0 << std::endl;
    const unsigned int N = N0;
    float y[N];
    float signal[N];
    int i = 0;

    for( unsigned int n = 0; n < N+1; n++ ) {
        signal[n] =  sig->Integral( 0, n ) / sig->Integral( 0, N );
        y[n]     = ( bkg->Integral( 0, n ) / bkg->Integral( 0, N ) );

        if( y[n] == 0 ) {
            i++;
            continue;
        }
       // std::cout << "x " << sig->GetXaxis()->GetBinCenter( n ) << ", sg " <<1- signal[n] << ", bkg " << y[n]<<  std::endl;
        if( signal[n] > 0.1 && printFlag == 0 ) {
                      cut90= sig->GetXaxis()->GetBinCenter( n );
            std::cout << "Better than 90 pc eff.("
                      << 1-signal[n] << "). Bin "
                      << n << ", x value "
                      << cut90
                      << ", Backg rejection "
                      << y[n] 
                      << std::endl;

            printFlag = 1;
        }
        if( signal[n] > 0.05 && printFlag2 == 0 ) {
                      cut95= sig->GetXaxis()->GetBinCenter( n );
            std::cout << "Better than 95 pc eff.("
                      << signal[n] << "). Bin "
                      << n << ", x value "
                      << cut95 
                      << ", Backg rejection "
                      << y[n] 
                      << std::endl;

            printFlag2 = 1;
        }
        //g_e->SetPoint(n-i,1-signal[n], y[n]);
        g_e->SetPoint( n - i, 1-signal[n] , y[n]);

        //use 1-n for signal because I am integrating from bin 0 to bin N, which is techincally asking how to reject signal.
    }
    g_e->SetTitle( Form( "ROCs_%i", tc ) );


    //g_e->GetXaxis()->SetRangeUser(0,1);
    //g_e->GetXaxis()->SetTitle("background efficiency");
    //g_e->GetYaxis()->SetRangeUser(0,1);
    //g_e->GetYaxis()->SetTitle("signal efficiency");
    //
    //if (tc ==1){
    //  g_e->Draw("apl same");
    //} else {
    //  g_e->Draw("pl same");
    //}
    //g_g->Draw("pl same");:w
    //
    //g_eg->Draw("pl same");
    //std::string name2 = name.str().replace(name.str().find("_signal_h"), 6, "");

    if( tc == 1 ) {
        leg->AddEntry( g_e, "current FLASHgg MVA", "lep" );
    } else {
        leg->AddEntry( g_e, "retrained MVA", "lep" );
    }
    std::pair<float, float> values = std::make_pair(cut90,cut95);
    std::pair < TGraph*, std::pair < float, float> > result = std::make_pair (g_e, values);
return result;
    //leg->Draw();
}

std::pair<double,double> workingPoint(TGraph *gr)
{
	int nPoints = gr->GetN();
	int l=0;

	TGraph *dr = new TGraph();
	double *x= gr->GetX();
	double *y= gr->GetY();
	std::pair<double,double> corner;
	corner.first=0;
	corner.second=1;
	std::pair<double,double> wp;
	double maxR =9999.;
	for( int i =1 ; i< nPoints-1 ; i++){

		double dx = corner.first-x[i];
		double dy = corner.second-y[i];
		dx= dx*dx;
		dy= dy*dy;
		double dr = sqrt(dx+dy);

		if (dr< maxR){
			maxR=dr;
			wp.first=x[i];
			wp.second=y[i];
		}
	}
	return wp;
}

TGraph *makeROC( TH1 *sig, TH1 *bkg, TLegend *leg, int tc )
{

    gStyle->SetOptStat( 0 );
    sig->Sumw2( 0 );
    bkg->Sumw2( 0 );
    TGraph *g_e = new TGraph();

    g_e->SetMarkerColor( tc );
    g_e->SetLineColor( tc );
    g_e->SetMarkerSize( 1 );
    //g_e->SetMarkerStyle(21);

    bool printFlag = 0;
    int N0 = bkg->GetNbinsX();
    const unsigned int N = N0;
    float y[N];
    float signal[N];
    int i = 0;
    for( unsigned int n = 0; n < N; n++ ) {


        signal[n] = sig->Integral( 0, n ) / sig->Integral( 0, N );

        y[n] = ( bkg->Integral( 0, n ) / bkg->Integral( 0, N ) );

        if( y[n] == 0 ) {
            i++;
            continue;
        }
        if( signal[n] > 0.98 && printFlag == 0 ) {

            std::cout << "Better than 97 pc eff.(" << signal[n] << "). Bin " << n << ", x value " << n *sig->GetBinWidth( n ) << std::endl;
            printFlag = 1;
        }
        g_e->SetPoint( n - i,  signal[n], 1-y[n] );
        //use 1-n for signal because I am integrating from bin 0 to bin N, which is techincally asking how to reject signal.


    }
    g_e->SetTitle( "ROCs" );



    g_e->GetXaxis()->SetTitle( "signal efficiency" );
    g_e->GetXaxis()->SetRangeUser( -0.1, 1.1 );
    g_e->GetYaxis()->SetRangeUser( -0.1, 1.1 );
    g_e->GetYaxis()->SetTitle( "background rejection" );
    if( tc == 1 ) {
        g_e->Draw( "apl same" );
    } else {
        g_e->Draw( "pl same" );
    }
    //g_g->Draw("pl same");:w
    //
    //g_eg->Draw("pl same");
    //std::string name2 = name.str().replace(name.str().find("_signal_h"), 6, "");

    if( tc == 1 ) {
        leg->AddEntry( g_e, "current FLASHgg MVA", "lep" );
    } else {
        leg->AddEntry( g_e, "retrained MVA", "lep" );
    }

    //leg->Draw();
    return g_e;
    }


void makeROCs( TString Nevent = "10000", TString Level = "VBF" )
{

	int nEvents = std::atoi(Nevent.Data()); 
    TString path;
    //if(Level =="VBF") {path="${WORKSPACE}/test_vbfmva_compare/";}
    //else if(Level =="VBFDiPhoDiJet") {
    path = "${WORKSPACE}/test_diphodijet_compare/";
    //}
//	else {
//	std::cout << " Please specify option 2 as 'VBF' or 'VBFDiPhoDiJet' only" << std::endl;
//	return;
//	}

    //
    TFile *f_sig = TFile::Open( path + "output_VBF_HToGG_M-125_13TeV-powheg-pythia6_numEvent" + Nevent + "_histos.root" );

    //TTree *tree_old_sig = (TTree*) f_sig->Get("VBFMVADumper/trees/vbf_m125_13TeV_All");
    //TTree *tree_old_sig = (TTree*) f_sig->Get("VBFMVADumper/trees/vbf_m125_13TeV_All");

   /* TH1F *old_sig = ( TH1F * ) f_sig->Get( Level + "MVADumper/histograms/vbf_m125_13TeV_Good" + Level + "Leg" + Level + "MVAValue" );
    std::cout << Level + "MVADumper/histograms/vbf_m125_13TeV_Good" + Level + "Leg" + Level + "MVAValue" << std::endl ;
    TH1F *new_sig = ( TH1F * ) f_sig->Get( Level + "MVADumperNew/histograms/vbf_m125_13TeV_Good" + Level + "New" + Level + "MVAValue" );

    TFile *f_bkg = TFile::Open( path + "output_GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent" + Nevent + "_histos.root" );

    TH1F *old_bkg = ( TH1F * ) f_bkg->Get( Level + "MVADumper/histograms/gamJet_13TeV_Good" + Level + "Leg" + Level + "MVAValue" );
    TH1F *new_bkg = ( TH1F * ) f_bkg->Get( Level + "MVADumperNew/histograms/gamJet_13TeV_Good" + Level + "New" + Level + "MVAValue" );*/
	
    
    TFile *f_bkg1 = TFile::Open(path + "output_DYJetsToLL_M-50_13TeV-madgraph-pythia8_numEvent"+Nevent+"_histos.root                  ");
	TFile *f_bkg2 = TFile::Open(path + "output_GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent"+Nevent+"_histos.root     ");
	TFile *f_bkg3 = TFile::Open(path + "output_GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent"+Nevent+"_histos.root ");
	TFile *f_sig1 = TFile::Open(path + "output_GluGluToHToGG_M-125_13TeV-powheg-pythia6_numEvent"+Nevent+"_histos.root                ");
	// TFile *f_sig2 = TFile::Open(path + "output_TTbarH_HToGG_M-125_13TeV_amcatnlo-pythia8-tauola_numEvent"+Nevent+"_histos.root        ");
	TFile *f_sig3 = TFile::Open(path + "output_VBF_HToGG_M-125_13TeV-powheg-pythia6_numEvent"+Nevent+"_histos.root                    ");
	//TFile *f_sig4 = TFile::Open(path + "output_WH_ZH_HToGG_M-125_13TeV_pythia6_numEvent"+Nevent+"_histos.root                         ");


	TH1F *old_bkg1 = (TH1F*)f_bkg1->Get(Level +"MVADumper/histograms/dyJets_13TeV_Good"+Level+"Leg"+Level+"MVAValue");
	TH1F *old_bkg2 = (TH1F*)f_bkg2->Get(Level +"MVADumper/histograms/gamJet_13TeV_Good"+Level+"Leg"+Level+"MVAValue");
	TH1F *old_bkg3 = (TH1F*)f_bkg3->Get(Level +"MVADumper/histograms/gamJet_13TeV_Good"+Level+"Leg"+Level+"MVAValue");
	TH1F *old_sig1 = (TH1F*)f_sig1->Get(Level +"MVADumper/histograms/ggh_m125_13TeV_Good"+Level+"Leg"+Level+"MVAValue");
	// TH1F *old_sig2 = (TH1F*)f_sig2->Get(Level +"MVADumper/histograms/tth_m125_13TeV_Good"+Level+"Leg"+Level+"MVAValue");
	TH1F *old_sig3 = (TH1F*)f_sig3->Get(Level +"MVADumper/histograms/vbf_m125_13TeV_Good"+Level+"Leg"+Level+"MVAValue");
	// TH1F *old_sig4 = (TH1F*)f_sig4->Get(Level +"MVADumper/histograms/wzh_m125_13TeV_Good"+Level+"Leg"+Level+"MVAValue");

	TH1F *new_bkg1 = (TH1F*)f_bkg1->Get(Level +"MVADumperNew/histograms/dyJets_13TeV_Good"+Level+"New"+Level+"MVAValue");
	TH1F *new_bkg2 = (TH1F*)f_bkg2->Get(Level +"MVADumperNew/histograms/gamJet_13TeV_Good"+Level+"New"+Level+"MVAValue");
	TH1F *new_bkg3 = (TH1F*)f_bkg3->Get(Level +"MVADumperNew/histograms/gamJet_13TeV_Good"+Level+"New"+Level+"MVAValue");
	TH1F *new_sig1 = (TH1F*)f_sig1->Get(Level +"MVADumperNew/histograms/ggh_m125_13TeV_Good"+Level+"New"+Level+"MVAValue");
	// TH1F *new_sig2 = (TH1F*)f_sig2->Get(Level +"MVADumperNew/histograms/tth_m125_13TeV_Good"+Level+"New"+Level+"MVAValue");
	TH1F *new_sig3 = (TH1F*)f_sig3->Get(Level +"MVADumperNew/histograms/vbf_m125_13TeV_Good"+Level+"New"+Level+"MVAValue");
	// TH1F *new_sig4 = (TH1F*)f_sig4->Get(Level +"MVADumperNew/histograms/wzh_m125_13TeV_Good"+Level+"New"+Level+"MVAValue");

	TH1F *old_sig = new TH1F("old_sig","old_sig",old_sig1->GetNbinsX(),-1.2,1.2);
	TH1F *old_bkg = new TH1F("old_bkg","old_bkg",old_bkg1->GetNbinsX(),-1.2,1.2);
	TH1F *new_sig = new TH1F("new_sig","new_sig",new_sig1->GetNbinsX(),-1.2,1.2);
	TH1F *new_bkg = new TH1F("new_bkg","new_bkg",new_bkg1->GetNbinsX(),-1.2,1.2);
    
	float wS[5];	
	float wB[4];

	wB[0]=0;
	wB[1]=4746.0/(nEvents);// DyJetsToLL
	wB[2]=17180.0*0.0379/(nEvents);//gamJets pt>40
	wB[3]=145400.0*0.001776/(nEvents);//gamJets pt  in 20->40

	wS[0]=0;
	wS[1]=43.92*2.28e-3/(nEvents); //ggH
	wS[2]=3.748*2.28e-3 /(nEvents); //ttH
	wS[3]=2.2496*2.28e-3/(nEvents); //VBH
	wS[4]=0.5608*2.28e-3/(nEvents); //WZH

	/*	wB[0]=0;
			wB[1]=1;
			wB[2]=1;
			wB[3]=1;

			wS[0]=1;
			wS[1]=1;
			wS[2]=1;
			wS[3]=1;
			wS[4]=1;*/

	new_bkg->Add(new_bkg1,wB[1]);
	new_bkg->Add(new_bkg2,wB[2]);
	new_bkg->Add(new_bkg3,wB[3]);
	new_bkg->Add(new_sig1,wS[1]); //ggH counts as background here

	//	new_sig->Add(new_sig2,wS[2]);
	new_sig->Add(new_sig3,wS[3]);
	//	new_sig->Add(new_sig4,wS[4]);

	old_bkg->Add(old_bkg1,wB[1]);
	old_bkg->Add(old_bkg2,wB[2]);
	old_bkg->Add(old_bkg3,wB[3]);
	old_bkg->Add(old_sig1,wS[1]);//ggH counts as background for VBF

	//	old_sig->Add(old_sig2,wS[2]);
	old_sig->Add(old_sig3,wS[3]);
	//	old_sig->Add(old_sig4,wS[4]);
    
    
	std::cout << "plotting" << std::endl;
	std::cout << " --> obj::old_bkg::" << old_bkg << " entries " << old_bkg->GetEntries() <<std::endl;
	std::cout << " --> obj::new_bkg::" << new_bkg << " entries " << new_bkg->GetEntries()   <<std::endl;
	std::cout << " --> obj::old_sig::" << old_sig << " entries " << old_sig->GetEntries()   <<std::endl;
	std::cout << " --> obj::new_sig::" << new_sig << " entries " << new_sig->GetEntries()   <<std::endl;
    
    TCanvas *c1 = new TCanvas( "c1", "", 600, 600 );
    c1->cd();
//	gStyle->SetNdivisions(810,"X");
//	gStyle->SetNdivisions(810,"Y");
    c1->SetGridx();
    c1->SetGridy();
    c1->Modified();
    c1->Update();

// old_bkg->Draw();
    std::vector<std::string> xxx;//={"signal","hgg"};
    xxx.push_back( "current MVA" );
    xxx.push_back( "retrained MVA" );

    int itc = 1;
    TLegend *leg = new TLegend( 0.3, 0.3, 0.65, 0.5 );

    TMultiGraph *mul  = new TMultiGraph();

	TGraph *oldg = makeROC2( old_sig, old_bkg,leg, itc);


    mul->Add( oldg );
	int itc0 = itc;
    itc++;
	
    TGraph *newg = makeROC2( new_sig, new_bkg,leg, itc);
	
    mul->Add(newg);
	int itc1 = itc;
    mul->SetTitle( ";background efficiency;signal efficiency" );
    mul->Draw( "APL*" );
	
    TLine *tl = new TLine(0,0,1,1);
	std::pair<double,double> wp = workingPoint(oldg);
	tl->SetLineColor(itc0);
	tl->DrawLine(0.,1.,wp.first,wp.second);
	std::pair<double,double> wp1 = workingPoint(newg);
	tl->SetLineColor(itc1);
	tl->DrawLine(0.,1.,wp1.first,wp1.second);
	
	tl->SetLineStyle(8);
	tl->SetLineColor(kBlue);
//	tl->DrawLine(wp.first,wp.second,wp1.first,wp1.second);

	TH1F *dummy = new TH1F("","",100,0,1);
	dummy->SetMarkerColor(kWhite);
	dummy->SetLineColor(kBlue);
	dummy->SetLineStyle(8);
	
	double r =sqrt((wp.first - wp1.first)*(wp.first - wp1.first) + (wp.second - wp1.second)*(wp.second - wp1.second)); 

	ostringstream text;
	text << " Improvement : " << r ;
//	leg->AddEntry(dummy,text.str().c_str(),"lp");
    leg->Draw();


	c1->SaveAs("plots/"+Level+"_ROCs_"+Nevent+".pdf");

	std::cout << " DIST " << r << std::endl;

    //TCanvas *c2 = new TCanvas("c2","",600,600);
    //c2->cd();
    //new_bkg->SetLineColor(kBlue);
    //new_sig->SetLineColor(kRed);
    //new_bkg->Draw();
    //new_sig->Draw("same");


}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

