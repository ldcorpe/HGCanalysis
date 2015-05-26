
#include <memory>
#include "TTree.h"
#include "TFile.h"
#include "TMVA/Reader.h"

void makeTrainTrees(int iO=0){
std::cout << "========== Eta-Reweighing ==========" << std::endl;
std::cout <<"[INFO] option specifed = " << iO << " - " ;
if (iO==0) {
std::cout << " please use: " << std::endl 
 << " 1 - PH1_A0_PU50 " << std::endl
 << " 2 - PH1_A1k_PU140 " << std::endl
 << " 3 - HGC_PU140 " << std::endl;
 return;
}
if (iO==1) std::cout << "PH1_A0_PU50 " << std::endl;
if (iO==2) std::cout << "PH1_A1k_PU140 " << std::endl;
if (iO==3) std::cout << "HGC_PU140 " << std::endl;
if (iO==-1) std::cout <<  " ALL SAMPLES " << std::endl;

if (iO==1 || iO==-1){

 std::cout<<"============== PH1_A0_PU50 ==============" <<std::endl;
	//------- Files to update ----------//
	TFile *_file1 = TFile::Open("HoverE_PH1_A0_PU50.root","update"); 
	//TFile *_file2 = TFile::Open("HoverE_PH1_A1k_PU140.root","udpate"); 
	//TFile *_file3 = TFile::Open("HoverE_HGC.root","update");

	//------- New Files --------------//
	TFile *Louie = new TFile ("HoverE_PH1_A0_PU50_weight.root","RECREATE");
	//TFile *L2 = new TFile ("HoverE_PH1_A1k_PU14_weight.root","RECREATE");
	//TFile *L3 = new TFile ("HoverE_HGC_weight.root","RECREATE");

  //--------- Variables of interest -------------//
	//Old tree
	float pt;
	float etaSC;
	float etaEcal;
	float sigmaIetaIeta;
	float hoe;
	int matchIndex;
	float weight;
	//New tree
	



	// ------- File 1, signal tree ----------//
	TTree *ts1 = (TTree*) _file1->Get("hgg/treePhoton");


/*int N = spt->GetNbinsX();
	for( int ih =1 ; ih<N+1 ; ih++){
	std::cout <<" BIN " << ih << " abs(eta) from " << spt->GetBinLowEdge(ih) << " to " << spt->GetBinLowEdge(ih)+0.1 << ", weight " <<1/( spt->GetBinContent(ih)/spt->Integral(0,N)) << std::endl;
}*/
	ts1->SetBranchAddress("pt"           ,&pt           );
	ts1->SetBranchAddress("etaSC"        ,&etaSC        );
	ts1->SetBranchAddress("etaEcal"      ,&etaEcal      );
	ts1->SetBranchAddress("sigmaIetaIeta",&sigmaIetaIeta);
	ts1->SetBranchAddress("hoe"          ,&hoe          );
	ts1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	/*TTree *ts2 = (TTree*) _file2->Get("hgg/treePhoton");
		TTree *tb2 = (TTree*) _file2->Get("hgg/treeBackground");
		TTree *ts3 = (TTree*) _file3->Get("hgg/treePhoton");
		TTree *tb3 = (TTree*) _file3->Get("hgg/treeBackground");*/
	


	//-----------Signal -----------------
	TTree *tPho = new TTree("treePhoton","example");
	tPho->Branch("weight", &weight, "weight/F");
	tPho->Branch("pt"           ,&pt           ,"pt/F"           );
	tPho->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tPho->Branch("etaEcal"        ,&etaEcal        ,"etaEcal/F"        );
	tPho->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tPho->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tPho->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );


	Long64_t nentries = ts1->GetEntries();
	std::cout << " Process Signal , with " << nentries << " entries "<< std::endl;
	ts1->Draw("abs(etaSC):pt>>spt(10,30,100,10,1.6,2.5)","abs(etaSC)>1.6 && abs(etaSC)<2.5 && pt>30 && matchIndex>-1","colz");

	for (Long64_t i = 0; i < nentries; i++){
		ts1->GetEntry(i);
		int bin = spt->FindBin(pt,fabs(etaSC));
		if (pt<30 || fabs(etaSC)<1.6 || fabs(etaSC)>2.5 || pt>10 ) {weight =0;}
		else {weight=(1./(float) spt->GetBinContent(bin));}///spt->Integral(0,N)) ;}
		if (i %100000==0|| i==nentries-1) std::cout <<  " SIGNAL events processed " << i << ", eta " <<etaSC<<", pt " << pt  <<", bin " << bin << ", w " << weight<< std::endl;
		
		tPho->Fill();
	}

	//--------- File 1, background tres --------------//
	TTree *tb1 = (TTree*) _file1->Get("hgg/treeBackground");
	tb1->SetBranchAddress("pt"           ,&pt           );
	tb1->SetBranchAddress("etaSC"        ,&etaSC        );
	tb1->SetBranchAddress("sigmaIetaIeta",&sigmaIetaIeta);
	tb1->SetBranchAddress("hoe"          ,&hoe          );
	tb1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	//-----------Background -----------------
	TTree *tBkg = new TTree("treeBackground","example");
	tBkg->Branch("weight",          &weight, "weight/F");
	tBkg->Branch("pt"           ,&pt           ,"pt/F"           );
	tBkg->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tBkg->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tBkg->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tBkg->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );

	nentries = tb1->GetEntries();
	tb1->Draw("abs(etaSC):pt>>bpt(10,30,100,10,1.6,2.5)","abs(etaSC)>1.6 && abs(etaSC)<2.5 && pt>30");
	std::cout << " Process Background , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		tb1->GetEntry(i);
		int bin = bpt->FindBin(pt,fabs(etaSC));
		if (pt<30 || fabs(etaSC)<1.6 || fabs(etaSC)>2.5 || pt>100 ) {weight =0;}
		else {weight=(1./(float) bpt->GetBinContent(bin));}///spt->Integral(0,N)) ;}
		if (i %100000==0 || i==nentries-1) std::cout <<  " BACKGROUND events processed " << i << std::endl;
		matchIndex=0;
		tBkg->Fill();
	}


	Louie->Write();
	std::cout << "wrote tree " << "HoverE_PH1_A0_PU50_weight.root"  <<std::endl;
}





//////////////////////////////////////////////////////////////////////////


if (iO==2|| iO==-1){
 std::cout<<"============== PH1_A1k_PU140 ==============" <<std::endl;

	//------- Files to update ----------//
	TFile *_file1 = TFile::Open("HoverE_PH1_A1k_PU140.root","udpate"); 
	//TFile *_file3 = TFile::Open("HoverE_HGC.root","update");

	//------- New Files --------------//
	TFile *Louie = new TFile ("HoverE_PH1_A1k_PU140_weight.root","RECREATE");
	//TFile *L3 = new TFile ("HoverE_HGC_weight.root","RECREATE");

  //--------- Variables of interest -------------//
	//Old tree
	float pt;
	float etaSC;
	float etaEcal;
	float sigmaIetaIeta;
	float hoe;
	int matchIndex;
	//New tree
	float weight;
	

	// ------- File 1, signal tree ----------//
	TTree *ts1 = (TTree*) _file1->Get("hgg/treePhoton");
	ts1->SetBranchAddress("pt"           ,&pt           );
	ts1->SetBranchAddress("etaSC"        ,&etaSC        );
	ts1->SetBranchAddress("etaEcal"      ,&etaEcal      );
	ts1->SetBranchAddress("sigmaIetaIeta",&sigmaIetaIeta);
	ts1->SetBranchAddress("hoe"          ,&hoe          );
	ts1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	/*TTree *ts2 = (TTree*) _file2->Get("hgg/treePhoton");
		TTree *tb2 = (TTree*) _file2->Get("hgg/treeBackground");
		TTree *ts3 = (TTree*) _file3->Get("hgg/treePhoton");
		TTree *tb3 = (TTree*) _file3->Get("hgg/treeBackground");*/



	//-----------Signal -----------------
	TTree *tPho = new TTree("treePhoton","example");
	tPho->Branch("weight", &weight, "weight/F");
	tPho->Branch("pt"           ,&pt           ,"pt/F"           );
	tPho->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tPho->Branch("etaEcal"        ,&etaEcal        ,"etaEcal/F"        );
	tPho->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tPho->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tPho->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );


	Long64_t nentries = ts1->GetEntries();
	ts1->Draw("abs(etaSC):pt>>spt(10,30,100,10,1.6,2.5)","abs(etaSC)>1.6 && abs(etaSC)<2.5 && pt>30 && matchIndex>-1");
	std::cout << " Process Signal , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		ts1->GetEntry(i);
		int bin = spt->FindBin(pt,fabs(etaSC));
		if (pt<30 ||pt>100 || fabs(etaSC)<1.6 || fabs(etaSC)>2.5 ) {weight =0;}
		else {weight=(1./(float) spt->GetBinContent(bin));}///spt->Integral(0,N)) ;}
		if (i %100000==0|| i==nentries-1) std::cout <<  " SIGNAL events processed " << i << std::endl;
		tPho->Fill();
	}

	//--------- File 1, background tres --------------//
	TTree *tb1 = (TTree*) _file1->Get("hgg/treeBackground");
	tb1->SetBranchAddress("pt"           ,&pt           );
	tb1->SetBranchAddress("etaSC"        ,&etaSC        );
	tb1->SetBranchAddress("sigmaIetaIeta",&sigmaIetaIeta);
	tb1->SetBranchAddress("hoe"          ,&hoe          );
	tb1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	//-----------Background -----------------
	TTree *tBkg = new TTree("treeBackground","example");
	tBkg->Branch("weight",          &weight, "weight/F");
	tBkg->Branch("pt"           ,&pt           ,"pt/F"           );
	tBkg->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tBkg->Branch("etaEcal"      ,&etaEcal      ,"etaEcal/F"      );
	tBkg->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tBkg->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tBkg->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );

	nentries = tb1->GetEntries();
	std::cout << " Process Background , with " << nentries << " entries "<< std::endl;
	tb1->Draw("abs(etaSC):pt>>bpt(10,30,100,10,1.6,2.5)","abs(etaSC)>1.6 && abs(etaSC)<2.5 && pt>30");

	for (Long64_t i = 0; i < nentries; i++){
		tb1->GetEntry(i);
		int bin = bpt->FindBin(pt,fabs(etaSC));
		if (pt<30|| pt>100 || fabs(etaSC)<1.6 || fabs(etaSC)>2.5 ) {weight =0;}
		else {weight=(1./(float) bpt->GetBinContent(bin));}///spt->Integral(0,N)) ;}
		if (i %100000==0 || i==nentries-1) std::cout <<  " BACKGROUND events processed " << i << std::endl;
		matchIndex=0;
		tBkg->Fill();
	}


	Louie->Write();
	std::cout << "wrote tree " << "HoverE_PH1_A1k_PU140_weight.root"  <<std::endl;
}
//////////////////////////////////////////////////////////////////////////


if (iO==3 || iO==-1){
 std::cout<<"============== HGC_PU140 ==============" <<std::endl;

	//------- Files to update ----------//
	TFile *_file1 = TFile::Open("HoverE_HGC_015.root","update");

	//------- New Files --------------//
	TFile *Louie = new TFile ("HoverE_HGC_weight.root","RECREATE");

  //--------- Variables of interest -------------//
	//Old tree
	float pt;
	float etaSC;
	float etaEcal;
	float sigmaEtaEta;
	float hoe;
	float lengthCompatibility;
	int matchIndex;
	//New tree
	float weight;
	

	// ------- File 1, signal tree ----------//
	TTree *ts1 = (TTree*) _file1->Get("hgg/treePhoton");
	ts1->SetBranchAddress("pt"           ,&pt           );
	ts1->SetBranchAddress("etaSC"        ,&etaSC        );
	ts1->SetBranchAddress("etaEcal"      ,&etaEcal      );
	ts1->SetBranchAddress("sigmaEtaEta",&sigmaEtaEta);
	ts1->SetBranchAddress("hoe"          ,&hoe          );
	ts1->SetBranchAddress("lengthCompatibility"          ,&lengthCompatibility         );
	ts1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	/*TTree *ts2 = (TTree*) _file2->Get("hgg/treePhoton");
		TTree *tb2 = (TTree*) _file2->Get("hgg/treeBackground");
		TTree *ts3 = (TTree*) _file3->Get("hgg/treePhoton");
		TTree *tb3 = (TTree*) _file3->Get("hgg/treeBackground");*/



	//-----------Signal -----------------
	TTree *tPho = new TTree("treePhoton","example");
	tPho->Branch("weight", &weight, "weight/F");
	tPho->Branch("pt"           ,&pt           ,"pt/F"           );
	tPho->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tPho->Branch("etaEcal"        ,&etaEcal        ,"etaEcal/F"        );
	tPho->Branch("sigmaEtaEta",&sigmaEtaEta,"sigmaEtaEta/F");
	tPho->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tPho->Branch("lengthCompatibility"          ,&lengthCompatibility          ,"lengthCompatibility/F"          );
	tPho->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );


	Long64_t nentries = ts1->GetEntries();
	std::cout << " Process Signal , with " << nentries << " entries "<< std::endl;
	ts1->Draw("abs(etaSC):pt>>spt(10,30,100,10,1.6,3)","abs(etaSC)>1.6 && abs(etaSC)<3 && pt>30 && matchIndex>-1 && lengthCompatibility>-10");

	for (Long64_t i = 0; i < nentries; i++){
		ts1->GetEntry(i);
		int bin = spt->FindBin(pt,fabs(etaSC));
		if (pt<30 ||pt>100 || fabs(etaSC)<1.6 || fabs(etaSC)>3 ) {weight =0;}
		else {weight=(1./(float) spt->GetBinContent(bin));}///spt->Integral(0,N)) ;}
		if (i %100000==0|| i==nentries-1) std::cout <<  " SIGNAL events processed " << i << std::endl;
		tPho->Fill();
	}

	//--------- File 1, background tres --------------//
	TTree *tb1 = (TTree*) _file1->Get("hgg/treeBackground");
	tb1->SetBranchAddress("pt"           ,&pt           );
	tb1->SetBranchAddress("etaSC"        ,&etaSC        );
	tb1->SetBranchAddress("sigmaEtaEta",&sigmaEtaEta);
	tb1->SetBranchAddress("hoe"          ,&hoe          );
	tb1->SetBranchAddress("lengthCompatibility"          ,&lengthCompatibility         );
	tb1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	//-----------Background -----------------
	TTree *tBkg = new TTree("treeBackground","example");
	tBkg->Branch("weight",          &weight, "weight/F");
	tBkg->Branch("pt"           ,&pt           ,"pt/F"           );
	tBkg->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tBkg->Branch("etaEcal"      ,&etaEcal      ,"etaEcal/F"      );
	tBkg->Branch("sigmaEtaEta",&sigmaEtaEta,"sigmaEtaEta/F");
	tBkg->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tBkg->Branch("lengthCompatibility"          ,&lengthCompatibility          ,"lengthCompatibility/F"          );
	tBkg->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );

	nentries = tb1->GetEntries();
	std::cout << " Process Background , with " << nentries << " entries "<< std::endl;
	tb1->Draw("abs(etaSC):pt>>bpt(10,30,100,10,1.6,3)","abs(etaSC)>1.6 && abs(etaSC)<3 && pt>30 && lengthCompatibility>-10 ");

	for (Long64_t i = 0; i < nentries; i++){
		tb1->GetEntry(i);
		int bin = bpt->FindBin(pt,fabs(etaSC));
		if (pt<30 || pt>100 || fabs(etaSC)<1.6 || fabs(etaSC)>3 ) {weight =0;}
		else {weight=(1./(float) bpt->GetBinContent(bin));}///spt->Integral(0,N)) ;}
		if (i %100000==0 || i==nentries-1) std::cout <<  " BACKGROUND events processed " << i << std::endl;
		matchIndex=0;
		tBkg->Fill();
	}


	Louie->Write();
	std::cout << "wrote tree " << "HoverE_HGC_weight.root"  <<std::endl;
}
}
