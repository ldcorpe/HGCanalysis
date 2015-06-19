
#include <memory>
#include "TTree.h"
#include "TFile.h"
#include "TMVA/Reader.h"

void makeMVATrees(int iO=0){
std::cout << "========== Reprocess MVA ==========" << std::endl;
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

	//------- Files to update ----------//
	TFile *_file1 = TFile::Open("HoverE_PH1_A0_PU50.root","update"); 
	//TFile *_file2 = TFile::Open("HoverE_PH1_A1k_PU140.root","udpate"); 
	//TFile *_file3 = TFile::Open("HoverE_HGC.root","update");

	//------- New Files --------------//
	TFile *Louie = new TFile ("HoverE_PH1_A0_PU50_newMVA.root","RECREATE");
	//TFile *L2 = new TFile ("HoverE_PH1_A1k_PU14_newMVA.root","RECREATE");
	//TFile *L3 = new TFile ("HoverE_HGC_newMVA.root","RECREATE");

  //--------- Variables of interest -------------//
	//Old tree
	float pt;
	float etaSC;
	float etaEcal;
	float sigmaIetaIeta;
	float hoe;
	int matchIndex;
	//New tree
	float newMVA;
	
	//-------MVA reader---------//
	std::unique_ptr<TMVA::Reader> Mva_;
	Mva_.reset( new TMVA::Reader( "!Color:Silent" ) );
	Mva_->AddVariable( "sigmaIetaIeta", &sigmaIetaIeta );
	Mva_->AddVariable( "hoe", &hoe );
	Mva_->AddVariable( "etaSC", &etaSC );
	Mva_->AddVariable( "pt", &pt );
	Mva_->BookMVA( "BDTG","weights/PH1_A0_PU50_PhotonID_BDTG.weights.xml");

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
	tPho->Branch("MVA", &newMVA, "MVA/F");
	tPho->Branch("pt"           ,&pt           ,"pt/F"           );
	tPho->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tPho->Branch("etaEcal"        ,&etaEcal        ,"etaEcal/F"        );
	tPho->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tPho->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tPho->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );


	Long64_t nentries = ts1->GetEntries();
	std::cout << " Process Signal , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		ts1->GetEntry(i);
		if (i %100000==0|| i==nentries-1) std::cout <<  " SIGNAL events processed " << i << std::endl;
		newMVA=  Mva_->EvaluateMVA( "BDTG" );;
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
	tBkg->Branch("MVA",          &newMVA, "MVA/F");
	tBkg->Branch("pt"           ,&pt           ,"pt/F"           );
	tBkg->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tBkg->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tBkg->Branch("hoe"          ,&hoe          ,"hoe/F"          );

	nentries = tb1->GetEntries();
	std::cout << " Process Background , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		tb1->GetEntry(i);
		if (i %100000==0 || i==nentries-1) std::cout <<  " BACKGROUND events processed " << i << std::endl;
		newMVA= Mva_->EvaluateMVA( "BDTG" );;
		tBkg->Fill();
	}


	Louie->Write();
	std::cout << "wrote tree " << "HoverE_PH1_A0_PU50_newMVA.root"  <<std::endl;
}





//////////////////////////////////////////////////////////////////////////


if (iO==2|| iO==-1){

	//------- Files to update ----------//
	TFile *_file1 = TFile::Open("HoverE_PH1_A1k_PU140.root","udpate"); 
	//TFile *_file3 = TFile::Open("HoverE_HGC.root","update");

	//------- New Files --------------//
	TFile *Louie = new TFile ("HoverE_PH1_A1k_PU140_newMVA.root","RECREATE");
	//TFile *L3 = new TFile ("HoverE_HGC_newMVA.root","RECREATE");

  //--------- Variables of interest -------------//
	//Old tree
	float pt;
	float etaSC;
	float etaEcal;
	float sigmaIetaIeta;
	float hoe;
	int matchIndex;
	//New tree
	float newMVA;
	
	//-------MVA reader---------//
	std::unique_ptr<TMVA::Reader> Mva_;
	Mva_.reset( new TMVA::Reader( "!Color:Silent" ) );
	Mva_->AddVariable( "sigmaIetaIeta", &sigmaIetaIeta );
	Mva_->AddVariable( "hoe", &hoe );
	Mva_->AddVariable( "etaSC", &etaSC );
	Mva_->AddVariable( "pt", &pt );
	Mva_->BookMVA( "BDTG","weights/PH1_A1k_PU140_PhotonID_BDTG.weights.xml");

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
	tPho->Branch("MVA", &newMVA, "MVA/F");
	tPho->Branch("pt"           ,&pt           ,"pt/F"           );
	tPho->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tPho->Branch("etaEcal"        ,&etaEcal        ,"etaEcal/F"        );
	tPho->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tPho->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tPho->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );


	Long64_t nentries = ts1->GetEntries();
	std::cout << " Process Signal , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		ts1->GetEntry(i);
		if (i %100000==0|| i==nentries-1) std::cout <<  " SIGNAL events processed " << i << std::endl;
		newMVA=  Mva_->EvaluateMVA( "BDTG" );;
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
	tBkg->Branch("MVA",          &newMVA, "MVA/F");
	tBkg->Branch("pt"           ,&pt           ,"pt/F"           );
	tBkg->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tBkg->Branch("etaEcal"      ,&etaEcal      ,"etaEcal/F"      );
	tBkg->Branch("sigmaIetaIeta",&sigmaIetaIeta,"sigmaIetaIeta/F");
	tBkg->Branch("hoe"          ,&hoe          ,"hoe/F"          );

	nentries = tb1->GetEntries();
	std::cout << " Process Background , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		tb1->GetEntry(i);
		if (i %100000==0 || i==nentries-1) std::cout <<  " BACKGROUND events processed " << i << std::endl;
		newMVA= Mva_->EvaluateMVA( "BDTG" );;
		tBkg->Fill();
	}


	Louie->Write();
	std::cout << "wrote tree " << "HoverE_PH1_A1k_PU140_newMVA.root"  <<std::endl;
}
//////////////////////////////////////////////////////////////////////////


if (iO==3 || iO==-1){

	//------- Files to update ----------//
	TFile *_file1 = TFile::Open("HoverE_HGC_015.root","update");
	//TFile *_file1 = TFile::Open("HoverEHggPU140_gamJet_ssz3_Isolation.root","update");

	//------- New Files --------------//
	TFile *Louie = new TFile ("HoverE_HGC_newMVA.root","RECREATE");

  //--------- Variables of interest -------------//
	//Old tree
	float pt;
	float etaSC;
	float etaEcal;
	float sigmaEtaEta;
	float hoe;
	float lengthCompatibility;
	float trkIso;
	float ecalIso;
	float hcalIso;
	int matchIndex;
	//New tree
	float newMVA;
	
	//-------MVA reader---------//
	std::unique_ptr<TMVA::Reader> Mva_;
	Mva_.reset( new TMVA::Reader( "!Color:Silent" ) );
	Mva_->AddVariable( "sigmaEtaEta", &sigmaEtaEta );
	Mva_->AddVariable( "hoe", &hoe );
	Mva_->AddVariable( "lengthCompatibility", &lengthCompatibility );
	//Mva_->AddVariable( "trkIso", &trkIso );
	//Mva_->AddVariable( "ecalIso", &ecalIso );
	//Mva_->AddVariable( "hcalIso", &hcalIso );
	Mva_->AddVariable( "etaSC", &etaSC );
	Mva_->AddVariable( "pt", &pt );
	Mva_->BookMVA( "BDTG","weights/HGC_PhotonID_BDTG.weights.xml");

	// ------- File 1, signal tree ----------//
	TTree *ts1 = (TTree*) _file1->Get("hgg/treePhoton");
	ts1->SetBranchAddress("pt"           ,&pt           );
	ts1->SetBranchAddress("etaSC"        ,&etaSC        );
	ts1->SetBranchAddress("etaEcal"      ,&etaEcal      );
	ts1->SetBranchAddress("sigmaEtaEta",&sigmaEtaEta);
	ts1->SetBranchAddress("hoe"          ,&hoe          );
	ts1->SetBranchAddress("lengthCompatibility"          ,&lengthCompatibility         );
//	ts1->SetBranchAddress("trkIso"          ,&trkIso         );
//	ts1->SetBranchAddress("ecalIso"          ,&ecalIso         );
//	ts1->SetBranchAddress("hcalIso"          ,&hcalIso         );
	ts1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	/*TTree *ts2 = (TTree*) _file2->Get("hgg/treePhoton");
		TTree *tb2 = (TTree*) _file2->Get("hgg/treeBackground");
		TTree *ts3 = (TTree*) _file3->Get("hgg/treePhoton");
		TTree *tb3 = (TTree*) _file3->Get("hgg/treeBackground");*/



	//-----------Signal -----------------
	TTree *tPho = new TTree("treePhoton","example");
	tPho->Branch("MVA", &newMVA, "MVA/F");
	tPho->Branch("pt"           ,&pt           ,"pt/F"           );
	tPho->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tPho->Branch("etaEcal"        ,&etaEcal        ,"etaEcal/F"        );
	tPho->Branch("sigmaEtaEta",&sigmaEtaEta,"sigmaEtaEta/F");
	tPho->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tPho->Branch("lengthCompatibility"          ,&lengthCompatibility          ,"lengthCompatibility/F"          );
//	tPho->Branch("trkIso"          ,&trkIso          ,"trkIso/F"          );
//	tPho->Branch("ecalIso"          ,&ecalIso          ,"ecalIso/F"          );
//	tPho->Branch("hcalIso"          ,&hcalIso          ,"hcalIso/F"          );
	tPho->Branch("matchIndex"   ,&matchIndex   ,"matchIndex/I"   );


	Long64_t nentries = ts1->GetEntries();
	std::cout << " Process Signal , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		ts1->GetEntry(i);
		if (i %100000==0|| i==nentries-1) std::cout <<  " SIGNAL events processed " << i << std::endl;
		newMVA=  Mva_->EvaluateMVA( "BDTG" );;
	//	std::cout << " trkIso " << trkIso << ", ecalIso " << ecalIso << ", kl " << lengthCompatibility<<  std::endl;
		tPho->Fill();
	}

	//--------- File 1, background tres --------------//
	TTree *tb1 = (TTree*) _file1->Get("hgg/treeBackground");
	tb1->SetBranchAddress("pt"           ,&pt           );
	tb1->SetBranchAddress("etaSC"        ,&etaSC        );
	tb1->SetBranchAddress("sigmaEtaEta",&sigmaEtaEta);
	tb1->SetBranchAddress("hoe"          ,&hoe          );
	tb1->SetBranchAddress("lengthCompatibility"          ,&lengthCompatibility         );
//	tb1->SetBranchAddress("trkIso"          ,&trkIso         );
//	tb1->SetBranchAddress("ecalIso"          ,&ecalIso         );
//	tb1->SetBranchAddress("hcalIso"          ,&hcalIso         );
	tb1->SetBranchAddress("matchIndex"   ,&matchIndex   );

	//-----------Background -----------------
	TTree *tBkg = new TTree("treeBackground","example");
	tBkg->Branch("MVA",          &newMVA, "MVA/F");
	tBkg->Branch("pt"           ,&pt           ,"pt/F"           );
	tBkg->Branch("etaSC"        ,&etaSC        ,"etaSC/F"        );
	tBkg->Branch("etaEcal"      ,&etaEcal      ,"etaEcal/F"      );
	tBkg->Branch("sigmaEtaEta",&sigmaEtaEta,"sigmaEtaEta/F");
	tBkg->Branch("hoe"          ,&hoe          ,"hoe/F"          );
	tBkg->Branch("lengthCompatibility"          ,&lengthCompatibility          ,"lengthCompatibility/F"          );
//	tBkg->Branch("trkIso"          ,&trkIso          ,"trkIso/F"          );
//	tBkg->Branch("ecalIso"          ,&ecalIso          ,"ecalIso/F"          );
//	tBkg->Branch("hcalIso"          ,&hcalIso          ,"hcalIso/F"          );

	nentries = tb1->GetEntries();
	std::cout << " Process Background , with " << nentries << " entries "<< std::endl;

	for (Long64_t i = 0; i < nentries; i++){
		tb1->GetEntry(i);
		if (i %100000==0 || i==nentries-1) std::cout <<  " BACKGROUND events processed " << i << std::endl;
		newMVA= Mva_->EvaluateMVA( "BDTG" );;
		tBkg->Fill();
	}


	Louie->Write();
	std::cout << "wrote tree " << "HoverE_HGC_newMVA.root"  <<std::endl;
}
}
