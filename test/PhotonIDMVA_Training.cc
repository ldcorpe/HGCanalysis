///==== include ====
#include "TFile.h"
#include "TChain.h"
#include "TMinuit.h"
#include <sstream>
#include <iostream>
#include "TMVA/Factory.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//#include "TMVAGui.C"

#endif
using namespace std;
// --------- MAIN -------------------
void PhotonIDMVA_Training(TString Nevent="10000", TString Level="PhotonID", bool skipBadEvents=true)
{
  // you must define $WORKSPACE first
  TString path;
	path="";
  
	bool useDiphotonPt = 0;
  bool usePhotonsPt = true;
  int nEvents= std::atoi(Nevent.Data());







  
  //TFile *inputB1 = TFile::Open(path + "HoverEHggPU140_QCD_ssz3_8.root ");
  //TFile *inputS1 = TFile::Open(path + "HoverE_EMPrefilter.root ");
  TFile *inputS1 = TFile::Open(path + "HoverE_HGC_weight.root");
 // TFile *inputS1 = TFile::Open(path + "HoverE_HGC.root ");
 /* TFile *inputB2 = TFile::Open(path + "output_GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent"+Nevent+"_histos.root     ");
  TFile *inputB3 = TFile::Open(path + "output_GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent"+Nevent+"_histos.root ");
  TFile *inputS1 = TFile::Open(path + "output_GluGluToHToGG_M-125_13TeV-powheg-pythia6_numEvent"+Nevent+"_histos.root                ");
  TFile *inputS2 = TFile::Open(path + "output_TTbarH_HToGG_M-125_13TeV_amcatnlo-pythia8-tauola_numEvent"+Nevent+"_histos.root        ");
  TFile *inputS3 = TFile::Open(path + "output_VBF_HToGG_M-125_13TeV-powheg-pythia6_numEvent"+Nevent+"_histos.root                    ");
  TFile *inputS4 = TFile::Open(path + "output_WH_ZH_HToGG_M-125_13TeV_pythia6_numEvent"+Nevent+"_histos.root                         ");*/
  
  //TTree *treeB1 = (TTree*)inputB1->Get("hgg/treeBackground");
 // TTree *treeB1 = (TTree*)inputS1->Get("hgg/treeBackground");
 // TTree *treeS1 = (TTree*)inputS1->Get("hgg/treePhoton");
  TTree *treeB1 = (TTree*)inputS1->Get("treeBackground");
  TTree *treeS1 = (TTree*)inputS1->Get("treePhoton");
 /* TTree *treeB2 = (TTree*)inputB2->Get(Level +"MVADumperNew/trees/gamJet_13TeV_GoodPhotonIDNew");
  TTree *treeB3 = (TTree*)inputB3->Get(Level +"MVADumperNew/trees/gamJet_13TeV_GoodPhotonIDNew");
  TTree *treeS1 = (TTree*)inputS1->Get(Level +"MVADumperNew/trees/ggh_m125_13TeV_GoodPhotonIDNew");
  TTree *treeS2 = (TTree*)inputS2->Get(Level +"MVADumperNew/trees/tth_m125_13TeV_GoodPhotonIDNew");
  TTree *treeS3 = (TTree*)inputS3->Get(Level +"MVADumperNew/trees/vbf_m125_13TeV_GoodPhotonIDNew");
  TTree *treeS4 = (TTree*)inputS4->Get(Level +"MVADumperNew/trees/wzh_m125_13TeV_GoodPhotonIDNew");*/


	// Declaration of leaf types
	float   weightS[5];
	float   weightB[4];
/*	float   sigmarv;
	float   sigmawv;
	float   vtxprob;
	float   leadptom;
	float   subleadptom;
	float   leadeta;
	float   subleadeta;
	float   CosPhi;
	float   leadmva;
	float   subleadmva;



if(Level =="PhotonID") {
		treeS1->SetBranchAddress("sigmarv   "      , &sigmarv          );
		treeS1->SetBranchAddress("sigmawv   "      , &sigmawv          );
		treeS1->SetBranchAddress("vtxprob   "      , &vtxprob          );
		treeS1->SetBranchAddress("leadptom  "      , &leadptom         );
		treeS1->SetBranchAddress("subleadptom"     , &subleadptom      );
		treeS1->SetBranchAddress("leadeta   "      , &leadeta          );
		treeS1->SetBranchAddress("subleadeta"      , &subleadeta       );
		treeS1->SetBranchAddress("CosPhi    "      , &CosPhi           );
		treeS1->SetBranchAddress("leadmva   "      , &leadmva          );
		treeS1->SetBranchAddress("subleadmva"      , &subleadmva       );
	}*/


   /* weightB[0]=0;
    weightB[1]=4746.0/nEvents;// DyJetsToLL
    weightB[2]=17180.0*0.0379/nEvents;//gamJets pt>40
    weightB[3]=145400.0*0.001776/nEvents;//gamJets pt  in 20->40

    weightS[0]=0;
    weightS[1]=43.92*2.28e-3/nEvents; //ggH
    weightS[2]=3.748*2.28e-3/nEvents; //ttH
    weightS[3]=2.2496*2.28e-3/nEvents; //VBH
    weightS[4]=0.5608*2.28e-3/nEvents; //WZH*/

		
		std::cout << "DEBUG treeS1 " << treeS1->GetEntries()<< ", B1 " << treeB1->GetEntries() << std::endl;

    weightS[1]=1;
    weightB[1]=1;



	// Create a new root output file.
	string outputFileName;
	if(Level =="PhotonID"){	
		outputFileName = "HGC_PhotonID";
	}

	// -- reader
	TFile* outputFile = TFile::Open((outputFileName+".root").c_str(), "RECREATE" );
	TMVA::Factory* factory = new TMVA::Factory(outputFileName.c_str(), outputFile,
			"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
	// -- variables
	if(Level =="PhotonID") {
    factory->AddVariable("sigmaEtaEta" );
    factory->AddVariable("hoe" );
    factory->AddVariable("lengthCompatibility" );
   // factory->AddVariable("trkIso" );
    //factory->AddVariable("ecalIso" );
    //factory->AddVariable("hcalIso" );
    factory->AddVariable("etaSC" );
    factory->AddVariable("pt" );
	}


	//event weights per tree (see below for setting event-wise weights)
	Double_t signalWeight = 1.0;
	Double_t backgroundWeight = 1.0;


	// ====== register trees ====================================================
	factory->AddSignalTree    ( treeS1,  weightS[1] );
//	factory->AddSignalTree    ( treeS2,  weightS[2] );
//	factory->AddSignalTree    ( treeS3,  weightS[3] );
//	factory->AddSignalTree    ( treeS4,  weightS[4] );
	factory->AddBackgroundTree( treeB1,  weightB[1] );
//	factory->AddBackgroundTree( treeB2,  weightB[2] );
//	factory->AddBackgroundTree( treeB3,  weightB[3] );


	// == supress the the negative points on the input variables
	// == this high correlation between variables
	TCut mycuts ="";// " leadPho_PToM > (60./120.) && sublPho_PToM> (30./120.)";
	TCut mycutb ="";// " leadPho_PToM> (60./120.) && sublPho_PToM> (30./120.)";
		mycuts ="pt>30 && pt<100 && abs(etaSC)>1.6 && abs(etaSC)<3 && matchIndex>-1 && lengthCompatibility>-20 && hoe<0.5";// Might need a better discriminant LC
		mycutb ="pt>30 && pt <100 &&  abs(etaSC)>1.6 && abs(etaSC)<3 && lengthCompatibility>-20 && hoe<0.5";// Might need a better discriminant LC
	//	mycutb ="pt>30 &&  abs(etaSC)>1.6 && abs(etaSC)<3 && lengthCompatibility>-10";// 
	//	//mycuts ="pt>30 &&  abs(etaSC)>1.6 && abs(etaSC)<2.8 && sigmaEtaEta<0.006";// Might need a better discriminant LC
	//	mycutb ="pt>30 &&  abs(etaSC)>1.6 && abs(etaSC)<2.8 && sigmaEtaEta<0.006";// 
	factory->SetSignalWeightExpression("weight");
	factory->SetBackgroundWeightExpression("weight");

	// tell the factory to use all remaining events in the trees after training for testing:
	factory->PrepareTrainingAndTestTree( mycuts, mycutb,
			"SplitMode=Random:NormMode=NumEvents:!V" );
	// Boosted Decision Trees: use BDTG ( Gradient Boost )
	factory->BookMethod( TMVA::Types::kBDT, "BDTG",
			"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5:MaxDepth=3:NegWeightTreatment=IgnoreNegWeights" );
//	factory->BookMethod( TMVA::Types::kBDT, "BDT",
	//	"!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
//	"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=15:MaxDepth=5" );
	// book Cuts
	//factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
	// "H:!V:FitMethod=GA:CutRangeMin[0]=20:CutRangeMax[0]=500:CutRangeMin[1]=20:CutRangeMax[1]=500:VarProp=FSmart:VarProp[4]=FMin:EffSel:Steps=30:Cycles=3:PopSize=500:SC_steps=10:SC_rate=5:SC_factor=0.95" );
	// ---- Now you can tell the factory to train, test, and evaluate the MVAs
	// Train MVAs using the set of training events
	factory->TrainAllMethods();
	// ---- Evaluate all MVAs using the set of test events
	factory->TestAllMethods();
	// ----- Evaluate and compare performance of all configured MVAs
	factory->EvaluateAllMethods();
	// --------------------------------------------------------------
	// Save the output
	outputFile->Close();
	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;
	delete factory;

	//	if (!gROOT->IsBatch()) TMVAGui( (outputFileName+".root").c_str() );
}
