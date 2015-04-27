// -*- C++ -*-
//
// Package:    BasicHggAnalyser
// Class:      BasicHggAnalyser
// 
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

//
// class declaration
//
// limit dR allowed for geometrical matching of reco/gen particles
float dRLim =0.5;

// information to be loaded into TTree
struct infoTruth_t {

	float  eReco_over_eTrue;
	float eta;
	float etaEcal;
	float etaSC;
	float etaSeed;
	float phi;
	float phiEcal;
	float phiSC;
	float phiSeed;
	int matchIndex;
	float pt;
	float ptReco;
	float eTrue;
	float eReco;
	float eSeed;
	float eClust1;
	float eClust2;
	float eClust3;
	float dEtaClust2;
	float dPhiClust2;
	float dEtaClust3;
	float dPhiClust3;
	float eSeed_over_eReco;
	float phiWidth;
	float etaWidth;
	int clustersSize;
	int nClusters09;

};

struct info_t {
	float pt;
	float eta;
	float phi;
	float eReco;
	float eTrue;
	int  matchIndex;
	float  eReco_over_eTrue;
	float mass;
	int clustersSize;
};

// .h class info
class BasicHggAnalyser : public edm::EDAnalyzer {
	public:
		explicit BasicHggAnalyser(const edm::ParameterSet&);
		~BasicHggAnalyser();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
		float resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc, const edm::PtrVector<reco::PFCluster>& clusters);
		float clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c, const edm::PtrVector<reco::PFCluster>& clusters);

	private:

		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		edm::EDGetTokenT<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >endcapRecHitCollection_      ; 
		edm::EDGetTokenT<edm::View<reco::SuperCluster> >endcapSuperClusterCollection_;
		edm::EDGetTokenT<edm::View<reco::PFCluster> >endcapClusterCollection_     ;
		edm::EDGetTokenT<edm::View<reco::GenParticle> >genParticlesCollection_     ;

		TTree *tree;
		TTree *treeTruth;
		info_t info;
		infoTruth_t infoTruth;
		TH1F *eta_h; 
		TH1F *phi_h; 
		TH1F *dEta_h; 
		TH2F *phiW_v_eRoT_h; 
		TH1F *dPhi_h;
		TH1F *nSC_h;
		TH1F *eRoT_OLD_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT_h;  //energy true over reco ie etrue/ ereco
		TH2F *dPhi_v_eRoT_h    ; 
		TH2F *nClust_v_eRoT_h  ; 
		TH2F *nClust90_v_eRoT_h;


};

// constructor
BasicHggAnalyser::BasicHggAnalyser(const edm::ParameterSet& iConfig):
	endcapRecHitCollection_(consumes <edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles"))))
{
	edm::Service<TFileService> fs_;
	eta_h         = fs_->make<TH1F>("eta_h","eta_h",100,-5,5);
	phi_h         = fs_->make<TH1F>("phi_h","phi_h",100,-5,5);
	dEta_h        = fs_->make<TH1F>("dEta_h","dEta_h",1000,-1,1);
	dPhi_h        = fs_->make<TH1F>("dPhi_h","dPhi_h",1000,3,3.3);
	eRoT_OLD_h        = fs_->make<TH1F>("eRoT_OLD_h","eRoT_OLD_h",1000,-2,2);
	eRoT_h        = fs_->make<TH1F>("eRoT_h","eRoT_h",1000,-2,2);
	phiW_v_eRoT_h        = fs_->make<TH2F>("phiW_v_eRoT_h","phiW_v_eRoT_h",100,0,0.1,100,0,1.3);
	dPhi_v_eRoT_h        = fs_->make<TH2F>("dPhi_v_eRoT_h","dPhi_v_eRoT_h",100,-0.1,0.1,100,0,1.3);
	nClust_v_eRoT_h        = fs_->make<TH2F>("nClust_v_eRoT_h","nCLust_v_eRoT_h",100,0,10,100,0,1.3);
	nClust90_v_eRoT_h        = fs_->make<TH2F>("nClust90_v_eRoT_h","nClust90_v_eRoT_h",100,0,10,100,0,1.3);
	nSC_h        = fs_->make<TH1F>("nSC","nSC",100,0,100);


	tree = fs_->make<TTree>("tree","");
	tree->Branch("eReco_over_eTrue"              ,&info.eReco_over_eTrue             ,"eReco_over_eTrue/F");
	tree->Branch("pt"              ,&info.pt             ,"pt/F");
	tree->Branch("eta"              ,&info.eta             ,"eta/F");
	tree->Branch("phi"              ,&info.phi             ,"phi/F");
	tree->Branch("eReco"              ,&info.eReco            ,"eReco/F");
	tree->Branch("eTrue"              ,&info.eTrue            ,"eTrue/F");
	tree->Branch("matchIndex"              ,&info.matchIndex            ,"matchIndex/I");
	tree->Branch("clustersSize"              ,&info.clustersSize            ,"clustersSize/I");

	treeTruth = fs_->make<TTree>("treeTruth","");
	treeTruth->Branch("pt"              ,&infoTruth.pt            ,"pt/F");
	treeTruth->Branch("ptReco"              ,&infoTruth.pt            ,"ptReco/F");
	treeTruth->Branch("eReco"              ,&infoTruth.eReco            ,"eReco/F");
	treeTruth->Branch("eTrue"              ,&infoTruth.eTrue            ,"eTrue/F");
	treeTruth->Branch("eClust1"              ,&infoTruth.eClust1            ,"eClust1/F");
	treeTruth->Branch("eClust2"              ,&infoTruth.eClust2            ,"eClust2/F");
	treeTruth->Branch("eClust3"              ,&infoTruth.eClust3            ,"eClust3/F");
	treeTruth->Branch("dEtaClust2"              ,&infoTruth.dEtaClust2            ,"dEtaClust2/F");
	treeTruth->Branch("dEtaClust3"              ,&infoTruth.dEtaClust3            ,"dEtaClust3/F");
	treeTruth->Branch("dPhiClust2"              ,&infoTruth.dPhiClust2            ,"dPhiClust2/F");
	treeTruth->Branch("dPhiClust3"              ,&infoTruth.dPhiClust3            ,"dPhiClust3/F");
	treeTruth->Branch("eSeed"              ,&infoTruth.eSeed            ,"eSeed/F");
	treeTruth->Branch("eSeed_over_eReco"   ,&infoTruth.eSeed_over_eReco             ,"eSeed_over_eReco/F");
	treeTruth->Branch("eReco_over_eTrue"              ,&infoTruth.eReco_over_eTrue             ,"eReco_over_eTrue/F");
	treeTruth->Branch("eta"              ,&infoTruth.eta             ,"eta/F");
	treeTruth->Branch("etaEcal"              ,&infoTruth.etaEcal             ,"etaEcal/F");
	treeTruth->Branch("etaSC"              ,&infoTruth.etaSC             ,"etaSC/F");
	treeTruth->Branch("etaSeed"              ,&infoTruth.etaSeed             ,"etaSeed/F");
	treeTruth->Branch("phi"              ,&infoTruth.phi             ,"phi/F");
	treeTruth->Branch("phiEcal"              ,&infoTruth.phiEcal             ,"phiEcal/F");
	treeTruth->Branch("phiSC"              ,&infoTruth.phiSC             ,"phiSC/F");
	treeTruth->Branch("phiSeed"              ,&infoTruth.phiSeed             ,"phiSeed/F");
	treeTruth->Branch("matchIndex"              ,&infoTruth.matchIndex            ,"matchIndex/I");
	treeTruth->Branch("clustersSize"              ,&infoTruth.clustersSize           ,"clustersSize/I");
	treeTruth->Branch("nClusters09"              ,&infoTruth.nClusters09           ,"nClusters09/I");
	treeTruth->Branch("etaWidth"              ,&infoTruth.etaWidth             ,"etaWidth/F");
	treeTruth->Branch("phiWidth"              ,&infoTruth.phiWidth             ,"phiWidth/F");
}

// destructor
BasicHggAnalyser::~BasicHggAnalyser()
{

}


//
// member functions
//

// ------------ method called for each event  ------------

float BasicHggAnalyser::resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){
	//float seedEta = sc->seed()->eta();
	float seedPhi = sc->seed()->phi();
	float total=0;
	for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
		//	std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;

	// if (seedEta - sc->clusters()[ic]->eta() >0.01) continue;
	 if (deltaPhi(seedPhi, sc->clusters()[ic]->phi()) >0.2) continue;

		for (unsigned int j =0 ; j < clusters.size() ; j++){

			if (clusters[j]->position()==(sc->clusters())[ic]->position()) {
				//		std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
				total = total +(clusters[j]->emEnergy());
				break;
			}
		}

	}

	return total;
}

float BasicHggAnalyser::clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c,const edm::PtrVector<reco::PFCluster>& clusters){

	float emEnergy=0;
	//for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
	//	std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
	for (unsigned int j =0 ; j < clusters.size() ; j++){

			//	std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
		if (clusters[j]->position()==(c->position())) {
			emEnergy=(clusters[j]->emEnergy());
			break;
		}
	}

	return emEnergy;
}

	void
BasicHggAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;


	Handle<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >  > HGCEERechits;
	iEvent.getByToken(endcapRecHitCollection_,HGCEERechits);
	//const PtrVector<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >& rechits = HGCEERechits->ptrVector();

	Handle<edm::View<reco::SuperCluster> > HGCEESCs;
	iEvent.getByToken(endcapSuperClusterCollection_,HGCEESCs);
	const PtrVector<reco::SuperCluster>& sclusters = HGCEESCs->ptrVector();

	Handle<edm::View<reco::PFCluster> > HGCEEClusters;
	iEvent.getByToken(endcapClusterCollection_,HGCEEClusters);
	const PtrVector<reco::PFCluster>& clusters = HGCEEClusters->ptrVector();

	Handle<edm::View<reco::GenParticle> > genParts;
	iEvent.getByToken(genParticlesCollection_,genParts);
	const PtrVector<reco::GenParticle>& gens = genParts->ptrVector();

	//	steco_over<< "[debug] number of rechits " << HGCEERechits->size() <<", SCs " << sclusters.size() << ", clusters " << clusters.size() << " gens " << gens.size() <<   std::endl;

	std::cout << "[DEBUG] 0" << std::endl;

	// initialise tree entries
	info.pt=-999.;
	info.eta=-999.;
	info.phi=-999.;
	info.eReco=-999.;
	info.eTrue=-999.;
	info.matchIndex=-999;
	infoTruth.eta=-999.;
	infoTruth.etaEcal=-999.;
	infoTruth.etaSC=-999.;
	infoTruth.etaSeed=-999.;
	infoTruth.phi=-999.;
	infoTruth.phiEcal=-999.;
	infoTruth.phiSC=-999.;
	infoTruth.phiSeed=-999.;
	infoTruth.etaWidth=-999.;
	infoTruth.phiWidth=-999.;
	infoTruth.matchIndex=-999;
	infoTruth.clustersSize=-999;
	infoTruth.nClusters09=-999;
	infoTruth.eReco  =-999.;         
	infoTruth.eTrue         =-999.;  
	infoTruth.eSeed           =-999.;
	infoTruth.eSeed_over_eReco=-999.;



	nSC_h->Fill(sclusters.size());

	for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree

		infoTruth.eta=-999.;
		infoTruth.etaEcal=-999.;
		infoTruth.etaSC=-999.;
		infoTruth.etaSeed=-999.;
		infoTruth.phi=-999.;
		infoTruth.phiSC=-999.;
		infoTruth.phiSeed=-999.;
		infoTruth.eReco_over_eTrue = -999.;
		infoTruth.matchIndex=-999;
		infoTruth.pt=-999.;
		infoTruth.ptReco=-999.;
		infoTruth.clustersSize=-999;
		infoTruth.nClusters09=-999;
		infoTruth.eReco  =-999.;         
		infoTruth.eTrue         =-999.;  
		infoTruth.eSeed           =-999.;
		infoTruth.eClust1           =-999.;
		infoTruth.eClust2           =-999.;
		infoTruth.eClust3           =-999.;
		infoTruth.dEtaClust2           =-999.;
		infoTruth.dEtaClust3           =-999.;
		infoTruth.dPhiClust2           =-999.;
		infoTruth.dPhiClust3           =-999.;
		infoTruth.eSeed_over_eReco=-999.;
		infoTruth.etaWidth=-999.;
		infoTruth.phiWidth=-999.;

		if (gens[igp]->pdgId() != 22 || gens[igp]->status() != 3) continue;
		assert(gens.size() >0); // only the case for the electron gun sample
		infoTruth.eta=gens[igp]->eta();
		infoTruth.phi      = gens[igp]->phi();
	 

		std::cout << "[debug] gen pdgid " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << ", e " << gens[igp]->energy()<<", et " <<gens[igp]->energy()/std::cosh(gens[igp]->eta()) << ", eta:phi" << gens[igp]->eta() << ":" << gens[igp]->phi() << std::endl;  // (gens[igp]->mother()->pdgId()) <<  std::endl;

		float dRBest =999;
		infoTruth.pt = gens[igp]->pt();
		infoTruth.eTrue = gens[igp]->energy();
		
		float zHGC =320.38;
		float rHGC = 32.0228;
		float RHGC = 152.153;
		if (infoTruth.eta <0) zHGC = -1*zHGC;
		math::XYZPoint vPos = gens[igp]->vertex();
		float theta = 2*std::atan(std::exp(-infoTruth.eta));
		float h = std::tan(theta)*(zHGC-vPos.z());
		if (h>RHGC || h<rHGC) continue;
		float eta_ECAL = -std::log(std::tan(std::asin(h/std::sqrt(h*h+zHGC*zHGC))/2)) ;
		if  (infoTruth.eta <0) eta_ECAL = -1*eta_ECAL;

	//	std::cout << "[photon] vtx z= " << vPos.z() << ", genEta = " << infoTruth.eta << ", ecalEta = " << eta_ECAL << std::endl;

		infoTruth.etaEcal =eta_ECAL;
		infoTruth.phiEcal = infoTruth.phi;

		for (unsigned int isc =0; isc < sclusters.size() ; isc++){ //subloop over sc's to find matches
			// calculate dR... dE = dEta, dP = dPhi
			float dE = fabs(sclusters[isc]->eta() - eta_ECAL);
			dE =dE*dE;
			float dP = deltaPhi(sclusters[isc]->phi(), gens[igp]->phi());
			dP =dP*dP;
			float dR = sqrt(dE +dP);
			std::cout << "[debug ] SC e " << (sclusters[isc])->energy() << ", et " <<  sclusters[isc]->energy()/std::cosh(eta_ECAL)<< ", eta:phi " << (sclusters[isc])->eta() << ":" <<(sclusters[isc])->phi() << ", dR " << dR << std::endl;

			if (dR < dRLim) { // only true if dR is both below limit value and smaller than previous best 
			//	std::cout << "check : ereco : " << resumEmEnergy(sclusters[isc], clusters) << ", etrue : " << infoTruth.eTrue << " dR = "<< dR << std::endl; 
				std::cout << "MATCH" << std::endl;
				// so, if we have a dR match, then store the corresponding match index.
			if (dR < dRBest){
				dRBest = dR;
				infoTruth.matchIndex = isc;
				}
			}

		}
		
		if (infoTruth.matchIndex <-1) { std::cout << "NO MATCH " << std::endl;};
		//float recoPt = -999.;
		if (infoTruth.matchIndex >-1) {
			//	recoPt = sclusters[infoTruth.matchIndex]->energy() /std::cosh(sclusters[infoTruth.matchIndex]->eta());
			//infoTruth.eReco = (sclusters[infoTruth.matchIndex]->energy());
			infoTruth.eReco = resumEmEnergy(sclusters[infoTruth.matchIndex], clusters);
			infoTruth.ptReco = resumEmEnergy(sclusters[infoTruth.matchIndex], clusters);
			infoTruth.eSeed = (clusterEmEnergy(sclusters[infoTruth.matchIndex]->seed(), clusters));
			if(sclusters[infoTruth.matchIndex]->clusters().size()>0){
			infoTruth.eClust1 = (clusterEmEnergy(sclusters[infoTruth.matchIndex]->clusters()[0], clusters));
	//		infotruth.dEtaClust1 = fabs(sclusters[infoTruth.matchIndex]->clusters()[0]->eta() - sclusters[infoTruth.matchIndex]->seed()->eta());
		//	infotruth.dPhiClust1 = fabs(sclusters[infoTruth.matchIndex]->clusters()[0]->phi() - sclusters[infoTruth.matchIndex]->seed()->phi());
			}
			if(sclusters[infoTruth.matchIndex]->clusters().size()>1){
			infoTruth.eClust2 = (clusterEmEnergy(sclusters[infoTruth.matchIndex]->clusters()[1], clusters));
			infoTruth.dEtaClust2 = fabs(sclusters[infoTruth.matchIndex]->clusters()[1]->eta() - sclusters[infoTruth.matchIndex]->seed()->eta());
			infoTruth.dPhiClust2 = fabs(sclusters[infoTruth.matchIndex]->clusters()[1]->phi() - sclusters[infoTruth.matchIndex]->seed()->phi());
			}
			if(sclusters[infoTruth.matchIndex]->clusters().size()>2){
			infoTruth.eClust3 = (clusterEmEnergy(sclusters[infoTruth.matchIndex]->clusters()[2], clusters));
			infoTruth.dEtaClust3 = fabs(sclusters[infoTruth.matchIndex]->clusters()[2]->eta() - sclusters[infoTruth.matchIndex]->seed()->eta());
			infoTruth.dPhiClust3 = fabs(sclusters[infoTruth.matchIndex]->clusters()[2]->phi() - sclusters[infoTruth.matchIndex]->seed()->phi());
			}
			
			infoTruth.eSeed_over_eReco = (infoTruth.eSeed/infoTruth.eReco);
			infoTruth.eReco_over_eTrue = (infoTruth.eReco/gens[igp]->energy());
			infoTruth.clustersSize = (int) sclusters[infoTruth.matchIndex]->clustersSize();
			infoTruth.etaWidth=sclusters[infoTruth.matchIndex]->etaWidth();
			infoTruth.phiWidth=sclusters[infoTruth.matchIndex]->phiWidth();
			infoTruth.etaSC    = sclusters[infoTruth.matchIndex]->eta();
			infoTruth.etaSeed  =sclusters[infoTruth.matchIndex]->seed()->eta();
			infoTruth.phiSC    = sclusters[infoTruth.matchIndex]->phi();
			infoTruth.phiSeed  =sclusters[infoTruth.matchIndex]->seed()->phi();


float etot =0;
			for (unsigned int ic =0 ; ic < sclusters[infoTruth.matchIndex]->clusters().size() ; ic++){
				infoTruth.nClusters09 = ic+1;
				etot +=(sclusters[infoTruth.matchIndex]->clusters())[ic]->energy();
			std::cout << "cluster "<< ic << ", energy " << (sclusters[infoTruth.matchIndex]->clusters())[ic]->energy() << ", percent "<< etot /infoTruth.eReco << std::endl;
				if((etot/infoTruth.eReco >0.9)) break;
			}
		}

		//if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  recoPt >40)
		if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  gens[igp]->pt() >40)
		{
			eRoT_OLD_h->Fill(sclusters[infoTruth.matchIndex]->energy()/gens[igp]->energy());
			eRoT_h->Fill(infoTruth.eReco/gens[igp]->energy());
			phiW_v_eRoT_h->Fill(infoTruth.phiWidth,infoTruth.eReco_over_eTrue);
			dPhi_v_eRoT_h ->Fill(infoTruth.phiSeed - infoTruth.phi,infoTruth.eReco_over_eTrue);
			nClust_v_eRoT_h ->Fill(infoTruth.clustersSize,infoTruth.eReco_over_eTrue);
			nClust90_v_eRoT_h->Fill(infoTruth.nClusters09,infoTruth.eReco_over_eTrue);


			if(infoTruth.eReco_over_eTrue <0.5){
				std::cout << "MATCH with eReco/eTrue <0.5" << std::endl;
				std::cout << "_____|  RECO    |  TRUE    |" <<  std::endl;
				std::cout << " eta | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->eta() <<" | "<< gens[igp]->eta()<< std::endl;
				std::cout << " phi | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->phi() <<" | "<< gens[igp]->phi()<< std::endl;
				std::cout << " e   | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->energy() <<" | "<< gens[igp]->energy()<< std::endl;
				std::cout << " e   | " << std::setprecision(5) << infoTruth.eReco <<" | "<< infoTruth.eTrue<< std::endl;
			}
		}

		treeTruth->Fill();
	}

	//--------------> End per-genPhoton tree <------------------


	//--------------> Begin per-SC tree <---------------------- 

	// loop over superclusters (eg reco particles). 
	for( unsigned int isc =0; isc< sclusters.size() ; isc++){ // isc = index_super_cluster

		info.eReco_over_eTrue=-999.;
		info.pt=-999.;
		info.eta=-999.;
		info.phi=-999.;
		info.eReco=-999.;
		info.eTrue=-999.;
		info.matchIndex=-999;
		info.clustersSize=-999;

		//std::cout << " sc " << isc<< " eta " << sclusters[isc]->eta() << ", phi " << sclusters[isc]->phi()<< std::endl;

		info.pt= sclusters[isc]->energy() /std::cosh(sclusters[isc]->eta());
		info.eta=sclusters[isc]->eta();
		info.phi=sclusters[isc]->phi();
		info.eReco = sclusters[isc]->energy();
		info.clustersSize = (int) sclusters[isc]->clustersSize();;

		// fill histograms with eta/phi info
		eta_h->Fill(info.eta);
		phi_h->Fill(info.phi);

		float dRBest = 999.; // dR best is used to find the gen-reco match with smallest dR.

		// now loop over gen particles
		for (unsigned int igp =0; igp < gens.size() ; igp++){ // igp = index_gen_particles

			float dE = sclusters[isc]->eta() - gens[igp]->eta();
			dE =dE*dE;
			float dP = sclusters[isc]->phi() - gens[igp]->phi();
			dP =dP*dP;
			float dR = sqrt(dE +dP);

			if (dR < dRLim && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
				dRBest = dR;

				// so, if we have a dR match, then store the corresponding match index.
				info.matchIndex = igp;
				info.eTrue = gens[igp]->energy();

			}
		}

		if(info.matchIndex >-1){	
			info.eReco_over_eTrue = info.eReco/info.eTrue; 
			//	eRoT_h->Fill(info.eReco/info.eTrue);
		}
		tree->Fill();

	}
	/*
		 if (eTrue[0] <0) {
		 eTrue[0] = gens[igp]->energy();
		 } else {
		 eTrue[1] = eTrue[0];
		 eTrue[0] = gens[igp]->energy();
		 }*/
	/*		}
				}
				std::cout << "[DEBUG] 3" << std::endl;*/
/*	if(sclusters.size()){
//	tree->Fill();
//	dEta_h->Fill(info.eta[0]+info.eta[1]);
//	dPhi_h->Fill(fabs(info.phi[0]-info.phi[1]));

if( info.matchIndex[0] > -1){
//	eRoT_h->Fill(info.eReco[0]/info.eTrue[info.matchIndex[0]]);
std::cout << "RECO 0 matched TRUE " << info.matchIndex[0] << std::endl;
std::cout << "eta_reco " << info.eta[0] << ", eta_true" << gens[info.matchIndex[0]]->eta() << std::endl;
std::cout << "phi_reco " << info.phi[0] << ", phi_true" << gens[info.matchIndex[0]]->phi() << std::endl;
std::cout << "e_reco " << info.eReco[0] << ", e_true" << gens[info.matchIndex[0]]->energy() << " (" << info.eTrue[info.matchIndex[0]] << ")" <<  std::endl;
}
if( info.matchIndex[1] > -1){
//	eRoT_h->Fill(info.eReco[1]/info.eTrue[info.matchIndex[1]]);
std::cout << "RECO 1 matched TRUE " << info.matchIndex[1] << std::endl;
std::cout << "eta_reco " << info.eta[1] << ", eta_true" << gens[info.matchIndex[1]]->eta() << std::endl;
std::cout << "phi_reco " << info.phi[1] << ", phi_true" << gens[info.matchIndex[1]]->phi() << std::endl;
std::cout << "e_reco " << info.eReco[1] << ", e_true" << gens[info.matchIndex[1]]->energy() << " (" << info.eTrue[info.matchIndex[1]] << ")" <<  std::endl;
}



}


*/
return ;
}




// ------------ method called once each job just before starting event loop  ------------
	void 
BasicHggAnalyser::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
BasicHggAnalyser::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
BasicHggAnalyser::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
BasicHggAnalyser::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
BasicHggAnalyser::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
BasicHggAnalyser::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BasicHggAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicHggAnalyser);
