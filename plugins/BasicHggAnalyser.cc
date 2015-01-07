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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
using namespace std;

//
// class declaration
//
// limit dR allowed for geometrical matching of reco/gen particles
float dRLim =0.1;

// information to be loaded into TTree
struct infoTruth_t {

float eta[2];
int matchIndex[2];

};

struct info_t {
	float pt;
	float eta[2];
	float phi[2];
	float eReco[2];
	float eTrue[2];
	int  matchIndex[2];
	float mass;
};

// .h class info
class BasicHggAnalyser : public edm::EDAnalyzer {
	public:
		explicit BasicHggAnalyser(const edm::ParameterSet&);
		~BasicHggAnalyser();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		edm::Service<TFileService> fs_;

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
		TH1F *dPhi_h; 
		TH1F *eToR_h;  //energy true over reco ie etrue/ ereco

};

// constructor
BasicHggAnalyser::BasicHggAnalyser(const edm::ParameterSet& iConfig):
	endcapRecHitCollection_(consumes <edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles"))))
{
}

 // destructor
BasicHggAnalyser::~BasicHggAnalyser()
{


}


//
// member functions
//

// ------------ method called for each event  ------------
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

	std::cout << "[debug] number of rechits " << HGCEERechits->size() <<", SCs " << sclusters.size() << ", clusters " << clusters.size() <<   std::endl;
	
	// initialise tree entries
	info.eta[0]=-999.;
	info.eta[1]=-999.;
	info.phi[0]=-999.;
	info.phi[1]=-999.;
	info.eReco[0]=-999.;
	info.eReco[1]=-999.;
	info.eTrue[0]=-999.;
	info.eTrue[1]=-999.;
	info.matchIndex[0]=-999; //matchIndex stores the relevant index of the matched particle, and is -1 if no match is found. 
	info.matchIndex[1]=-999;

	infoTruth.eta[0]=-999.;
	infoTruth.eta[1]=-999.;
	infoTruth.matchIndex[0]=-999;
	infoTruth.matchIndex[1]=-999;


	for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree

		assert(gens.size() ==2);
		infoTruth.eta[igp]=gens[igp]->eta();

		float dRBest =999;

		for (unsigned int isc =0; isc < sclusters.size() ; isc++){ //suloop over sc's to find matches

			// calculate dR... dE = dEta, dP = dPhi
			float dE = sclusters[isc]->eta() - gens[igp]->eta();
			dE =dE*dE;
			float dP = sclusters[isc]->phi() - gens[igp]->phi();
			dP =dP*dP;
			float dR = sqrt(dE +dP);

			if (dR < dRLim && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
				dRBest = dR;

				// so, if we have a dR match, then store the corresponding match index.
				infoTruth.matchIndex[igp] = isc;
			}

		}

	}

	treeTruth->Fill();

	// buffers for eReco and eTrue so they can be stored in the right order
	//float eReco[2]={-999.,-999.};
	//float eTrue[2]={-999.,-999.};

	// loop over superclusters (eg reco particles). 
	for( unsigned int isc =0; isc< sclusters.size() ; isc++){ // isc = index_super_cluster
		std::cout << " sc " << isc<< " eta " << sclusters[isc]->eta() << ", phi " << sclusters[isc]->phi()<< std::endl;
		assert(isc<3);
		// fill tree info with  eta/phi/ energy info
		info.eta[isc]=sclusters[isc]->eta();
		info.phi[isc]=sclusters[isc]->phi();
		info.eReco[isc] = sclusters[isc]->energy();

		// in principle it would be nice to store the leading photon energy first...
		/*
			 if (eReco[0] <0) {
			 eReco[0] = sclusters[isc]->energy();
			 } else {
			 eReco[1] = eReco[0];
			 eReco[0] = sclusters[isc]->energy();
			 }*/

		// fill histograms with eta/phi info
		eta_h->Fill(info.eta[isc]);
		phi_h->Fill(info.phi[isc]);

		float dRBest = 999.; // dR best is used to find the gen-reco match with smallest dR.

		// now loop over gen particles
		for (unsigned int igp =0; igp < gens.size() ; igp++){ // igp = index_gen_particles

			info.eTrue[igp] = gens[igp]->energy();

			// calculate dR... dE = dEta, dP = dPhi
			float dE = sclusters[isc]->eta() - gens[igp]->eta();
			dE =dE*dE;
			float dP = sclusters[isc]->phi() - gens[igp]->phi();
			dP =dP*dP;
			float dR = sqrt(dE +dP);

			if (dR < dRLim && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
				dRBest = dR;

				// so, if we have a dR match, then store the corresponding match index.
				info.matchIndex[isc] = igp;
			}

			//std::cout << "[debug] pdg id " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << std::endl;


			/*
				 if (eTrue[0] <0) {
				 eTrue[0] = gens[igp]->energy();
				 } else {
				 eTrue[1] = eTrue[0];
				 eTrue[0] = gens[igp]->energy();
				 }*/
		}
	}
	if(sclusters.size()){
		tree->Fill();
		dEta_h->Fill(info.eta[0]+info.eta[1]);
		dPhi_h->Fill(fabs(info.phi[0]-info.phi[1]));

		if( info.matchIndex[0] > -1){
			eToR_h->Fill(info.eReco[0]/info.eTrue[info.matchIndex[0]]);
			std::cout << "RECO 0 matched TRUE " << info.matchIndex[0] << std::endl;
			std::cout << "eta_reco " << info.eta[0] << ", eta_true" << gens[info.matchIndex[0]]->eta() << std::endl;
			std::cout << "phi_reco " << info.phi[0] << ", phi_true" << gens[info.matchIndex[0]]->phi() << std::endl;
			std::cout << "e_reco " << info.eReco[0] << ", e_true" << gens[info.matchIndex[0]]->energy() << " (" << info.eTrue[info.matchIndex[0]] << ")" <<  std::endl;
		}
		if( info.matchIndex[1] > -1){
			eToR_h->Fill(info.eReco[1]/info.eTrue[info.matchIndex[1]]);
			std::cout << "RECO 1 matched TRUE " << info.matchIndex[1] << std::endl;
			std::cout << "eta_reco " << info.eta[1] << ", eta_true" << gens[info.matchIndex[1]]->eta() << std::endl;
			std::cout << "phi_reco " << info.phi[1] << ", phi_true" << gens[info.matchIndex[1]]->phi() << std::endl;
			std::cout << "e_reco " << info.eReco[1] << ", e_true" << gens[info.matchIndex[1]]->energy() << " (" << info.eTrue[info.matchIndex[1]] << ")" <<  std::endl;
		}



	}



	return ;
}


// ------------ method called once each job just before starting event loop  ------------
	void 
BasicHggAnalyser::beginJob()
{

	eta_h         = fs_->make<TH1F>("eta_h","eta_h",100,-5,5);
	phi_h         = fs_->make<TH1F>("phi_h","phi_h",100,-5,5);
	dEta_h        = fs_->make<TH1F>("dEta_h","dEta_h",1000,-1,1);
	dPhi_h        = fs_->make<TH1F>("dPhi_h","dPhi_h",1000,3,3.3);
	eToR_h        = fs_->make<TH1F>("eToR_h","eToR_h",1000,-2,2);


	tree = fs_->make<TTree>("tree","");
	tree->Branch("eta"              ,info.eta             ,"eta[2]/F");
	tree->Branch("phi"              ,info.phi             ,"phi[2]/F");
	tree->Branch("eReco"              ,info.eReco            ,"eReco[2]/F");
	tree->Branch("eTrue"              ,info.eTrue            ,"eTrue[2]/F");
	tree->Branch("matchIndex"              ,info.matchIndex            ,"matchIndex[2]/I");
	
	treeTruth = fs_->make<TTree>("treeTruth","");
	treeTruth->Branch("eta"              ,infoTruth.eta             ,"eta[2]/F");
	treeTruth->Branch("matchIndex"              ,infoTruth.matchIndex            ,"matchIndex[2]/I");
	return ;
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
