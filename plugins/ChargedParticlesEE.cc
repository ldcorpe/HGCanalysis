// -*- C++ -*-
//
// Package:    ChargedParticlesEE
// Class:      ChargedParticlesEE
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
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
using namespace std;

//
// class declaration
//
// limit dR allowed for geometrical matching of reco/gen particles
//float dRLim =0.05;

// information to be loaded into TTree
struct infoTruth_t {

	float  eReco_over_eTrue;
int nCharged2 ;
int nCharged2_5 ;
int nCharged3 ;
float hPt;
};


// .h class info
class ChargedParticlesEE : public edm::EDAnalyzer {
	public:
		explicit ChargedParticlesEE(const edm::ParameterSet&);
		~ChargedParticlesEE();

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
	//	info_t info;
		infoTruth_t infoTruth;

		TProfile *n2_p;
		TProfile *n2_5_p;
		TProfile *n3_p;
};

// constructor
ChargedParticlesEE::ChargedParticlesEE(const edm::ParameterSet& iConfig):
//	endcapRecHitCollection_(consumes <edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
//	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
//	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles"))))
{
	edm::Service<TFileService> fs_;
	
	n2_p = fs_->make<TProfile>("n2_p","n2_p",100,0,1000,0,20);
	n2_5_p = fs_->make<TProfile>("n2_5_p","n2_5_p",100,0,1000,0,20);
	n3_p = fs_->make<TProfile>("n3_p","n3_p",100,0,1000,0,20);


	treeTruth = fs_->make<TTree>("tree","");
	treeTruth->Branch("nCharged2"             ,&infoTruth.nCharged2           ,"nCharged2/I");
	treeTruth->Branch("nCharged2_5"           ,&infoTruth.nCharged2_5           ,"nCharged2_5/I");
	treeTruth->Branch("nCharged3"             ,&infoTruth.nCharged3           ,"nCharged3/I");
	treeTruth->Branch("hPt"             ,&infoTruth.hPt          ,"hPt/F");

}

// destructor
ChargedParticlesEE::~ChargedParticlesEE()
{

}


//
// member functions
//

// ------------ method called for each event  ------------

/*float ChargedParticlesEE::resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){

	float total=0;
	for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
		//	std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
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

float ChargedParticlesEE::clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c,const edm::PtrVector<reco::PFCluster>& clusters){

	float emEnergy=0;
	//for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
	//	std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
	for (unsigned int j =0 ; j < clusters.size() ; j++){

		if (clusters[j]->position()==(c->position())) {
			//		std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
			emEnergy=(clusters[j]->emEnergy());
			break;
		}
	}

	return emEnergy;
}*/

	void
ChargedParticlesEE::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

/*
	Handle<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >  > HGCEERechits;
	iEvent.getByToken(endcapRecHitCollection_,HGCEERechits);
	//const PtrVector<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >& rechits = HGCEERechits->ptrVector();

	Handle<edm::View<reco::SuperCluster> > HGCEESCs;
	iEvent.getByToken(endcapSuperClusterCollection_,HGCEESCs);
	const PtrVector<reco::SuperCluster>& sclusters = HGCEESCs->ptrVector();

	Handle<edm::View<reco::PFCluster> > HGCEEClusters;
	iEvent.getByToken(endcapClusterCollection_,HGCEEClusters);
	const PtrVector<reco::PFCluster>& clusters = HGCEEClusters->ptrVector();*/

	Handle<edm::View<reco::GenParticle> > genParts;
	iEvent.getByToken(genParticlesCollection_,genParts);
	const PtrVector<reco::GenParticle>& gens = genParts->ptrVector();

	//	steco_over<< "[debug] number of rechits " << HGCEERechits->size() <<", SCs " << sclusters.size() << ", clusters " << clusters.size() << " gens " << gens.size() <<   std::endl;

	std::cout << "[DEBUG] 0" << std::endl;

	// initialise tree entries
	infoTruth.nCharged2  = 0;       
	infoTruth.nCharged2_5 =0; 
	infoTruth.nCharged3 = 0;
	infoTruth.hPt = -999.;



	//nSC_h->Fill(sclusters.size());

	for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree


		if (gens[igp]->pdgId() == 25 ) {infoTruth.hPt = gens[igp]->pt();}
		
	if(gens[igp]->charge() != 0 && fabs(gens[igp]->eta())>1.5 && fabs(gens[igp]->eta())<3. && gens[igp]->pt() >2.)	{
	infoTruth.nCharged2++;
	}
	if(gens[igp]->charge() != 0 && fabs(gens[igp]->eta())>1.5 && fabs(gens[igp]->eta())<3. && gens[igp]->pt() >2.5)	{
	infoTruth.nCharged2_5++;	
	}
	if(gens[igp]->charge() != 0 && fabs(gens[igp]->eta())>1.5 && fabs(gens[igp]->eta())<3. && gens[igp]->pt() >3.)	{
	infoTruth.nCharged3++;	
	}

	}

	n2_p->Fill(infoTruth.hPt,infoTruth.nCharged2);
	n2_5_p->Fill(infoTruth.hPt,infoTruth.nCharged2_5);
	n3_p->Fill(infoTruth.hPt,infoTruth.nCharged3);


	treeTruth->Fill();
	//--------------> End per-genPhoton tree <------------------


	//--------------> Begin per-SC tree <---------------------- 

	// loop over superclusters (eg reco particles). 
/*	for( unsigned int isc =0; isc< sclusters.size() ; isc++){ // isc = index_super_cluster

		info.eReco_over_eTrue=-999.;
		info.pt=-999.;
		info.eta=-999.;
		info.phi=-999.;
		info.eReco=-999.;
		info.eTrue=-999.;
		info.matchIndex=-999;
		info.clustersSize=-999;

		std::cout << " sc " << isc<< " eta " << sclusters[isc]->eta() << ", phi " << sclusters[isc]->phi()<< std::endl;

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
	
	//	 if (eTrue[0] <0) {
	//	 eTrue[0] = gens[igp]->energy();
	//	 } else {
	//	 eTrue[1] = eTrue[0];
	//	 eTrue[0] = gens[igp]->energy();
	//	 }
			}
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
ChargedParticlesEE::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
ChargedParticlesEE::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
ChargedParticlesEE::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
ChargedParticlesEE::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
ChargedParticlesEE::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
ChargedParticlesEE::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ChargedParticlesEE::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ChargedParticlesEE);
