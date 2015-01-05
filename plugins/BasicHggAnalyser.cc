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


#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/Utilities/interface/InputTag.h"

using namespace std;

//
// class declaration
//

class BasicHggAnalyser : public edm::EDAnalyzer {
public:
  explicit BasicHggAnalyser(const edm::ParameterSet&);
  ~BasicHggAnalyser();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
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
  
  // ----------member data ---------------------------
  
  //HepMC::IO_GenEvent ascii_out;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BasicHggAnalyser::BasicHggAnalyser(const edm::ParameterSet& iConfig):
endcapRecHitCollection_(consumes <edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE"))))
{
  // actually, pset is NOT in use - we keep it here just for illustratory putposes
 //	 HepMC::IO_GenEvent ascii_out("IO_GenEvent.dat",std::ios::out);
}


BasicHggAnalyser::~BasicHggAnalyser()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

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
	
	std::cout << "debug: rechits " << HGCEERechits->size() <<", SCs " << sclusters.size() << ", clusters " << clusters.size() <<   std::endl;
	return ;
}


// ------------ method called once each job just before starting event loop  ------------
	void 
BasicHggAnalyser::beginJob()
{
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
