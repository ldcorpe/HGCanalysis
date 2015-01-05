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
BasicHggAnalyser::BasicHggAnalyser(const edm::ParameterSet& iConfig)
{
  // actually, pset is NOT in use - we keep it here just for illustratory putposes
  HepMC::IO_GenEvent ascii_out("IO_GenEvent.dat",std::ios::out);
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

   // here's an example of accessing GenEventInfoProduct
   /*
     Handle< GenEventInfoProduct > GenInfoHandle;
     iEvent.getByLabel( "generator", GenInfoHandle );
     double qScale = GenInfoHandle->qScale();
     double pthat = ( GenInfoHandle->hasBinningValues() ? 
     (GenInfoHandle->binningValues())[0] : 0.0);
     cout << " qScale = " << qScale << " pthat = " << pthat << endl;
   */

   //Handle< HepMCProduct > EvtHandle ;
   //iEvent.getByLabel( "generator", EvtHandle ) ;
  
   //const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;

   //write GenEvent into HepMC file
   //ascii_out << Evt;
   std::cout << "test" << std::endl;
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
