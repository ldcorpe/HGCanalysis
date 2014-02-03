#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"


#include <iostream>

using namespace std;

//
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig ) : geometryDefined_(false)
{
  ddViewName_  = iConfig.getUntrackedParameter<std::string>("ddViewName", "");
  eeHits_      = iConfig.getUntrackedParameter<std::string>("eeHits",     "HGCHitsEE");
  heHits_      = iConfig.getUntrackedParameter<std::string>("heHits",     "HGCHitsHE");
  genSource_   = iConfig.getUntrackedParameter<std::string>("genSource",  "genParticles");

  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("run",      &simEvt_.run,      "run/I");
  t_->Branch("lumi",     &simEvt_.lumi,     "lumi/I");
  t_->Branch("event",    &simEvt_.event,    "event/I");
  t_->Branch("nee",      &simEvt_.nee,      "nee/I");
  t_->Branch("ee_layer",  simEvt_.ee_layer, "ee_layer[nee]/I");
  t_->Branch("ee_edep",   simEvt_.ee_edep,  "ee_edep[nee]/F");
  t_->Branch("ee_t",      simEvt_.ee_t,     "ee_t[nee]/F");
}

//
HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer()
{
}

//
void HGCSimHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  cout << "[HGCSimHitsAnalyzer::analyze] start" << endl;

  if(!geometryDefined_){
    //get handle to the DD
    edm::ESTransientHandle<DDCompactView> ddViewH;
    iSetup.get<IdealGeometryRecord>().get( ddViewName_, ddViewH );
    
    //safety check
    if (!ddViewH.isValid() || !ddViewH.description()) {
      cout << "[HGCSimHitsAnalyzer::analyze] Handle for DD is not valid or does not contain any valid description" << endl;
      return;
    }    

    geometryDefined_=defineGeometry(ddViewH);
    if(geometryDefined_) { std::cout << "[HGCSimHitsAnalyzer::analyze] geometry has been defined" << endl;                      }
    else                 { std::cout << "[HGCSimHitsAnalyzer::analyze] could not define geometry from DD view" << endl; return; }
  }

  //get hits and gen level information
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);

  edm::Handle<edm::PCaloHitContainer> caloHitsEE;
  iEvent.getByLabel(edm::InputTag("g4SimHits",eeHits_),caloHitsEE); 

  edm::Handle<edm::PCaloHitContainer> caloHitsHE;
  iEvent.getByLabel(edm::InputTag("g4SimHits",heHits_),caloHitsHE); 

  simEvt_.run    = iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();
  simEvt_.event  = iEvent.id().event();

  simEvt_.nee=0;
  analyzeEEHits(caloHitsEE,genParticles);
  
  if(simEvt_.nee>0)  t_->Fill();
}

//
bool HGCSimHitsAnalyzer::defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH)
{
   //const DDCompactView &pToDDView=*ddViewH;
  return true;
}

//
void HGCSimHitsAnalyzer::analyzeEEHits(edm::Handle<edm::PCaloHitContainer> &caloHits, edm::Handle<edm::View<reco::Candidate> > &gen)
{
  if(!caloHits.isValid()) return;
  for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it) 
    {
      HGCEEDetId detId(hit_it->id());

      simEvt_.ee_layer[simEvt_.nee] = detId.layer();
      simEvt_.ee_edep[simEvt_.nee]  = hit_it->energy();
      simEvt_.ee_t[simEvt_.nee]     = hit_it->time();
      simEvt_.nee++;

//       std::cout << hex << "0x" << uint32_t(detId) << dec 
// 		<< " | isFwd=" << detId.isForward()
// 		<< " isEE=" << detId.isEE() 
// 		<< " zside=" << detId.zside() 
// 		<< " layer=" << detId.layer() 
// 		<< " subsec=" << detId.subsector() 
// 		<< " module=" << detId.module()
// 		<< " cell=" << detId.cell() 
// 		<< " | den=" << den
// 		<< std::endl;
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);