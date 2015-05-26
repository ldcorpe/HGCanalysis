// -*- C++ -*-
//
// Package:    HoverEAnalyzer_MVA
// Class:      HoverEAnalyzer_MVA
// 
//
//

class ElectronSeedGenerator ;
class SeedFilter ;
class EgammaHcalIsolation ;
class ElectronHcalHelper ;
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
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "RecoCaloTools/Selectors/interface/CaloDualConeSelector.h"
#include "RecoEcal/EgammaClusterAlgos/interface/HGCALShowerBasedEmIdentificationLC2.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "PCAShowerAnalysis.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TMVA/Reader.h"

using namespace std;

int eventIndex=0;
//
// class declaration
//

// information to be loaded into TTree
struct infoPhoton_t {

	float eta;
	float etaSC;
	float etaEcal;
	float x;
	float y;
	float phi;
	float phiSC;
	float phiEcal;
	int matchIndex;
	float pt;
	float recopt;
	float et;
	float dRBest;
	float eTrue;
	float eReco;
	float hoe;
	float sigmaEtaEta;
	float sigmaRR;
	float lengthCompatibility;
  int converted;
	int iProcess;
	float phiCrackDistance;
	int invalidNeighbour;
	int eventIndex;
	float MVA;

};

struct infoBackground_t {
	float eta;
	float etaSC;
	float etaEcal;
	float x;
	float y;
	float phi;
	float phiSC;
	float phiEcal;
	int matchIndex;
	float pt;
	float recopt;
	float et;
	float eTrue;
	float eReco;
	float hoe;
	float dRPhoton;
	float sigmaEtaEta;
	float sigmaRR;
	float lengthCompatibility;
	int iProcess;
	float phiCrackDistance;
	int invalidNeighbour;
	int eventIndex;
	float MVA;
	float dRBest;
};

// .h class info
class HoverEAnalyzer_MVA : public edm::EDAnalyzer {
	public:
		explicit HoverEAnalyzer_MVA(const edm::ParameterSet&);
		~HoverEAnalyzer_MVA();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		double fillEmIdVars(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters,const  edm::PtrVector<HGCRecHit>& rcs, const HGCalGeometry *geom , const edm::Event& iEvent, const edm::EventSetup& iSetup);
float annemarieEnergy( const edm::PtrVector<HGCRecHit>& rechitvec, const edm::Ptr<reco::SuperCluster>& sc, 	std::vector<std::string> geometrySource_, const edm::Event& iEvent, const edm::EventSetup& iSetup, float eta_ECAL, int sig=1);
		float phiCrackDistance (float x, float y);
		void getMaximumCell(const edm::PtrVector<HGCRecHit>& rechitvec,const double & phimax,const    double & etamax,std::vector<HGCEEDetId> & detidmax);
		double absWeight(const unsigned layer, const bool dedx=false);

		void fillDetIdStack(std::vector<HGCEEDetId> & detidmax_vec, std::vector<HGCEEDetId> & detidstack);

	private:

		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

 
	edm::FileInPath MVAweightfile_;
	 unique_ptr<TMVA::Reader>Mva_;

		ElectronHcalHelper * hcalHelperEndcap_ ;
		edm::ESHandle<CaloGeometry> caloGeom_ ;
		unsigned long long caloGeomCacheId_ ;
		edm::ESHandle<CaloTopology> caloTopo_;
		unsigned long long caloTopoCacheId_;
		edm::EDGetTokenT<edm::View<HGCRecHit> >endcapRecHitCollection_      ; 
		edm::EDGetTokenT<edm::View<reco::SuperCluster> >endcapSuperClusterCollection_;
		edm::EDGetTokenT<edm::View<reco::PFCluster> >endcapClusterCollection_     ;
		edm::EDGetTokenT<edm::View<reco::GenParticle> >genParticlesCollection_     ;
		edm::EDGetTokenT<edm::View<reco::PFRecHit> > eeRecHitCollection_     ;

  std::string g4TracksSource_, g4VerticesSource_; 
		std::vector<std::string> geometrySource_;
  	const HGCalGeometry * hgcEEGeom_;
    
		ROOT::Math::XYZPoint truthVtx_;
		edm::Handle<HGCRecHitCollection> recHits_;	
		int PU;
		int process;
		int genMatchPU;
  double cellSize_;
  bool doLogWeight_;
  double mipE_;
  unsigned nSR_;
  unsigned debug_;
  bool singleGamma_;


		TTree *treePhoton;
		TTree *treeBackground;
		infoBackground_t infoBackground;
		infoPhoton_t infoPhoton;

		float sigmaEtaEta_;
		float hoe_;
		float lengthCompatibility_;

};

// constructor
HoverEAnalyzer_MVA::HoverEAnalyzer_MVA(const edm::ParameterSet& iConfig):
	hcalHelperEndcap_(0), caloGeom_(0), caloGeomCacheId_(0), caloTopo_(0), caloTopoCacheId_(0),
	endcapRecHitCollection_(consumes <edm::View<HGCRecHit> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles")))),
	eeRecHitCollection_(consumes <edm::View<reco::PFRecHit> > (iConfig.getUntrackedParameter<edm::InputTag>("eeRecHitCollection",edm::InputTag("particleFlowRecHitHGCEE:Cleaned"))))

{ 
  MVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "MVAweightfile" );

	sigmaEtaEta_=0;
	hoe_ =0;
	lengthCompatibility_=0;

	Mva_.reset( new TMVA::Reader( "!Color:Silent" ) );
	Mva_->AddVariable( "sigmaEtaEta", &sigmaEtaEta_ );
	Mva_->AddVariable( "hoe", &hoe_ );
	Mva_->AddVariable( "lengthCompatibility", &lengthCompatibility_ );

	Mva_->BookMVA( "BDTG",MVAweightfile_.fullPath() );
	g4TracksSource_           = iConfig.getUntrackedParameter<std::string>("g4TracksSource");
	g4VerticesSource_         = iConfig.getUntrackedParameter<std::string>("g4VerticesSource");

	geometrySource_ = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
	ElectronHcalHelper::Configuration hcalCfgEndcap ;
	hcalCfgEndcap.hOverEConeSize = iConfig.getParameter<double>("hOverEConeSize") ;
	hcalCfgEndcap.hOverEMethod = iConfig.getParameter<int>("hOverEMethodEndcap") ;
	PU  = iConfig.getParameter<int>("PU") ;
	process  = iConfig.getParameter<int>("process_") ;
	genMatchPU  = iConfig.getParameter<int>("genMatchPU_") ;
	if (hcalCfgEndcap.hOverEConeSize>0)
	{
		hcalCfgEndcap.useTowers = true ;
		hcalCfgEndcap.hcalTowers = iConfig.getParameter<edm::InputTag>("hcalTowers") ;
		//here the HCAL clusters
		if (hcalCfgEndcap.hOverEMethod==3)
		{ 
			hcalCfgEndcap.hcalClusters = iConfig.getParameter<edm::InputTag>("endcapHCALClusters") ;
		}
	}
	hcalCfgEndcap.hOverEPtMin = iConfig.getParameter<double>("hOverEPtMin") ;
	hcalHelperEndcap_ = new ElectronHcalHelper(hcalCfgEndcap) ;


	edm::Service<TFileService> fs_;
	//	hoe_Sig_h         = fs_->make<TH1F>("hoe_Sig_h","hoe_Sig__h",1000,0,1);
	//	hoe_PU_h         = fs_->make<TH1F>("hoe_PU_h","hoe_PU__h",1000,0,1);


	treePhoton = fs_->make<TTree>("treePhoton","");
	treePhoton->Branch("pt"                 ,&infoPhoton.pt                 ,"pt/F");
	treePhoton->Branch("hoe"                ,&infoPhoton.hoe                ,"hoe/F");
	treePhoton->Branch("eReco"              ,&infoPhoton.eReco              ,"eReco/F");
	treePhoton->Branch("eTrue"              ,&infoPhoton.eTrue              ,"eTrue/F");
	treePhoton->Branch("eta"                ,&infoPhoton.eta                ,"eta/F");
	treePhoton->Branch("etaSC"              ,&infoPhoton.etaSC              ,"etaSC/F");
	treePhoton->Branch("etaEcal"            ,&infoPhoton.etaEcal            ,"etaEcal/F");
	treePhoton->Branch("phi"                ,&infoPhoton.phi                ,"phi/F");
	treePhoton->Branch("x"                  ,&infoPhoton.x                  ,"x/F");
	treePhoton->Branch("y"                  ,&infoPhoton.y                  ,"y/F");
	treePhoton->Branch("phiSC"              ,&infoPhoton.phiSC              ,"phiSC/F");
	treePhoton->Branch("phiEcal"            ,&infoPhoton.phiEcal            ,"phiEcal/F");
	treePhoton->Branch("matchIndex"         ,&infoPhoton.matchIndex         ,"matchIndex/I");
	treePhoton->Branch("sigmaEtaEta"        ,&infoPhoton.sigmaEtaEta        ,"sigmaEtaEta/F");
	treePhoton->Branch("sigmaRR"            ,&infoPhoton.sigmaRR            ,"sigmaRR/F");
	treePhoton->Branch("dRBest"            ,&infoPhoton.dRBest            ,"dRBest/F");
	treePhoton->Branch("lengthCompatibility",&infoPhoton.lengthCompatibility,"lengthCompatibility/F");
	treePhoton->Branch("phiCrackDistance",&infoPhoton.phiCrackDistance,"phiCrackDistance/F");
	treePhoton->Branch("iProcess"           ,&infoPhoton.iProcess           ,"iProcess/I");
	treePhoton->Branch("converted"           ,&infoPhoton.converted           ,"converted/I");
	treePhoton->Branch("invalidNeighbour"           ,&infoPhoton.invalidNeighbour           ,"invalidNeighbour/I");
	treePhoton->Branch("eventIndex"           ,&infoPhoton.eventIndex           ,"eventIndex/I");
	treePhoton->Branch("MVA"           ,&infoPhoton.MVA          ,"MVA/F");

	treeBackground = fs_->make<TTree>("treeBackground","");
	treeBackground->Branch("pt"                 ,&infoBackground.pt                 ,"pt/F");
	treeBackground->Branch("hoe"                ,&infoBackground.hoe                ,"hoe/F");
	treeBackground->Branch("dRPhoton"                ,&infoBackground.dRPhoton                ,"dRPhoton/F");
	treeBackground->Branch("eReco"              ,&infoBackground.eReco              ,"eReco/F");
	treeBackground->Branch("eTrue"              ,&infoBackground.eTrue              ,"eTrue/F");
	treeBackground->Branch("eta"                ,&infoBackground.eta                ,"eta/F");
	treeBackground->Branch("etaSC"              ,&infoBackground.etaSC              ,"etaSC/F");
	treeBackground->Branch("etaEcal"            ,&infoBackground.etaEcal            ,"etaEcal/F");
	treeBackground->Branch("phi"                ,&infoBackground.phi                ,"phi/F");
	treeBackground->Branch("phiSC"              ,&infoBackground.phiSC              ,"phiSC/F");
	treeBackground->Branch("phiEcal"            ,&infoBackground.phiEcal            ,"phiEcal/F");
	treeBackground->Branch("x"                  ,&infoBackground.x                  ,"x/F");
	treeBackground->Branch("y"                  ,&infoBackground.y                  ,"y/F");
	treeBackground->Branch("matchIndex"         ,&infoBackground.matchIndex         ,"matchIndex/I");
	treeBackground->Branch("sigmaEtaEta"        ,&infoBackground.sigmaEtaEta        ,"sigmaEtaEta/F");
	treeBackground->Branch("sigmaRR"            ,&infoBackground.sigmaRR       ,"sigmaRR/F");
	treeBackground->Branch("lengthCompatibility",&infoBackground.lengthCompatibility,"lengthCompatibility/F");
	treeBackground->Branch("phiCrackDistance",&infoBackground.phiCrackDistance,"phiCrackDistance/F");
	treeBackground->Branch("iProcess"           ,&infoBackground.iProcess           ,"iProcess/I");
	treeBackground->Branch("invalidNeighbour"           ,&infoBackground.invalidNeighbour           ,"invalidNeighbour/I");
	treeBackground->Branch("eventIndex"           ,&infoBackground.eventIndex           ,"eventIndex/I");
	treeBackground->Branch("MVA"           ,&infoBackground.MVA          ,"MVA/F");
	treeBackground->Branch("dRBest"            ,&infoBackground.dRBest            ,"dRBest/F");
	cellSize_ = 1;
	doLogWeight_ = true;
	mipE_ = 0.0000551;
	nSR_ = 5;
	truthVtx_ = ROOT::Math::XYZPoint(0,0,0);
}



// destructor
HoverEAnalyzer_MVA::~HoverEAnalyzer_MVA()
{

}


double HoverEAnalyzer_MVA::fillEmIdVars(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters,const edm::PtrVector<HGCRecHit> &  rcs, 	const HGCalGeometry *geom,const edm::Event& iEvent, const edm::EventSetup& iSetup ){

	HGCALShowerBasedEmIdentificationLC2 test(1, geom);
	PCAShowerAnalysis pcaShowerAnalysis(iEvent,iSetup);

	for (unsigned int j =0 ; j < clusters.size() ; j++){

		if (clusters[j]->position()==(sc->seed()->position())) {
			GlobalPoint pcaShowerPos;
			GlobalVector pcaShowerDir;
			pcaShowerAnalysis.showerParameters((clusters[j].get()),pcaShowerPos,pcaShowerDir);

			//test.setShowerPosition(clusters[j]->position());
			//test.setShowerDirection(clusters[j]->axis());
			//
			math::XYZPoint showerPos (pcaShowerPos.x(),pcaShowerPos.y(),pcaShowerPos.z());
			math::XYZVector showerDir (pcaShowerDir.x(),pcaShowerDir.y(),pcaShowerDir.z());

			test.setShowerPosition(showerPos);
			test.setShowerDirection(showerDir);

			std::cout << " DEBUG clusters[j]->position " << clusters[j]->position() << ", clusters[j]->axis() " << clusters[j]->axis() << std::endl;
			std::cout << " DEBUG  pcaShowerPos " << pcaShowerPos << ", pcaShowerDir " << pcaShowerDir << std::endl;

			infoBackground.sigmaEtaEta = test.HGCALShowerBasedEmIdentificationLC2::sigmaetaeta( *(clusters[j].get()), rcs);
			infoPhoton.sigmaEtaEta = infoBackground.sigmaEtaEta; 
			sigmaEtaEta_ = infoPhoton.sigmaEtaEta;
			infoBackground.lengthCompatibility = test.HGCALShowerBasedEmIdentificationLC2::lengthCompatibility( *(clusters[j].get()), rcs);
			//std::cout << "DEBUG shower Start Claude " << test.HGCALShowerBasedEmIdentificationLC2::startPosition0( *(clusters[j].get()))<< ", ME " << test.HGCALShowerBasedEmIdentificationLC2::startPosition( *(clusters[j].get()), rcs) << std::endl; 
			infoPhoton.lengthCompatibility = infoBackground.lengthCompatibility;
			lengthCompatibility_=infoPhoton.lengthCompatibility;

			break;
		}
	}

	return infoPhoton.sigmaEtaEta;
}


	void
HoverEAnalyzer_MVA::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	


	eventIndex++;

	std::cout << "=======================================================" << std::endl;
	using namespace edm;
	iEvent.getByLabel(edm::InputTag("HGCalRecHit:HGCEERecHits"),recHits_);
	if (hcalHelperEndcap_)
	{
		hcalHelperEndcap_->checkSetup(iSetup) ;
		hcalHelperEndcap_->readEvent(iEvent) ;
	}
	// get calo geometry

	edm::ESHandle<HGCalGeometry> geomH;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],geomH);
	const HGCalGeometry *geom=geomH.product();


	if (caloGeomCacheId_!=iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
		iSetup.get<CaloGeometryRecord>().get(caloGeom_);
		caloGeomCacheId_=iSetup.get<CaloGeometryRecord>().cacheIdentifier();
	}
	if (caloTopoCacheId_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
		caloTopoCacheId_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
		iSetup.get<CaloTopologyRecord>().get(caloTopo_);
	}

	//	Handle<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >  > HGCEERechits_;
	//	iEvent.getByToken(endcapRecHitCollection_,HGCEERechits_);
	//	auto  HGCEERechits =HGCEERechits_.product();
	//const edm::SortedCollection<HGCRecHit>& rechits = HGCEERechits->product();

	Handle<edm::View<reco::SuperCluster> > HGCEESCs;
	iEvent.getByToken(endcapSuperClusterCollection_,HGCEESCs);
	const PtrVector<reco::SuperCluster>& sclusters = HGCEESCs->ptrVector();

	Handle<edm::View<reco::PFCluster> > HGCEEClusters;
	iEvent.getByToken(endcapClusterCollection_,HGCEEClusters);
	const PtrVector<reco::PFCluster>& clusters = HGCEEClusters->ptrVector();

	Handle<edm::View<reco::GenParticle> > genParts;
	iEvent.getByToken(genParticlesCollection_,genParts);
	const PtrVector<reco::GenParticle>& gens = genParts->ptrVector();


	Handle<edm::View<HGCRecHit> > eeRecHits;
	iEvent.getByToken(endcapRecHitCollection_, eeRecHits);
	const edm::PtrVector<HGCRecHit>& rechitvec = eeRecHits->ptrVector();

	edm::Handle<edm::SimTrackContainer> SimTk;
	iEvent.getByLabel(g4TracksSource_,SimTk);
	edm::Handle<edm::SimVertexContainer> SimVtx;
	iEvent.getByLabel(g4VerticesSource_,SimVtx);


	int photonFound =0;
	int photonFoundIndex =-999;
	float photonFoundSCEta = -999.;
	float photonFoundSCPhi = -999.;


	for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree

		infoPhoton.eta=-999.;
		infoPhoton.etaSC=-999.;
		infoPhoton.etaEcal=-999.;
		infoPhoton.phi=-999.;
		infoPhoton.phiSC=-999.;
		infoPhoton.phiEcal=-999.;
		infoPhoton.matchIndex=-999;
		infoPhoton.pt=-999.;
		infoPhoton.et=-999.;
		infoPhoton.recopt=-999.;
		infoPhoton.hoe  =-999.;         
		infoPhoton.eReco  =-999.;         
		infoPhoton.eTrue         =-999.;  
		infoPhoton.sigmaEtaEta=-999.;
		infoPhoton.lengthCompatibility=-999.;
		infoPhoton.converted=-999.;
		infoPhoton.invalidNeighbour=-999.;
		infoPhoton.eventIndex=eventIndex;

		float dRBest =999;
	//	float dRMax= 0.01;
		float dRMax= 0.1;

		if (fabs(gens[igp]->pdgId()) != 22) continue;
		if(gens[igp]->status() != 1) continue;
		if(gens[igp]->mother()->status() != 44) continue;
		infoPhoton.eta   = gens[igp]->eta();
		infoPhoton.phi   = gens[igp]->phi();
		infoPhoton.pt    = gens[igp]->pt();
		infoPhoton.et    = (gens[igp]->energy())/std::cosh(infoPhoton.eta);
		infoPhoton.eTrue = gens[igp]->energy();

		float zHGC =320.38;
		float rHGC = 32.0228;
		float RHGC = 152.153;
		if (infoPhoton.eta <0) zHGC = -1*zHGC;
		math::XYZPoint vPos = gens[igp]->vertex();

		float theta = 2*std::atan(std::exp(-infoPhoton.eta));
		float h = std::tan(theta)*(zHGC-vPos.z());
		if (h>RHGC || h<rHGC);// {
		//		std::cout << "debug skip because h>RHGC || h<rHGC gives 1" << std::endl; 
		//		continue;
		//		}
		float eta_ECAL = -std::log(std::tan(std::asin(h/std::sqrt(h*h+zHGC*zHGC))/2)) ;
		if  (infoPhoton.eta <0) eta_ECAL = -1*eta_ECAL;

		infoPhoton.etaEcal =eta_ECAL;
		infoPhoton.phiEcal = infoPhoton.phi;
		infoPhoton.x = (zHGC/std::cos(theta))* std::sin(theta) *std::cos(infoPhoton.phi);
		infoPhoton.y = (zHGC/std::cos(theta))* std::sin(theta) *std::sin(infoPhoton.phi);
		infoPhoton.phiCrackDistance = phiCrackDistance(infoPhoton.x, infoPhoton.y);
		if(fabs(gens[igp]->eta()) <1.5 ||fabs(gens[igp]->eta())  > 3) continue;
		truthVtx_=gens[igp]->vertex();
		math::XYZVectorD hitPos=getInteractionPositionLC(SimTk.product(),SimVtx.product(), gens[igp]->pt()).pos;
		const double z = std::fabs(hitPos.z());
		infoPhoton.converted = (unsigned)(z < 317 && z > 1e-3);


		for (unsigned int isc =0; isc < sclusters.size() ; isc++){ //subloop over sc's to find matches

			int detector = sclusters[isc]->seed()->hitsAndFractions()[0].first.subdetId() ;

			if( detector!=HGCEE) continue;

			float dE = (sclusters[isc]->seed()->eta() - eta_ECAL);
			dE =dE*dE;
			float dP = deltaPhi(sclusters[isc]->seed()->phi() , gens[igp]->phi());
			dP =dP*dP;
			float dR = sqrt(dE +dP);

			//	std::cout << "DEBUG SC " << isc << ", gen e, " << gens[igp]->eta() << ", " << gens[igp]->phi() << ", sc " << sclusters[isc]->eta() << " " <<  sclusters[isc]->phi() << ", dR " << dR <<", dRBest " << dRBest << " , detector!=HGCEE " << (detector!=HGCEE) << std::endl;

			if (dR < dRBest && dR<dRMax) { // only true if dR is both below limit value and smaller than previous best dR.
				dRBest = dR;
				infoPhoton.matchIndex = isc;
			  infoPhoton.dRBest =dRBest;
			}
		}

		//float recoPt = -999.;
		if (infoPhoton.matchIndex >-1) {
			infoPhoton.eReco = annemarieEnergy(rechitvec,sclusters[infoPhoton.matchIndex],geometrySource_,iEvent, iSetup, eta_ECAL ) ;


			photonFound =1;
			photonFoundIndex= infoPhoton.matchIndex;
			infoPhoton.recopt =  infoPhoton.eReco/std::cosh(infoPhoton.etaSC);
			infoPhoton.etaSC    = sclusters[infoPhoton.matchIndex]->eta();
			infoPhoton.phiSC    = sclusters[infoPhoton.matchIndex]->phi();
			photonFoundSCEta= infoPhoton.etaSC;
			photonFoundSCPhi= infoPhoton.phiSC;

			std::cout << "debug - found photon. gen e,p "<< infoPhoton.eta <<", " << infoPhoton.phi << ", e,p ecal "<< infoPhoton.etaEcal << ", " << infoPhoton.phiEcal << ", e,p SC " << infoPhoton.etaSC<< ", "<< infoPhoton.phiSC << " ereco " << infoPhoton.eReco <<", eTrue " << infoPhoton.eTrue <<" CONVERTED " << infoPhoton.converted<< std::endl;
			std::cout << "debug - photonFoundIndex"<< photonFoundIndex <<", photonFoundSCEta"  << photonFoundSCEta<< ", photonFoundSCPhi " <<  photonFoundSCPhi<<std::endl;


			double sieie =  fillEmIdVars(sclusters[infoPhoton.matchIndex], clusters,  rechitvec, geom,iEvent, iSetup);
			if (sieie);


			const auto & scl = (*HGCEESCs)[infoPhoton.matchIndex] ;
			float had =hcalHelperEndcap_->HCALClustersBehindSC(scl);

			float  scle =  sclusters[infoPhoton.matchIndex]->energy();
			float hoe =had/scle;
			infoPhoton.hoe = hoe;
			hoe_ = hoe;
			infoPhoton.MVA = Mva_->EvaluateMVA( "BDTG" );

		}
			treePhoton->Fill();


	}

	//--------------> End per-genPhoton tree <------------------

	if(photonFound){ //only want SCs which are NOT the photon.
		//--------------> Begin per-Background tree <---------------------- 
		int nSCAnalyzed =0;

		for( unsigned int isc =0; isc< sclusters.size() ; isc++){ // isc = index_super_cluster

			infoBackground.pt=-999.;
			infoBackground.hoe=-999.;
			infoBackground.dRPhoton=-999.;
			infoBackground.eta=-999.;
			infoBackground.phi=-999.;
			infoBackground.eReco=-999.;
			infoBackground.eTrue=-999.;
			infoBackground.matchIndex=-999;
			infoBackground.invalidNeighbour=-999;
			infoBackground.eventIndex= eventIndex;

			float dEPho = (sclusters[isc]->eta() - photonFoundSCEta);
			dEPho =dEPho*dEPho;
			float dPPho = deltaPhi(sclusters[isc]->phi() , photonFoundSCPhi);
			dPPho =dPPho*dPPho;
			float dRPho = sqrt(dEPho +dPPho);
			infoBackground.dRPhoton = dRPho;

			std::cout << "debug - photonFoundIndex"<< photonFoundIndex <<", photonFoundSCEta"  << photonFoundSCEta<< ", photonFoundSCPhi " <<  photonFoundSCPhi<< "  sclusters[isc]->eta() " << sclusters[isc]->eta() <<", sclusters[isc]->phi()" << sclusters[isc]->phi()<< std::endl;
			std::cout << " debug dRPho " << dRPho << ", dEPho " << dEPho << ", dPPho " << dPPho << std::endl;	
			if ( dRPho<1. ) continue;
			//infoBackground.pt= sclusters[isc]->energy()/std::cosh(sclusters[isc]->eta());

			//	if (infoBackground.pt <4) continue;
			infoBackground.pt= sclusters[isc]->energy() /std::cosh(sclusters[isc]->eta());

			if( infoBackground.pt <20.) continue;
			//	infoBackground.eta=sclusters[isc]->eta();
			infoBackground.etaSC=sclusters[isc]->eta();
			infoBackground.eReco = annemarieEnergy(rechitvec,sclusters[isc],geometrySource_,iEvent, iSetup, infoBackground.eta ,0) ;
			//	infoBackground.phi=sclusters[isc]->phi();
			infoBackground.phiSC=sclusters[isc]->phi();
			infoBackground.dRBest = infoPhoton.dRBest;
			//	std::cout << " SC " << isc <<"  added, pt " <<infoBackground.pt <<  std::endl;
			//	std::cout << " jet match truth " << infoBackground.eta << ", " << infoBackground.phi << ", reco " << infoBackground.etaSC << ", " << infoBackground.phiSC << std::endl;



			double sieie =  fillEmIdVars(sclusters[isc], clusters,  rechitvec,  geom,iEvent, iSetup  );
			if (sieie);

			const auto & scl = (*HGCEESCs)[isc] ;
			float had =hcalHelperEndcap_->HCALClustersBehindSC(scl);

			float  scle =  sclusters[isc]->energy();
			float hoe =had/scle;
			infoBackground.hoe = hoe;
			hoe_ =hoe;
			//	std::cout << " debug quark/gluon  " <<gens[igp]->pdgId()<<", status " <<  gens[igp]->status() <<"  pass " << gens[igp]->eta() << ", " << gens[igp]->phi() <<", status " << gens[igp]->status() << ", mother status " << gens[igp]->mother()->status() << ", dRbest " << dRBest << ", eReco " <<  infoBackground.eReco << ", eSC " << scle << ", etrue " << infoBackground.eTrue  << ", eReco/eTrue >0.7 ?" << (infoBackground.eReco/infoBackground.eTrue >0.7)<< ", eSC/eTrue >0.7 ?" << (scle/infoBackground.eTrue >0.7)<< std::endl;
			infoBackground.MVA = Mva_->EvaluateMVA( "BDTG" );

			treeBackground->Fill();


		}
		std::cout << "nSCAnalyzed " << nSCAnalyzed << std::endl;
	}


	/*
		 if(photonFound){ //only want SCs which are NOT the photon.
//--------------> Begin per-Background tree <---------------------- 
    actory->AddVariable("sigmaEtaEta" );
int nSCAnalyzed =0;

for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree

infoBackground.pt=-999.;
infoBackground.et=-999.;
infoBackground.hoe=-999.;
infoBackground.dRPhoton=-999.;
infoBackground.eta=-999.;
infoBackground.phi=-999.;
infoBackground.eReco=-999.;
infoBackground.eTrue=-999.;
infoBackground.matchIndex=-999;
infoBackground.invalidNeighbour=-999;
infoBackground.eventIndex= eventIndex;

float dRBest =999;

if (fabs(gens[igp]->pdgId()) > 8 && gens[igp]->pdgId()!=21  ) continue;
if(gens[igp]->status() != 23) {continue; }
//	std::cout << " debug --- status 23 pass " << gens[igp]->eta() << ", " << gens[igp]->phi() << std::endl;
//	if(gens[igp]->mother()->status() != 44) {continue; }
//	std::cout << " debug --- --- mother status 44  pass " << gens[igp]->eta() << ", " << gens[igp]->phi() << std::endl;

infoBackground.eta   = gens[igp]->eta();
infoBackground.phi   = gens[igp]->phi();
infoBackground.pt    = gens[igp]->pt();
infoBackground.et    = (gens[igp]->energy())/std::cosh(infoBackground.eta);
infoBackground.eTrue = gens[igp]->energy();

float zHGC =320.38;
float rHGC = 32.0228;
float RHGC = 152.153;
if (infoBackground.eta <0) zHGC = -1*zHGC;
math::XYZPoint vPos = gens[igp]->vertex();

float theta = 2*std::atan(std::exp(-infoBackground.eta));
float h = std::tan(theta)*(zHGC-vPos.z());
if (h>RHGC || h<rHGC);// {
//		std::cout << "debug skip because h>RHGC || h<rHGC gives 1" << std::endl; 
//		continue;
//		}
float eta_ECAL = -std::log(std::tan(std::asin(h/std::sqrt(h*h+zHGC*zHGC))/2)) ;
if  (infoBackground.eta <0) eta_ECAL = -1*eta_ECAL;

infoBackground.etaEcal =eta_ECAL;
infoBackground.phiEcal = infoBackground.phi;
infoBackground.x = (zHGC/std::cos(theta))* std::sin(theta) *std::cos(infoBackground.phi);
infoBackground.y = (zHGC/std::cos(theta))* std::sin(theta) *std::sin(infoBackground.phi);
infoBackground.phiCrackDistance = phiCrackDistance(infoBackground.x, infoBackground.y);
if(fabs(eta_ECAL) <1.5 ||fabs(eta_ECAL)  > 3) continue;
truthVtx_=gens[igp]->vertex();

for( unsigned int isc =0; isc< sclusters.size() ; isc++){ // isc = index_super_cluster

if ((int) isc == photonFoundIndex) continue;
int detectorSC = sclusters[isc]->seed()->hitsAndFractions()[0].first.subdetId() ;
if ( detectorSC!=HGCEE) { continue;}

infoBackground.pt= sclusters[isc]->energy() /std::cosh(sclusters[isc]->eta());
if(infoBackground.pt <20.) continue;
float dE = (sclusters[isc]->eta() - gens[igp]->eta());
dE =dE*dE;
float dP = deltaPhi(sclusters[isc]->phi() , gens[igp]->phi());
dP =dP*dP;
float dR = sqrt(dE +dP);




if (dR < dRBest && dR<0.5) { // only true if dR is both below limit value and smaller than previous best dR.
	dRBest = dR;
	infoBackground.matchIndex = isc;
}
}


nSCAnalyzed++;


if (infoBackground.matchIndex>-1){
	float dEPho = (sclusters[infoBackground.matchIndex]->eta() - photonFoundSCEta);
	dEPho =dEPho*dEPho;
	float dPPho = deltaPhi(sclusters[infoBackground.matchIndex]->phi() , photonFoundSCPhi);
	dPPho =dPPho*dPPho;
	float dRPho = sqrt(dEPho +dPPho);
	infoBackground.dRPhoton = dRPho;
	std::cout << "debug - photonFoundIndex"<< photonFoundIndex <<", photonFoundSCEta"  << photonFoundSCEta<< ", photonFoundSCPhi " <<  photonFoundSCPhi<< "  sclusters[infoBackground.matchIndex]->eta() " << sclusters[infoBackground.matchIndex]->eta() <<", sclusters[infoBackground.matchIndex]->phi()" << sclusters[infoBackground.matchIndex]->phi()<< std::endl;
	std::cout << " debug dRPho " << dRPho << ", dEPho " << dEPho << ", dPPho " << dPPho << std::endl;	
	if ( dRPho<0.5 ) continue;
	//infoBackground.pt= sclusters[infoBackground.matchIndex]->energy()/std::cosh(sclusters[infoBackground.matchIndex]->eta());

	//	if (infoBackground.pt <4) continue;
	infoBackground.et= sclusters[infoBackground.matchIndex]->energy() /std::cosh(sclusters[infoBackground.matchIndex]->eta());
	//	infoBackground.eta=sclusters[infoBackground.matchIndex]->eta();
	infoBackground.etaSC=sclusters[infoBackground.matchIndex]->eta();
	infoBackground.eReco = annemarieEnergy(rechitvec,sclusters[infoBackground.matchIndex],geometrySource_,iEvent, iSetup, infoBackground.eta ,0) ;
	//	infoBackground.phi=sclusters[infoBackground.matchIndex]->phi();
	infoBackground.phiSC=sclusters[infoBackground.matchIndex]->phi();
	//	std::cout << " SC " << infoBackground.matchIndex <<"  added, pt " <<infoBackground.pt <<  std::endl;
	//	std::cout << " jet match truth " << infoBackground.eta << ", " << infoBackground.phi << ", reco " << infoBackground.etaSC << ", " << infoBackground.phiSC << std::endl;



	double sieie =  fillEmIdVars(sclusters[infoBackground.matchIndex], clusters,  rechitvec,  geom, iEvent, iSetup);
	if (sieie);

	const auto & scl = (*HGCEESCs)[infoBackground.matchIndex] ;
	float had =hcalHelperEndcap_->HCALClustersBehindSC(scl);

	float  scle =  sclusters[infoBackground.matchIndex]->energy();
	float hoe =had/scle;
	infoBackground.hoe = hoe;
	std::cout <<" DEBUG SC " << infoBackground.lengthCompatibility << std::endl;
	//	std::cout << " debug quark/gluon  " <<gens[igp]->pdgId()<<", status " <<  gens[igp]->status() <<"  pass " << gens[igp]->eta() << ", " << gens[igp]->phi() <<", status " << gens[igp]->status() << ", mother status " << gens[igp]->mother()->status() << ", dRbest " << dRBest << ", eReco " <<  infoBackground.eReco << ", eSC " << scle << ", etrue " << infoBackground.eTrue  << ", eReco/eTrue >0.7 ?" << (infoBackground.eReco/infoBackground.eTrue >0.7)<< ", eSC/eTrue >0.7 ?" << (scle/infoBackground.eTrue >0.7)<< std::endl;

	treeBackground->Fill();
}

}
}
*/

return ;
}

float HoverEAnalyzer_MVA::annemarieEnergy( const edm::PtrVector<HGCRecHit>& rechitvec, const edm::Ptr<reco::SuperCluster>& sc, 	std::vector<std::string> geometrySource_, const edm::Event& iEvent, const edm::EventSetup& iSetup, float eta_ECAL, int sig){

	std::map<int,const HGCalGeometry *> hgcGeometries;
	edm::ESHandle<HGCalGeometry> hgcGeo;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
	hgcGeometries[0]=hgcGeo.product();

	int nLayers_=30;
	std::vector<HGCEEDetId> detidmax;
	detidmax.resize(nLayers_,HGCEEDetId());

	std::vector<HGCEEDetId> detidstack;
	detidstack.resize(9*nLayers_,HGCEEDetId());
	double  phimax = infoPhoton.phi;
	double  etamax = infoPhoton.eta;//sc;
	if (sig){
		phimax = infoPhoton.phi;
		etamax = infoPhoton.eta;//sc;
	} else {
		phimax = infoBackground.phi;
		etamax = infoBackground.eta;//sc;
	}


	PCAShowerAnalysis pcaShowerAnalysis(iEvent,iSetup);

	GlobalPoint pcaShowerPos;
	GlobalVector pcaShowerDir;
	pcaShowerAnalysis.showerParameters((sc->seed().get()),pcaShowerPos,pcaShowerDir);


	double clus_eta = sc->seed()->eta();

	getMaximumCell(rechitvec,phimax,etamax,detidmax);

	fillDetIdStack(detidmax,detidstack);

	HGCALShowerBasedEmIdentificationLC2 test(1, hgcEEGeom_);

	const HGCalTopology& topology = hgcEEGeom_->topology();
	double sRR =  test.HGCALShowerBasedEmIdentificationLC2::sigmaRR(detidstack , recHits_);
	infoPhoton.sigmaRR =sRR;
	infoBackground.sigmaRR =sRR;

	float energyReco=0;

	const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
	double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
	double mip = 0.0000551;
	double weight[30] = {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};
	for (int iLayer =0; iLayer < 30 ; iLayer++){
		std::vector<double> Exy;
		Exy.resize(9,0);
		if (topology.valid(detidmax[iLayer])) {//fillNeighbours(detidmax[iL],Exy,info);
			auto detidmaxL = detidmax[iLayer]; 
			std::vector<HGCEEDetId> neighbours;
			neighbours.resize(9,HGCEEDetId());
			DetId tmp = topology.goSouth(detidmaxL);
			if (topology.valid(tmp)) {
				neighbours[1] = HGCEEDetId(tmp);
				tmp = topology.goWest(neighbours[1]);
				if (topology.valid(tmp)) neighbours[0] = HGCEEDetId(tmp);
				tmp = topology.goEast(neighbours[1]);
				if (topology.valid(tmp)) neighbours[2] = HGCEEDetId(tmp);
			}

			neighbours[4] = detidmaxL;
			tmp = topology.goWest(neighbours[4]);
			if (topology.valid(tmp)) neighbours[3] = HGCEEDetId(tmp);
			tmp = topology.goEast(neighbours[4]);
			if (topology.valid(tmp)) neighbours[5] = HGCEEDetId(tmp);

			tmp = topology.goNorth(detidmaxL);
			if (topology.valid(tmp)) {
				neighbours[7] = HGCEEDetId(tmp);
				tmp = topology.goWest(neighbours[7]);
				if (topology.valid(tmp)) neighbours[6] = HGCEEDetId(tmp);
				tmp = topology.goEast(neighbours[7]);
				if (topology.valid(tmp)) neighbours[8] = HGCEEDetId(tmp);
			}

			GlobalPoint center = hgcEEGeom_->getPosition(detidmaxL);


			for (unsigned idx(0);idx<9;++idx){
				if (!topology.valid(neighbours[idx])|| neighbours[idx].det()!=DetId::Forward || neighbours[idx].subdetId()!=HGCEE) {
					infoPhoton.invalidNeighbour = 1;
					infoBackground.invalidNeighbour = 1;
					continue;
				}
				HGCRecHitCollection::const_iterator theHit = recHits_->find(neighbours[idx]);
				if (theHit==recHits_->end()) continue;
				if (neighbours[idx].det()!=DetId::Forward || neighbours[idx].subdetId()!=HGCEE) continue;
				GlobalPoint cellPos = hgcEEGeom_->getPosition(neighbours[idx]);
				double posx = cellPos.x();
				double posy = cellPos.y();
				double posz = cellPos.z();
				if (fabs(posx-center.x())>sqrt(2) ||
						fabs(posy-center.y())>sqrt(2)){
					infoPhoton.invalidNeighbour =1;
					infoBackground.invalidNeighbour =1;
				}

				double costheta = fabs(posz)/sqrt(posz*posz+posx*posx+posy*posy);
				double energy = theHit->energy()/mipE_*costheta;//in MIP
				Exy[idx] = energy;
				double scale = mip*corr;
				if (scale && weight[0]);
				energyReco += energy*absWeight(iLayer);

			}

		} else {
			continue;
		}
	}
	energyReco=energyReco/fabs(tanh(eta_ECAL));
	double pars[3] = {69.5,4.5,-0.8};
	double paro[3] = {-50,0,0};  
	double offset = paro[0] + paro[1]*fabs(eta_ECAL) + paro[2]*eta_ECAL*eta_ECAL;
	double slope = pars[0] + pars[1]*fabs(eta_ECAL) + pars[2]*eta_ECAL*eta_ECAL;
	energyReco= (energyReco-offset)/slope;
	return energyReco;

}

float HoverEAnalyzer_MVA::phiCrackDistance (float x, float y){
	float d=999;
	float dmin=999;
	for (int N =0; N< 9; N++){

		float a=tan(N*0.349+0.1745);
		d = fabs( a*x -y)/sqrt(a*a + 1);

		if (d<dmin) {
			dmin=d;
		}
	}

	return dmin;
}

void HoverEAnalyzer_MVA::getMaximumCell(const edm::PtrVector<HGCRecHit>& rechitvec,const double & phimax,const double & etamax,std::vector<HGCEEDetId> & detidmax){
	unsigned int nLayers_=30;
	int debug_= 0;
	std::vector<double> dRmin;
	dRmin.resize(nLayers_,10);
	if (debug_) std::cout << " - Processing " << rechitvec.size() << " rechits, phimax "<< phimax << ", etamax " << etamax << "1 Nlayers " << nLayers_ <<std::endl;
	for (unsigned iH(0); iH<rechitvec.size(); ++iH){//loop on rechits
		const HGCRecHit & lHit = *(rechitvec[iH]);
		const HGCEEDetId & hgcid = lHit.detid();
		unsigned layer = hgcid.layer()-1;
		if (layer >= nLayers_) {
			std::cout << " -- Warning! Wrong layer number: " << layer << " max is set to " << nLayers_ << std::endl;
			continue;
		}
		GlobalPoint cellPos = hgcEEGeom_->getPosition(hgcid);
		double posx = cellPos.x();
		double posy = cellPos.y();
		double posz = cellPos.z();
		if (debug_>1) {
			std::cout << " --  RecoHit " << iH << "/" << rechitvec.size() << " -- layer " << layer << " id " << hgcid << std::endl
				<< " --  position x,y,z " << posx << "," << posy << "," << posz << std::endl;
		}

		ROOT::Math::XYZVector pos(posx-truthVtx_.x(),posy-truthVtx_.y(),posz-truthVtx_.z());
		double deta = fabs(pos.eta()-etamax);
		double dphi = deltaPhi(pos.phi(),phimax);

		double dR = sqrt(pow(deta,2)+pow(dphi,2));
		if (dR<dRmin[layer]) {
			dRmin[layer] = dR;
			detidmax[layer] = hgcid;

		}

	}//loop on rechits

}

double HoverEAnalyzer_MVA::absWeight(const unsigned layer, const bool dedx){
	if (dedx==false){
		if (layer == 0) return 0.08696;
		if (layer == 1) return 1;//0.92
		if (layer == 2) return 0.646989;//88.16/95.4=0.92
		if (layer == 3) return 0.617619;//51.245/95.4=0.537
		if (layer == 4) return 0.646989;
		if (layer == 5) return 0.617619;
		if (layer == 6) return 0.646989;
		if (layer == 7) return 0.617619;
		if (layer == 8) return 0.646989;
		if (layer == 9) return 0.617619;
		if (layer == 10) return 0.646989;
		if (layer == 11) return 0.942829;//74.45/95.4=0.78
		if (layer == 12) return 0.859702;//102.174/95.4=1.071
		if (layer == 13) return 0.942829;
		if (layer == 14) return 0.859702;
		if (layer == 15) return 0.942829;
		if (layer == 16) return 0.859702;
		if (layer == 17) return 0.942829;
		if (layer == 18) return 0.859702;
		if (layer == 19) return 0.942829;
		if (layer == 20) return 0.859702;
		if (layer == 21) return 1.37644;//105.39/95.4=1.1047
		if (layer == 22) return 1.30447;//131.476/95.4=1.378
		if (layer == 23) return 1.37644;
		if (layer == 24) return 1.30447;
		if (layer == 25) return 1.37644;
		if (layer == 26) return 1.30447;
		if (layer == 27) return 1.37644;
		if (layer == 28) return 1.30447;
		if (layer == 29) return 1.37644;//1.79662;//
	}
	else {
		if (layer == 0) return 0.06588;
		if (layer == 1) return 1;//95.4/95.4=1
		if (layer == 2) return 0.92;//88.16/95.4=0.92
		if (layer == 3) return 0.537;//51.245/95.4=0.537
		if (layer == 4) return 0.92;
		if (layer == 5) return 0.537;
		if (layer == 6) return 0.92;
		if (layer == 7) return 0.537;
		if (layer == 8) return 0.92;
		if (layer == 9) return 0.537;
		if (layer == 10) return 0.92;
		if (layer == 11) return 0.78;//74.45/95.4=0.78
		if (layer == 12) return 1.071;//102.174/95.4=1.071
		if (layer == 13) return 0.78;
		if (layer == 14) return 1.071;
		if (layer == 15) return 0.78;
		if (layer == 16) return 1.071;
		if (layer == 17) return 0.78;
		if (layer == 18) return 1.071;
		if (layer == 19) return 0.78;
		if (layer == 20) return 1.071;
		if (layer == 21) return 1.1047;//105.39/95.4=1.1047
		if (layer == 22) return 1.378;//131.476/95.4=1.378
		if (layer == 23) return 1.1047;
		if (layer == 24) return 1.378;
		if (layer == 25) return 1.1047;
		if (layer == 26) return 1.378;
		if (layer == 27) return 1.1047;
		if (layer == 28) return 1.378;
		if (layer == 29) return 1.1047;
	}
	return 1;
}


void HoverEAnalyzer_MVA::fillDetIdStack(std::vector<HGCEEDetId> & detidmax_vec, std::vector<HGCEEDetId> & detidstack){
	unsigned int nLayers_=30;
	for( unsigned int iL =0;iL< nLayers_ ; iL++){
		auto detidmax = detidmax_vec[iL];
		const HGCalTopology& topology = hgcEEGeom_->topology();
		std::vector<HGCEEDetId> neighbours;
		neighbours.resize(9,HGCEEDetId());
		DetId tmp = topology.goSouth(detidmax);
		if (topology.valid(tmp)) {
			neighbours[1] = HGCEEDetId(tmp);
			tmp = topology.goWest(neighbours[1]);
			if (topology.valid(tmp)) neighbours[0] = HGCEEDetId(tmp);
			tmp = topology.goEast(neighbours[1]);
			if (topology.valid(tmp)) neighbours[2] = HGCEEDetId(tmp);
		}

		neighbours[4] = detidmax;
		tmp = topology.goWest(neighbours[4]);
		if (topology.valid(tmp)) neighbours[3] = HGCEEDetId(tmp);
		tmp = topology.goEast(neighbours[4]);
		if (topology.valid(tmp)) neighbours[5] = HGCEEDetId(tmp);

		tmp = topology.goNorth(detidmax);
		if (topology.valid(tmp)) {
			neighbours[7] = HGCEEDetId(tmp);
			tmp = topology.goWest(neighbours[7]);
			if (topology.valid(tmp)) neighbours[6] = HGCEEDetId(tmp);
			tmp = topology.goEast(neighbours[7]);
			if (topology.valid(tmp)) neighbours[8] = HGCEEDetId(tmp);
		}

		for (unsigned idx(0);idx<9;++idx){

			detidstack[iL*9+idx]=neighbours[idx];

		}
	}
}


// ------------ method called once each job just before starting event loop  ------------
	void 
HoverEAnalyzer_MVA::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
HoverEAnalyzer_MVA::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
HoverEAnalyzer_MVA::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
	edm::ESHandle<HGCalGeometry> hgcGeo;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
	hgcEEGeom_=hgcGeo.product();
}

// ------------ method called when ending the processing of a run  ------------
	void 
HoverEAnalyzer_MVA::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
HoverEAnalyzer_MVA::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
HoverEAnalyzer_MVA::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoverEAnalyzer_MVA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoverEAnalyzer_MVA);
