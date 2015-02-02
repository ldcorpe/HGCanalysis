// -*- C++ -*-
//
// Package:    HoverEAnalyzer
// Class:      HoverEAnalyzer
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


using namespace std;

//
// class declaration
//
// limit dR allowed for geometrical matching of reco/gen particles
float dRlim2 =0.05;

// information to be loaded into TTree
struct infoTruth_t {

	float  eReco_over_eTrue;
	float eta;
	float etaSC;
	float etaSeed;
	float phi;
	float phiSC;
	float phiSeed;
	int matchIndex;
	float pt;
	float eTrue;
	float eReco;
	float hoe;
	float eSeed;
	float eSeed_over_eReco;
	float phiWidth;
	float etaWidth;
	int clustersSize;
	int nClusters09;
	float sigmaEtaEtaNoDir;
	float showerStartPosNoDir;
	float lengthCompatibilityNoDir;
	float sigmaEtaEta;
	float showerStartPos;
	float lengthCompatibility;
	int iProcess;

};

struct info_t {
	float pt;
	float eta;
	float hoe;
	float phi;
	float eReco;
	float eTrue;
	int  matchIndex;
	float  eReco_over_eTrue;
	float mass;
	int clustersSize;
	float sigmaEtaEta;
	float showerStartPos;
	float lengthCompatibility;
	int iProcess;
};

// .h class info
class HoverEAnalyzer : public edm::EDAnalyzer {
	public:
		explicit HoverEAnalyzer(const edm::ParameterSet&);
		~HoverEAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		double fillEmIdVarsNoDir(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters,const  edm::SortedCollection<HGCRecHit>& rcs, const HGCalGeometry *geom );
		double fillEmIdVars(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters,const  edm::SortedCollection<HGCRecHit>& rcs, const HGCalGeometry *geom );
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



		ElectronHcalHelper * hcalHelperEndcap_ ;
		edm::ESHandle<CaloGeometry> caloGeom_ ;
		unsigned long long caloGeomCacheId_ ;
		edm::ESHandle<CaloTopology> caloTopo_;
		unsigned long long caloTopoCacheId_;
		edm::EDGetTokenT<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >endcapRecHitCollection_      ; 
		edm::EDGetTokenT<edm::View<reco::SuperCluster> >endcapSuperClusterCollection_;
		edm::EDGetTokenT<edm::View<reco::PFCluster> >endcapClusterCollection_     ;
		edm::EDGetTokenT<edm::View<reco::GenParticle> >genParticlesCollection_     ;
		edm::EDGetTokenT<edm::View<reco::PFRecHit> > eeRecHitCollection_     ;
		std::string geometrySource_;

		int PU;
		int process;


		TH1F *hoe_hgg_h; 
		TH1F *hoe_zee_h; 
		TH1F *hoe_qcd_h; 

		TH1F *see_hgg_h; 
		TH1F *see_zee_h; 
		TH1F *see_qcd_h; 

		TH1F *ssz_hgg_h; 
		TH1F *ssz_zee_h; 
		TH1F *ssz_qcd_h; 

		TH1F *lcp_hgg_h; 
		TH1F *lcp_zee_h; 
		TH1F *lcp_qcd_h; 


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
HoverEAnalyzer::HoverEAnalyzer(const edm::ParameterSet& iConfig):
	hcalHelperEndcap_(0), caloGeom_(0), caloGeomCacheId_(0), caloTopo_(0), caloTopoCacheId_(0),
	endcapRecHitCollection_(consumes  <edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > > (iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles")))),
	eeRecHitCollection_(consumes <edm::View<reco::PFRecHit> > (iConfig.getUntrackedParameter<edm::InputTag>("eeRecHitCollection",edm::InputTag("particleFlowRecHitHGCEE:Cleaned"))))

{

	geometrySource_ = iConfig.getUntrackedParameter< std::string >("geometrySource");
	ElectronHcalHelper::Configuration hcalCfgEndcap ;
	hcalCfgEndcap.hOverEConeSize = iConfig.getParameter<double>("hOverEConeSize") ;
	hcalCfgEndcap.hOverEMethod = iConfig.getParameter<int>("hOverEMethodEndcap") ;
	PU  = iConfig.getParameter<int>("PU") ;
	process  = iConfig.getParameter<int>("process_") ;
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

	hoe_hgg_h = fs_->make<TH1F>("hoe_hgg_h","hoe_hgg_h",200,0,5);
	hoe_zee_h = fs_->make<TH1F>("hoe_zee_h","hoe_zee_h",200,0,5);
	hoe_qcd_h = fs_->make<TH1F>("hoe_qcd_h","hoe_qcd_h",200,0,5);

	see_hgg_h = fs_->make<TH1F>("see_hgg_h","see_hgg_h",200,0,0.025);
	see_zee_h = fs_->make<TH1F>("see_zee_h","see_zee_h",200,0,0.025);
	see_qcd_h = fs_->make<TH1F>("see_qcd_h","see_qcd_h",200,0,0.025);

	ssz_hgg_h = fs_->make<TH1F>("ssz_hgg_h","ssz_hgg_h",100,300,350);
	ssz_zee_h = fs_->make<TH1F>("ssz_zee_h","ssz_zee_h",100,300,350);
	ssz_qcd_h = fs_->make<TH1F>("ssz_qcd_h","ssz_qcd_h",100,300,350);

	lcp_hgg_h = fs_->make<TH1F>("lcp_hgg_h","lcp_hgg_h",200,-5,5);
	lcp_zee_h = fs_->make<TH1F>("lcp_zee_h","lcp_zee_h",200,-5,5);
	lcp_qcd_h = fs_->make<TH1F>("lcp_qcd_h","lcp_qcd_h",200,-5,5);

	hoe_hgg_h->SetLineColor(kRed);
	see_hgg_h->SetLineColor(kRed);
	ssz_hgg_h->SetLineColor(kRed);
	lcp_hgg_h->SetLineColor(kRed);

	hoe_zee_h->SetLineColor(kGreen);
	see_zee_h->SetLineColor(kGreen);
	ssz_zee_h->SetLineColor(kGreen);
	lcp_zee_h->SetLineColor(kGreen);

	hoe_qcd_h->SetLineColor(kBlue);
	see_qcd_h->SetLineColor(kBlue);
	ssz_qcd_h->SetLineColor(kBlue);
	lcp_qcd_h->SetLineColor(kBlue);

	hoe_hgg_h->Sumw2();
	see_hgg_h->Sumw2();
	ssz_hgg_h->Sumw2();
	lcp_hgg_h->Sumw2();

	hoe_zee_h->Sumw2();
	see_zee_h->Sumw2();
	ssz_zee_h->Sumw2();
	lcp_zee_h->Sumw2();

	hoe_qcd_h->Sumw2();
	see_qcd_h->Sumw2();
	ssz_qcd_h->Sumw2();
	lcp_qcd_h->Sumw2();


	tree = fs_->make<TTree>("tree","");
	tree->Branch("eReco_over_eTrue"              ,&info.eReco_over_eTrue             ,"eReco_over_eTrue/F");
	tree->Branch("pt"              ,&info.pt             ,"pt/F");
	tree->Branch("eta"              ,&info.eta             ,"eta/F");
	tree->Branch("hoe"              ,&info.hoe            ,"hoe/F");
	tree->Branch("phi"              ,&info.phi             ,"phi/F");
	tree->Branch("eReco"              ,&info.eReco            ,"eReco/F");
	tree->Branch("eTrue"              ,&info.eTrue            ,"eTrue/F");
	tree->Branch("matchIndex"              ,&info.matchIndex            ,"matchIndex/I");
	tree->Branch("clustersSize"              ,&info.clustersSize            ,"clustersSize/I");
	tree->Branch("iProcess"              ,&info.iProcess         ,"iProcess/I");
	tree->Branch("sigmaEtaEta"              ,&info.sigmaEtaEta            ,"sigmaEtaEta/F");
	tree->Branch("showerStartPos"              ,&info.showerStartPos        ,"showerStartPos/F");
	tree->Branch("lengthCompatibility"              ,&info.lengthCompatibility       ,"lengthCompatibility/F");

	treeTruth = fs_->make<TTree>("treeTruth","");
	treeTruth->Branch("pt"              ,&infoTruth.pt            ,"pt/F");
	treeTruth->Branch("hoe"              ,&infoTruth.hoe            ,"hoe/F");
	treeTruth->Branch("eReco"              ,&infoTruth.eReco            ,"eReco/F");
	treeTruth->Branch("eTrue"              ,&infoTruth.eTrue            ,"eTrue/F");
	treeTruth->Branch("eSeed"              ,&infoTruth.eSeed            ,"eSeed/F");
	treeTruth->Branch("eSeed_over_eReco"   ,&infoTruth.eSeed_over_eReco             ,"eSeed_over_eReco/F");
	treeTruth->Branch("eReco_over_eTrue"              ,&infoTruth.eReco_over_eTrue             ,"eReco_over_eTrue/F");
	treeTruth->Branch("eta"              ,&infoTruth.eta             ,"eta/F");
	treeTruth->Branch("etaSC"              ,&infoTruth.etaSC             ,"etaSC/F");
	treeTruth->Branch("etaSeed"              ,&infoTruth.etaSeed             ,"etaSeed/F");
	treeTruth->Branch("phi"              ,&infoTruth.phi             ,"phi/F");
	treeTruth->Branch("phiSC"              ,&infoTruth.phiSC             ,"phiSC/F");
	treeTruth->Branch("phiSeed"              ,&infoTruth.phiSeed             ,"phiSeed/F");
	treeTruth->Branch("matchIndex"              ,&infoTruth.matchIndex            ,"matchIndex/I");
	treeTruth->Branch("clustersSize"              ,&infoTruth.clustersSize           ,"clustersSize/I");
	treeTruth->Branch("nClusters09"              ,&infoTruth.nClusters09           ,"nClusters09/I");
	treeTruth->Branch("etaWidth"              ,&infoTruth.etaWidth             ,"etaWidth/F");
	treeTruth->Branch("phiWidth"              ,&infoTruth.phiWidth             ,"phiWidth/F");
	treeTruth->Branch("sigmaEtaEta"              ,&infoTruth.sigmaEtaEta            ,"sigmaEtaEta/F");
	treeTruth->Branch("showerStartPos"              ,&infoTruth.showerStartPos        ,"showerStartPos/F");
	treeTruth->Branch("lengthCompatibility"              ,&infoTruth.lengthCompatibility       ,"lengthCompatibility/F");
	treeTruth->Branch("sigmaEtaEtaNoDir"              ,&infoTruth.sigmaEtaEtaNoDir            ,"sigmaEtaEtaNoDir/F");
	treeTruth->Branch("showerStartPosNoDir"              ,&infoTruth.showerStartPosNoDir        ,"showerStartPosNoDir/F");
	treeTruth->Branch("lengthCompatibilityNoDir"              ,&infoTruth.lengthCompatibilityNoDir       ,"lengthCompatibilityNoDir/F");
	treeTruth->Branch("iProcess"              ,&infoTruth.iProcess         ,"iProcess/I");
}

// destructor
HoverEAnalyzer::~HoverEAnalyzer()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
double HoverEAnalyzer::fillEmIdVarsNoDir(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters,const edm::SortedCollection<HGCRecHit>&  rcs, 	const HGCalGeometry *geom ){

	HGCALShowerBasedEmIdentificationLC2 test(0, geom);

	double see = 0;
	for (unsigned int j =0 ; j < clusters.size() ; j++){

		if (clusters[j]->position()==(sc->seed()->position())) {
			//	std::cout << " HOVERE rcs size " << rcs.size() << std::endl;
			//		std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
			// test.setShowerPosition(clusters[j]->position());
			//test.setShowerDirection(clusters[j]->axis());
			see =  test.HGCALShowerBasedEmIdentificationLC2::sigmaetaeta( *(clusters[j].get()), rcs);
			infoTruth.sigmaEtaEtaNoDir = see;
			infoTruth.sigmaEtaEtaNoDir = test.HGCALShowerBasedEmIdentificationLC2::sigmaetaeta( *(clusters[j].get()), rcs);
			infoTruth.lengthCompatibilityNoDir = test.HGCALShowerBasedEmIdentificationLC2::lengthCompatibility( *(clusters[j].get()), rcs);
			infoTruth.showerStartPosNoDir = fabs(test.HGCALShowerBasedEmIdentificationLC2::startPosition( *(clusters[j].get()), rcs).z());

			//	std::cout << " sieie " <<infoTruth.sigmaEtaEta << ", lC " <<infoTruth.lengthCompatibility  << " ,ssp " <<	infoTruth.showerStartPos << std::endl;

			break;
		}
	}

	return see;
}
double HoverEAnalyzer::fillEmIdVars(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters,const edm::SortedCollection<HGCRecHit>&  rcs, 	const HGCalGeometry *geom ){

	HGCALShowerBasedEmIdentificationLC2 test(PU, geom);
	//test.setShowerPosition(sc->seed()->position());
	//test.setShowerDirection(sc->seed()->axis());

	double see = 0;
	for (unsigned int j =0 ; j < clusters.size() ; j++){

		if (clusters[j]->position()==(sc->seed()->position())) {
			//	std::cout << " HOVERE rcs size " << rcs.size() << std::endl;
			//		std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
			test.setShowerPosition(clusters[j]->position());
			test.setShowerDirection(clusters[j]->axis());
			see =  test.HGCALShowerBasedEmIdentificationLC2::sigmaetaeta( *(clusters[j].get()), rcs);
			info.sigmaEtaEta = see;
			info.sigmaEtaEta = test.HGCALShowerBasedEmIdentificationLC2::sigmaetaeta( *(clusters[j].get()), rcs);
			info.lengthCompatibility = test.HGCALShowerBasedEmIdentificationLC2::lengthCompatibility( *(clusters[j].get()), rcs);
			info.showerStartPos = fabs(test.HGCALShowerBasedEmIdentificationLC2::startPosition( *(clusters[j].get()), rcs).z());
			infoTruth.sigmaEtaEta = test.HGCALShowerBasedEmIdentificationLC2::sigmaetaeta( *(clusters[j].get()), rcs);
			infoTruth.lengthCompatibility = test.HGCALShowerBasedEmIdentificationLC2::lengthCompatibility( *(clusters[j].get()), rcs);
			infoTruth.showerStartPos = fabs(test.HGCALShowerBasedEmIdentificationLC2::startPosition( *(clusters[j].get()), rcs).z());
			//		std::cout << " sieie " <<infoTruth.sigmaEtaEta << ", lC " <<infoTruth.lengthCompatibility  << " ,ssp " <<	infoTruth.showerStartPos << std::endl;

			break;
		}
	}

	return see;
}

float HoverEAnalyzer::resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){

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

float HoverEAnalyzer::clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c,const edm::PtrVector<reco::PFCluster>& clusters){

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
}

	void
HoverEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	if (hcalHelperEndcap_)
	{
		hcalHelperEndcap_->checkSetup(iSetup) ;
		hcalHelperEndcap_->readEvent(iEvent) ;
	}
	// get calo geometry

	edm::ESHandle<HGCalGeometry> geomH;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_,geomH);
	const HGCalGeometry *geom=geomH.product();


	if (caloGeomCacheId_!=iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
		iSetup.get<CaloGeometryRecord>().get(caloGeom_);
		caloGeomCacheId_=iSetup.get<CaloGeometryRecord>().cacheIdentifier();
	}
	if (caloTopoCacheId_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
		caloTopoCacheId_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
		iSetup.get<CaloTopologyRecord>().get(caloTopo_);
	}

	Handle<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >  > HGCEERechits_;
	iEvent.getByToken(endcapRecHitCollection_,HGCEERechits_);
	auto  HGCEERechits =HGCEERechits_.product();
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


	Handle<edm::View<reco::PFRecHit> > eeRecHits;
	iEvent.getByToken(eeRecHitCollection_, eeRecHits);
	const PtrVector<reco::PFRecHit>& rcs = eeRecHits->ptrVector();
	if(rcs.size());


	//	steco_over<< "[debug] number of rechits " << HGCEERechits->size() <<", SCs " << sclusters.size() << ", clusters " << clusters.size() << " gens " << gens.size() <<   std::endl;
	//
	std::cout << "PU " << PU << std::endl;


	int iProcess =process ;

	if(PU >40){
	//	std::cout << "high PU, only look at process " << process << std::endl;
	}
	if (process <2){
	//	std::cout << "[DEBUG] Process " << iProcess;
		infoTruth.iProcess =iProcess;
		info.iProcess =iProcess;

		if (iProcess==0) std::cout << " -> Hgg Photons" << std::endl;
		if (iProcess==1)std::cout << " -> Zee" << std::endl;
		if (iProcess==2 )std::cout << " -> QCD" << std::endl;

		// initialise tree entries
		info.pt=-999.;
		info.eta=-999.;
		info.phi=-999.;
		info.eReco=-999.;
		info.eTrue=-999.;
		info.matchIndex=-999;
		infoTruth.eta=-999.;
		infoTruth.etaSC=-999.;
		infoTruth.etaSeed=-999.;
		infoTruth.phi=-999.;
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
			infoTruth.etaSC=-999.;
			infoTruth.etaSeed=-999.;
			infoTruth.phi=-999.;
			infoTruth.phiSC=-999.;
			infoTruth.phiSeed=-999.;
			infoTruth.eReco_over_eTrue = -999.;
			infoTruth.matchIndex=-999;
			infoTruth.pt=-999.;
			infoTruth.clustersSize=-999;
			infoTruth.nClusters09=-999;
			infoTruth.hoe  =-999.;         
			infoTruth.eReco  =-999.;         
			infoTruth.eTrue         =-999.;  
			infoTruth.eSeed           =-999.;
			infoTruth.eSeed_over_eReco=-999.;
			infoTruth.etaWidth=-999.;
			infoTruth.phiWidth=-999.;

			infoTruth.sigmaEtaEta=-999.;
			infoTruth.showerStartPos=-999.;
			infoTruth.lengthCompatibility=-999.;
			infoTruth.sigmaEtaEtaNoDir=-999.;
			infoTruth.showerStartPosNoDir=-999.;
			infoTruth.lengthCompatibilityNoDir=-999.;

			if (iProcess==0 && (fabs(gens[igp]->pdgId()) != 22 || gens[igp]->status() != 3)) {
				//std::cout <<  "process 0 - pass "<< (iProcess==0 && (fabs(gens[igp]->pdgId()) != 22 || gens[igp]->status() != 3)) << std::endl;
				//std::cout <<  "process 0 - pass "<< (fabs(gens[igp]->pdgId()) != 22 || gens[igp]->status() != 3) << std::endl;
				continue; }
			if (iProcess==1 && (fabs(gens[igp]->pdgId()) != 11 || gens[igp]->status() != 3)) continue;
			if (iProcess==2 && ( !(fabs(gens[igp]->pdgId()) <8 || gens[igp]->pdgId() == 21) || gens[igp]->status() != 3)) continue;

			if (gens[igp]->status() == 3) {
		//		std::cout << "[debug] gen pdgid " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << ", e " << gens[igp]->energy() << ", mother " << std::endl;// (gens[igp]->mother()->pdgId()) <<  std::endl
			};


			assert(gens.size() >0); // only the case for the electron gun sample
			infoTruth.eta=gens[igp]->eta();
			infoTruth.phi      = gens[igp]->phi();

			//		std::cout << "[debug] gen pdgid " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << ", e " << gens[igp]->energy() << ", mother " << std::endl;// (gens[igp]->mother()->pdgId()) <<  std::endl;

			float dRBest =999;
			infoTruth.pt = gens[igp]->pt();
			infoTruth.eTrue = gens[igp]->energy();

			for (unsigned int isc =0; isc < sclusters.size() ; isc++){ //subloop over sc's to find matches

				// calculate dR... dE = dEta, dP = dPhi
				float dE = sclusters[isc]->eta() - gens[igp]->eta();
				dE =dE*dE;
				float dP = sclusters[isc]->phi() - gens[igp]->phi();
				dP =dP*dP;
				float dR = sqrt(dE +dP);

				if (dR < dRlim2 && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
					dRBest = dR;

					// so, if we have a dR match, then store the corresponding match index.
					infoTruth.matchIndex = isc;
				}

			}

			//float recoPt = -999.;
			if (infoTruth.matchIndex >-1) {
				//	recoPt = sclusters[infoTruth.matchIndex]->energy() /std::cosh(sclusters[infoTruth.matchIndex]->eta());
				//infoTruth.eReco = (sclusters[infoTruth.matchIndex]->energy());
				infoTruth.eReco = resumEmEnergy(sclusters[infoTruth.matchIndex], clusters);
				infoTruth.eSeed = (clusterEmEnergy(sclusters[infoTruth.matchIndex]->seed(), clusters));
				infoTruth.eSeed_over_eReco = (infoTruth.eSeed/infoTruth.eReco);
				infoTruth.eReco_over_eTrue = (infoTruth.eReco/gens[igp]->energy());
				infoTruth.clustersSize = (int) sclusters[infoTruth.matchIndex]->clustersSize();
				infoTruth.etaWidth=sclusters[infoTruth.matchIndex]->etaWidth();
				infoTruth.phiWidth=sclusters[infoTruth.matchIndex]->phiWidth();
				infoTruth.etaSC    = sclusters[infoTruth.matchIndex]->eta();
				infoTruth.etaSeed  =sclusters[infoTruth.matchIndex]->seed()->eta();
				infoTruth.phiSC    = sclusters[infoTruth.matchIndex]->phi();
				infoTruth.phiSeed  =sclusters[infoTruth.matchIndex]->seed()->phi();

				int detector = sclusters[infoTruth.matchIndex]->seed()->hitsAndFractions()[0].first.subdetId() ;

				if ( detector==HGCEE) {
					//std::vector<float> locCov = EcalClusterTools::localCovariances( *(sclusters[infoTruth.matchIndex]->seed().get()), &(*eeRecHits), &(*caloTopo_));	
					//	 std::vector<float> locCov = EcalClusterTools::localCovariances( *(sclusters[infoTruth.matchIndex]->seed().get()), (const edm::SortedCollection<EcalRecHit>*) &(*HGCEERechits), &(*caloTopo_));	


					double sieie =  fillEmIdVars(sclusters[infoTruth.matchIndex], clusters,  *HGCEERechits, geom);
					if (sieie);
					// sieie =  fillEmIdVarsNoDir(sclusters[infoTruth.matchIndex], clusters,  *HGCEERechits, geom);	
					//		std::cout << "sigmaIeIE " << sieie  << std::endl ;

					//	std::cout << "test -1" << std::endl;
					const auto & scl = (*HGCEESCs)[infoTruth.matchIndex] ;
					//	std::cout << "test 0" << std::endl;
					//float had =hcalHelperEndcap_->HCALClustersBehindSC((sclusters)[infoTruth.matchIndex]);
					float had =hcalHelperEndcap_->HCALClustersBehindSC(scl);
					//		float had =hcalHelperEndcap_->HCALClustersBehindSC(scl)
					//		std::cout << "test 1" << std::endl;

					float  scle =  sclusters[infoTruth.matchIndex]->energy();
					float hoe =had/scle;
					infoTruth.hoe = hoe;


					for (unsigned int ic =0 ; ic < sclusters[infoTruth.matchIndex]->clusters().size() ; ic++){
						infoTruth.nClusters09 = ic+1; 
						if((sclusters[infoTruth.matchIndex]->clusters())[ic]->energy()/infoTruth.eReco >0.9) break;
					}
				}

				//if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  recoPt >40)
				if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  gens[igp]->pt() >10)
				{
					eRoT_OLD_h->Fill(sclusters[infoTruth.matchIndex]->energy()/gens[igp]->energy());
					eRoT_h->Fill(infoTruth.eReco/gens[igp]->energy());
					phiW_v_eRoT_h->Fill(infoTruth.phiWidth,infoTruth.eReco_over_eTrue);
					dPhi_v_eRoT_h ->Fill(infoTruth.phiSeed - infoTruth.phi,infoTruth.eReco_over_eTrue);
					nClust_v_eRoT_h ->Fill(infoTruth.clustersSize,infoTruth.eReco_over_eTrue);
					nClust90_v_eRoT_h->Fill(infoTruth.nClusters09,infoTruth.eReco_over_eTrue);

					if(iProcess ==0){
						hoe_hgg_h->Fill(infoTruth.hoe);
						see_hgg_h->Fill(infoTruth.sigmaEtaEta);
						ssz_hgg_h->Fill(infoTruth.showerStartPos);
						lcp_hgg_h->Fill(infoTruth.lengthCompatibility);

					}


					if(iProcess ==1){
						hoe_zee_h->Fill(infoTruth.hoe);
						see_zee_h->Fill(infoTruth.sigmaEtaEta);
						ssz_zee_h->Fill(infoTruth.showerStartPos);
						lcp_zee_h->Fill(infoTruth.lengthCompatibility); 			
					}

					if(iProcess ==2){
						hoe_qcd_h->Fill(infoTruth.hoe);
						see_qcd_h->Fill(infoTruth.sigmaEtaEta);
						ssz_qcd_h->Fill(infoTruth.showerStartPos);
						lcp_qcd_h->Fill(infoTruth.lengthCompatibility);

					}
				}


			}

			treeTruth->Fill();
		}
	}

	//--------------> End per-genPhoton tree <------------------

	if ( process ==2){
		//--------------> Begin per-SC tree <---------------------- 

		// loop over superclusters (eg reco particles). 
		int nSCAnalyzed =0;
		for( unsigned int isc =0; isc< sclusters.size() ; isc++){ // isc = index_super_cluster
			
			 int detectorSC = sclusters[isc]->seed()->hitsAndFractions()[0].first.subdetId() ;

				if ( detectorSC!=HGCEE) { continue;}
			nSCAnalyzed++;
			info.eReco_over_eTrue=-999.;
			info.pt=-999.;
			info.hoe=-999.;
			info.eta=-999.;
			info.phi=-999.;
			info.eReco=-999.;
			info.eTrue=-999.;
			info.matchIndex=-999;
			info.clustersSize=-999;

			//	std::cout << " sc " << isc<< " eta " << sclusters[isc]->eta() << ", phi " << sclusters[isc]->phi()<< std::endl;

			info.pt= sclusters[isc]->energy() /std::cosh(sclusters[isc]->eta());
			info.eta=sclusters[isc]->eta();
			info.phi=sclusters[isc]->phi();
			info.eReco = sclusters[isc]->energy();
			info.clustersSize = (int) sclusters[isc]->clustersSize();

			// fill histograms with eta/phi info
			eta_h->Fill(info.eta);
			phi_h->Fill(info.phi);

		//	float dRBest = 999.; // dR best is used to find the gen-reco match with smallest dR.


			int detector = sclusters[isc]->seed()->hitsAndFractions()[0].first.subdetId() ;

			if ( detector==HGCEE) {


				double sieie =  fillEmIdVars(sclusters[isc], clusters,  *HGCEERechits, geom);	
				if (sieie);
				// sieie =  fillEmIdVarsNoDir(sclusters[infoTruth.matchIndex], clusters,  *HGCEERechits, geom);	
				//		std::cout << "test -1"<< sieie << std::endl;
				const auto & scl = (*HGCEESCs)[isc] ;
				//		std::cout << "test 0" << std::endl;
				//	float had =hcalHelperEndcap_->HCALClustersBehindSC((*sclusters)[isc]);
				float had =hcalHelperEndcap_->HCALClustersBehindSC(scl);
				//		float had =hcalHelperEndcap_->HCALClustersBehindSC(scl)
				//		std::cout << "test 1" << std::endl;

				float  scle =  sclusters[isc]->energy();
				float hoe =had/scle;
				info.hoe = hoe;
			}


			if(iProcess ==2){
				hoe_qcd_h->Fill(info.hoe);
				see_qcd_h->Fill(info.sigmaEtaEta);
				ssz_qcd_h->Fill(info.showerStartPos);
				lcp_qcd_h->Fill(info.lengthCompatibility);
			}

			tree->Fill();

		}

		std::cout << " nSCAnalyzed = " << nSCAnalyzed << std::endl;
	}
	return ;
}




// ------------ method called once each job just before starting event loop  ------------
	void 
HoverEAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
HoverEAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
HoverEAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
HoverEAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
HoverEAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
HoverEAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoverEAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoverEAnalyzer);
