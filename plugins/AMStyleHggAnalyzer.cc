// -*- C++ -*-
//
// Package:    AMStyleHggAnalyser
// Class:      AMStyleHggAnalyser
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"

#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "RecoEcal/EgammaClusterAlgos/interface/HGCALShowerBasedEmIdentificationLC2.h"

#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include "PCAShowerAnalysis.h"
using namespace std;


namespace {
template<typename DETID>
std::pair<int,int> getlayer(const unsigned rawid) {
DETID id(rawid);
return std::make_pair(id.zside(),id.layer());
}
}

//
// class declaration
//
// limit dR allowed for geometrical matching of reco/gen particles
float maxDr =1;

// information to be loaded into TTree
struct infoTruth_t {

	float eReco_over_eTrue;
	float eta;
	float etaSC;
	float etaSeed;
	float etaEcal;
	float phi;
	float phiSC;
	float phiSeed;
	float phiEcal;
	int   matchIndex;
	float pt;
	float eTrue;
	float eReco;
	float eRecoCleaned;
	float eReco3x3;
	float eReco3x3Claude;
	float eReco3x3AnneMarie;
	float eReco3x3RecHit;
	float eReco5x5Claude;
	float eReco5x5ClaudeSC;
	float eSeed;
	float eSeed_over_eReco;
	float phiWidth;
	float etaWidth;
	float phiWidthNew;
	float etaWidthNew;
	int clustersSize;
	int nClusters09;
	int converted;
	float x;
	float y;
	float phiCrackDistance;
	float sRR;

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
class AMStyleHggAnalyser : public edm::EDAnalyzer {
	public:
		explicit AMStyleHggAnalyser(const edm::ParameterSet&);
		~AMStyleHggAnalyser();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		double DeltaPhi(const double & phi1, const double & phi2);
		double correctRecHitEnergy(const edm::Ptr<reco::PFRecHit>& rh);
		double correctRecHitEnergy(const edm::Ref<std::vector<reco::PFRecHit> >& rh);
		float resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc, const edm::PtrVector<reco::PFCluster>& clusters);
		float resumEmEnergyTest( const edm::Ptr<reco::SuperCluster>& sc, const edm::PtrVector<reco::PFCluster>& clusters);
		float clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c, const edm::PtrVector<reco::PFCluster>& clusters);
		float claudeEnergy( const edm::Ptr<reco::SuperCluster>& sc, 	std::vector<std::string> geometrySource_, const edm::Event& iEvent, const edm::EventSetup& iSetup, float window);
		float annemarieEnergy( const edm::PtrVector<HGCRecHit>& rechitvec, const edm::Ptr<reco::SuperCluster>& sc, 	std::vector<std::string> geometrySource_, const edm::Event& iEvent, const edm::EventSetup& iSetup, float window);
		float get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rcs, const edm::PtrVector<reco::PFCluster>& clusters,  const edm::Ptr<reco::SuperCluster>& sc, const edm::EventSetup& iSetup);
		float phiCrackDistance (float x, float y);
void getMaximumCell(const edm::PtrVector<HGCRecHit>& rechitvec,const double & phimax,const double & etamax,std::vector<HGCEEDetId> & detidmax);
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

		edm::EDGetTokenT<edm::View<HGCRecHit> >endcapRecHitCollection_      ; 
		edm::EDGetTokenT<edm::View<reco::SuperCluster> >endcapSuperClusterCollection_;
		edm::EDGetTokenT<edm::View<reco::PFCluster> >endcapClusterCollection_     ;
		edm::EDGetTokenT<edm::View<reco::GenParticle> >genParticlesCollection_     ;
		edm::EDGetTokenT<edm::View<int> >genParticlesInts_     ;
		edm::EDGetTokenT<edm::View<reco::PFRecHit> > eeRecHitCollection_     ;
		
		std::string g4TracksSource_, g4VerticesSource_;
		std::vector<std::string> geometrySource_;
  	const HGCalGeometry * hgcEEGeom_;
    
		ROOT::Math::XYZPoint truthVtx_;

		edm::Handle<HGCRecHitCollection> recHits_;	
		TTree *tree;
		TTree *treeTruth;
		info_t info;
		infoTruth_t infoTruth;
		TH1F *rechit_h; 
		TH1F *phi_h; 
		TH1F *dEta_h; 
		TH1F *dPhi_Unconv_h;
		TH1F *dEta_Unconv_h; 
		TH1F *dPhi_h;
		TH2F *phiW_v_eRoT_h; 
		TH1F *nSC_h;
		TH1F *eRoT_OLD_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT3x3_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT3x3Claude_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT3x3ClaudeSC_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT5x5Claude_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT5x5ClaudeSC_h;  //energy true over reco ie etrue/ ereco
		TH2F *dPhi_v_eRoT_h    ; 
		TH2F *nClust_v_eRoT_h  ; 
		TH2F *nClust90_v_eRoT_h;
		const TGraphErrors * _hgcOverburdenParam;
		const TGraphErrors * _hgcLambdaOverburdenParam;
		const std::vector<double> _weights_ee;
  double cellSize_;
  bool doLogWeight_;
  double mipE_;
  unsigned nSR_;
  unsigned debug_;
  bool singleGamma_;



};

// constructor
AMStyleHggAnalyser::AMStyleHggAnalyser(const edm::ParameterSet& iConfig):
	endcapRecHitCollection_(consumes <edm::View<HGCRecHit> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles")))),
	genParticlesInts_(consumes <edm::View<int> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles")))),
	eeRecHitCollection_(consumes <edm::View<reco::PFRecHit> > (iConfig.getUntrackedParameter<edm::InputTag>("eeRecHitCollection",edm::InputTag("particleFlowRecHitHGCEE:Cleaned")))),
	_hgcOverburdenParam(nullptr),
	_hgcLambdaOverburdenParam(nullptr),
	_weights_ee(iConfig.getParameter<std::vector<double> >("weights_ee"))
{	

	//SIM TRACK
  g4TracksSource_           = iConfig.getUntrackedParameter<std::string>("g4TracksSource");
  g4VerticesSource_         = iConfig.getUntrackedParameter<std::string>("g4VerticesSource");
	// SIM TRACK
	geometrySource_ = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
	edm::Service<TFileService> fs_;
	rechit_h         = fs_->make<TH1F>("rechit_h","rechit_h",1000000,0,1);
	phi_h         = fs_->make<TH1F>("phi_h","phi_h",100,-5,5);
	dEta_h        = fs_->make<TH1F>("dEta_h","dEta_h",1000,-1,1);
	dPhi_h        = fs_->make<TH1F>("dPhi_h","dPhi_h",1000,-3.3,3.3);
	dEta_Unconv_h        = fs_->make<TH1F>("dEta_Unconv_h","dEta_Unconv_h",1000,-1,1);
	dPhi_Unconv_h        = fs_->make<TH1F>("dPhi_Unconv_h","dPhi_Unconv_h",1000,-3.3,3.3);
	eRoT_OLD_h        = fs_->make<TH1F>("eRoT_OLD_h","eRoT_OLD_h",1000,-2,2);
	eRoT_h        = fs_->make<TH1F>("eRoT_h","eRoT_h",1000,-2,2);
	eRoT3x3_h        = fs_->make<TH1F>("eRoT3x3_h","eRoT3x3_h",1000,-2,2);
	eRoT3x3Claude_h        = fs_->make<TH1F>("eRoT3x3Claude_h","eRoT3x3Claude_h",1000,-2,2);
	eRoT3x3ClaudeSC_h        = fs_->make<TH1F>("eRoT3x3ClaudeSC_h","eRoT3x3ClaudeSC_h",1000,-2,2);
	eRoT5x5Claude_h        = fs_->make<TH1F>("eRoT5x5Claude_h","eRoT5x5Claude_h",1000,-2,2);
	eRoT5x5ClaudeSC_h        = fs_->make<TH1F>("eRoT5x5ClaudeSC_h","eRoT5x5ClaudeSC_h",1000,-2,2);
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
	treeTruth->Branch("eReco"              ,&infoTruth.eReco            ,"eReco/F");
	treeTruth->Branch("eRecoCleaned"              ,&infoTruth.eRecoCleaned            ,"eRecoCleaned/F");
	treeTruth->Branch("eReco3x3"              ,&infoTruth.eReco3x3            ,"eReco3x3/F");
	treeTruth->Branch("eReco3x3Claude"              ,&infoTruth.eReco3x3Claude            ,"eReco3x3Claude/F");
	treeTruth->Branch("eReco3x3AnneMarie"              ,&infoTruth.eReco3x3AnneMarie            ,"eReco3x3AnneMarie/F");
	treeTruth->Branch("eReco3x3RecHit"              ,&infoTruth.eReco3x3RecHit            ,"eReco3x3RecHit/F");
	treeTruth->Branch("eReco5x5Claude"              ,&infoTruth.eReco5x5Claude            ,"eReco5x5Claude/F");
	treeTruth->Branch("eReco5x5ClaudeSC"              ,&infoTruth.eReco5x5ClaudeSC            ,"eReco5x5ClaudeSC/F");
	treeTruth->Branch("eTrue"              ,&infoTruth.eTrue            ,"eTrue/F");
	treeTruth->Branch("eSeed"              ,&infoTruth.eSeed            ,"eSeed/F");
	treeTruth->Branch("eSeed_over_eReco"   ,&infoTruth.eSeed_over_eReco             ,"eSeed_over_eReco/F");
	treeTruth->Branch("eReco_over_eTrue"              ,&infoTruth.eReco_over_eTrue             ,"eReco_over_eTrue/F");
	treeTruth->Branch("eta"              ,&infoTruth.eta             ,"eta/F");
	treeTruth->Branch("etaEcal"              ,&infoTruth.etaEcal             ,"etaEcal/F");
	treeTruth->Branch("etaSC"              ,&infoTruth.etaSC             ,"etaSC/F");
	treeTruth->Branch("etaSeed"              ,&infoTruth.etaSeed             ,"etaSeed/F");
	treeTruth->Branch("phi"              ,&infoTruth.phi             ,"phi/F");
	treeTruth->Branch("phiSC"              ,&infoTruth.phiSC             ,"phiSC/F");
	treeTruth->Branch("x"              ,&infoTruth.x             ,"x/F");
	treeTruth->Branch("y"              ,&infoTruth.y             ,"y/F");
	treeTruth->Branch("phiSC"              ,&infoTruth.phiSC             ,"phiSC/F");
	treeTruth->Branch("phiEcal"              ,&infoTruth.phiEcal            ,"phiEcal/F");
	treeTruth->Branch("phiSeed"              ,&infoTruth.phiSeed             ,"phiSeed/F");
	treeTruth->Branch("matchIndex"              ,&infoTruth.matchIndex            ,"matchIndex/I");
	treeTruth->Branch("clustersSize"              ,&infoTruth.clustersSize           ,"clustersSize/I");
	treeTruth->Branch("nClusters09"              ,&infoTruth.nClusters09           ,"nClusters09/I");
	treeTruth->Branch("etaWidth"              ,&infoTruth.etaWidth             ,"etaWidth/F");
	treeTruth->Branch("phiWidth"              ,&infoTruth.phiWidth             ,"phiWidth/F");
	treeTruth->Branch("etaWidthNew"              ,&infoTruth.etaWidthNew             ,"etaWidthNew/F");
	treeTruth->Branch("phiWidthNew"              ,&infoTruth.phiWidthNew             ,"phiWidthNew/F");
	treeTruth->Branch("converted"              ,&infoTruth.converted           ,"converted/I");
	treeTruth->Branch("phiCrackDistance"              ,&infoTruth.phiCrackDistance           ,"phiCrackDistance/F");
	treeTruth->Branch("sRR"              ,&infoTruth.sRR       ,"sRR/F");

	if(iConfig.exists("hgcOverburdenParamFile"))
	{
		edm::FileInPath fp = iConfig.getParameter<edm::FileInPath>("hgcOverburdenParamFile");
		TFile *fIn=TFile::Open(fp.fullPath().c_str());
		if(fIn)
		{
			_hgcOverburdenParam=(const TGraphErrors *) fIn->Get("x0Overburden");
			_hgcLambdaOverburdenParam = (const TGraphErrors *) fIn->Get("lambdaOverburden");
			fIn->Close();
		}
	}
	cellSize_ = 1;
	doLogWeight_ = true;
	mipE_ = 0.0000551;
	nSR_ = 5;
	truthVtx_ = ROOT::Math::XYZPoint(0,0,0);
}

// destructor
AMStyleHggAnalyser::~AMStyleHggAnalyser()
{

}

void AMStyleHggAnalyser::getMaximumCell(const edm::PtrVector<HGCRecHit>& rechitvec,const double & phimax,const double & etamax,std::vector<HGCEEDetId> & detidmax){
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
    //double energy = lHit.energy();
    
    ROOT::Math::XYZVector pos(posx-truthVtx_.x(),posy-truthVtx_.y(),posz-truthVtx_.z());
    double deta = fabs(pos.eta()-etamax);
    double dphi = DeltaPhi(pos.phi(),phimax);
    
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    if (dR<dRmin[layer]) {
      dRmin[layer] = dR;
      detidmax[layer] = hgcid;

    }
    
  }//loop on rechits
  
}

//
// member functions
//

// ------------ method called for each event  ------------

void AMStyleHggAnalyser::fillDetIdStack(std::vector<HGCEEDetId> & detidmax_vec,
				   std::vector<HGCEEDetId> & detidstack){
	unsigned int nLayers_=30;
	for( unsigned int iL =0;iL< nLayers_ ; iL++){
	auto detidmax = detidmax_vec[iL];
  const HGCalTopology& topology = hgcEEGeom_->topology();
//	if topology.valid(detidmax){	}
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

float AMStyleHggAnalyser::resumEmEnergyTest( const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){

	float total=0;
	//	float newtotal =0;
	for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
		//	//std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
		for (unsigned int j =0 ; j < clusters.size() ; j++){
		
		double clus_eta = clusters[j]->eta();
		double clus_phi = clusters[j]->phi();
	
    if (fabs(clus_eta - sc->seed()->eta())>0.025) continue;
		if (fabs(deltaPhi(clus_phi , sc->seed()->phi()))>0.11) continue;

			if (clusters[j]->position()==(sc->clusters())[ic]->position()) {
				//		//std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
				total = total +(clusters[j]->emEnergy());
				break;
			}
		}

	}




	return total;
}
float AMStyleHggAnalyser::resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){

	float total=0;
	for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
		//	//std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
		for (unsigned int j =0 ; j < clusters.size() ; j++){

			if (clusters[j]->position()==(sc->clusters())[ic]->position()) {
				//		//std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
				total = total +(clusters[j]->emEnergy());

				std::vector< reco::PFRecHitFraction > recs =  clusters[j]->recHitFractions();

}
				////std::cout<< "SIZE OF RECHITS (alt)" << recs.size() <<", layer " << clusters[j]->layer() <<  std::endl;
		}

	}

	return total;
}
double AMStyleHggAnalyser::DeltaPhi(const double & phi1, const double & phi2){
	double dphi = phi1 - phi2;
	if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
	if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
	return dphi;
}


//float AMStyleHggAnalyser::get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rechits, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc){
double AMStyleHggAnalyser::correctRecHitEnergy(const edm::Ptr<reco::PFRecHit>& rh){

	/*_coef_a(conf.getParameter<double>("effMip_to_InverseGeV_a")),
		_coef_b(conf.getParameter<double>("effMip_to_InverseGeV_b")),
		_coef_c(conf.getParameter<double>("effMip_to_InverseGeV_c")),
		_coef_d(conf.exists("effMip_to_InverseGeV_d") ? conf.getParameter<double>("effMip_to_InverseGeV_d") : 0),
		_coef_e(conf.exists("effMip_to_InverseGeV_e") ? conf.getParameter<double>("effMip_to_InverseGeV_e") : 1.0),*/

	double _coef_a =80.0837; 
	//	double _coef_b =-107.229; 
	double _coef_c =0.0472817; 
	double _coef_d =-0.266294; 
	double _coef_e =0.34684;
	double _mipValueInGeV_ee = 55.1*1e-6;

	bool _isEMCalibration = 1;

	double eCorr = 0.0;
	//	double e_ee(0.0);
	const double rh_eta = rh->positionREP().eta();
	const double abs_eta = std::abs(rh_eta);
	const double effMIP_to_InvGeV = _coef_a*fabs(std::tanh(rh_eta))/(1.-(_coef_c*pow(rh_eta,2)+_coef_d*abs_eta+_coef_e));
	if (effMIP_to_InvGeV) ;

	const std::vector<double>* weights = nullptr;
	weights = &_weights_ee;

	//DetId theid = rh->detId();
	std::pair<int,int> zside_layer ;
	double mip_value = 0.0;
	double mip2gev = 0.0;
	double hgcOverburdenWeight = 0.0;

	zside_layer = getlayer<HGCEEDetId>(rh->detId());
	mip_value = _mipValueInGeV_ee;
	mip2gev = effMIP_to_InvGeV;

	if(_isEMCalibration) {
		if(zside_layer.second==1 && _hgcOverburdenParam) {
			hgcOverburdenWeight = _hgcOverburdenParam->Eval(abs_eta);
		}
	} /*else {
			if(zside_layer.second==1 && _hgcOverburdenParam) {
			hgcOverburdenWeight = _hgcLambdaOverburdenParam->Eval(abs_eta);
			}
			} */
	//	e_ee += ((*weights)[zside_layer.second-1] + hgcOverburdenWeight)*rh->energy()/(mip_value*std::tanh(abs_eta));

	const int layer = zside_layer.second;
	const double energy_MIP = rh->energy()/mip_value;
	const double added_energy = ((*weights)[layer-1]+hgcOverburdenWeight)*energy_MIP/mip2gev;
	eCorr += added_energy;

	return eCorr;
}
//double AMStyleHggAnalyser::correctRecHitEnergy(const edm::Ptr<reco::PFRecHit>& rh){
double AMStyleHggAnalyser::correctRecHitEnergy(const edm::Ref<std::vector<reco::PFRecHit> >& rh){

	/*_coef_a(conf.getParameter<double>("effMip_to_InverseGeV_a")),
		_coef_b(conf.getParameter<double>("effMip_to_InverseGeV_b")),
		_coef_c(conf.getParameter<double>("effMip_to_InverseGeV_c")),
		_coef_d(conf.exists("effMip_to_InverseGeV_d") ? conf.getParameter<double>("effMip_to_InverseGeV_d") : 0),
		_coef_e(conf.exists("effMip_to_InverseGeV_e") ? conf.getParameter<double>("effMip_to_InverseGeV_e") : 1.0),*/

	double _coef_a =80.0837; 
	//	double _coef_b =-107.229; 
	double _coef_c =0.0472817; 
	double _coef_d =-0.266294; 
	double _coef_e =0.34684;
	double _mipValueInGeV_ee = 55.1*1e-6;

	bool _isEMCalibration = 1;

	double eCorr = 0.0;
	//	double e_ee(0.0);
	const double rh_eta = rh->positionREP().eta();
	const double abs_eta = std::abs(rh_eta);
	const double effMIP_to_InvGeV = _coef_a*fabs(std::tanh(rh_eta))/(1.-(_coef_c*pow(rh_eta,2)+_coef_d*abs_eta+_coef_e));
	if (effMIP_to_InvGeV) ;

	const std::vector<double>* weights = nullptr;
	weights = &_weights_ee;

	//DetId theid = rh->detId();
	std::pair<int,int> zside_layer ;
	double mip_value = 0.0;
	double mip2gev = 0.0;
	double hgcOverburdenWeight = 0.0;

	zside_layer = getlayer<HGCEEDetId>(rh->detId());
	mip_value = _mipValueInGeV_ee;
	mip2gev = effMIP_to_InvGeV;

	if(_isEMCalibration) {
		if(zside_layer.second==1 && _hgcOverburdenParam) {
			hgcOverburdenWeight = _hgcOverburdenParam->Eval(abs_eta);
		}
	} /*else {
			if(zside_layer.second==1 && _hgcOverburdenParam) {
			hgcOverburdenWeight = _hgcLambdaOverburdenParam->Eval(abs_eta);
			}
			} */
	//	e_ee += ((*weights)[zside_layer.second-1] + hgcOverburdenWeight)*rh->energy()/(mip_value*std::tanh(abs_eta));

	const int layer = zside_layer.second;
	const double energy_MIP = rh->energy()/mip_value;
	const double added_energy = ((*weights)[layer-1]+hgcOverburdenWeight)*energy_MIP/mip2gev;
	eCorr += added_energy;

	return eCorr;
}


//float AMStyleHggAnalyser::get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rechits, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc){
float AMStyleHggAnalyser::annemarieEnergy( const edm::PtrVector<HGCRecHit>& rechitvec, const edm::Ptr<reco::SuperCluster>& sc, 	std::vector<std::string> geometrySource_, const edm::Event& iEvent, const edm::EventSetup& iSetup, float eta_ECAL){
	//////////////// TEST//////////////
	
//	std::cout << "[debug anne marie] 0 " << std::endl;
	std::map<int,const HGCalGeometry *> hgcGeometries;
	edm::ESHandle<HGCalGeometry> hgcGeo;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
	hgcGeometries[0]=hgcGeo.product();

  int nLayers_=30;
	std::vector<HGCEEDetId> detidmax;
  detidmax.resize(nLayers_,HGCEEDetId());

  std::vector<HGCEEDetId> detidstack;
  detidstack.resize(9*nLayers_,HGCEEDetId());

const double & phimax = infoTruth.phi;
const double & etamax = infoTruth.eta;//sc;
double clus_eta = sc->seed()->eta();

  getMaximumCell(rechitvec,phimax,etamax,detidmax);
	
	
	fillDetIdStack(detidmax,detidstack);

	HGCALShowerBasedEmIdentificationLC2 test(1, hgcEEGeom_);
	//std::cout << " -- True eta-phi " << etamax << " " << phimax << std::endl    ;

  const HGCalTopology& topology = hgcEEGeom_->topology();
	double sRR =  test.HGCALShowerBasedEmIdentificationLC2::sigmaRR(detidstack , recHits_);
	infoTruth.sRR =sRR;
	
	std::cout << " DEBUG SRR " << sRR << std::endl;

  //float maxE=0;
	float energyReco=0;
	float energyRH=0;
		
		const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;

		double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
		double mip = 0.0000551;
		double weight[30] = {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};
		for (int iLayer =0; iLayer < 30 ; iLayer++){
			//			std::cout << "LAYER "<< iLayer <<", clu "<< clu << ", dRmin "  << std::endl;
			//		std::cout << "POS " << pcaShowerPos << ", AXIS " << pcaShowerDir << std::endl;
    std::vector<double> Exy;
    Exy.resize(9,0);
//	std::cout << "[debug anne marie] 2 " << detidmax.size()  << std::endl;
    if (topology.valid(detidmax[iLayer])) {//fillNeighbours(detidmax[iL],Exy,info);
//	std::cout << "[debug anne marie] 2 .1" << std::endl;
auto detidmaxL = detidmax[iLayer]; 
  std::vector<HGCEEDetId> neighbours;
  neighbours.resize(9,HGCEEDetId());
//	std::cout << "[debug anne marie] 2 .2" << std::endl;
  DetId tmp = topology.goSouth(detidmaxL);
  if (topology.valid(tmp)) {
    neighbours[1] = HGCEEDetId(tmp);
//	std::cout << "[debug anne marie] 2 .3" << std::endl;
    tmp = topology.goWest(neighbours[1]);
    if (topology.valid(tmp)) neighbours[0] = HGCEEDetId(tmp);
    tmp = topology.goEast(neighbours[1]);
//	std::cout << "[debug anne marie] 2 .4" << std::endl;
    if (topology.valid(tmp)) neighbours[2] = HGCEEDetId(tmp);
  }

//	std::cout << "[debug anne marie] 3.1 " << std::endl;
  neighbours[4] = detidmaxL;
  tmp = topology.goWest(neighbours[4]);
  if (topology.valid(tmp)) neighbours[3] = HGCEEDetId(tmp);
///	std::cout << "[debug anne marie] 3 .2" << std::endl;
  tmp = topology.goEast(neighbours[4]);
  if (topology.valid(tmp)) neighbours[5] = HGCEEDetId(tmp);
//	std::cout << "[debug anne marie] 3 .3" << std::endl;

  tmp = topology.goNorth(detidmaxL);
  if (topology.valid(tmp)) {
//	std::cout << "[debug anne marie] 3 .4" << std::endl;
    neighbours[7] = HGCEEDetId(tmp);
    tmp = topology.goWest(neighbours[7]);
    if (topology.valid(tmp)) neighbours[6] = HGCEEDetId(tmp);
    tmp = topology.goEast(neighbours[7]);
//	std::cout << "[debug anne marie] 3 .5" << std::endl;
    if (topology.valid(tmp)) neighbours[8] = HGCEEDetId(tmp);
  }
//	std::cout << "[debug anne marie] 4 " << std::endl;
  for (unsigned idx(0);idx<9;++idx){
    if (!topology.valid(neighbours[idx])) {
     // info.noPhiCrack = false;
      continue;
    }
    HGCRecHitCollection::const_iterator theHit = recHits_->find(neighbours[idx]);
    if (theHit==recHits_->end()) continue;
    if (neighbours[idx].det()!=DetId::Forward || neighbours[idx].subdetId()!=HGCEE) continue;
    GlobalPoint cellPos = hgcEEGeom_->getPosition(neighbours[idx]);
   // if ((neighbours[idx].sector()*neighbours[idx].subsector()) != (detidmaxL.sector()*detidmaxL.subsector())) info.noPhiCrack = false;
    double posx = cellPos.x();
    double posy = cellPos.y();
    double posz = cellPos.z();
    double costheta = fabs(posz)/sqrt(posz*posz+posx*posx+posy*posy);
    double energy = theHit->energy()/mipE_*costheta;//in MIP
    Exy[idx] = energy;
		double scale = mip*corr;
		if (scale && weight[0]);
		energyReco += energy*absWeight(iLayer);
		energyRH += theHit->energy();
		rechit_h->Fill(theHit->energy());

	//	std::cout <<" layer " << iLayer<< ", neigh " << idx <<" energyReco " <<  energyReco<< ", e add " <<energy*absWeight(iLayer) << std::endl;

  }

		} else {
continue;
		}
}
//std::cout << "[debug anne marie] 5 " << energyReco << " tanh thing " << fabs(tanh(eta_ECAL)) << "eta ECAL " << eta_ECAL << std::endl;
  energyReco=energyReco/fabs(tanh(eta_ECAL));
//std::cout << "[debug anne marie] 5.5 " << energyReco <<  std::endl;
  //calibration for signal region 2: 3*3 cm^2
	  double pars[3] = {69.5,4.5,-0.8};
	  double paro[3] = {-50,0,0};  
		 double offset = paro[0] + paro[1]*fabs(eta_ECAL) + paro[2]*eta_ECAL*eta_ECAL;
		 double slope = pars[0] + pars[1]*fabs(eta_ECAL) + pars[2]*eta_ECAL*eta_ECAL;
		 energyReco= (energyReco-offset)/slope;
//std::cout << "[debug anne marie] 6 " << energyReco <<  std::endl;
	infoTruth.eReco3x3RecHit = energyRH;
	return energyReco;

}
float AMStyleHggAnalyser::claudeEnergy( const edm::Ptr<reco::SuperCluster>& sc, 	std::vector<std::string> geometrySource_, const edm::Event& iEvent, const edm::EventSetup& iSetup, float window){

	//////////////// TEST//////////////

	std::map<int,const HGCalGeometry *> hgcGeometries;
	edm::ESHandle<HGCalGeometry> hgcGeo;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
	hgcGeometries[0]=hgcGeo.product();

	float energy=0;
	float etaWidthNum= 0;
	float phiWidthNum= 0;
	assert (infoTruth.eReco>-999);
	float denom = infoTruth.eReco;
	PCAShowerAnalysis pcaShowerAnalysis(iEvent,iSetup);
	//auto sc = sclusters[infoTruth.matchIndex]	;
	int clu =0; 
	for (auto itcl=(sc)->clustersBegin(); itcl!=(sc)->clustersEnd(); itcl++) {
		GlobalPoint pcaShowerPos;
		GlobalVector pcaShowerDir;
		pcaShowerAnalysis.showerParameters(&(**itcl),pcaShowerPos,pcaShowerDir);

		const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
		double clus_eta = (*itcl)->eta();
		double clus_phi = (*itcl)->phi();
		
		if (abs(clus_eta - sc->seed()->eta())>0.025) continue;
		if (abs(deltaPhi(clus_phi , sc->seed()->phi()))>0.11) continue;
		

		double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
		double mip = 0.0000551;
		double weight[30] = {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};
		for (int iLayer =0; iLayer < 30 ; iLayer++){
			//			std::cout << "LAYER "<< iLayer <<", clu "<< clu << ", dRmin "  << std::endl;
			//		std::cout << "POS " << pcaShowerPos << ", AXIS " << pcaShowerDir << std::endl;
			float dRmin = 9999;
			float dRmax = 0.5;
			if (dRmax ); // stop compilation errros for unsed variables.

			float bestX;
			float bestY;
			float bestZ;

			float baryX = 0;
			float baryY = 0;
			float baryZ = 0;
			float totWeight =0;

			//find position of cell nearest to axis in that layer.		
			for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
				const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
				HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);
				if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
					const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
					GlobalPoint cellPos = hgcGeometries[0]->getPosition(hgcid_);
					const int layer = hgcid_.layer();
					if (layer != iLayer  ) continue;
					//				std::cout << "available rechits "  << " X " << cellPos.x() << ", Y " << cellPos.y() << ", Z " << cellPos.z() << std::endl; 
					GlobalPoint intercept(pcaShowerPos);
					intercept += ((cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z())*pcaShowerDir;
					//			std::cout << "INTERCEPT " << intercept << std::endl; 
					float dEta = (intercept.eta() - cellPos.eta() );
					float dPhi = DeltaPhi(intercept.phi(), cellPos.phi());
					float dR = sqrt (dPhi*dPhi + dEta*dEta);
					if ((dR < dRmin))
					{
						//	if (dR> dRmax) continue;
						dRmin = dR;
						bestX = cellPos.x();
						bestY = cellPos.y();
						bestZ = cellPos.z();
					}
					// recompute calibrated SC energy
					//	std::cout << "n Neih " << theHit->neighbours8().size() << std::endl;
					double scale = mip*corr;
					float e = theHit->energy()*weight[layer-1]/scale;
					totWeight += e;
					baryX += e*cellPos.x();
					baryY += e*cellPos.y();
					baryZ += e*cellPos.z();

				}
			}
			//		std::cout << "nearest rechit " << dRmin << " best X " << bestX << ", bestY " << bestY << ", best Z " << bestZ << std::endl; 
			baryX = baryX/totWeight;	
			baryY = baryY/totWeight;	
			baryZ = baryZ/totWeight;	
			//	std::cout << "barycenter " << baryX << ", " << baryY << ", " << baryZ << std::endl;
			//GlobalPoint bestCellPos(bestX, bestY, bestZ);
			if( bestX && bestY && bestZ) ; //avoid compialtion err for unused var.
			GlobalPoint bestCellPos(baryX, baryY, baryZ);
			//	GlobalPoint baryPos(baryX, baryY, baryZ);

			//		std::cout << "LAYER "<< iLayer <<", clu "<< clu << std::endl;
			//Now sum over the nearest cells around that best cell...
			int nRHadded =0;
			for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
				const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
				HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);
				if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
					const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
					GlobalPoint cellPos = hgcGeometries[0]->getPosition(hgcid_);
					const int layer = hgcid_.layer();
					if (layer != iLayer  ) continue;
					//	GlobalPoint intercept(pcaShowerPos);
					//	intercept += (-1*(cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z())*pcaShowerDir;
					float dx = fabs(bestCellPos.x() - cellPos.x() );
					float dy = fabs(bestCellPos.y() - cellPos.y() );
					//	float dR = sqrt (dx*dx + dy*dy);
					//std::cout << "dx "<< dx<< ", dy " << dy << ", dr " << dR  << "true? "<< (dR < (1.2*sqrt(2))) << std::endl;
					//	if ((dR < (1.2*sqrt(2))))
					if (dx <window && dy <window)
					{
						// recompute calibrated SC energy
						double scale = mip*corr;
						//	std::cout << "n Neih " << theHit->neighbours8().size() << std::endl;
						energy += theHit->energy()*weight[layer-1]/scale;
						float dEta = cellPos.eta() - sc->seed()->eta();
						float dPhi = DeltaPhi(cellPos.phi() , sc->seed()->phi());

						etaWidthNum +=( theHit->energy()*weight[layer-1]/scale)*dEta*dEta;
						phiWidthNum +=( theHit->energy()*weight[layer-1]/scale)*dPhi*dPhi;


						nRHadded++;
						//		std::cout << "add rechit "<< nRHadded << "  with dx " << dx <<", dy " << dy << ", energy " << theHit->energy()*weight[layer-1]/scale << ", tot energy " << energy << std::endl; 
					}
				}
			}

		}
		clu++;
	}
	//////////////// TEST//////////////

	infoTruth.etaWidthNew = sqrt(etaWidthNum/ denom);
	infoTruth.phiWidthNew = sqrt(phiWidthNum/ denom);
	return energy;
}
float AMStyleHggAnalyser::get3x3EmEnergy(const  edm::PtrVector<reco::PFRecHit>& rcs, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc,const edm::EventSetup& iSetup ){


	double _coef_a =80.0837; 
	double _coef_b =-107.229; 

	float totalE =0;
	float newtotalE =0;
	float eta = sc->seed()->eta();
	float phi = sc->seed()->phi();
	//	float eta = infoTruth.eta;
	//	float phi = infoTruth.phi;
	float debug_max =0;
	//	float dRmin =9999.;
	int match =-1;
	DetId  id_; 

	for (unsigned int iLayer =0; iLayer<30 ; iLayer++){  // iLayer = iLayer

		float dRmin =9999.;
		match =-1;

		//std::cout << " DEBUG DEBUG nRH " << rcs.size() << std::endl;
		for (unsigned int k=0; k < rcs.size() ; k++) {
			int hgc_layer = HGCEEDetId(rcs[k]->detId()).layer();

			if(hgc_layer != (int) iLayer) continue;
			//		//std::cout << "TEST, rechits rho,eta, phi " << rcs[k]->positionREP() << std::endl;
			float rcEta = rcs[k]->positionREP().eta(); 
			float rcPhi = rcs[k]->positionREP().phi(); 
			//	refPos= geom->getPosition(hit_it->id()) ;
			float dEta =0;
			float dPhi =0;

			dEta = rcEta-eta;
			dPhi = DeltaPhi(rcPhi ,phi);
			if(rcs[k]->energy() > debug_max){ debug_max =rcs[k]->energy();}

			float dR = sqrt (dEta*dEta + dPhi*dPhi);

			if(dR < dRmin  ){
				dRmin = dR;
				match =k;
				// id_ = hit_it->id();
				//§	//std::cout << "DEBUG dRmin " << dRmin << ", match " << match << std::endl; 
			}
			//k++;
		}


		//std::cout << "[info] eta,phi match for layer "<< iLayer <<" has index " << match << ", DRMIN " << dRmin <<  std::endl; 
		//std::cout << "[DEBUG]i max rechit index "<< match<< std::endl; 
		if (match==-1) continue;	
		//std::cout << ", and "<< rcs[match]->neighbours8().size() << "neighbours "<<  std::endl; 

		float maxE =rcs[match]->energy();
		//	//std::cout << "DEBUG central: " <<  maxE << std::endl;
		int maxEIndex =-1;
		for( unsigned int iNeighbour=0; iNeighbour< rcs[match]->neighbours8().size(); iNeighbour++ ){


			//	float e =1;
			//	float e = (rcs[match]->neighbours8()[iNeighbour])->energy();
			float e2 = correctRecHitEnergy((rcs[match]->neighbours8()[iNeighbour]));
			//std::cout << "DEBUG " <<  e2 << ", nNeigh " << (rcs[match]->neighbours8()[iNeighbour])->neighbours8().size()<< std::endl;
			if (e2 > maxE) {

				maxE=e2;
				maxEIndex=iNeighbour;
				//eta =(rcs[match]->neighbours8()[iNeighbour])->positionREP().eta();
				//phi =(rcs[match]->neighbours8()[iNeighbour])->positionREP().phi();

			}

		}
		//std::cout << "[info] maxE has index " << maxEIndex << " and energy " << maxE << std::endl; 
		//	//std::cout << "[info]  **CORR* maxE has index " << maxEIndex << " and energy " << maxE << std::endl; 
		if (maxEIndex>-1){

			totalE = totalE + rcs[match]->neighbours8()[maxEIndex]->energy();
			newtotalE = newtotalE + correctRecHitEnergy((rcs[match]->neighbours8()[maxEIndex]));
			//		//std::cout << "[INFO]  best match is around	"<< std::endl;
			//		//std::cout << "[INFO] 	**CORR** 	Addding Rechit 0 " << "gives "<<correctRecHitEnergy((rcs[match]->neighbours8()[maxEIndex])) << std::endl;

			for( unsigned int iNeighbour2=0; iNeighbour2< ((rcs[match]->neighbours8()[maxEIndex])->neighbours8()).size(); iNeighbour2++ ){


				totalE = totalE + ((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2]->energy();
				newtotalE = newtotalE + correctRecHitEnergy((((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2] ));
				//		//std::cout << "[INFO] 		Addding Rechit  "<< iNeighbour2 << "gives "<< totalE << std::endl;
				//std::cout << "[INFO] **CORR** 	Addding Rechit  "<< iNeighbour2 << "gives "<< correctRecHitEnergy((((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2] )) << ", w nerighbours " <<  (((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2] )->neighbours8().size()  << std::endl ;


			}



		}  else {

			//	//std::cout << "[INFO]  best match is central	"<< std::endl;
			totalE = totalE + rcs[match]->energy();
			newtotalE = newtotalE + correctRecHitEnergy( rcs[match]);
			//	//std::cout << "[INFO] 		Addding Rechit 0 " << "gives "<< totalE << std::endl;
			//	//std::cout << "[INFO] **CORR** 	Addding Rechit 0 " << "gives "<<correctRecHitEnergy( rcs[match])
			//		<< std::endl;

			for( unsigned int iNeighbour2=0; iNeighbour2< (rcs[match]->neighbours8()).size(); iNeighbour2++ ){


				totalE = totalE + (rcs[match]->neighbours8())[iNeighbour2]->energy();
				newtotalE = newtotalE + correctRecHitEnergy(((rcs[match]->neighbours8())[iNeighbour2]));

				//		//std::cout << "[INFO] 		Addding Rechit  "<< iNeighbour2 << "gives "<< totalE << std::endl;
				//		//std::cout << "[INFO] **CORR**	Addding Rechit  "<< iNeighbour2 << "gives "<< correctRecHitEnergy(((rcs[match]->neighbours8())[iNeighbour2]))
				//		<< std::endl;
			}

		}

		//	//std::cout << "[INFO] Layer "<< iLayer  << " subtotal "<< totalE << std::endl;
		//	//std::cout << "[INFO] **CORR** Layer "<< iLayer  << " subtotal "<< newtotalE << std::endl;

	}
	//decode position
	//uint32_t recoDetId(recHits->find(id_)->id());

	//const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	//int layer( ((rcs[match]->detId() >> 19) & 0x1f) );// + layerCtrOffset-1 );

	//int hgc_layer = HGCEEDetId(rcs[match]->detId()).layer();
	//refPos= geom->getPosition(id_) ;
	totalE = std::max(0., totalE-_coef_b/_coef_a);
	newtotalE = std::max(0., newtotalE-_coef_b/_coef_a);

	////std::cout << "NEIGHOURS " << recHits->find(id_)->neighbours4()->size()<< std::endl;
	if(match>-1){
		//		//std::cout << "SC E " << sc->energy()  << ", calculated E " << totalE  <<std::endl;
		//		//std::cout << "**CORR** SC E " << sc->energy()  << ", calculated E " << newtotalE  <<std::endl;
		//	//std::cout << "match rechit eta, phi, layer " << rcs[match]->positionREP().eta() << ", " <<rcs[match]->positionREP().eta()  << ", " << layer <<  std::endl;
		//	//std::cout << "match rechit eta, phi, layer (alt) " << rcs[match]->positionREP().eta() << ", " <<rcs[match]->positionREP().eta()  << ", " << hgc_layer <<  std::endl;
	}
	//	totalE = std::max(0., totalE-_coef_b/_coef_a);
	//return (totalE);				
	return (newtotalE);				
}


float AMStyleHggAnalyser::phiCrackDistance (float x, float y){
	float d=999;
	float dmin=999;
	//there are six lines to consider, which go through the origin to the point (cos(N*0.349+0.1745), sin(N*0.349+0.1745));
	//0.349 is 20 deg in rad, 0.1745 is teh 10 deg offset
	// these are the lines y = tan(N*pi/6)x., or tan(N*pi/6)* -y =0
	// Distance( ax+by+c=0, (x0,y0)) = fabs(a0 + by0 +c)/sqrt(a*a+b*b)

	for (int N =0; N< 9; N++){

		float a=tan(N*0.349+0.1745);
		//d = fabs( a*x -y)/sqrt(a*a + 1);
		d = fabs( a*x -y)/sqrt(a*a + 1);

		if (d<dmin) {
		
			dmin=d;

		}

//std::cout << "debug N" << N <<", x" << x << ", y " << y << ", a " << a<< ", d " << d << ", dmin " << dmin << std::endl;
	}
	
	return dmin;
}




float AMStyleHggAnalyser::clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c,const edm::PtrVector<reco::PFCluster>& clusters){

	float emEnergy=0;
	//for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
	//	//std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
	for (unsigned int j =0 ; j < clusters.size() ; j++){

		if (clusters[j]->position()==(c->position())) {
			//		//std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
			emEnergy=(clusters[j]->emEnergy());
			break;
		}
	}

	return emEnergy;
}

	void
AMStyleHggAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	using namespace edm;
	iEvent.getByLabel(edm::InputTag("HGCalRecHit:HGCEERecHits"),recHits_);

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

	//SIM TRACK
	//generator level particles
	edm::Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByLabel("genParticles", genParticles);
	//	auto gens2 = genParticles.product();
	//	size_t maxGenParts(genParticles->size());
	//	if(maxGenParts>2) maxGenParts=2;
	//	if(genParticles->size()>maxGenParts) std::cout << "[Warning] found more than " << maxGenParts << " gen particles, will save only first " << maxGenParts << std::endl;


	//Geant4 collections
	edm::Handle<edm::SimTrackContainer> SimTk;
	iEvent.getByLabel(g4TracksSource_,SimTk);
	edm::Handle<edm::SimVertexContainer> SimVtx;
	iEvent.getByLabel(g4VerticesSource_,SimVtx); 
	edm::Handle<edm::View<int> > genBarcodes;
	//	iEvent.getByLabel("genParticlesTag",genBarcodes);
	iEvent.getByToken(genParticlesInts_,genBarcodes);
	auto gBarcodes =  genBarcodes->ptrVector();
	//SIM TRACK




	nSC_h->Fill(sclusters.size());
	int count =0;

	//std::cout << "++++++++++++++++ " << gens.size() <<", " << gBarcodes.size() << std::endl;

	for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree

		infoTruth.eta=-999.;
		infoTruth.etaSC=-999.;
		infoTruth.etaSeed=-999.;
		infoTruth.phi=-999.;
		infoTruth.y=-999.;
		infoTruth.x=-999.;
		infoTruth.phiSC=-999.;
		infoTruth.phiSeed=-999.;
		infoTruth.eReco_over_eTrue = -999.;
		infoTruth.matchIndex=-999;
		infoTruth.pt=-999.;
		infoTruth.clustersSize=-999;
		infoTruth.nClusters09=-999;
		infoTruth.eReco  =-999.;         
		infoTruth.eRecoCleaned  =-999.;         
		infoTruth.eReco3x3Claude  =-999.;         
		infoTruth.eReco3x3AnneMarie  =-999.;         
		infoTruth.eReco5x5Claude  =-999.;         
		infoTruth.eReco5x5ClaudeSC  =-999.;         
		infoTruth.eTrue         =-999.;  
		infoTruth.eSeed           =-999.;
		infoTruth.eSeed_over_eReco=-999.;
		infoTruth.etaWidth=-999.;
		infoTruth.phiWidth=-999.;
		infoTruth.etaWidthNew=-999.;
		infoTruth.phiWidthNew=-999.;
		infoTruth.converted =-999.; //SIM TRACK
		infoTruth.phiCrackDistance =-999.; //SIM TRACK
	// initialise tree entries
	info.pt=-999.;
	info.eta=-999.;
	info.phi=-999.;
	info.eReco=-999.;
	////info.eRecoCleaned=-999.;
	info.eTrue=-999.;
	info.matchIndex=-999;
	infoTruth.eta=-999.;
	infoTruth.etaSC=-999.;
	infoTruth.etaSeed=-999.;
	infoTruth.phi=-999.;
	infoTruth.x=-999.;
	infoTruth.y=-999.;
	infoTruth.phiSC=-999.;
	infoTruth.phiSeed=-999.;
	infoTruth.etaWidth=-999.;
	infoTruth.phiWidth=-999.;
	infoTruth.etaWidthNew=-999.;
	infoTruth.phiWidthNew=-999.;
	infoTruth.matchIndex=-999;
	infoTruth.clustersSize=-999;
	infoTruth.nClusters09=-999;
	infoTruth.eReco  =-999.;         
	infoTruth.eRecoCleaned  =-999.;         
	infoTruth.eTrue         =-999.;  
	infoTruth.eSeed           =-999.;
	infoTruth.eSeed_over_eReco=-999.;
	infoTruth.sRR=-999.;

		//std::cout << "debug " << genBarcodes->at(igp) << std::endl; 

	//	if (gens[igp]->pdgId() != 22 ) continue;
	//	if( gens[igp]->status() != 3) continue;
		//if( gens[igp]->pt() <10.) continue;


		
		
    //if ((gens.size()>2 && (gens[igp]->pdgId() != 22 || gens[igp]->status() != 44)) || (gens.size()==2 && (gens[igp]->pdgId() != 22 || gens[igp]->status() != 1)) ) continue;
		//if (igp >4 ){
		//if( gens[igp]->mother()->status() != 23) continue;
		//}
		if (gens[igp]->pdgId() != 22) continue;
		if (gens[igp]->status() != 1) continue;
		if (gens[igp]->mother()->status() != 44) continue;
		std::cout << 	"Gen particle " << igp << ", status "<< gens[igp]->status() << ", id " <<  gens[igp]->pdgId() << ", eta " << gens[igp]->eta() << ", phi " << gens[igp]->phi() <<" vertex " << gens[igp]->vertex()<< std::endl;
	if(igp>3)	std::cout << " -- mother "<< ", status "<< gens[igp]->mother()->status() << ", id " <<  gens[igp]->mother()->pdgId() << ", eta " << gens[igp]->mother()->eta() << ", phi " << gens[igp]->mother()->phi() <<" vertex " << gens[igp]->mother()->vertex()<< std::endl;
		assert(gens.size() >0); // only the case for the electron gun sample
		infoTruth.eta      = gens[igp]->eta();
		infoTruth.phi      = gens[igp]->phi();
		truthVtx_=gens[igp]->vertex();

		size_t nHitsBeforeHGC(0);
		//math::XYZVectorD hitPos=getInteractionPositionLC(SimTk.product(),SimVtx.product(),*(gens[igp].get())).pos;
		math::XYZVectorD hitPos=getInteractionPositionLC(SimTk.product(),SimVtx.product(), gens[igp]->pt()).pos;//,*(gens[igp].get())).pos;
		const double z = std::abs(hitPos.z());
		//double z=0;
		nHitsBeforeHGC += (unsigned)(z < 317 && z > 1e-3);



		//	std::cout << "[debug] gen "<< count << "  pdgid " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << ", z " << z << ", converted " << nHitsBeforeHGC << ", pt " << gens[igp]->pt() << ", barcode " << *((gBarcodes[igp]).get()) <<  std::endl;// (gens[igp]->mother()->pdgId()) <<  std::endl;

		count++;
		infoTruth.converted = nHitsBeforeHGC;




		float dRBest =999;
		infoTruth.pt = gens[igp]->pt();
		infoTruth.eTrue = gens[igp]->energy();
		std::cout << " debug " << infoTruth.eta <<", "<< infoTruth.phi << std::endl;

		float zHGC =320.38;
		float rHGC = 32.0228;
		float RHGC = 152.153;
		if (infoTruth.eta <0) zHGC = -1*zHGC;
		math::XYZPoint vPos = gens[igp]->vertex();

		float theta = 2*std::atan(std::exp(-infoTruth.eta));
		float h = std::tan(theta)*(zHGC-vPos.z());
		if (h>RHGC || h<rHGC);// {
//		std::cout << "debug skip because h>RHGC || h<rHGC gives 1" << std::endl; 
//		continue;
//		}
		float eta_ECAL = -std::log(std::tan(std::asin(h/std::sqrt(h*h+zHGC*zHGC))/2)) ;
		if  (infoTruth.eta <0) eta_ECAL = -1*eta_ECAL;

		infoTruth.etaEcal =eta_ECAL;
		infoTruth.phiEcal = infoTruth.phi;
		infoTruth.x = (zHGC/std::cos(theta))* std::sin(theta) *std::cos(infoTruth.phi);
		infoTruth.y = (zHGC/std::cos(theta))* std::sin(theta) *std::sin(infoTruth.phi);

	//	std::cout << "detector face hit x" << infoTruth.x << ", y " << infoTruth.y << ", z "<< (zHGC/std::cos(theta))* std::sin(theta) * sinh(eta_ECAL) << std::endl;
		infoTruth.phiCrackDistance = phiCrackDistance(infoTruth.x, infoTruth.y);
	//	std::cout << "debug "<<igp <<"  out of " << gens.size() <<  std::endl;

		//	if (fabs(eta_ECAL) >2.7 || fabs(eta_ECAL)<1.6) continue; //FIXME might want to remove this later
	//	std::cout << "[debug] gen "<< count << "  pdgid " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << ", z " << z << ", converted " << nHitsBeforeHGC << ", pt " << gens[igp]->pt() << ", barcode " << *((gBarcodes[igp]).get()) <<  std::endl;// (gens[igp]->mother()->pdgId()) <<  std::endl;

		for (unsigned int isc =0; isc < sclusters.size() ; isc++){ //subloop over sc's to find matches
			// calculate dR... dE = dEta, dP = dPhi
			float dE = fabs(sclusters[isc]->eta() - eta_ECAL);
			dE =dE*dE;
			float dP = DeltaPhi(sclusters[isc]->phi(), gens[igp]->phi());
			dP =dP*dP;
			float dR = sqrt(dE +dP);
			//std::cout << "[debug ] SC e " << (sclusters[isc])->energy() << ", et " <<  sclusters[isc]->energy()/std::cosh(eta_ECAL)<< ", eta:phi " << (sclusters[isc])->eta() << ":" <<(sclusters[isc])->phi() << ", dR " << dR << std::endl;

			if (dR < maxDr) { // only true if dR is both below limit value and smaller than previous best 
				//	std::cout << "check : ereco : " << resumEmEnergy(sclusters[isc], clusters) << ", etrue : " << infoTruth.eTrue << " dR = "<< dR << std::endl; 
			//	std::cout << "MATCH" << std::endl;
				// so, if we have a dR match, then store the corresponding match index.
				if (dR < dRBest){
					dRBest = dR;
					infoTruth.matchIndex = isc;
	//	std::cout << " debug dR " << dR << std::endl;
				}
			}

		}

		if (infoTruth.matchIndex <-1) { std::cout << "NO MATCH " << std::endl;};

		//float recoPt = -999.;
		if (infoTruth.matchIndex >-1) {
			//	recoPt = sclusters[infoTruth.matchIndex]->energy() /std::cosh(sclusters[infoTruth.matchIndex]->eta());
			//infoTruth.eReco = (sclusters[infoTruth.matchIndex]->energy());
			infoTruth.eReco = resumEmEnergy(sclusters[infoTruth.matchIndex], clusters);
			infoTruth.eRecoCleaned = resumEmEnergyTest(sclusters[infoTruth.matchIndex], clusters);
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

			std::cout << " gen " << infoTruth.eta << ", " << infoTruth.phi << ", SC " << infoTruth.etaSC << ", " << infoTruth.phiSC <<", idx " <<  infoTruth.matchIndex <<  ", drbest " << dRBest << std::endl;

			auto 	sc = sclusters[infoTruth.matchIndex];
			for(unsigned int l =0 ; l<sc->clusters().size(); l++){

				if (sc->clusters()[l]->eta() == sc->seed()->eta()) continue;

				dEta_h->Fill(sc->seed()->eta() - (sc->clusters())[l]->eta(), sc->clusters()[l]->energy()/infoTruth.eReco);
				dPhi_h->Fill(DeltaPhi(sc->seed()->phi() , (sc->clusters()[l])->phi()),sc->clusters()[l]->energy()/infoTruth.eReco);

				if(infoTruth.phiWidth<0.012){
					dEta_Unconv_h->Fill(sc->seed()->eta() - (sc->clusters())[l]->eta());
					dPhi_Unconv_h->Fill(DeltaPhi(sc->seed()->phi() , (sc->clusters()[l])->phi()));
				}


			}


 //std::cout <<  "debug debug  " << std::endl;
			infoTruth.eReco3x3AnneMarie  =annemarieEnergy(rechitvec,sclusters[infoTruth.matchIndex],geometrySource_,iEvent, iSetup, eta_ECAL );

		std::cout << " photon " << infoTruth.eta << ", " << infoTruth.phi << " ereco " << infoTruth.eReco3x3AnneMarie << ", etrue " << infoTruth.eTrue << ", CONVERTED " << infoTruth.converted <<" ,eReco/eTrue > 0.7 ?"<< (infoTruth.eReco3x3AnneMarie/infoTruth.eTrue>0.7) <<  std::endl;
	//		std::cout << "AM energy " << infoTruth.eReco3x3AnneMarie << ", etrue " << infoTruth.eTrue <<", eta ECAL " << eta_ECAL <<  std::endl;
			infoTruth.eReco5x5Claude  =claudeEnergy(sclusters[infoTruth.matchIndex],geometrySource_,iEvent, iSetup, 1. ) ;
			//	infoTruth.eReco5x5ClaudeSC  =annemarieEnergy(sclusters[infoTruth.matchIndex],geometrySource_,iEvent, iSetup, 2. ) ;
			//		float test = get3x3EmEnergy(  rcs, clusters, sclusters[infoTruth.matchIndex], iSetup);
			//		infoTruth.eReco3x3  = test;
			//std::cout << "*********** E TRUE " << infoTruth.eTrue << ", E RECO (OLD) "<< infoTruth.eReco  <<", E RECO (NEW) " << test  << " E CLAUDE " << infoTruth.eReco3x3Claude<< std::endl;
			//	if (test) ;

			for (unsigned int ic =0 ; ic < sclusters[infoTruth.matchIndex]->clusters().size() ; ic++){
				infoTruth.nClusters09 = ic+1; 
				if((sclusters[infoTruth.matchIndex]->clusters())[ic]->energy()/infoTruth.eReco >0.9) break;
			}
		}

		//if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  recoPt >40)
		if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  gens[igp]->pt() >40)
		{
			//	eRoT_OLD_h->Fill(sclusters[infoTruth.matchIndex]->energy()/gens[igp]->energy());
			eRoT_h->Fill(infoTruth.eReco/gens[igp]->energy());
			//eRoT3x3_h->Fill(infoTruth.eReco3x3/gens[igp]->energy());
			eRoT3x3Claude_h->Fill(infoTruth.eReco3x3Claude/infoTruth.eTrue);
			eRoT5x5Claude_h->Fill(infoTruth.eReco5x5Claude/infoTruth.eTrue);
			eRoT5x5ClaudeSC_h->Fill(infoTruth.eReco5x5ClaudeSC/infoTruth.eTrue);
			phiW_v_eRoT_h->Fill(infoTruth.phiWidth,infoTruth.eReco_over_eTrue);
			dPhi_v_eRoT_h ->Fill(infoTruth.phiSeed - infoTruth.phi,infoTruth.eReco_over_eTrue);
			nClust_v_eRoT_h ->Fill(infoTruth.clustersSize,infoTruth.eReco_over_eTrue);
			nClust90_v_eRoT_h->Fill(infoTruth.nClusters09,infoTruth.eReco_over_eTrue);


			if(infoTruth.eReco_over_eTrue <0.7){
				/*		//std::cout << "MATCH with eReco/eTrue <0.7" << std::endl;
				//std::cout << "_____|  RECO    |  TRUE    |" <<  std::endl;
				//std::cout << " eta | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->eta() <<" | "<< gens[igp]->eta()<< std::endl;
				//std::cout << " phi | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->phi() <<" | "<< gens[igp]->phi()<< std::endl;
				//std::cout << " e   | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->energy() <<" | "<< gens[igp]->energy()<< std::endl;
				//std::cout << " e   | " << std::setprecision(5) << infoTruth.eReco <<" | "<< infoTruth.eTrue<< std::endl;*/
			}
		}

		treeTruth->Fill();
	}

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

	//	//std::cout << " sc " << isc<< " eta " << sclusters[isc]->eta() << ", phi " << sclusters[isc]->phi()<< std::endl;

	info.pt= sclusters[isc]->energy() /std::cosh(sclusters[isc]->eta());
	info.eta=sclusters[isc]->eta();
	info.phi=sclusters[isc]->phi();
	info.eReco = sclusters[isc]->energy();
	info.clustersSize = (int) sclusters[isc]->clustersSize();;

	// fill histograms with eta/phi info
	rechit_h->Fill(info.eta);
	phi_h->Fill(info.phi);

	float dRBest = 999.; // dR best is used to find the gen-reco match with smallest dR.

	// now loop over gen particles
	for (unsigned int igp =0; igp < gens.size() ; igp++){ // igp = index_gen_particles

	float dE = sclusters[isc]->eta() - gens[igp]->eta();
	dE =dE*dE;
	float dP = sclusters[isc]->phi() - gens[igp]->phi();
	dP =dP*dP;
	float dR = sqrt(dE +dP);

	if (dR < maxDr && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
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

	}*/
	/*
		 if (eTrue[0] <0) {
		 eTrue[0] = gens[igp]->energy();
		 } else {
		 eTrue[1] = eTrue[0];
		 eTrue[0] = gens[igp]->energy();
		 }*/
	/*		}
				}
				if(sclusters.size()){
//	tree->Fill();
//	dEta_h->Fill(info.eta[0]+info.eta[1]);
//	dPhi_h->Fill(fabs(info.phi[0]-info.phi[1]));

if( info.matchIndex[0] > -1){
//	eRoT_h->Fill(info.eReco[0]/info.eTrue[info.matchIndex[0]]);
//std::cout << "RECO 0 matched TRUE " << info.matchIndex[0] << std::endl;
//std::cout << "eta_reco " << info.eta[0] << ", eta_true" << gens[info.matchIndex[0]]->eta() << std::endl;
//std::cout << "phi_reco " << info.phi[0] << ", phi_true" << gens[info.matchIndex[0]]->phi() << std::endl;
//std::cout << "e_reco " << info.eReco[0] << ", e_true" << gens[info.matchIndex[0]]->energy() << " (" << info.eTrue[info.matchIndex[0]] << ")" <<  std::endl;
}
if( info.matchIndex[1] > -1){
//	eRoT_h->Fill(info.eReco[1]/info.eTrue[info.matchIndex[1]]);
//std::cout << "RECO 1 matched TRUE " << info.matchIndex[1] << std::endl;
//std::cout << "eta_reco " << info.eta[1] << ", eta_true" << gens[info.matchIndex[1]]->eta() << std::endl;
//std::cout << "phi_reco " << info.phi[1] << ", phi_true" << gens[info.matchIndex[1]]->phi() << std::endl;
//std::cout << "e_reco " << info.eReco[1] << ", e_true" << gens[info.matchIndex[1]]->energy() << " (" << info.eTrue[info.matchIndex[1]] << ")" <<  std::endl;
}



}


*/
return ;
}

double AMStyleHggAnalyser::absWeight(const unsigned layer, const bool dedx){
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



// ------------ method called once each job just before starting event loop  ------------
	void 
AMStyleHggAnalyser::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
AMStyleHggAnalyser::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
AMStyleHggAnalyser::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
	edm::ESHandle<HGCalGeometry> hgcGeo;
	iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
	hgcEEGeom_=hgcGeo.product();
}

// ------------ method called when ending the processing of a run  ------------
	void 
AMStyleHggAnalyser::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
AMStyleHggAnalyser::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
AMStyleHggAnalyser::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AMStyleHggAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AMStyleHggAnalyser);
