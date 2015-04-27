// -*- C++ -*-
// //
// // Package:    Test3x3Analyzer
// // Class:      Test3x3Analyzer
// // 
// //
// //
//
//
// // system include files
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
#include "DataFormats/Math/interface/Point3D.h"

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
#include "TFile.h"
#include "TGraphErrors.h"

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
float maxDr3 =0.05;

// information to be loaded into TTree
struct infoTruth_t {

	float eReco_over_eTrue;
	float eta;
	float etaSC;
	float etaSeed;
	float phi;
	float phiSC;
	float phiSeed;
	int   matchIndex;
	float pt;
	float eTrue;
	float eReco;
	float eReco3x3;
	float eReco3x3Simple;
	float eSeed;
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
class Test3x3Analyzer : public edm::EDAnalyzer {
	public:
		explicit Test3x3Analyzer(const edm::ParameterSet&);
		~Test3x3Analyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		double DeltaPhi(const double & phi1, const double & phi2);
		std::pair<double,double>  axisEtaPhi( const edm::Ptr<reco::SuperCluster>& sc,  	const edm::Ptr<reco::GenParticle>& gen, double z);
		double correctRecHitEnergy(const edm::Ptr<reco::PFRecHit>& rh);
		double correctRecHitEnergy(const edm::Ref<std::vector<reco::PFRecHit> >& rh);
		float resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc, const edm::PtrVector<reco::PFCluster>& clusters);
		float resumEmEnergyTest( const edm::Ptr<reco::SuperCluster>& sc, const edm::PtrVector<reco::PFCluster>& clusters,const  edm::PtrVector<reco::PFRecHit>& rcs );
		float clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c, const edm::PtrVector<reco::PFCluster>& clusters);
		//	float get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rechits, const edm::PtrVector<reco::PFCluster>& clusters,  const edm::Ptr<reco::SuperCluster>& sc);
		float get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rcs, const edm::PtrVector<reco::PFCluster>& clusters,  const edm::Ptr<reco::SuperCluster>& sc, const edm::Ptr<reco::GenParticle>& gens, const edm::EventSetup& iSetup);
		float get3x3EmEnergySimple(const edm::PtrVector<reco::PFRecHit>& rcs, const edm::PtrVector<reco::PFCluster>& clusters,  const edm::Ptr<reco::SuperCluster>& sc, const edm::Ptr<reco::GenParticle>& gens, const edm::EventSetup& iSetup);

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
		edm::EDGetTokenT<edm::View<reco::PFRecHit> > eeRecHitCollection_     ;

		std::vector<std::string> geometrySource_;

		edm::Handle<HGCRecHitCollection> recHits_;	
		TTree *tree;
		TTree *treeTruth;
		info_t info;
		infoTruth_t infoTruth;
		TH1F *C_e_h; 
		TH1F *C_e_resum_h; 
		TH1F *C_e_frac_h; 
		TH1F *eta_h; 
		TH1F *phi_h; 
		TH1F *dEta_h; 
		TH2F *phiW_v_eRoT_h; 
		TH1F *dPhi_h;
		TH1F *nSC_h;
		TH1F *eRoT_OLD_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT3x3_h;  //energy true over reco ie etrue/ ereco
		TH1F *eRoT3x3Simple_h;  //energy true over reco ie etrue/ ereco
		TH2F *dPhi_v_eRoT_h    ; 
		TH2F *nClust_v_eRoT_h  ; 
		TH2F *nClust90_v_eRoT_h;
		const TGraphErrors * _hgcOverburdenParam;
		const TGraphErrors * _hgcLambdaOverburdenParam;
		const std::vector<double> _weights_ee;


};

// constructor
Test3x3Analyzer::Test3x3Analyzer(const edm::ParameterSet& iConfig):
	endcapRecHitCollection_(consumes <edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
	endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
	endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
	genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles")))),
	eeRecHitCollection_(consumes <edm::View<reco::PFRecHit> > (iConfig.getUntrackedParameter<edm::InputTag>("eeRecHitCollection",edm::InputTag("particleFlowRecHitHGCEE:Cleaned")))),
	_hgcOverburdenParam(nullptr),
	_hgcLambdaOverburdenParam(nullptr),
	_weights_ee(iConfig.getParameter<std::vector<double> >("weights_ee"))
{

	geometrySource_ = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
	edm::Service<TFileService> fs_;
	C_e_h         = fs_->make<TH1F>("C_e_h","C_e_h",200,0,200);
	C_e_resum_h   = fs_->make<TH1F>("C_e_resum_h","C_e_resum_h",200,0,200);
	C_e_frac_h   = fs_->make<TH1F>("C_e_frac_h","C_e_frac_h",2000,-1,1);
	eta_h         = fs_->make<TH1F>("eta_h","eta_h",100,-5,5);
	phi_h         = fs_->make<TH1F>("phi_h","phi_h",100,-5,5);
	dEta_h        = fs_->make<TH1F>("dEta_h","dEta_h",1000,-1,1);
	dPhi_h        = fs_->make<TH1F>("dPhi_h","dPhi_h",1000,3,3.3);
	eRoT_OLD_h        = fs_->make<TH1F>("eRoT_OLD_h","eRoT_OLD_h",1000,-2,2);
	eRoT_h        = fs_->make<TH1F>("eRoT_h","eRoT_h",1000,-2,2);
	eRoT3x3_h        = fs_->make<TH1F>("eRoT3x3_h","eRoT3x3_h",1000,-2,2);
	eRoT3x3Simple_h        = fs_->make<TH1F>("eRoT3x3Simple_h","eRoT3x3Simple_h",1000,-2,2);
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
	treeTruth->Branch("eReco3x3"              ,&infoTruth.eReco3x3            ,"eReco3x3/F");
	treeTruth->Branch("eReco3x3Simple"              ,&infoTruth.eReco3x3Simple            ,"eReco3xSimple3/F");
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
}

// destructor
Test3x3Analyzer::~Test3x3Analyzer()
{

}


//
// member functions
//

// ------------ method called for each event  ------------

std::pair<double,double>  Test3x3Analyzer::axisEtaPhi( const edm::Ptr<reco::SuperCluster>& sc,  	const edm::Ptr<reco::GenParticle>& gen, double z){

	math::XYZPoint sPos = sc->seed()->position();
	math::XYZPoint vPos = gen->vertex();
	double t = ( z- sPos.z())/(vPos.z()- sPos.z());
	double linex = ((-sPos.x() + vPos.x())*t + sPos.x());
	double liney = ((-sPos.y() + vPos.y())*t + sPos.y());

	math::XYZPoint p (linex, liney, z);
	//p.setX()= linex;
	//p.setY()= liney;
	//p.setZ() =z;

	double eta = p.eta();
	double phi = p.phi();

	std::pair<double, double> res;

	res.first =eta;
	res.second =phi;

	return res;
}

float Test3x3Analyzer::resumEmEnergyTest( const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters, const  edm::PtrVector<reco::PFRecHit>& rcs){

	double _coef_a =80.0837; 
	double _coef_b =-107.229; 
	float total=0;
	float totalnew=0;
	float subtotal=0;
	//	float newtotal =0;
	for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
		//	std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
		for (unsigned int j =0 ; j < clusters.size() ; j++){

			if (clusters[j]->position()==(sc->clusters())[ic]->position()) {
				//		std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
				total = total +(clusters[j]->emEnergy());
				subtotal=0;	

				for (unsigned int i =0 ; i <clusters[j]->recHitFractions().size() ; i++){

					DetId  recid =  (sc->clusters())[ic]->hitsAndFractions()[i].first;

					for (unsigned int irh=0; irh < rcs.size() ; irh++){

						if (rcs[irh]->detId() == recid)
						{

							subtotal = subtotal+ correctRecHitEnergy(rcs[irh]);
							//	std::cout << " add " << i <<"th rechit, subtotal "<< subtotal<< std::endl;
							break;

						}

					}
				}

				subtotal = std::max(0., subtotal-_coef_b/_coef_a);
				totalnew = totalnew + subtotal;


				//		std::cout<< "RECHITS ENERGY emEnergy" << clusters[j]->emEnergy() <<", resum " << subtotal <<  std::endl;


				break;
			}
		}
	}

	C_e_h->Fill(total);
	C_e_resum_h->Fill(totalnew);
	if ((sc->energy() /std::cosh(sc->eta()))>40.) {
		C_e_frac_h->Fill((totalnew-total)/total);
	}
	return total;
}


float Test3x3Analyzer::resumEmEnergy(const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){

	float total=0;
	for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
		//	std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
		for (unsigned int j =0 ; j < clusters.size() ; j++){

			if (clusters[j]->position()==(sc->clusters())[ic]->position()) {
				//		std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
				total = total +(clusters[j]->emEnergy());

				std::vector< reco::PFRecHitFraction > recs =  clusters[j]->recHitFractions();


				//std::cout<< "SIZE OF RECHITS (alt)" << recs.size() <<", layer " << clusters[j]->layer() <<  std::endl;

				break;
			}
		}

	}

	return total;
}
double Test3x3Analyzer::DeltaPhi(const double & phi1, const double & phi2){
	double dphi = phi1 - phi2;
	if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
	if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
	return dphi;
}


//float Test3x3Analyzer::get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rechits, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc){
double Test3x3Analyzer::correctRecHitEnergy(const edm::Ptr<reco::PFRecHit>& rh){

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
//double Test3x3Analyzer::correctRecHitEnergy(const edm::Ptr<reco::PFRecHit>& rh){
double Test3x3Analyzer::correctRecHitEnergy(const edm::Ref<std::vector<reco::PFRecHit> >& rh){

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

float Test3x3Analyzer::get3x3EmEnergySimple(const  edm::PtrVector<reco::PFRecHit>& rcs, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc,const edm::Ptr<reco::GenParticle>& gen, const edm::EventSetup& iSetup ){


	double _coef_a =80.0837; 
	double _coef_b =-107.229; 

	float totalE =0;
	float newtotalE =0;
	//	float eta = sc->seed()->eta();
	//	float phi = sc->seed()->phi();
	//	float eta = infoTruth.eta;
	//	float phi = infoTruth.phi;
	float debug_max =0;
	//	float dRmin =9999.;
	int match =-1;
	DetId  id_; 

	for (unsigned int iLayer =0; iLayer<30 ; iLayer++){  // iLayer = iLayer

		bool out =1;
		float dRmin =9999.;
		match =-1;
		for (unsigned int k=0; k < rcs.size() ; k++) {
			int hgc_layer = HGCEEDetId(rcs[k]->detId()).layer();

			if(hgc_layer != (int) iLayer) continue;
			//		std::cout << "TEST, rechits rho,eta, phi " << rcs[k]->positionREP() << std::endl;
			float layer_z =  rcs[k]->position().z() ; //this isa strange way to get z_layer, but since only rechits on that layer are processed at this step, we know we can ask their z easily.
			float eta = axisEtaPhi(sc, gen,  layer_z).first; 
			float phi = axisEtaPhi(sc, gen,  layer_z).second;
			if (eta * sc->seed()->eta() <0 ) continue; // we want to only consider rechits on the right side. if the SC we are interested in is in EE+, then rechits with -eta will cause continue; Likewise, if SC is in EE-, rechits in EE+ will cause continue;
			if(out) {	
				std::cout << "Layer " << iLayer << ", z " << layer_z << " eta,phi "<< eta << ", " << phi << " | cf sc z, eta, phi " <<  sc->seed()->z() << ", "  <<  sc->seed()->eta() << ", " <<  sc->seed()->phi() << "| vtz z = " << gen->vertex().z()<< std::endl;
				out =0;
			}

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
				//§	std::cout << "DEBUG dRmin " << dRmin << ", match " << match << std::endl; 
			}
			//k++;
		}

		//	std::cout << "[info] eta,phi match for layer "<< iLayer <<" has index " << match << ", DRMIN " << dRmin <<  std::endl; 
		//	std::cout << "[DEBUG]i max rechit index "<< match<< std::endl; 
		if (match==-1) continue;	
		//		std::cout << ", and "<< rcs[match]->neighbours8().size() << "neighbours "<<  std::endl; 

		//		float maxE =rcs[match]->energy();
		//	std::cout << "DEBUG central: " <<  maxE << std::endl;
		//	std::cout << "[info] maxE has index " << maxEIndex << " and energy " << maxE << std::endl; 
		//	std::cout << "[info]  **CORR* maxE has index " << maxEIndex << " and energy " << maxE << std::endl; 

		//	std::cout << "[INFO]  best match is central	"<< std::endl;
		totalE = totalE + rcs[match]->energy();
		newtotalE = newtotalE + correctRecHitEnergy( rcs[match]);
		//	std::cout << "[INFO] 		Addding Rechit 0 " << "gives "<< totalE << std::endl;
		//	std::cout << "[INFO] **CORR** 	Addding Rechit 0 " << "gives "<<correctRecHitEnergy( rcs[match])
		//		<< std::endl;

		for( unsigned int iNeighbour2=0; iNeighbour2< (rcs[match]->neighbours8()).size(); iNeighbour2++ ){


			totalE = totalE + (rcs[match]->neighbours8())[iNeighbour2]->energy();
			newtotalE = newtotalE + correctRecHitEnergy(((rcs[match]->neighbours8())[iNeighbour2]));

			//		std::cout << "[INFO] 		Addding Rechit  "<< iNeighbour2 << "gives "<< totalE << std::endl;
			//		std::cout << "[INFO] **CORR**	Addding Rechit  "<< iNeighbour2 << "gives "<< correctRecHitEnergy(((rcs[match]->neighbours8())[iNeighbour2]))
			//		<< std::endl;
		}



		//	std::cout << "[INFO] Layer "<< iLayer  << " subtotal "<< totalE << std::endl;
		//	std::cout << "[INFO] **CORR** Layer "<< iLayer  << " subtotal "<< newtotalE << std::endl;

	}
	//decode position
	//uint32_t recoDetId(recHits->find(id_)->id());

	//const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	//int layer( ((rcs[match]->detId() >> 19) & 0x1f) );// + layerCtrOffset-1 );

	//int hgc_layer = HGCEEDetId(rcs[match]->detId()).layer();
	//refPos= geom->getPosition(id_) ;
	totalE = std::max(0., totalE-_coef_b/_coef_a);
	newtotalE = std::max(0., newtotalE-_coef_b/_coef_a);

	//std::cout << "NEIGHOURS " << recHits->find(id_)->neighbours4()->size()<< std::endl;
	if(match>-1){
		//		std::cout << "SC E " << sc->energy()  << ", calculated E " << totalE  <<std::endl;
		//		std::cout << "**CORR** SC E " << sc->energy()  << ", calculated E " << newtotalE  <<std::endl;
		//	std::cout << "match rechit eta, phi, layer " << rcs[match]->positionREP().eta() << ", " <<rcs[match]->positionREP().eta()  << ", " << layer <<  std::endl;
		//	std::cout << "match rechit eta, phi, layer (alt) " << rcs[match]->positionREP().eta() << ", " <<rcs[match]->positionREP().eta()  << ", " << hgc_layer <<  std::endl;
	}
	//	totalE = std::max(0., totalE-_coef_b/_coef_a);
	//return (totalE);				
	return (newtotalE);				
}


//float Test3x3Analyzer::get3x3EmEnergy(const edm::PtrVector<reco::PFRecHit>& rechits, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc){
float Test3x3Analyzer::get3x3EmEnergy(const  edm::PtrVector<reco::PFRecHit>& rcs, const edm::PtrVector<reco::PFCluster>& clusters, const edm::Ptr<reco::SuperCluster>& sc,const edm::Ptr<reco::GenParticle>& gen, const edm::EventSetup& iSetup ){


	double _coef_a =80.0837; 
	double _coef_b =-107.229; 

	float totalE =0;
	float newtotalE =0;
	//	float eta = sc->seed()->eta();
	//	float phi = sc->seed()->phi();
	//	float eta = infoTruth.eta;
	//	float phi = infoTruth.phi;
	float debug_max =0;
	//	float dRmin =9999.;
	int match =-1;
	DetId  id_; 

	for (unsigned int iLayer =0; iLayer<30 ; iLayer++){  // iLayer = iLayer

		bool out =1;
		float dRmin =9999.;
		match =-1;
		for (unsigned int k=0; k < rcs.size() ; k++) {
			int hgc_layer = HGCEEDetId(rcs[k]->detId()).layer();

			if(hgc_layer != (int) iLayer) continue;
			//		std::cout << "TEST, rechits rho,eta, phi " << rcs[k]->positionREP() << std::endl;
			float layer_z =  rcs[k]->position().z() ; //this isa strange way to get z_layer, but since only rechits on that layer are processed at this step, we know we can ask their z easily.
			float eta = axisEtaPhi(sc, gen,  layer_z).first; 
			float phi = axisEtaPhi(sc, gen,  layer_z).second;
			if (eta * sc->seed()->eta() <0 ) continue; // we want to only consider rechits on the right side. if the SC we are interested in is in EE+, then rechits with -eta will cause continue; Likewise, if SC is in EE-, rechits in EE+ will cause continue;
			if(out) {	
				std::cout << "Layer " << iLayer << ", z " << layer_z << " eta,phi "<< eta << ", " << phi << " | cf sc z, eta, phi " <<  sc->seed()->z() << ", "  <<  sc->seed()->eta() << ", " <<  sc->seed()->phi() << "| vtz z = " << gen->vertex().z()<< std::endl;
				out =0;
			}

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
				//§	std::cout << "DEBUG dRmin " << dRmin << ", match " << match << std::endl; 
			}
			//k++;
		}

		//	std::cout << "[info] eta,phi match for layer "<< iLayer <<" has index " << match << ", DRMIN " << dRmin <<  std::endl; 
		//	std::cout << "[DEBUG]i max rechit index "<< match<< std::endl; 
		if (match==-1) continue;	
		//		std::cout << ", and "<< rcs[match]->neighbours8().size() << "neighbours "<<  std::endl; 

		float maxE =rcs[match]->energy();
		//	std::cout << "DEBUG central: " <<  maxE << std::endl;
		int maxEIndex =-1;
		for( unsigned int iNeighbour=0; iNeighbour< rcs[match]->neighbours8().size(); iNeighbour++ ){


			//	float e =1;
			//	float e = (rcs[match]->neighbours8()[iNeighbour])->energy();
			float e2 = correctRecHitEnergy((rcs[match]->neighbours8()[iNeighbour]));
			//		std::cout << "DEBUG " <<  e2 << ", nNeigh " << (rcs[match]->neighbours8()[iNeighbour])->neighbours8().size()<< std::endl;
			if (e2 > maxE) {

				maxE=e2;
				maxEIndex=iNeighbour;
				//eta =(rcs[match]->neighbours8()[iNeighbour])->positionREP().eta();
				//phi =(rcs[match]->neighbours8()[iNeighbour])->positionREP().phi();

			}

		}
		//	std::cout << "[info] maxE has index " << maxEIndex << " and energy " << maxE << std::endl; 
		//	std::cout << "[info]  **CORR* maxE has index " << maxEIndex << " and energy " << maxE << std::endl; 
		if (maxEIndex>-1){

			totalE = totalE + rcs[match]->neighbours8()[maxEIndex]->energy();
			newtotalE = newtotalE + correctRecHitEnergy((rcs[match]->neighbours8()[maxEIndex]));
			//		std::cout << "[INFO]  best match is around	"<< std::endl;
			//		std::cout << "[INFO] 	**CORR** 	Addding Rechit 0 " << "gives "<<correctRecHitEnergy((rcs[match]->neighbours8()[maxEIndex])) << std::endl;

			for( unsigned int iNeighbour2=0; iNeighbour2< ((rcs[match]->neighbours8()[maxEIndex])->neighbours8()).size(); iNeighbour2++ ){


				totalE = totalE + ((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2]->energy();
				newtotalE = newtotalE + correctRecHitEnergy((((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2] ));
				//		std::cout << "[INFO] 		Addding Rechit  "<< iNeighbour2 << "gives "<< totalE << std::endl;
				//			std::cout << "[INFO] **CORR** 	Addding Rechit  "<< iNeighbour2 << "gives "<< correctRecHitEnergy((((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2] )) << ", w nerighbours " <<  (((rcs[match]->neighbours8()[maxEIndex])->neighbours8())[iNeighbour2] )->neighbours8().size()  << std::endl ;


			}



		}  else {

			//	std::cout << "[INFO]  best match is central	"<< std::endl;
			totalE = totalE + rcs[match]->energy();
			newtotalE = newtotalE + correctRecHitEnergy( rcs[match]);
			//	std::cout << "[INFO] 		Addding Rechit 0 " << "gives "<< totalE << std::endl;
			//	std::cout << "[INFO] **CORR** 	Addding Rechit 0 " << "gives "<<correctRecHitEnergy( rcs[match])
			//		<< std::endl;

			for( unsigned int iNeighbour2=0; iNeighbour2< (rcs[match]->neighbours8()).size(); iNeighbour2++ ){


				totalE = totalE + (rcs[match]->neighbours8())[iNeighbour2]->energy();
				newtotalE = newtotalE + correctRecHitEnergy(((rcs[match]->neighbours8())[iNeighbour2]));

				//		std::cout << "[INFO] 		Addding Rechit  "<< iNeighbour2 << "gives "<< totalE << std::endl;
				//		std::cout << "[INFO] **CORR**	Addding Rechit  "<< iNeighbour2 << "gives "<< correctRecHitEnergy(((rcs[match]->neighbours8())[iNeighbour2]))
				//		<< std::endl;
			}

		}

		//	std::cout << "[INFO] Layer "<< iLayer  << " subtotal "<< totalE << std::endl;
		//	std::cout << "[INFO] **CORR** Layer "<< iLayer  << " subtotal "<< newtotalE << std::endl;

	}
	//decode position
	//uint32_t recoDetId(recHits->find(id_)->id());

	//const GlobalPoint refPos( std::move( geom->getPosition(recoDetId) ) );
	//int layer( ((rcs[match]->detId() >> 19) & 0x1f) );// + layerCtrOffset-1 );

	//int hgc_layer = HGCEEDetId(rcs[match]->detId()).layer();
	//refPos= geom->getPosition(id_) ;
	totalE = std::max(0., totalE-_coef_b/_coef_a);
	newtotalE = std::max(0., newtotalE-_coef_b/_coef_a);

	//std::cout << "NEIGHOURS " << recHits->find(id_)->neighbours4()->size()<< std::endl;
	if(match>-1){
		//		std::cout << "SC E " << sc->energy()  << ", calculated E " << totalE  <<std::endl;
		//		std::cout << "**CORR** SC E " << sc->energy()  << ", calculated E " << newtotalE  <<std::endl;
		//	std::cout << "match rechit eta, phi, layer " << rcs[match]->positionREP().eta() << ", " <<rcs[match]->positionREP().eta()  << ", " << layer <<  std::endl;
		//	std::cout << "match rechit eta, phi, layer (alt) " << rcs[match]->positionREP().eta() << ", " <<rcs[match]->positionREP().eta()  << ", " << hgc_layer <<  std::endl;
	}
	//	totalE = std::max(0., totalE-_coef_b/_coef_a);
	//return (totalE);				
	return (newtotalE);				
}




float Test3x3Analyzer::clusterEmEnergy(const edm::Ptr<reco::CaloCluster>& c,const edm::PtrVector<reco::PFCluster>& clusters){

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
Test3x3Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	iEvent.getByLabel(edm::InputTag("HGCalRecHit:HGCEERecHits"),recHits_);

	/*	Handle<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >  > HGCEERechits;
			iEvent.getByToken(endcapRecHitCollection_,HGCEERechits);
			iEvent.getByLabel(edm::InputTag("HGCalRecHit",hitCollections_[i]),recHits);*/
	//	const PtrVector<edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> > >& rechits = recHits_->ptrVector();
	//	std::cout << "RECHIT SIZE " << recHits_->size() << std::endl;

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
		infoTruth.eReco  =-999.;         
		infoTruth.eTrue         =-999.;  
		infoTruth.eSeed           =-999.;
		infoTruth.eSeed_over_eReco=-999.;
		infoTruth.etaWidth=-999.;
		infoTruth.phiWidth=-999.;

		if (gens[igp]->pdgId() != 22 || gens[igp]->status() != 3) continue;
		assert(gens.size() >0); // only the case for the electron gun sample
		infoTruth.eta=gens[igp]->eta();
		infoTruth.phi      = gens[igp]->phi();

		std::cout << "[debug] gen pdgid " << gens[igp]->pdgId() << ", vertex " << gens[igp]->vertex() << ", e " << gens[igp]->energy() << ", mother " << std::endl;// (gens[igp]->mother()->pdgId()) <<  std::endl;

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

			if (dR < maxDr3 && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
				dRBest = dR;

				// so, if we have a dR match, then store the corresponding match index.
				infoTruth.matchIndex = isc;
			}

		}

		//float recoPt = -999.;
		if (infoTruth.matchIndex >-1) {
			//	recoPt = sclusters[infoTruth.matchIndex]->energy() /std::cosh(sclusters[infoTruth.matchIndex]->eta());
			//infoTruth.eReco = (sclusters[infoTruth.matchIndex]->energy());
			infoTruth.eReco = resumEmEnergyTest(sclusters[infoTruth.matchIndex], clusters, rcs);
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

			float test = get3x3EmEnergy(  rcs, clusters, sclusters[infoTruth.matchIndex], gens[igp], iSetup);
			float testSimple = get3x3EmEnergySimple(  rcs, clusters, sclusters[infoTruth.matchIndex], gens[igp], iSetup);
			std::cout << "*********** E TRUE " << infoTruth.eTrue << ", E RECO (OLD) "<< infoTruth.eReco  <<", E RECO (NEW) " << test << ", E Reco (SIMPLE) " << testSimple << std::endl;
			//	if (test) ;
			infoTruth.eReco3x3  = test;
			infoTruth.eReco3x3Simple  = testSimple;

			for (unsigned int ic =0 ; ic < sclusters[infoTruth.matchIndex]->clusters().size() ; ic++){
				infoTruth.nClusters09 = ic+1; 
				if((sclusters[infoTruth.matchIndex]->clusters())[ic]->energy()/infoTruth.eReco >0.9) break;
			}
		}

		//if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  recoPt >40)
		if(fabs(infoTruth.eta) > 1.6 && fabs(infoTruth.eta) <2.8 && infoTruth.matchIndex > -1 &&  gens[igp]->pt() >40)
		{
			eRoT_OLD_h->Fill(sclusters[infoTruth.matchIndex]->energy()/gens[igp]->energy());
			eRoT_h->Fill(infoTruth.eReco/gens[igp]->energy());
			eRoT3x3_h->Fill(infoTruth.eReco3x3/gens[igp]->energy());
			eRoT3x3Simple_h->Fill(infoTruth.eReco3x3Simple/gens[igp]->energy());
			phiW_v_eRoT_h->Fill(infoTruth.phiWidth,infoTruth.eReco_over_eTrue);
			dPhi_v_eRoT_h ->Fill(infoTruth.phiSeed - infoTruth.phi,infoTruth.eReco_over_eTrue);
			nClust_v_eRoT_h ->Fill(infoTruth.clustersSize,infoTruth.eReco_over_eTrue);
			nClust90_v_eRoT_h->Fill(infoTruth.nClusters09,infoTruth.eReco_over_eTrue);


			if(infoTruth.eReco_over_eTrue <0.7){
				/*		std::cout << "MATCH with eReco/eTrue <0.7" << std::endl;
							std::cout << "_____|  RECO    |  TRUE    |" <<  std::endl;
							std::cout << " eta | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->eta() <<" | "<< gens[igp]->eta()<< std::endl;
							std::cout << " phi | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->phi() <<" | "<< gens[igp]->phi()<< std::endl;
							std::cout << " e   | " << std::setprecision(5) << sclusters[infoTruth.matchIndex]->energy() <<" | "<< gens[igp]->energy()<< std::endl;
							std::cout << " e   | " << std::setprecision(5) << infoTruth.eReco <<" | "<< infoTruth.eTrue<< std::endl;*/
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

		//	std::cout << " sc " << isc<< " eta " << sclusters[isc]->eta() << ", phi " << sclusters[isc]->phi()<< std::endl;

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

			if (dR < maxDr3 && dR < dRBest) { // only true if dR is both below limit value and smaller than previous best dR.
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
				if(sclusters.size()){
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
Test3x3Analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
Test3x3Analyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
Test3x3Analyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
Test3x3Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
Test3x3Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
Test3x3Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Test3x3Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Test3x3Analyzer);

