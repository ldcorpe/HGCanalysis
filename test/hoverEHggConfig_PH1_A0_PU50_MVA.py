import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

process.load('FWCore.MessageService.MessageLogger_cfi')
#v6 geometry
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV6MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV6Muon_cff')
#v5 geometry
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
#v4 geometry
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV4MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV4Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGCEE_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')


process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False)
                                        #SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

ffile=0
step=-1
preFix='Single13_CMSSW_6_2_0_SLHC18'
doFullAnalysis=True
import os,sys
from UserCode.HGCanalysis.storeTools_cff import fillFromStore

process_ =0
subsample =0;

import os,sys
if(len(sys.argv)>2):
	process_= int(sys.argv[2])
	print 'index %d'%(process_)
if (len(sys.argv) ==0):
	print 'no index! default is 0'


processList = ["gamJet","ZEE","QCD"]
#fileNames
if(process_ ==0): 
	with open('sample/140PU/gamJet_age140PU.txt','r') as f:
		fileNames = f.readlines()

if(process_ ==1): 
	with open('sample/140PU/Zee.txt','r') as f:
		fileNames = f.readlines()


if(process_ ==2): 
	with open('sample/140PU/Qcd.txt','r') as f:
		fileNames = f.readlines()

tot = len(fileNames)
N= tot/3


#process.GlobalTag.globaltag = 'auto:upgradePLS3'

print 'processing %s files'%(processList[process_])
#print '%d, %s'%(process_*N, ((process_+1)*(N)-1))
for x in range(0,tot):
	print fileNames[x]

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(fileNames),
                            #fileNames=cms.untracked.vstring("/store/mc/GEM2019Upg14DR/GJet_Pt-15to3000_Tune4C_14TeV_pythia8/GEN-SIM-RECO/Phase1age1kfixJan23_PU140BX25_PH1_1K_FB_V2-v2/70000/DA561C36-B8BF-E411-88FE-0025902008A4.root"),
                            skipEvents=cms.untracked.uint32(subsample*1000))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30000) )


#load the analyzer
import getpass
whoami=getpass.getuser()
outputTag=preFix.replace('/','_')
process.TFileService = cms.Service("TFileService", fileName = cms.string('HoverEHggPU140_%s_ssz3_%d_.root'%(processList[process_],subsample)))

process.hgg = cms.EDAnalyzer("HoverEAnalyzer_Phase1",
											  geometrySource = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive', 'HGCalHEScintillatorSensitive'),
												eeRecHitCollection = cms.untracked.InputTag("particleFlowRecHitECAL:Cleaned"),
												endcapPhotonCollection = cms.untracked.InputTag("photons"),
												endcapClusterCollection = cms.untracked.InputTag("particleFlowClusterECAL"),
												genParticlesTag =  cms.untracked.InputTag("genParticles"),
                             g4VerticesSource  = cms.untracked.string('g4SimHits'),
                             g4TracksSource    = cms.untracked.string('g4SimHits'),
												endcapHCALClusters= cms.InputTag("particleFlowClusterHCAL"),
												hcalTowers = cms.InputTag("towerMaker"),
												hOverEPtMin = cms.double(2.),
												hOverEMethodEndcap = cms.int32(3),
												hOverEConeSize = cms.double(0.05),
												PU = cms.int32(140),
												process_ = cms.int32(1),
												genMatchPU_ = cms.int32(0),
												MVAweightfile = cms.FileInPath("UserCode/HGCanalysis/data/PH1_A0_PU50_PhotonID_BDTG.weights.xml",)
		)


#from Configuration.Generator.Pythia8CommonSettings_cfi import *
#from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

#source = cms.Source("EmptySource")
#generator = cms.EDFilter('Pythia8GeneratorFilter',
#		comEnergy = cms.double(13000.0),
#		crossSection = cms.untracked.double(1.0),
#		filterEfficiency = cms.untracked.double(1.0),
#		maxEventsToPrint = cms.untracked.int32(0),
#		pythiaHepMCVerbosity = cms.untracked.bool(False),
#		pythiaPylistVerbosity = cms.untracked.int32(0),
#
#		PythiaParameters = cms.PSet(
#		pythia8CommonSettingsBlock,
#			pythia8CUEP8M1SettingsBlock,
#			processParameters = cms.vstring(
#			'PromptPhoton:qg2qgamma = on       ! prompt photon production',
#				'PromptPhoton:qqbar2ggamma = on    ! prompt photon production',
#				'PromptPhoton:gg2ggamma = on       ! prompt photon production',
#				'PhaseSpace:pTHatMin = 40.         ! minimum pt hat for hard interactions', 
#				'PhaseSpace:pTHatMax = -1          ! maximum pt hat for hard interactions'),
#			parameterSets = cms.vstring('pythia8CommonSettings',
#				'pythia8CUEP8M1Settings',
#				'processParameters')
#			)
#		)

process.gj_filter = cms.EDFilter("PythiaFilterGammaGamma",
				PtSeedThr = cms.untracked.double(5.0),
				EtaSeedThr = cms.untracked.double(2.8),
				PtGammaThr = cms.untracked.double(0.0),
				EtaGammaThr = cms.untracked.double(2.8),
				PtElThr = cms.untracked.double(2.0),
				EtaElThr = cms.untracked.double(2.8),
				dRSeedMax = cms.untracked.double(0.0),
				dPhiSeedMax = cms.untracked.double(0.2),
				dEtaSeedMax = cms.untracked.double(0.12),
				dRNarrowCone = cms.untracked.double(0.02),
				PtTkThr = cms.untracked.double(1.6),
				EtaTkThr = cms.untracked.double(2.2),
				dRTkMax = cms.untracked.double(0.2),
				PtMinCandidate1 = cms.untracked.double(15.),
				PtMinCandidate2 = cms.untracked.double(15.),
				EtaMaxCandidate = cms.untracked.double(3.0),
				NTkConeMax = cms.untracked.int32(2),
				NTkConeSum = cms.untracked.int32(4),
				InvMassWide = cms.untracked.double(80.0),
				InvMassNarrow = cms.untracked.double(14000.0),
				EnergyCut = cms.untracked.double(1.0),
				AcceptPrompts = cms.untracked.bool(True),
				PromptPtThreshold = cms.untracked.double(15.0)   

	)

#ProductionFilterSequence = cms.Sequence(generator*gj_filter)


#run it
process.p = cms.Path(#process.analysis
			# process.particleFlowRecHitHGCEELC*
			#process.gj_filter*
			process.hgg
			)

