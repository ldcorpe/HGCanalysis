from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'phoID_100515_PH1_A1k_PU140_v2'
config.General.workArea = 'crab_projects'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hoverEHggConfig_PH1_A1k_PU140_MVA.py'
#config.JobType.pyCfgParams = ['output=alcareco.root', 'skim=none', 'type=ALCARECO','doTree=0', 'jsonFile=json.txt', 'secondaryOutput=ntuple.root', 'isCrab=1']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles = ['PH1_A0_PU50_PhotonID_BDTG.weights.xml']
config.Data.inputDataset = '/GJet_Pt-15to3000_Tune4C_14TeV_pythia8/GEM2019Upg14DR-Phase1age1kfixJan23_PU140BX25_PH1_1K_FB_V2-v2/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outLFNDirBase = '' #add a path to your area
config.Data.publication = False
config.Site.storageSite = "T2_CH_CERN"
