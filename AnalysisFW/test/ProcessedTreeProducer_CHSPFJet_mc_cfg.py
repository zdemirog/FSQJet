import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.GlobalTag.globaltag = "PHYS14_25_V2::All"

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"root://eoscms//eos/cms/store/user/ilknur/rootfiles_for_testing/2015_MC/QCD_Pt_15to3000_Tune4C_Flat_13TeV_pythia8_AODSIM_NoPU_PHYS14_25_V1_v1.root"
)
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_CHSJet_mc.root'))

process.ak8 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    isCHS           = cms.untracked.bool(True),
    pfjets          = cms.InputTag('ak8PFJets'),
    pfjetschs       = cms.InputTag('ak8PFJetsCHS'),
    #calojets        = cms.InputTag('ak7CaloJets'),
    genjets         = cms.untracked.InputTag('ak8GenJets'),                         
    ## database entry for the uncertainties ######
    PFPayloadName        = cms.string(''),
    PFPayloadNameCHS     = cms.string(''),                         
    #CaloPayloadName = cms.string(''),
    jecUncSrc         = cms.string(''),
    jecUncSrcCHS       = cms.string(''),                         
    jecUncSrcNames  = cms.vstring(''),
    jecL1FastFile          = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L1FastJet_AK8PF.txt'),
    jecL2RelativeFile      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L2Relative_AK8PF.txt'),
    jecL3AbsoluteFile      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L3Absolute_AK8PF.txt'),
    jecL2L3ResidualFile    = cms.string(''),
    jecL1FastFileCHS          = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L1FastJet_AK8PFchs.txt'),
    jecL2RelativeFileCHS      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L2Relative_AK8PFchs.txt'),
    jecL3AbsoluteFileCHS      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L3Absolute_AK8PFchs.txt'),
    jecL2L3ResidualFileCHS    = cms.string(''),

    ## calojet ID and extender for the JTA #######
    #calojetID       = cms.InputTag('ak7JetID'),
    #calojetExtender = cms.InputTag('ak7JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## rho #######################################
    srcCaloRho      = cms.InputTag('fixedGridRhoFastjetAllCalo'),
    srcPFRho        = cms.InputTag('fixedGridRhoFastjetAll'),
    #srcPFRho        = cms.InputTag('fixedGridRhoFastjetAllCalo'),                         
    ## MC $ Generator flags ######################
    isMCarlo        = cms.untracked.bool(True),
    useGenInfo      = cms.untracked.bool(True),
    ## simulated PU ##############################
    srcPU           = cms.untracked.InputTag('addPileupInfo'),
    ## preselection cuts #########################
    maxY            = cms.double(5.0), 
    minPFPt         = cms.double(5),
    minPFFatPt      = cms.double(10),
    maxPFFatEta     = cms.double(2.5),
    #minCaloPt       = cms.double(5),
    minGenPt        = cms.untracked.double(5),
    minNPFJets      = cms.int32(1),
    #minNCaloJets    = cms.int32(1), 
    minJJMass       = cms.double(-1),
    ##minTwrEta       = cms.double(0.0),
    ##maxTwrEta       = cms.double(5.0),                         
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_ZeroBias_v1'
                                  ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT")
    ## jec services ##############################
    ##pfjecService    = cms.string('ak7PFL1FastL2L3'),
    ##calojecService  = cms.string('ak7CaloL1L2L3')
)

process.ak5 = process.ak8.clone(
    pfjets           = 'ak5PFJets',
    pfjetschs       =  'ak5PFJetsCHS',
    #calojets         = 'ak5CaloJets',
    genjets          = 'ak5GenJets',
    # rho #######################################                                                                                                                                      
    #srcCaloRho      = cms.InputTag('ak5CaloJets','rho'),
    #srcPFRho        = cms.InputTag('ak5PFJets','rho'),
    #calojetID        = 'ak5JetID',
    #calojetExtender  = 'ak5JetExtender',
    jecL1FastFile          = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L1FastJet_AK5PF.txt'),
    jecL2RelativeFile      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L2Relative_AK5PF.txt'),
    jecL3AbsoluteFile      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L3Absolute_AK5PF.txt'),
    jecL2L3ResidualFile    = cms.string(''),
    jecL1FastFileCHS          = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L1FastJet_AK5PFchs.txt'),
    jecL2RelativeFileCHS      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L2Relative_AK5PFchs.txt'),
    jecL3AbsoluteFileCHS      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L3Absolute_AK5PFchs.txt'),
    jecL2L3ResidualFileCHS    = cms.string(''),

    #pfjecService     = 'ak5PFL1FastL2L3',
    #calojecService   = 'ak5CaloL1L2L3',
    printTriggerMenu = False 
)
process.ak4 = process.ak8.clone(
    pfjets           = 'ak4PFJets',
    #calojets         = 'ak5CaloJets',
    pfjetschs       =  'ak4PFJetsCHS',
    genjets          = 'ak4GenJets',
    #calojetID        = 'ak4JetID',                                                                                                                                                     
    #calojetExtender  = 'ak4JetExtender',                                                                                                                                                  # rho #######################################                                                                                                                                      
    #srcCaloRho      = cms.InputTag('ak4CaloJets','rho'),
    #srcPFRho        = cms.InputTag('ak4PFJets','rho'),
    jecL1FastFile          = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L1FastJet_AK4PF.txt'),
    jecL2RelativeFile      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L2Relative_AK4PF.txt'),
    jecL3AbsoluteFile      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L3Absolute_AK4PF.txt'),
    jecL2L3ResidualFile    = cms.string(''),
    jecL1FastFileCHS          = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L1FastJet_AK4PFchs.txt'),
    jecL2RelativeFileCHS      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L2Relative_AK4PFchs.txt'),
    jecL3AbsoluteFileCHS      = cms.string('PHYS14_V4_MC/PHYS14_V4_MC_L3Absolute_AK4PFchs.txt'),
    jecL2L3ResidualFileCHS    = cms.string(''),

    #pfjecService     = 'ak5PFL1FastL2L3',                                                                                                                                              
    #calojecService   = 'ak5CaloL1L2L3',                                                                                                                                                
    printTriggerMenu = False
)

process.path = cms.Path(process.ak8 * process.ak5 * process.ak4)
#process.path = cms.Path(process.ak5)


