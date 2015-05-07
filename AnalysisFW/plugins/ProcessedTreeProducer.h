#ifndef ProcessedTreeProducer_h
#define ProcessedTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "FSQJet/AnalysisFW/interface/QCDJet.h"
#include "FSQJet/AnalysisFW/interface/QCDEvent.h"
#include "FSQJet/AnalysisFW/interface/QCDEventHdr.h"
#include "FSQJet/AnalysisFW/interface/QCDCaloJet.h"
#include "FSQJet/AnalysisFW/interface/QCDPFJet.h"
#include "FSQJet/AnalysisFW/interface/QCDTower.h"
#include "FSQJet/AnalysisFW/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;

class ProcessedTreeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit ProcessedTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~ProcessedTreeProducer();
  private:  
    void buildTree();


    static bool sort_calojets(QCDCaloJet j1, QCDCaloJet j2) {
      return j1.ptCor() > j2.ptCor();
    }
    static bool sort_pfjets(QCDPFJet j1, QCDPFJet j2) {
      return j1.ptCor() > j2.ptCor();
    }
    
    static bool sort_tower(QCDTower t1, QCDTower t2) {
      return t1.twr_et() > t2.twr_et();
    }
    //---- configurable parameters --------  
    bool   mIsMCarlo;
    bool   mIsCHS;
    bool   mUseGenInfo;
    bool   mPrintTriggerMenu;
    bool   isPFJecUncSet_,isPFJecUncSetCHS_;
    int    mGoodVtxNdof,mMinNCaloJets,mMinNPFJets;
    double mGoodVtxZ; 
    double mMinCaloPt,mMinPFPt,mMinPFFatPt,mMaxPFFatEta,mMinGenPt,mMaxY,mMinJJMass,mMinTwrEt,mMinTwrEta,mMaxTwrEta;
    //std::string mCaloJECservice;
    std::string mPFJECservice;
    std::string mPFPayloadName;
    //std::string mCaloPayloadName;
    std::string mPFJECUncSrc;
    std::string mPFPayloadNameNoCHS;
    std::string mPFPayloadNameCHS;
    std::string mPFJECUncSrcNoCHS;
    std::string mPFJECUncSrcCHS;
    // add JEC 53x//
    std::string mJECL1FastFile, mJECL2RelativeFile,mJECL3AbsoluteFile, mJECL2L3ResidualFile,
	mJECL1FastFileCHS, mJECL2RelativeFileCHS,mJECL3AbsoluteFileCHS, mJECL2L3ResidualFileCHS;

    // ---- Jet Corrector Parameter --//                                                                                                                                                  
    JetCorrectorParameters *L1Fast, *L2Relative, *L3Absolute, *L2L3Residual;
    vector<JetCorrectorParameters> vecL1Fast, vecL2Relative, vecL3Absolute, vecL2L3Residual;
    FactorizedJetCorrector *jecL1Fast, *jecL2Relative, *jecL3Absolute, *jecL2L3Residual;

    //******************************//

    
    std::vector<std::string> mPFJECUncSrcNames;
    edm::InputTag mCaloJetsName;
    edm::InputTag mPFJetsName;
    edm::InputTag mPFJetsNameCHS;
    edm::InputTag mGenJetsName;
    edm::InputTag mCaloJetID;
    edm::InputTag mCaloJetExtender;
    edm::InputTag mOfflineVertices;
    edm::InputTag mSrcCaloRho;
    edm::InputTag mSrcPFRho;
    edm::InputTag mSrcPU;
    //---- TRIGGER -------------------------
    std::string   processName_;
    std::vector<std::string> triggerNames_;
    std::vector<unsigned int> triggerIndex_;
    edm::InputTag triggerResultsTag_;
    edm::InputTag triggerEventTag_;
    edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
    edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
    HLTConfigProvider hltConfig_;
    //---- CORRECTORS ----------------------
    const JetCorrector *mPFJEC;
    const JetCorrector *mCALOJEC;
    JetCorrectionUncertainty *mCALOUnc;
    JetCorrectionUncertainty *mPFUnc;
    std::vector<JetCorrectionUncertainty*> mPFUncSrc;
    
    edm::Service<TFileService> fs;
    TTree *mTree;
    TTree *mTree2;

    TH1F *mTriggerPassHisto,*mTriggerNamesHisto; 
    //---- TREE variables --------
    QCDEvent *mEvent;
    QCDEvent *mEvent2;
};

#endif
