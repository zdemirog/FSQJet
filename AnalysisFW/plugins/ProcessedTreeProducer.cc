#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include "TTree.h"
#include <vector>
#include <cassert>
#include <TLorentzVector.h>

#include "FSQJet/AnalysisFW/plugins/ProcessedTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

ProcessedTreeProducer::ProcessedTreeProducer(edm::ParameterSet const& cfg) 
{
    //mPFJECservice      = cfg.getParameter<std::string>               ("pfjecService");
  mGoodVtxNdof       = cfg.getParameter<double>                    ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>                    ("goodVtxZ");

  mMinPFPt           = cfg.getParameter<double>                    ("minPFPt");
  //mMinTwrEt          = cfg.getParameter<double>                    ("minTwrEt");
  mMinPFFatPt        = cfg.getParameter<double>                    ("minPFFatPt");
  mMaxPFFatEta       = cfg.getParameter<double>                    ("maxPFFatEta");
  mMinJJMass         = cfg.getParameter<double>                    ("minJJMass");
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  //mMinTwrEta         = cfg.getParameter<double>                    ("minTwrEta");
  //mMaxTwrEta         = cfg.getParameter<double>                    ("maxTwrEta");
  mMinNPFJets        = cfg.getParameter<int>                       ("minNPFJets");
 
  mOfflineVertices   = cfg.getParameter<edm::InputTag>             ("offlineVertices");
  mPFJetsName        = cfg.getParameter<edm::InputTag>             ("pfjets");
  mPFJetsNameCHS     = cfg.getParameter<edm::InputTag>             ("pfjetschs");
  mSrcCaloRho        = cfg.getParameter<edm::InputTag>             ("srcCaloRho");
  mSrcPFRho          = cfg.getParameter<edm::InputTag>             ("srcPFRho");
  mSrcPU             = cfg.getUntrackedParameter<edm::InputTag>    ("srcPU",edm::InputTag("addPileupInfo"));
  mGenJetsName       = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  mPrintTriggerMenu  = cfg.getUntrackedParameter<bool>             ("printTriggerMenu",false);
  mIsMCarlo          = cfg.getUntrackedParameter<bool>             ("isMCarlo",false);
  mIsCHS             = cfg.getUntrackedParameter<bool>             ("isCHS",false);
  mUseGenInfo        = cfg.getUntrackedParameter<bool>             ("useGenInfo",false);
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",5);
  processName_       = cfg.getParameter<std::string>               ("processName");
  triggerNames_      = cfg.getParameter<std::vector<std::string> > ("triggerName");
  triggerResultsTag_ = cfg.getParameter<edm::InputTag>             ("triggerResults");
  triggerEventTag_   = cfg.getParameter<edm::InputTag>             ("triggerEvent");

  mPFPayloadNameNoCHS= cfg.getParameter<std::string>               ("PFPayloadName");
  mPFPayloadNameCHS  = cfg.getParameter<std::string>               ("PFPayloadNameCHS");
  mPFJECUncSrcNoCHS  = cfg.getParameter<std::string>               ("jecUncSrc");
  mPFJECUncSrcCHS    = cfg.getParameter<std::string>               ("jecUncSrcCHS");
  mPFJECUncSrcNames  = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");

  mJECL1FastFile          = cfg.getParameter<std::string>          ("jecL1FastFile");
  mJECL1FastFileCHS       = cfg.getParameter<std::string>          ("jecL1FastFileCHS");
  mJECL2RelativeFile      = cfg.getParameter<std::string>          ("jecL2RelativeFile");
  mJECL2RelativeFileCHS   = cfg.getParameter<std::string>          ("jecL2RelativeFileCHS");
  mJECL3AbsoluteFile      = cfg.getParameter<std::string>          ("jecL3AbsoluteFile");
  mJECL3AbsoluteFileCHS   = cfg.getParameter<std::string>          ("jecL3AbsoluteFileCHS");
  mJECL2L3ResidualFile    = cfg.getParameter<std::string>          ("jecL2L3ResidualFile");
  mJECL2L3ResidualFileCHS = cfg.getParameter<std::string>          ("jecL2L3ResidualFileCHS");
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::beginJob() 
{
  mTree = fs->make<TTree>("ProcessedTree","ProcessedTree");
  mEvent = new QCDEvent();
  mTree->Branch("events","QCDEvent",&mEvent);
  mTriggerNamesHisto = fs->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  mTriggerNamesHisto->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<triggerNames_.size();i++)
    mTriggerNamesHisto->Fill(triggerNames_[i].c_str(),1);
  mTriggerPassHisto = fs->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  mTriggerPassHisto->SetBit(TH1::kCanRebin);
  isPFJecUncSet_ = false;
  isPFJecUncSetCHS_ = false;
  
  // --- Initializing the jet correctors ---- //                                                                     
  if(mIsCHS) {
      L1Fast       = new JetCorrectorParameters(mJECL1FastFileCHS.c_str());
      L2Relative   = new JetCorrectorParameters(mJECL2RelativeFileCHS.c_str());
      L3Absolute   = new JetCorrectorParameters(mJECL3AbsoluteFileCHS.c_str());
      if(!mIsMCarlo)
	  L2L3Residual = new JetCorrectorParameters(mJECL2L3ResidualFileCHS.c_str());
  } // if(mIsCHS)                                                                                                    
  else {
      L1Fast       = new JetCorrectorParameters(mJECL1FastFile.c_str());
      L2Relative   = new JetCorrectorParameters(mJECL2RelativeFile.c_str());
      L3Absolute   = new JetCorrectorParameters(mJECL3AbsoluteFile.c_str());
      if(!mIsMCarlo)
	  L2L3Residual = new JetCorrectorParameters(mJECL2L3ResidualFile.c_str());
  } // else -- non-chs    


  vecL1Fast.push_back(*L1Fast);
  vecL2Relative.push_back(*L2Relative);
  vecL3Absolute.push_back(*L3Absolute);
  if(!mIsMCarlo)
      vecL2L3Residual.push_back(*L2L3Residual);

  jecL1Fast       = new FactorizedJetCorrector(vecL1Fast);
  jecL2Relative   = new FactorizedJetCorrector(vecL2Relative);
  jecL3Absolute   = new FactorizedJetCorrector(vecL3Absolute);
  if(!mIsMCarlo)
      jecL2L3Residual = new FactorizedJetCorrector(vecL2L3Residual);

} 
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::endJob() 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger names in (new) config
      cout<<"New trigger menu found !!!"<<endl;
      triggerIndex_.clear(); 
      const unsigned int n(hltConfig_.size());
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
        cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" ";  
        if (triggerIndex_[itrig] >= n)
          cout<<"does not exist in the current menu"<<endl;
        else
          cout<<"exists"<<endl;
      }
      cout << "Available TriggerNames are: " << endl;
      if (mPrintTriggerMenu)
        hltConfig_.dump("Triggers");
    }
  } 
  else {
    cout << "ProcessedTreeProducer::analyze:"
         << " config extraction failure with process name "
         << processName_ << endl;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::analyze(edm::Event const& event, edm::EventSetup const& iSetup) 
{ 
  
  vector<QCDPFJet>      mPFJets;
  vector<QCDJet>        mPFFatJets;
  vector<QCDPFJet>      tmpPFJets;
  vector<QCDPFJet>      mPFJetsCHS;
  vector<LorentzVector> mGenJets;
  QCDEventHdr mEvtHdr; 
  QCDMET mPFMet;
  //-------------- Basic Event Info ------------------------------
  mEvtHdr.setRun(event.id().run());
  mEvtHdr.setEvt(event.id().event());
  mEvtHdr.setLumi(event.luminosityBlock());
  mEvtHdr.setBunch(event.bunchCrossing());
  //-------------- Beam Spot --------------------------------------
  Handle<reco::BeamSpot> beamSpot;
  event.getByLabel("offlineBeamSpot", beamSpot);
  if (beamSpot.isValid())
    mEvtHdr.setBS(beamSpot->x0(),beamSpot->y0(),beamSpot->z0());
  else
    mEvtHdr.setBS(-999,-999,-999);

  //-------------- HCAL Noise Summary -----------------------------
  Handle<bool> noiseSummary; 	 
  /*if (!mIsMCarlo) {
    event.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), noiseSummary);         
    mEvtHdr.setHCALNoise(*noiseSummary);
  }
  else*/
    mEvtHdr.setHCALNoise(true);

  //-------------- Trigger Info -----------------------------------
  event.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "ProcessedTreeProducer::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  event.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "ProcessedTreeProducer::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }
  vector<int> L1Prescales,HLTPrescales,Fired;
  vector<vector<LorentzVector> > mL1Objects,mHLTObjects;
  // sanity check
  assert(triggerResultsHandle_->size() == hltConfig_.size());
  //------ loop over all trigger names ---------
  for(unsigned itrig=0;itrig<triggerNames_.size() && !mIsMCarlo;itrig++) {
    bool accept(false);
    int preL1(-1);
    int preHLT(-1);
    int tmpFired(-1); 
    vector<LorentzVector> vvL1,vvHLT; 
    if (triggerIndex_[itrig] < hltConfig_.size()) {
      accept = triggerResultsHandle_->accept(triggerIndex_[itrig]);
      const std::pair<int,int> prescales(hltConfig_.prescaleValues(event,iSetup,triggerNames_[itrig]));
      preL1    = prescales.first;
      preHLT   = prescales.second;
      if (!accept)
        tmpFired = 0;
      else {
        mTriggerPassHisto->Fill(triggerNames_[itrig].c_str(),1);
        tmpFired = 1;
      }
      //--------- modules on this trigger path--------------
      const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_[itrig]));
      const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex_[itrig]));
      bool foundL1(false);
      for(unsigned int j=0; j<=moduleIndex; ++j) {
        const string& moduleLabel(moduleLabels[j]);
        const string  moduleType(hltConfig_.moduleType(moduleLabel));
        //--------check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
        if (filterIndex<triggerEventHandle_->sizeFilters()) {
          const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
          const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
          const size_type nI(VIDS.size());
          const size_type nK(KEYS.size());
          assert(nI==nK);
          const size_type n(max(nI,nK));
          const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
          if (foundL1) {
            for(size_type i=0; i!=n; ++i) {
              const TriggerObject& TO(TOC[KEYS[i]]);
              TLorentzVector P4;
              P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
              LorentzVector qcdhltobj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
              vvHLT.push_back(qcdhltobj);
              //cout<<TO.pt()<<endl;
            }
          }
          else { 
            for(size_type i=0; i!=n; ++i) {
              const TriggerObject& TO(TOC[KEYS[i]]);
              TLorentzVector P4;
              P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
              LorentzVector qcdl1obj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
              vvL1.push_back(qcdl1obj);
              //cout<<TO.pt()<<endl;  
            }
            foundL1 = true; 
          }
        }
      }// loop over modules
    }// if the trigger exists in the menu
    //cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" "<<accept<<" "<<tmpFired<<endl;
    Fired.push_back(tmpFired);
    L1Prescales.push_back(preL1);
    HLTPrescales.push_back(preHLT);
    mL1Objects.push_back(vvL1);
    mHLTObjects.push_back(vvHLT);
  }// loop over trigger names  
  mEvent->setTrigDecision(Fired);
  mEvent->setPrescales(L1Prescales,HLTPrescales);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  //-------------- Vertex Info -----------------------------------
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(mOfflineVertices,recVtxs);
  //------------- reject events without reco vertices ------------
  int VtxGood(0);
  bool isPVgood(false);
  float PVx(0),PVy(0),PVz(0),PVndof(0);
  for(VertexCollection::const_iterator i_vtx = recVtxs->begin(); i_vtx != recVtxs->end(); i_vtx++) {
    int index = i_vtx-recVtxs->begin();
    if (index == 0) {
      PVx    = i_vtx->x();
      PVy    = i_vtx->y();
      PVz    = i_vtx->z();
      PVndof = i_vtx->ndof();
    }
    if (!(i_vtx->isFake()) && i_vtx->ndof() >= mGoodVtxNdof && fabs(i_vtx->z()) <= mGoodVtxZ) {
      if (index == 0) {
        isPVgood = true;
      }
      VtxGood++;
    }
  }
  mEvtHdr.setVertices(recVtxs->size(),VtxGood);
  mEvtHdr.setPV(isPVgood,PVndof,PVx,PVy,PVz);
  //-------------- Rho ------------------------------------------------
  Handle<double> rhoCalo;
  event.getByLabel(mSrcCaloRho,rhoCalo);
  Handle<double> rhoPF;
  event.getByLabel(mSrcPFRho,rhoPF);
  mEvtHdr.setRho(*rhoCalo,*rhoPF);
  
  //cout<<" rho pf "<<*rhoPF<<endl;

  //-------------- Generator Info -------------------------------------
  Handle<GenEventInfoProduct> hEventInfo;
  //-------------- Simulated PU Info ----------------------------------
  Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (mIsMCarlo && mUseGenInfo) { 
    event.getByLabel("generator", hEventInfo);
    mEvtHdr.setPthat(hEventInfo->binningValues()[0]);
    mEvtHdr.setWeight(hEventInfo->weight());
    event.getByLabel(mSrcPU, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PUI;
    int nbx = PupInfo->size();
    int ootpuEarly(0),ootpuLate(0),intpu(0);
    float Tnpv = -1.; // new variable for computing pileup weight factor for the event
    for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() < 0)
        ootpuEarly += PUI->getPU_NumInteractions();
      else if (PUI->getBunchCrossing() > 0)
        ootpuLate += PUI->getPU_NumInteractions();
      else {
        intpu += PUI->getPU_NumInteractions(); 
        Tnpv = PUI->getTrueNumInteractions();
       } 
    } 
     
    mEvtHdr.setPU(nbx,ootpuEarly,ootpuLate,intpu);
    mEvtHdr.setTrPu(Tnpv);
  } 
  else {
    mEvtHdr.setPthat(0);
    mEvtHdr.setWeight(0); 
    mEvtHdr.setPU(0,0,0,0);
    mEvtHdr.setTrPu(0);
  }
  //---------------- Jets ---------------------------------------------

  //mPFJEC   = JetCorrector::getJetCorrector(mPFJECservice,iSetup);
  // ----- Initializing for JEC uncertainty sources -------------- //                                                                                                                    
  std::string mPFPayloadName, mPFJECUncSrc ;
  if(mIsCHS) {
      mPFPayloadName = mPFPayloadNameCHS;
      mPFJECUncSrc   = mPFJECUncSrcCHS;
  }
  else {
      mPFPayloadName = mPFPayloadNameNoCHS;
      mPFJECUncSrc   = mPFJECUncSrcNoCHS;
  } // else                                                                                                                                                                           


  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParColl;
  if (mPFPayloadName != "" && !isPFJecUncSet_){
    iSetup.get<JetCorrectionsRecord>().get(mPFPayloadName,PFJetCorParColl); 
    JetCorrectorParameters const& PFJetCorPar = (*PFJetCorParColl)["Uncertainty"];
    mPFUnc = new JetCorrectionUncertainty(PFJetCorPar);
    if (mPFJECUncSrc != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        JetCorrectorParameters *par = new JetCorrectorParameters(mPFJECUncSrc,mPFJECUncSrcNames[isrc]); 
        JetCorrectionUncertainty *tmpUnc = new JetCorrectionUncertainty(*par);
        mPFUncSrc.push_back(tmpUnc);
      }
    }
    isPFJecUncSet_ = true;
  }
  
  Handle<GenJetCollection>  genjets;  
  if (mIsMCarlo) {
    event.getByLabel(mGenJetsName,genjets);
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
        mGenJets.push_back(i_gen->p4());
      }
    }
  }
  //----------- PFJets -------------------------
  Handle<PFJetCollection>   pfjets;
  if(mIsCHS) {event.getByLabel(mPFJetsNameCHS,pfjets);}
  else event.getByLabel(mPFJetsName,pfjets);

  for(PFJetCollection::const_iterator i_pfjet = pfjets->begin(); i_pfjet != pfjets->end(); i_pfjet++) {
    QCDPFJet qcdpfjet;
    // *********************JEC********************//
    vector<double> JecFactors; JecFactors.clear();
    // ---- Evaluating the L1Fast correction factor ---- //
    LorentzVector JetP4 = i_pfjet->p4(); // ---- accessing the 4-vector of the jet ---- //
    TLorentzVector tmpJet;
    tmpJet.SetPxPyPzE(JetP4.px(),JetP4.py(),JetP4.pz(), JetP4.energy());
    TLorentzVector UnCorrectedJet = tmpJet;
                                                                                                                   
    jecL1Fast->setJetPt(UnCorrectedJet.Pt());
    jecL1Fast->setJetA(i_pfjet->jetArea());
    jecL1Fast->setRho(*rhoPF);
    jecL1Fast->setJetEta(UnCorrectedJet.Eta());

    double corFactorL1Fast = jecL1Fast->getCorrection();
    //cout<<"L1Fast Cor Factor = "<<corFactorL1Fast<<endl;                                                                                                                        
    TLorentzVector tmpJetL1FastCorrected =
		   UnCorrectedJet * corFactorL1Fast; // ---- getting the jet corrected at L1Fast level ----- //    

      // ---- Evaluating the L2Relative correction factor ---- //                                                                                                                   
      jecL2Relative->setJetPt(tmpJetL1FastCorrected.Pt());
      jecL2Relative->setJetEta(tmpJetL1FastCorrected.Eta());
      double corFactorL2Relative = jecL2Relative->getCorrection();
      // cout<<"L2Relative Cor Factor"<<corFactorL2Relative<<endl;                                                                                                                   
      TLorentzVector tmpJetL1FastL2RelativeCorrected =
             	    tmpJetL1FastCorrected * corFactorL2Relative; //  ---- getting the jet corrected at L1Fast*L2Relative level ----- //                                            

      // ---- Evaluating the L3Absolute correction factor ---- //                                                                                                                   
     jecL3Absolute->setJetPt(tmpJetL1FastL2RelativeCorrected.Pt());
     jecL3Absolute->setJetEta(tmpJetL1FastL2RelativeCorrected.Eta());
     double corFactorL3Absolute = jecL3Absolute->getCorrection();
     ///cout<<"L3Absolute Cor Factor"<<corFactorL3Absolute<<endl;                                                                                                                   
     TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteCorrected =
            	    tmpJetL1FastL2RelativeCorrected * corFactorL3Absolute; // -- getting the jet corrected at L1Fast*L2Relative*L3Absolute level -- //       
 
     // ---- Evaluating the L2L3Rsidual correction factor ---- //                                                                                                                  
     TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected;
     double corFactorL2L3Residual(-1.);
     if(!mIsMCarlo) {
	 jecL2L3Residual->setJetPt(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pt());
	 jecL2L3Residual->setJetEta(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Eta());

	 corFactorL2L3Residual = jecL2L3Residual->getCorrection();
	 //cout<<"L2L3Rsidual Cor Factor"<<corFactorL2L3Residual<<endl;                                                                                                                
         tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected =
				  tmpJetL1FastL2RelativeL3AbsoluteCorrected * corFactorL2L3Residual; //  -- getting the jet corrected at L1Fast*L2Relative*L3Absolute*L2L3Residual level -- /    
     } // if(!mIsMCarlo)

     LorentzVector correctedJetP4; double CorFactor;
     if(!mIsMCarlo) { // -- if data, then take the full L1FastL2RelativeL3AbsoluteL2L3Residual corrected jet - //                                                                  
	 correctedJetP4 = LorentzVector(tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Px(),
					tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Py(),
					tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Pz(),
					tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.E());

	 CorFactor = tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.E()/UnCorrectedJet.E();
     }
     else { // -- if mc, the take the L1FastL2RelativeL3Absolute corrected jet --//                                                                                                
	 correctedJetP4 = LorentzVector(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Px(),
					tmpJetL1FastL2RelativeL3AbsoluteCorrected.Py(),
					tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pz(),
					tmpJetL1FastL2RelativeL3AbsoluteCorrected.E());


	 CorFactor = tmpJetL1FastL2RelativeL3AbsoluteCorrected.E()/UnCorrectedJet.E();
     }


     // ----- Storing the JEC correction factor in each label at a vector --- //                                                                                                   
     JecFactors.push_back(corFactorL1Fast); JecFactors.push_back(corFactorL2Relative);
     JecFactors.push_back(corFactorL3Absolute); JecFactors.push_back(corFactorL2L3Residual);
     JecFactors.push_back(CorFactor);
     //cout <<" Tot Correction Factor "<<CorFactor<<endl;

     //double scale = mPFJEC->correction(*i_pfjet,event,iSetup);

     //*************************************************************************************

    //---- preselection -----------------
    if (fabs(i_pfjet->y()) > mMaxY) continue;
    //---- vertex association -----------
    //---- get the vector of tracks -----
    reco::TrackRefVector vTrks(i_pfjet->getTrackRefs());
    float sumTrkPt(0.0),sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0),beta(0.0),betaStar(0.0);
    //---- loop over the tracks of the jet ----
    for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
      if (recVtxs->size() == 0) break;
      sumTrkPt += (*i_trk)->pt();
      //---- loop over all vertices ----------------------------
      for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
        //---- loop over the tracks associated with the vertex ---
        if (!((*recVtxs)[ivtx].isFake()) && (*recVtxs)[ivtx].ndof() >= mGoodVtxNdof && fabs((*recVtxs)[ivtx].z()) <= mGoodVtxZ) {
          for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk) {
            //---- match the jet track to the track from the vertex ----
            reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
            //---- check if the tracks match -------------------------
            if (trkRef == (*i_trk)) {
              if (ivtx == 0) {
                sumTrkPtBeta += (*i_trk)->pt();
              }
              else {
                sumTrkPtBetaStar += (*i_trk)->pt();
              }   
              break;
            }
          }
        } 
      }
    }
    if (sumTrkPt > 0) {
      beta     = sumTrkPtBeta/sumTrkPt;
      betaStar = sumTrkPtBetaStar/sumTrkPt;
    }
    qcdpfjet.setBeta(beta);
    qcdpfjet.setBetaStar(betaStar);
    //---- jec uncertainty --------------
    double unc(0.0);
    vector<float> uncSrc(0);
    if (mPFPayloadName != "") {
      mPFUnc->setJetEta(i_pfjet->eta());
      mPFUnc->setJetPt(CorFactor * i_pfjet->pt());
      //mPFUnc->setJetPt(Scale * i_pfjet->pt());
      unc = mPFUnc->getUncertainty(true);
    }
    if (mPFJECUncSrc != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        mPFUncSrc[isrc]->setJetEta(i_pfjet->eta());
        mPFUncSrc[isrc]->setJetPt(CorFactor * i_pfjet->pt());
	//mPFUncSrc[isrc]->setJetPt(Scale * i_pfjet->pt());
        float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
        uncSrc.push_back(unc1);
      }
    }
   
    // --- declaring the new jet with corrected P4 and unaltered energy fractions ---//                                                                                                  
    qcdpfjet.setP4(correctedJetP4);	       
    qcdpfjet.setCor(CorFactor); 
    //qcdpfjet.setP4(i_pfjet->p4());
    //qcdpfjet.setCor(scale);
    qcdpfjet.setUnc(unc);
    qcdpfjet.setUncSrc(uncSrc);
    qcdpfjet.setArea(i_pfjet->jetArea());
    qcdpfjet.setJecLabels(JecFactors);
    double chf   = i_pfjet->chargedHadronEnergyFraction();
    double nhf   = (i_pfjet->neutralHadronEnergy() + i_pfjet->HFHadronEnergy())/i_pfjet->energy();
    double phf   = i_pfjet->photonEnergyFraction();
    double elf   = i_pfjet->electronEnergyFraction();
    double muf   = i_pfjet->muonEnergyFraction();
    double hf_hf = i_pfjet->HFHadronEnergyFraction();
    double hf_phf= i_pfjet->HFEMEnergyFraction();
    int hf_hm    = i_pfjet->HFHadronMultiplicity();
    int hf_phm   = i_pfjet->HFEMMultiplicity();
    int chm      = i_pfjet->chargedHadronMultiplicity();
    int nhm      = i_pfjet->neutralHadronMultiplicity();
    int phm      = i_pfjet->photonMultiplicity();
    int elm      = i_pfjet->electronMultiplicity();
    int mum      = i_pfjet->muonMultiplicity();
    int npr      = i_pfjet->chargedMultiplicity() + i_pfjet->neutralMultiplicity();
    bool looseID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4));
    bool tightID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4));
    qcdpfjet.setLooseID(looseID);
    qcdpfjet.setTightID(tightID);
    qcdpfjet.setFrac(chf,nhf,phf,elf,muf);
    qcdpfjet.setMulti(npr,chm,nhm,phm,elm,mum);
    qcdpfjet.setHFFrac(hf_hf,hf_phf);
    qcdpfjet.setHFMulti(hf_hm,hf_phm);
    
    
    if (mIsMCarlo) {
      GenJetCollection::const_iterator i_matched;
      float rmin(999);
      for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
        double deltaR = reco::deltaR(*i_pfjet,*i_gen);
        if (deltaR < rmin) {
          rmin = deltaR;
          i_matched = i_gen;
        }
      }
      if (genjets->size() == 0) {
        LorentzVector tmpP4(0.0,0.0,0.0,0.0);
        qcdpfjet.setGen(tmpP4,0);
      }
      else
        qcdpfjet.setGen(i_matched->p4(),rmin);
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0); 
      qcdpfjet.setGen(tmpP4,0);
    }
    
    if (qcdpfjet.ptCor() >= mMinPFPt)
      mPFJets.push_back(qcdpfjet);
    if (qcdpfjet.ptCor() >= mMinPFFatPt && fabs(qcdpfjet.eta()) < mMaxPFFatEta && qcdpfjet.looseID())
      tmpPFJets.push_back(qcdpfjet);
  }// PF Jet
   
  sort(tmpPFJets.begin(),tmpPFJets.end(),sort_pfjets);

  if (tmpPFJets.size()>1) {
    LorentzVector lead[2], fat[2]; 
    float sumPt[2],sumPtUnc[2];
    for(unsigned i = 0; i<2; i++) {
      lead[i]     = tmpPFJets[i].p4()*tmpPFJets[i].cor();
      fat[i]      = tmpPFJets[i].p4()*tmpPFJets[i].cor();
      sumPt[i]    = tmpPFJets[i].ptCor();
      sumPtUnc[i] = tmpPFJets[i].ptCor() * tmpPFJets[i].unc();
    }
    double rmax = 1.1;
    for(unsigned i = 2; i<tmpPFJets.size(); i++) {
      LorentzVector cand = tmpPFJets[i].p4();
      double dR1 = deltaR(lead[0],cand);
      double dR2 = deltaR(lead[1],cand);
      int index(-1);
      if (dR1 < dR2 && dR1 < rmax) 
        index = 0;
      if (dR1 > dR2 && dR2 < rmax)
        index = 1;
      if (index > -1) {
        fat[index]      += cand * tmpPFJets[i].cor();
        sumPt[index]    += tmpPFJets[i].ptCor();
        sumPtUnc[index] += tmpPFJets[i].ptCor()*tmpPFJets[i].unc();
      } 
    }

    QCDJet fatJet[2];
    vector<float> uncSrc(0);
    for(unsigned i = 0; i<2; i++) { 
      fatJet[i].setP4(fat[i]);
      fatJet[i].setLooseID(tmpPFJets[i].looseID());
      fatJet[i].setTightID(tmpPFJets[i].tightID());
      fatJet[i].setCor(1.0);
      fatJet[i].setArea(0.0);
      fatJet[i].setUncSrc(uncSrc); 
      if (sumPt[i] > 0)
        fatJet[i].setUnc(sumPtUnc[i]/sumPt[i]);
      else
        fatJet[i].setUnc(0.0); 
      fatJet[i].setGen(tmpPFJets[i].genp4(),tmpPFJets[i].genR());
    }
    if (fatJet[0].pt()>fatJet[1].pt()) {
      mPFFatJets.push_back(fatJet[0]); 
      mPFFatJets.push_back(fatJet[1]);
    }
    else {
      mPFFatJets.push_back(fatJet[1]); 
      mPFFatJets.push_back(fatJet[0]);
    }
  }
    
  //---------------- met ---------------------------------------------
  Handle<PFMETCollection> pfmet;
  event.getByLabel("pfMet",pfmet);
  mPFMet.setVar((*pfmet)[0].et(),(*pfmet)[0].sumEt(),(*pfmet)[0].phi());
  //-------------- fill the tree -------------------------------------  
  sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
  mEvent->setEvtHdr(mEvtHdr);
  mEvent->setPFJets(mPFJets);
  mEvent->setFatJets(mPFFatJets);
  mEvent->setGenJets(mGenJets);
  mEvent->setPFMET(mPFMet);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  
  if ((mEvent->nPFJets() >= (unsigned)mMinNPFJets)){
      //if ((mEvent->pfmjjcor(0) >= mMinJJMass) ||(mEvent->fatmjjcor(0) >= mMinJJMass)) {
	  mTree->Fill();
	  //}
  }
  
}
//////////////////////////////////////////////////////////////////////////////////////////
ProcessedTreeProducer::~ProcessedTreeProducer() 
{
}

DEFINE_FWK_MODULE(ProcessedTreeProducer);
