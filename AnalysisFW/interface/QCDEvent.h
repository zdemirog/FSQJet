#ifndef QCDEvent_h
#define QCDEvent_h
#include "FSQJet/AnalysisFW/interface/QCDJet.h"
#include "FSQJet/AnalysisFW/interface/QCDMET.h"
#include "FSQJet/AnalysisFW/interface/QCDCaloJet.h"
#include "FSQJet/AnalysisFW/interface/QCDPFJet.h"
#include "FSQJet/AnalysisFW/interface/QCDTower.h"
#include "FSQJet/AnalysisFW/interface/QCDEventHdr.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include <vector>

class QCDEvent 
{
    public:
      typedef reco::Particle::LorentzVector LorentzVector;
      //------------ Constructor ------------------------------
      QCDEvent();
      //------------ Destructor -------------------------------
      ~QCDEvent();
      //------------ Set methods ------------------------------
      void setCaloMET(const QCDMET& fCaloMET)                     {CaloMet_ = fCaloMET;}
      void setPFMET(const QCDMET& fPFMET)                         {PFMet_ = fPFMET;}
      void setEvtHdr(const QCDEventHdr& fEvtHdr)                  {EvtHdr_ = fEvtHdr;}
      void setCaloJets(const std::vector<QCDCaloJet>& fCaloJets);
      void setPFJets(const std::vector<QCDPFJet>& fPFJets);
      void setCaloTowers(const std::vector<QCDTower>& fCaloTowers);
      void setFatJets(const std::vector<QCDJet>& fFatJets);
      void setGenJets(const std::vector<LorentzVector>& fGenJets);
      void setL1Obj(const std::vector<std::vector<LorentzVector> >& fL1Obj);
      void setHLTObj(const std::vector<std::vector<LorentzVector> >& fHLTObj);
      void setPrescales(const std::vector<int>& fPreL1, const std::vector<int>& fPreHLT) {L1Prescale_ = fPreL1; HLTPrescale_ = fPreHLT;}
      void setTrigDecision(const std::vector<int>& fTrigDecision) {TriggerDecision_ = fTrigDecision;}                           
      //------------ Get methods ------------------------------- 
      unsigned int nTriggers()                         const {return TriggerDecision_.size();}
      unsigned int nL1Obj(int i)                       const {return L1Obj_[i].size();}
      unsigned int nHLTObj(int i)                      const {return HLTObj_[i].size();}
      unsigned int nPFJets()                           const {return PFJets_.size();}
      unsigned int nTower()                            const {return CaloTowers_.size();}
      unsigned int nFatJets()                          const {return FatJets_.size();}
      unsigned int nCaloJets()                         const {return CaloJets_.size();}
      unsigned int nGenJets()                          const {return GenJets_.size();}
      int nGoodJets(int unc, int id, float ymax, float ptmin, std::vector<QCDJet> jets);
      int fired(int i)                                 const {return TriggerDecision_[i];}
      int preL1(int i)                                 const {return L1Prescale_[i];}
      int preHLT(int i)                                const {return HLTPrescale_[i];}
      float pfmjj();
      float calomjj();
      float genmjj(); 
      float pfmjjcor(int unc);
      float pfmjjcor(int unc,int src);
      float fatmjjcor(int unc);
      float calomjjcor(int unc);
      float pfmjjgen();
      float calomjjgen();
      const QCDMET&        calomet()                   const {return CaloMet_;}
      const QCDMET&        pfmet()                     const {return PFMet_;} 
      const LorentzVector& hltobj(int itrig, int iobj) const {return (HLTObj_[itrig])[iobj];}  
      const LorentzVector& l1obj(int itrig, int iobj)  const {return (L1Obj_[itrig])[iobj];}   
      const LorentzVector& genjet(int i)               const {return GenJets_[i];}
      const QCDPFJet&      pfjet(int i)                const {return PFJets_[i];}
      const QCDTower&      calotower(int i)            const {return CaloTowers_[i];}
      const QCDJet&        fatjet(int i)               const {return FatJets_[i];}
      const QCDCaloJet&    calojet(int i)              const {return CaloJets_[i];}
      const QCDEventHdr&   evtHdr()                    const {return EvtHdr_;}
 
    private:
      //---- event header (contains all the event info) --------------
      QCDEventHdr                              EvtHdr_;
      //---- CALO met object -----------------------------------------
      QCDMET                                   CaloMet_;
      //---- PF met object -------------------------------------------
      QCDMET                                   PFMet_; 
      //---- trigger decision vector --------------------------------- 
      std::vector<int>                         TriggerDecision_;
      //---- L1 prescale vector --------------------------------------
      std::vector<int>                         L1Prescale_;
      //---- HLT prescale vector -------------------------------------
      std::vector<int>                         HLTPrescale_;
      //---- HLT objects ---------------------------------------------  
      std::vector<std::vector<LorentzVector> > HLTObj_;
      //---- L1 objects ----------------------------------------------
      std::vector<std::vector<LorentzVector> > L1Obj_;
      //---- Genjets -------------------------------------------------
      std::vector<LorentzVector>               GenJets_;
      //---- CaloJets ------------------------------------------------ 
      std::vector<QCDCaloJet>                  CaloJets_;
      //---- PFJets --------------------------------------------------
      std::vector<QCDPFJet>                    PFJets_;
      //-------Tower----------------------------------------------------
      std::vector<QCDTower>                    CaloTowers_;
      //---- FatJets -------------------------------------------------
      std::vector<QCDJet>                      FatJets_;
};
#endif
