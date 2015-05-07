#ifndef QCDJet_h
#define QCDJet_h
#include "DataFormats/JetReco/interface/Jet.h"
//-------- Generic Jet class for QCD analyses ---------------
class QCDJet 
{
   public:
     typedef reco::Particle::LorentzVector LorentzVector;
     //------------ Constructor ------------------------------
     QCDJet() {}
     //------------ Destructor -------------------------------
     ~QCDJet() {}
     //------------ Sett methods -----------------------------
     void setP4(LorentzVector fP4) {P4_ = fP4;}
     void setGen(LorentzVector fP4, float fgenR) {genP4_ = fP4;genR_ = fgenR;}
     void setCor(float fCor)                     {cor_  = fCor;} 
     void setUnc(float fUnc)                     {unc_  = fUnc;} 
     void setUncSrc(std::vector<float> fUncSrc)  {uncSrc_ = fUncSrc;} 
     void setJecLabels(std::vector<double> fJecLabels) {jecLabels_ = fJecLabels;}

    //void setTower(std::vector<float> ftwr_Et,std::vector<float> ftwr_Eta,std::vector<float> ftwr_Phi)  {twr_et_ =ftwr_Et;
    //twr_eta_ = ftwr_Eta;twr_phi_=ftwr_Phi;    } 
    
     void setArea(float fArea)                   {area_ = fArea;}
     void setLooseID(bool fLooseID)              {looseID_ = fLooseID;} 
     void setTightID(bool fTightID)              {tightID_ = fTightID;}
     //------------ Get methods ------------------------------
     const LorentzVector& p4()    const {return P4_;}
     const LorentzVector& genp4() const {return genP4_;}
     float pt()                   const {return P4_.pt();}
     float genpt()                const {return genP4_.pt();}
     float geneta()               const {return genP4_.eta();} 
     float genR()                 const {return genR_;} 
     float ptCor()                const {return cor_ * P4_.pt();}
     float e()                    const {return P4_.energy();}
     float eCor()                 const {return cor_ * P4_.energy();}
     float eta()                  const {return P4_.eta();}
     float y()                    const {return P4_.Rapidity();}
     float phi()                  const {return P4_.phi();}
     float mass()                 const {return P4_.mass();}
     float cor()                  const {return cor_;}
     float unc()                  const {return unc_;} 
     float uncSrc(int i)          const {return uncSrc_[i];}
    //float twr_et(int i)          const {return twr_et_[i];}
    // float twr_eta(int i)          const {return twr_eta_[i];}
    // float twr_phi(int i)          const {return twr_phi_[i];}
     float area()                 const {return area_;} 
     bool  looseID()              const {return looseID_;}
     bool  tightID()              const {return tightID_;}

   private:
     //------ jet 4-momentum vector------------------
     LorentzVector P4_;
     //------ matched genjet 4-momentum vector-------
     LorentzVector genP4_;
     //------ matching radius -----------------------
     float genR_;
     //------ jec factor ----------------------------
     float cor_;
     //------ jec uncertainty -----------------------
     float unc_;
     //------ jec uncertainty sources ---------------
     std::vector<float> uncSrc_;
    // ----- All components of JEC Factor ----------                                                                                                 
    std::vector<double> jecLabels_;

    //----------calo tower---------------------------
    //std::vector<float> twr_et_;
    //std::vector<float> twr_eta_;
    //std::vector<float> twr_phi_;
     //------ jet area ------------------------------
     float area_;
     //------ loose ID flag -------------------------
     bool  looseID_;
     //------ tight ID flag -------------------------
     bool  tightID_;
};
#endif
