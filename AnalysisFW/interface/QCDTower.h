#ifndef QCDTower_h
#define QCDTower_h
#include "FSQJet/AnalysisFW/interface/QCDJet.h"                                                                                                                                                                                  
class QCDTower: public QCDJet  {
public:
     //------------ Constructor ------------------------------
     QCDTower() {twr_et_=0;twr_eta_=0;twr_phi_=0;}
     //------------ Destructor -------------------------------
     ~QCDTower() {}
     //------------ Set methods ------------------------------
    
     void setTower(float ftwr_et, float ftwr_eta, float ftwr_phi){twr_et_ = ftwr_et; twr_eta_ = ftwr_eta; twr_phi_ = ftwr_phi;}
    
     //------------ Get methods ------------------------------ 
     
     float twr_et()   const {return twr_et_; } 
     float twr_eta()  const {return twr_eta_;}
     float twr_phi()  const {return twr_phi_;} 
   private:
     //---- tower transverse energy ----
     float twr_et_;
     //---- tower eta ----
     float twr_eta_;
     //---- tower phi ------------
     float twr_phi_;
     
    };
#endif    
