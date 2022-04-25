#ifndef EpProManager_h
#define EpProManager_h

#include "TVector2.h"
#include "TString.h"
#include "StMessMgr.h"
#include "../StPicoEvent/StPicoTrack.h"

class TProfile2D;
class TProfile;
class TH1F;
class TH1D;
class TH2F;
class TH2D;

// sub event plane configurations
const int _numSubEvents = 3;
const char _subEventModes[_numSubEvents] = {'r','e','e'};
const float _subEventParams[_numSubEvents][2] =
{{9,12}, //EPD-C ring-range
 {-2.0,-1.1}, //TPC-A psuedorapidity
 {-1.0,0}}; //TPC-B psuedorapidity

class EpProManager
{
  public:
    EpProManager();
    ~EpProManager();


    // Recenter and  Shift Correction
    void InitEP();
    void WriteEP();
    // Event Plane method
    // Fill plots
    void FillEpdQa(Int_t iring, Double_t phi_epd_center, Double_t eta_epd_center,Double_t phi_epd_random, Double_t eta_epd_random); 
    void FillEpdMip(Double_t eta_epd_center, Int_t position, Double_t mip); 
    void FillEpdQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy); 
    void FillTpcAQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track);
    void FillTpcBQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track);
    void FillSubEpQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    
    void FillPsiRaw(Int_t isub, Double_t psi); 
    /*void FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k); // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
    void FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k); // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
    void FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k); // i = vertex pos/neg, k = ShiftOrder*/

    // Scalor Product method

  private:

    //raw EventPlane distribution (before recenter and shift)

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    //TH1F* href;
    //TH1F* hcent;
  
    TH2D *h2_RingIDvsEta_center;
    TH2D *h2_RingIDvsEta_random;
  
    TH1D *h_phi_center;
    TH1D *h_eta_center;
    TH1D *h_phi_random;
    TH1D *h_eta_random;
  
    TProfile2D *p_sector_eta;
  
    TH1F *h_eta[9];
    TH2D *h2_eta_pt[9];
  
    TProfile2D *p_EPDQx[7];
    TProfile2D *p_EPDQy[7];
    TProfile2D *p_TPCqx_A;
    TProfile2D *p_TPCqy_A;
    TProfile2D *p_TPCqx_B;
    TProfile2D *p_TPCqy_B;
    TString mout_recenter;
    
    TH1F *hPsiRawSub[_numSubEvents];
    TH1F *hPsiRecenteredSub[_numSubEvents];	
    TH1F *hPsiShiftedSub[_numSubEvents];	

    TH1F *hQxSub[_numSubEvents];
    TH1F *hQySub[_numSubEvents];
    TH2F *hQVectorSub[_numSubEvents];
    // Event Plane method
    TProfile2D *p_mq1x_Sub_EP[_numSubEvents];
    TProfile2D *p_mq1y_Sub_EP[_numSubEvents];

    /*TProfile2D *p_mq2y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq2y_Full_EP[2];    // 0 = vertex pos/neg

    TProfile2D *p_mq3x_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq3y_Full_EP[2];    // 0 = vertex pos/neg*/

    // Shift Correction | x axis is RunIndex, y axis is Centrality
    // Event Plane method
    /*TProfile2D *p_mcos2_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin2_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos2_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin2_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos2_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin2_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder

    TProfile2D *p_mcos3_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin3_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos3_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin3_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos3_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin3_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder*/


    static TString mVStr[2];

    ClassDef(EpProManager,1)
};

#endif
