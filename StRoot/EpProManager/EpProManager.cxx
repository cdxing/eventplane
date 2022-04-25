#include "EpProManager.h"
#include "../ConstManager/ConstManager.h"
//#include "StTriFlowConstants.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"

ClassImp(EpProManager)

    TString EpProManager::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

EpProManager::EpProManager()
{
}

//---------------------------------------------------------------------------------

EpProManager::~EpProManager()
{
}


//----------------------------------------------------------------------------

void EpProManager::InitEP()
{
    h2_RingIDvsEta_center = new TH2D("ringIDvsEta_center","rungIDvsEta_center",16,0,16,40,-6,-2);
    h2_RingIDvsEta_random = new TH2D("ringIDvsEta_random","rungIDvsEta_random",16,0,16,40,-6,-2);
    h_phi_center = new TH1D("h_phi_center","h_phi_center",630,-3.15,3.15);
    h_eta_center = new TH1D("h_eta_center","h_eta_center",300,-6,6);
    h_phi_random = new TH1D("h_phi_random","h_phi_random",630,-3.15,3.15);
    h_eta_random = new TH1D("h_eta_random","h_eta_random",300,-6,6);

    p_sector_eta=new TProfile2D("p_sector_eta", "", 300, -7.0, -2.0, 12, 0.5, 12.5);
    // recenter parameter for EPD EP
    /*for(int i=0; i<7; i++){
        p_EPDQx[i] = new TProfile2D(Form("EPDQx_ring%d_recen",i), "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
        p_EPDQy[i] = new TProfile2D(Form("EPDQy_ring%d_recen",i), "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    }
    // recenter parameter for TPC EP
    p_TPCqx_A = new TProfile2D("TPCqx_A_recen", "",  ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,500, -0.5, 499.5);
    p_TPCqy_A = new TProfile2D("TPCqy_A_recen", "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    p_TPCqx_B = new TProfile2D("TPCqx_B_recen", "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    p_TPCqy_B = new TProfile2D("TPCqy_B_recen", "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);*/
    

    for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
    {
    	//
    	float param1 = _subEventParams[sub][0];// - COMrapidity;
    	float param2 = _subEventParams[sub][1];// - COMrapidity;
    	
    	TString form = Form(" - #sub event #in %c (%1.1f, %1.1f)",_subEventModes[sub],param1 ,param2);
    	
    	hPsiRawSub[sub] = new TH1F(Form("hPsiRawSub_{%d}",sub + 1), "#Psi distribution (raw)" + form,1024,-7.0,7.0);
    	hPsiRecenteredSub[sub] = new TH1F(Form("hPsiRecenteredSub_{%d}",sub + 1), "#Psi distribution (recentered)" + form,1024,-7.0,7.0);	
    	hPsiShiftedSub[sub] = new TH1F(Form("hPsiShiftedSub_{%d}",sub + 1), "#Psi distribution (shifted)" + form,1024,-7.0,7.0);		
    	
    	hQxSub[sub] = new TH1F(Form("hQxSub_{%d}",sub + 1), "Q_{x} distribution",10000,-100,100);	
    	hQySub[sub] = new TH1F(Form("hQySub_{%d}",sub + 1), "Q_{x} distribution",10000,-100,100);	
    	hQVectorSub[sub] = new TH2F(Form("hqVectorSub_{%d}",sub + 1),"#vec{Q} distribution",200,-100,100,200,-100,100);
    	p_mq1x_Sub_EP[sub] = new TProfile2D(Form("p_mq1x_Sub_EP_%d",sub + 1),"mQx1 for recenter correction Qx" + form,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    	p_mq1y_Sub_EP[sub] = new TProfile2D(Form("p_mq1y_Sub_EP_%d",sub + 1),"mQy1 for recenter correction Qx" + form,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    } // event plane Psi histograms
    cout << "End of event plane Histograms" << endl;
}
/*void EpProManager::InitEP()
{
    //===============================                          
    //  Define Histograms                           
    //===============================                                            

    h2_RingIDvsEta_center = new TH2D("ringIDvsEta_center","rungIDvsEta_center",16,0,16,40,-6,-2);
    h2_RingIDvsEta_random = new TH2D("ringIDvsEta_random","rungIDvsEta_random",16,0,16,40,-6,-2);
    h_phi_center = new TH1D("h_phi_center","h_phi_center",630,-3.15,3.15);
    h_eta_center = new TH1D("h_eta_center","h_eta_center",300,-6,6);
    h_phi_random = new TH1D("h_phi_random","h_phi_random",630,-3.15,3.15);
    h_eta_random = new TH1D("h_eta_random","h_eta_random",300,-6,6);


    for(int i=0; i<9; i++){
        h_eta[i] = new TH1F(Form("eta_%d",i),Form("eta_%d",i),300,-2.5,0.5);
        h2_eta_pt[i] = new TH2D(Form("eta_pt_%d",i),Form("eta_pt_%d",i),300,-2.5,0.5,200,0.0,5.0);
    }


    cout << "End of vent plane Histograms" << endl;
 
    for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            for(Int_t k = 0; k < 5; k++) // Shift Order
            {
                TString ProName;
                // Event Plane method
                ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_mcos2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_msin2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_mcos2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_msin2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

                ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_mcos3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_msin3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_mcos3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_msin3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

                // Scalor Product method
            }
        }

        for(Int_t k = 0; k < 5; k++) // Shift Order
        {
            TString ProName_Full;
            // Event Plane method
            ProName_Full = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_mcos2_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
            ProName_Full = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_msin2_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

            ProName_Full = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_mcos3_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
            ProName_Full = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_msin3_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

        }
    }

}
*/

//----------------------------------------------------------------------------
// Event Plane method
void EpProManager::FillEpdQa(Int_t iring, Double_t phi_epd_center, Double_t eta_epd_center,Double_t phi_epd_random, Double_t eta_epd_random)
{
  h_phi_center -> Fill(phi_epd_center);
  h_eta_center -> Fill(eta_epd_center);
  h_phi_random -> Fill(phi_epd_random);
  h_eta_random -> Fill(eta_epd_random);
  h2_RingIDvsEta_center -> Fill(iring,eta_epd_center);
  h2_RingIDvsEta_random -> Fill(iring,eta_epd_random);

}
void EpProManager::FillEpdMip(Double_t eta_epd_center, Int_t position, Double_t mip)
{
  p_sector_eta->Fill(eta_epd_center, position, mip);
}
void EpProManager::FillEpdQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_EPDQx[isub]->Fill(centnumber, runindex, qx);
  p_EPDQy[isub]->Fill(centnumber, runindex, qy);
}

void EpProManager::FillTpcAQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track)
{
  Double_t d_pT = track->pPt();
  Double_t d_phi = track->pMom().Phi();
  Double_t qx_TPC_A = d_pT*TMath::Cos(1.0*d_phi);
  Double_t qy_TPC_A = d_pT*TMath::Sin(1.0*d_phi);
  p_TPCqx_A->Fill(centnumber, runindex, qx_TPC_A);
  p_TPCqy_A->Fill(centnumber, runindex, qy_TPC_A);

}
void EpProManager::FillTpcBQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track)
{
  Double_t d_pT = track->pPt();
  Double_t d_phi = track->pMom().Phi();
  Double_t qx_TPC_A = d_pT*TMath::Cos(1.0*d_phi);
  Double_t qy_TPC_A = d_pT*TMath::Sin(1.0*d_phi);
  p_TPCqx_B->Fill(centnumber, runindex, qx_TPC_A);
  p_TPCqy_B->Fill(centnumber, runindex, qy_TPC_A);

}
void EpProManager::FillSubEpQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq1x_Sub_EP[isub]->Fill(centnumber, runindex, qx);
  p_mq1y_Sub_EP[isub]->Fill(centnumber, runindex, qy);
  hQxSub[isub]->Fill(qx);
  hQySub[isub]->Fill(qy);
  hQVectorSub[isub]->Fill(qx,qy);
}

    
void EpProManager::FillPsiRaw(Int_t isub, Double_t psi)
{
  hPsiRawSub[isub]->Fill(psi);
}


/*void EpProManager::FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}*/

/*void EpProManager::FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}*/

/*void EpProManager::FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}*/

//----------------------------------------------------------------------------

void EpProManager::WriteEP()
{
  h_phi_center -> Write();
  h_eta_center -> Write();
  h_phi_random -> Write();
  h_eta_random -> Write();
  h2_RingIDvsEta_center -> Write();
  h2_RingIDvsEta_random -> Write();
  p_sector_eta -> Write();
  /*for(int i=0; i<7; i++){
      p_EPDQx[i]->Write();
      p_EPDQy[i]->Write();
  }
  p_TPCqx_A->Write();
  p_TPCqy_A->Write();
  p_TPCqx_B->Write();
  p_TPCqy_B->Write();*/
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    hPsiRawSub[sub] ->Write();
    p_mq1x_Sub_EP[sub]->Write();
    p_mq1y_Sub_EP[sub]->Write();
    hQxSub[sub]->Write();
    hQySub[sub]->Write();
    hQVectorSub[sub]->Write();
  }
    /*for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            for(Int_t k = 0; k < 5; k++) // Shift Order
            {
                // Event Plane method
                p_mcos2_East_EP[i][j][k]->Write();
                p_msin2_East_EP[i][j][k]->Write();
                p_mcos2_West_EP[i][j][k]->Write();
                p_msin2_West_EP[i][j][k]->Write();

                p_mcos3_East_EP[i][j][k]->Write();
                p_msin3_East_EP[i][j][k]->Write();
                p_mcos3_West_EP[i][j][k]->Write();
                p_msin3_West_EP[i][j][k]->Write();

            }
        }
        for(Int_t k = 0; k < 5; k++) // Shift Order
        {
            // Event Plane method
            p_mcos2_Full_EP[i][k]->Write();
            p_msin2_Full_EP[i][k]->Write();
            p_mcos3_Full_EP[i][k]->Write();
            p_msin3_Full_EP[i][k]->Write();

        }
    }*/
}

//----------------------------------------------------------------------------


