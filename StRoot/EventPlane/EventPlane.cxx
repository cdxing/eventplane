#include "EventPlane.h"
#include "StRoot/CutManager/CutManager.h"
#include "StRoot/ConstManager/ConstManager.h"
//#include "../HistManager/HistManager.h"
#include "StRoot/HistManager/HistManager.h" // Load STARLibrary header files
#include "StRoot/EpProManager/EpProManager.h" // Load STARLibrary header files
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/Run/run.h"
#include "StThreeVectorF.hh"
//IClass header files
#include "StRoot/IClasses/IEventPlane.h"
#include "StRoot/IClasses/IEvent.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
ClassImp(EventPlane)

    StRefMultCorr* EventPlane::mRefMultCorr = NULL;
    //-----------------------------------------------------------------------------
    EventPlane::EventPlane(const char* name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName)
: StMaker(name)
{
    configs.read(configFileName);
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mEnergy = 3;

    mOutPut_EP=Form("%s_recenterpar.root",jobid);
}

//----------------------------------------------------------------------------- 
EventPlane::~EventPlane()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t EventPlane::Init() 
{
    if(!mRefMultCorr)
    {
        mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    }
    // qapid
    mCutManager = new CutManager(configs);
    //mHistManager = new HistManager();
    //mHistManager->InitQAPID();

    // eventplane 
    // EPD
    mEpdGeom = new StEpdGeom();
    mEpProManager = new EpProManager();
    mEpProManager->InitEP();
    //const int numSubEvents = 3;
    //char subEventModes[_numSubEvents] = {'r','e','e'};

    mFile_EP = new TFile(mOutPut_EP.Data(),"RECREATE");
    mFile_EP->cd();
    return 0;

}

//----------------------------------------------------------------------------- 
Int_t EventPlane::Finish() 
{


    if(mOutPut_EP != "")
    {
        mFile_EP->cd();
        //mHistManager->WriteQAPID();
        mEpProManager->WriteEP();
        mFile_EP->Close();
    }
    return kStOK;
}

//----------------------------------------------------------------------------- 
void EventPlane::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t EventPlane::Make() 
{
    if(!mPicoDstMaker) 
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) 
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
    const Int_t nTracks = mPicoDst->numberOfTracks();
    Int_t TrkMult = 0;
    for(Int_t i = 0; i < nTracks; i++) // TrkMult
    {
        StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
        if(!track->isPrimary()) continue; // Only Primary Tracks
        //mHistManager->FillTrackQA(track,(TVector3)mPicoEvent->primaryVertex());
        //StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(mVertexPos);
	TrkMult ++;
    } // TrkMult

    // RefMult
    Int_t runId = mPicoEvent->runId();

    Int_t refMult = mPicoEvent->refMult();
    Float_t vz = mPicoEvent->primaryVertex().Z();
    //Float_t vx = mPicoEvent->primaryVertex().X();
    //Float_t vy = mPicoEvent->primaryVertex().Y();

    //Float_t vzvpd = mPicoEvent->vzVpd();
    Int_t TOF_Mul = mPicoEvent->btofTrayMultiplicity();
    //Int_t nMatchedToF = mPicoEvent->nBTOFMatch();
    Float_t zdcX = mPicoEvent->ZDCx();
    
    //mHistManager->FillEventQA(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
    //mHistManager->FillEventCut(0);

    // runIndex
    const int runIndex = GetRunIndex(runId);
    // event plane IClasses
    //Start with clean events
    IEvent * theEvent = new IEvent;   
    IEvent * theEvent_wt = new IEvent;   
    IEvent * theEvent_tpc = new IEvent;   
    IEvent subEvents[_numSubEvents];
    IEvent subEvents_wt[_numSubEvents];
    IEvent subEvents_tpc[_numSubEvents];
    theEvent->ClearEvent();
    theEvent_wt->ClearEvent();
    theEvent_tpc->ClearEvent();
    for (int i = 0; i < _numSubEvents; i++)
    {
    	subEvents[i].ClearEvent();
    	subEvents_wt[i].ClearEvent();
    	subEvents_tpc[i].ClearEvent();
    }	

    // Event Cut
    if(mCutManager->passEventCut(mPicoDst)) // event cut
    {
	
        //mHistManager->FillEventQaCut(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
        //mHistManager->FillEventCut(1);
        mRefMultCorr->init(runId);
        mRefMultCorr->initEvent(refMult, vz, zdcX);
        const Int_t cent16 = mRefMultCorr->getCentralityBin16();
        const Int_t cent9 = mRefMultCorr->getCentralityBin9();
        //const Double_t reweight = mRefMultCorr->getWeight();
	//std::cout << "refMult: " << refMult << std::endl;
	//std::cout << "TrkMult: " << TrkMult << std::endl;
        //const int cent16 = mCutManager->getCentrality(TrkMult);
	
	//std::cout << "cent16: " << cent16 << std::endl;
	//std::cout << "cent9: " << cent9 << std::endl;
        //const double reweight = 1.0;
        //if(cent16 >  15 || cent16 < 0) return 0;
        if(cent9 > 8  || cent9 < 0) return 0;
        //mHistManager->FillEventCent(cent9);
        //mHistManager->FillEventCut(2);
	
	/// qapid
        //const Int_t nToFMatched = mCutManager->getMatchedToF();
        TVector3 mVertexPos = mPicoDst->event()->primaryVertex();
        //float mField = mPicoEvent->bField();
        //Int_t N_pp = 0, N_pm = 0, N_kp = 0, N_km = 0, N_pr = 0;

	/// eventplane
	// TPC eventplane
        //double qx_TPC_A=0.0, qy_TPC_A=0.0;
        //double qx_TPC_B=0.0, qy_TPC_B=0.0;

        //cout << "nTracks = " << nTracks << endl;
        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
            StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
            if(!track->isPrimary()) continue; // Only Primary Tracks
            //mHistManager->FillTrackCut(0);
            if(!mCutManager->passTrackBasic(track)) continue;
            //mHistManager->FillTrackCut(1);
            //mHistManager->FillTrackPhysics(track );
		
            //StPicoPhysicalHelix helix = track->helix(mField);
            //Float_t dca = helix.geometricSignedDistance(mVertexPos);
            //Short_t  s_charge = track->charge();
            Float_t dca=track->gDCA(mVertexPos).Mag();
	    //if(mCutManager->isTofTrack(mPicoDst,track)) mHistManager->FillTrackTof(mPicoDst,track);
	    /*if(mCutManager->isProton(track))
	    {
	   	mHistManager->FillProton(mPicoDst,track,configs.y_mid); 
		N_pr++;
	    }
	    if(mCutManager->isKaon(mPicoDst,track))
	    {
	   	mHistManager->FillKaon(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_kp++;
		else N_km++;
	    }
	    if(mCutManager->isPion(mPicoDst,track))
	    {
	   	mHistManager->FillPion(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_pp++;
		else N_pm++;
	    }
	    */
	    //TPC EP
            if(!mCutManager->passTrackEP(track,dca)) continue;
	    Double_t eta = track->pMom().Eta();
	    Double_t phi = track->pMom().Phi();
	    Double_t pt = track->pMom().Perp();
	    IEventPlane eventPlane(phi, pt);
	    eventPlane.SetEta(eta);
	    theEvent_tpc->AddEPParticle(eventPlane);
            /*if(eta > -2.0 && eta < -1.1)    // sub event A eta(-2,-1.25)   // mid-rapidity is 1.06  -2.0-1.4, -1.05-0.65
            {
	   	mEpProManager->FillTpcAQvec(cent16,runIndex, track); 
            }

            if(eta > -1.0 && eta < 0.)    // sub event B eta(-1.15,0) // -1.3-0.685,  -0.65-0
            {
	   	mEpProManager->FillTpcBQvec(cent16,runIndex,track); 
            }*/
        } // track loop

	/// eventplane
	// EPD information
        Int_t nepdHits = mPicoDst->numberOfEpdHits();
        StPicoEpdHit *epdHit;
        TVector3 StraightLine_center;
        TVector3 StraightLine_random;
        double phi_epd_center = {0.0};
        double phi_epd_random = {0.0};
        double eta_epd_center = {0.0};
        double eta_epd_random = {0.0};
        double mip;
        double TileWeight           = {0};

        //double Qx_EPD[7]={0};
        //double Qy_EPD[7]={0};
        //double weight_EPD[7]={0};
        for(Int_t iHit=0; iHit<nepdHits; iHit++){ // EPD loop
            epdHit = mPicoDst->epdHit(iHit);
            mip = epdHit->nMIP();
            int iring = epdHit->row() -1;//(1~16)-1 -> 0-15

            if( !epdHit) continue;
            if( !epdHit->isGood())
	    { 
            	std::cout << "note good epd hit "  << std::endl;
		    continue;
	    }
            int position = epdHit->position();
            int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner mose
            //std::cout << "ringgroup = " << ringgroup << std::endl;
            if(ringgroup == -1) continue;
            //std::cout << "epdHit id = " << epdHit->id() << std::endl;
            //if(epdHit->id() > 0 ) continue; // only East side for FXT 

            StraightLine_center = mEpdGeom->TileCenter(epdHit->id())        - mPicoEvent->primaryVertex();  // the collision is from midille to side of the TPC
            StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - mPicoEvent->primaryVertex();
            //StraightLine_center = mEpdGeom->TileCenter(epdHit->id())       ;
            //StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id());

            phi_epd_center = StraightLine_center.Phi();
            eta_epd_center = StraightLine_center.Eta();
            phi_epd_random = StraightLine_random.Phi();
            eta_epd_random = StraightLine_random.Eta();
	    mEpProManager->FillEpdQa(iring, phi_epd_center, eta_epd_center, phi_epd_random, eta_epd_random);
            if(mip < configs.epd_threshold) continue;
            TileWeight = (mip > configs.epd_max_weight) ? configs.epd_max_weight : mip;
            //std::cout << "tileweight  = " << TileWeight << std::endl;

	    mEpProManager->FillEpdMip(eta_epd_center,position,TileWeight);
	    Double_t eta_wt = mEpProManager->GetEtaWeight(cent9,eta_epd_center);
	    //std::cout<< "centrality : "<< cent9 << "eta: " << eta_epd_center << "eta weight: " << eta_wt << std::endl;

	    IEventPlane eventPlane(phi_epd_center, TileWeight);
	    IEventPlane eventPlane_wt(phi_epd_center, TileWeight*eta_wt);
	    //eventPlane.SetTileID(epdHit->id());
	    //std::cout << eta_epd_center << std::endl;
	    eventPlane.SetEta(eta_epd_center);
	    eventPlane_wt.SetEta(eta_epd_center);
	    theEvent->AddEPParticle(eventPlane);
	    theEvent_wt->AddEPParticle(eventPlane_wt);
            /*for(int i=0; i<4; i++){
                if(i == ringgroup){
                    Qx_EPD[i] += TileWeight * cos(1.0*phi_epd_center);
                    Qy_EPD[i] += TileWeight * sin(1.0*phi_epd_center);
                    weight_EPD[i] += TileWeight;
                }
            }
            if(ringgroup == 0 || ringgroup == 1)
            {
                Qx_EPD[4] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[4] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[4] += TileWeight;
            }
            if(ringgroup == 2 || ringgroup == 3)
            {
                Qx_EPD[5] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[5] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[5] += TileWeight;
            }
            if(ringgroup == 0 || ringgroup == 1 || ringgroup == 2 || ringgroup == 3)
            {
                Qx_EPD[6] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[6] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[6] += TileWeight;
            }*/
        } // EPD loop
    	/*for(int i=0; i<7; i++){
    	    if(weight_EPD[i] == 0. || Qx_EPD[i]== 0.0 || Qy_EPD[i] == 0.0) return 0;
    	    Qx_EPD[i] /= weight_EPD[i];
    	    Qy_EPD[i] /= weight_EPD[i];
    	}
    	for(int i=0; i<7; i++){
	    mEpProManager->FillEpdQvec(i,cent16,runIndex,Qx_EPD[i],Qy_EPD[i]);
	}*/
	//double psi_ievt = theEvent->GetEventPsi(1);
	//std::cout << "ievent psi = " << psi_ievt << std::endl;
	//mHistManager->FillPIDMult(N_pp , N_pm , N_kp , N_km , N_pr ); /// qapid
	for (int sub = 0; sub < _numSubEvents; sub++)
	{
	  //subEvents[sub] = theEvent->GetSubEvent(subEventModes[sub], subEventParams[sub][0] - COMrapidity, subEventParams[sub][1] - COMrapidity);
	  subEvents[sub] = theEvent->GetSubEvent(_subEventModes[sub], _subEventParams[sub][0], _subEventParams[sub][1]);
	  subEvents_wt[sub] = theEvent_wt->GetSubEvent(_subEventModes[sub], _subEventParams[sub][0], _subEventParams[sub][1]);
	  subEvents_tpc[sub] = theEvent_tpc->GetSubEvent(_subEventModes[sub], _subEventParams_tpc[sub][0], _subEventParams_tpc[sub][1]);
	}
	//std::cout << "Size of wt TPC-East sub EP = " << subEvents_tpc[0].GetEPParticles().size() << std::endl;
	//std::cout << "Size of wt TPC-West sub EP = " << subEvents_tpc[1].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  EPD-C sub EP = " << subEvents[0].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  TPC-A sub EP = " << subEvents[1].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  TPC-B sub EP = " << subEvents[2].GetEPParticles().size() << std::endl << std::endl;
	if (subEvents[0].GetEPParticles().size() < 2
	//|| subEvents[1].GetEPParticles().size() < 2
	|| subEvents[1].GetEPParticles().size() < 2)
		{return 0;}	
	
	float subQx[_numSubEvents];
	float subQy[_numSubEvents];
	float subQx_wt[_numSubEvents];
	float subQy_wt[_numSubEvents];
	float subQx_tpc[_numSubEvents];
	float subQy_tpc[_numSubEvents];
	mEpProManager->FillPsiRawSubs(subEvents[0].GetEventPsi(1),subEvents[1].GetEventPsi(1));
	mEpProManager->FillPsiRawSubs_wt(subEvents_wt[0].GetEventPsi(1),subEvents_wt[1].GetEventPsi(1));
	mEpProManager->FillPsiRawSubs_tpc(subEvents_tpc[0].GetEventPsi(2),subEvents_tpc[1].GetEventPsi(2));
	mEpProManager->FillPsiRawFull(theEvent->GetEventPsi(1));
	mEpProManager->FillPsiRawFull_wt(theEvent_wt->GetEventPsi(1));
	mEpProManager->FillPsiRawFull_tpc(theEvent_tpc->GetEventPsi(2));
	for (int sub = 0; sub < _numSubEvents; sub++) // raw event plane, store recenter parameter
	{
	  mEpProManager->FillPsiRaw(sub, subEvents[sub].GetEventPsi(1));
	  mEpProManager->FillPsiRaw_wt(sub, subEvents_wt[sub].GetEventPsi(1));
	  mEpProManager->FillPsiRaw_tpc(sub, subEvents_wt[sub].GetEventPsi(2));
     	  subQx[sub] = subEvents[sub].GetQx(1);
     	  subQy[sub] = subEvents[sub].GetQy(1);
     	  subQx_wt[sub] = subEvents_wt[sub].GetQx(1);
     	  subQy_wt[sub] = subEvents_wt[sub].GetQy(1);
     	  subQx_wt[sub] = subEvents_tpc[sub].GetQx(2);
     	  subQy_wt[sub] = subEvents_tpc[sub].GetQy(2);
	  mEpProManager->FillSubEpQvec(sub, cent16, runIndex, subQx[sub], subQy[sub]);
	  mEpProManager->FillSubEpQvec_wt(sub, cent16, runIndex, subQx_wt[sub], subQy_wt[sub]);
	  mEpProManager->FillSubEpQvec_tpc(sub, cent16, runIndex, subQx_tpc[sub], subQy_tpc[sub]);
	
	} // raw event plane, store recenter parameter
	
	

    } // event cut

    return kStOK;
}

int EventPlane::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<nrun; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!" << endl;
    return runIndex;
}

/*int EventPlane::Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={4, 9,17,30,50,78, 116,170,205};
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}*/
