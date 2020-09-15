#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h>
#include <AsgTools/AnaToolHandle.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include "TMath.h"
#include <vector>
#include <EventLoop/Algorithm.h>
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>

#define EL_RETURN_CHECK( CONTEXT, EXP )                   \
    do                                                    \
    {                                                     \
        if( ! EXP.isSuccess() )                           \
        {            pk                                     \
            Error( CONTEXT,                               \
                XAOD_MESSAGE( "Failed to execute: %s" ),  \
                #EXP );                                   \
            return EL::StatusCode::FAILURE;               \
       }                                                  \
    } while( false )

//tools:

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "PathResolver/PathResolver.h"

#include <JetInterface/IJetSelector.h>
#include <JetResolution/IJERTool.h>
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetCalibTools/IJetCalibrationTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"


#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"
//#include "MuonEfficiencyCorrections/IMuonEfficiencyScaleFactors.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "AsgTools/ToolHandle.h"
#include "AsgTools/AnaToolHandle.h"

#include "IsolationSelection/IsolationSelectionTool.h"
#include "AsgAnalysisInterfaces/IPileupReweightingTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TriggerMatchingTool/MatchingTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"

using namespace std;
using namespace CP;
using namespace xAOD;
using namespace InDet;

using namespace Trig;
using namespace TrigConf;

class MyxAODAnalysis : public EL::AnaAlgorithm
{
public:
    // this is a standard algorithm constructor
    MyxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);
    // these are the functions inherited from Algorithm
    virtual StatusCode initialize () override;
    virtual StatusCode execute () override;
    virtual StatusCode finalize () override;

private:
	std::string _outputName; //!
	//asg::AnaToolHandle<IGoodRunsListSelectionTool> m_grl; //!
	GoodRunsListSelectionTool *m_grl_1; //!
    GoodRunsListSelectionTool *m_grl_2; //!
	GoodRunsListSelectionTool *m_grl_3; //!
	GoodRunsListSelectionTool *m_grl_4; //!

	// trigger tools
	Trig::TrigDecisionTool *m_trigDecisionTool; //!
    TrigConf::xAODConfigTool *m_trigConfigTool; //!
	Trig::MatchingTool *m_trigMatchingTool; //!
	
	// track selection tool
	InDet::InDetTrackSelectionTool *m_trackSelection[2]; //!
	// jet tools
	asg::AnaToolHandle<IJetSelector> m_jetCleaning_1; //! 
	asg::AnaToolHandle<IJetSelector> m_jetCleaning_2; //! 
	asg::AnaToolHandle<IJERTool> m_JERTool; //!
	asg::AnaToolHandle<IJetCalibrationTool> m_JetCalibration; //!
	//JetCalibrationTool *m_JetCalibration; //! 
   
	JetVertexTaggerTool* pjvtag = 0; //!
    ToolHandle<IJetUpdateJvt> hjvtagup; //!
	// muon tools
    CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //!
    CP::MuonSelectionTool *m_muonSelection[3]; //!
    IsolationSelectionTool *m_IsoSelection[2]; //!

	/////////////////////
	CP::MuonEfficiencyScaleFactors *m_effi_corr[3]; //!
	// systematics tools
    CP::SystematicSet statup; //!
    CP::SystematicSet statdown; //!
    CP::SystematicSet sysup; //!
    CP::SystematicSet sysdown; //!
    CP::SystematicSet loptstatup; //!
    CP::SystematicSet loptstatdown; //!
    CP::SystematicSet loptsysup; //!
    CP::SystematicSet loptsysdown; //!

    // PRW tool
    asg::AnaToolHandle<CP::IPileupReweightingTool> m_PRWTool; //!


	int m_numCleanEvents; //! 
    int m_eventCounter; //!
    int m_passeventCounter; //! 
	const float mumass = 0.1057; //!
	std::string outputName; //!
	TTree *pairs; //!
	float RunNumber; //!
    double EventNumber; //!
    float bcid; //!
    float lbn; //!
    float averageIntPerXing; //! 
    float actualIntPerXing; //!

	
	float PileupWeight;//!
    float actualmufromtool; //!
	float averagemufromtool; //!
	float ScaledActualMuFromTool; //!
	float ScaledAverageMuFromTool; //!

    float utps66; //!
    float utps64; //!
    float utps44; //!
	float utps44_tag; //!
	float utps66_tag;//!
	float utps10_10;//!
	float utps10_6;//!
	float utps11_6;//!
	float utps64_tag;//!
	float utps11_6_tag;//!

	float zbt20;//!
	float zbt26;//!
	float zbt50;//!
	float zbt10_10;//!
	float zbt_14_14;//!
    
    float Noftracks; //!
    float Nofbkgs; //!
    float Nvtxs; //!
    float primary_vx_z; //!
    float truth_primary_vx_z; //!
    float closest_pup; //!

	std::vector<float> *vxs = 0; //!
    std::vector<int> *vxtype = 0; //!
	std::vector<float> *vxsumpt = 0; //!

    std::vector<float> *recoeff1 = 0; //! 
    std::vector<float> *recoeff2 = 0; //!
    std::vector<float> *recotosub1 = 0; //!
    std::vector<float> *recotosub2 = 0; //!
    std::vector<float> *recopt1 = 0; //!
    std::vector<float> *recopt2 = 0; //!
    std::vector<float> *recod01 = 0; //!
    std::vector<float> *recod02 = 0; //!
    std::vector<float> *recoiso1 = 0; //!
    std::vector<float> *recoiso2 = 0; //!
    std::vector<float> *recod0s1 = 0; //!
    std::vector<float> *recod0s2 = 0; //!
    std::vector<float> *recotrig1 = 0; //!
    std::vector<float> *recotrig2 = 0; //!
    std::vector<float> *recoutrig1 = 0; //!
    std::vector<float> *recoutrig2 = 0; //!
    std::vector<float> *recoeta1 = 0; //!
    std::vector<float> *recoeta2 = 0; //!
    std::vector<float> *recophi1 = 0; //!
    std::vector<float> *recophi2 = 0; //!
    std::vector<float> *recoch1 = 0; //!
    std::vector<float> *recoch2 = 0; //!
    std::vector<float> *pairpt = 0; //!
    std::vector<float> *pairy = 0; //!
    std::vector<float> *paireta = 0; //!
    std::vector<float> *pairm = 0; //!
    std::vector<float> *pairphi = 0; //!

    std::vector<std::vector<float>> *recoerr1 = 0; //!
    std::vector<std::vector<float>> *recoerr2 = 0; //!



    std::vector<float> *jet4_pt = 0; //! 
    std::vector<float> *jet4_uncalibpt = 0; //! 
    std::vector<float> *jet4_eta = 0; //!
    std::vector<float> *jet4_phi = 0; //!
    std::vector<float> *jet4_ntrks = 0; //!
    std::vector<float> *jet4_ntrks5 = 0; //!
    std::vector<float> *jet4_ntrks6 = 0; //!
    std::vector<float> *jet4_jvf = 0; //!
	std::vector<int> *jet_cleaning_1 = 0;
	std::vector<int> *jet_cleaning_2 = 0;

    std::vector<float> *trks_true = 0; //!
    std::vector<float> *trks_eta = 0; //!
    std::vector<float> *trks_phi = 0; //!
    std::vector<float> *trks_pt = 0; //!
    std::vector<float> *trks_sinthetaz0 = 0; //! 
    std::vector<float> *trks_d0sig = 0; 



    //truth stuff
    float Upspt;
    float Upsy; 
    float Upsm; 
    float Upsphi; 
    float UpsID; 
    float UpsParentID;

    float truemupt1; 
    float truemupt2; 
    float truemueta1;
    float truemueta2;
    float truemuphi1;
    float truemuphi2;
    float truemuch1; 
    float truemuch2; 

    std::vector<float> *truesib_recoq = 0; 
    std::vector<float> *truesib_recoz0 = 0;
    std::vector<float> *truesib_eta = 0;
    std::vector<float> *truesib_phi = 0;
    std::vector<float> *truesib_pt = 0; 
    std::vector<float> *truesib_ID = 0; 

    std::vector<float> *truetrk_recoq = 0; 
    std::vector<float> *truetrk_recoz0 = 0;
    std::vector<float> *truetrk_eta = 0;
    std::vector<float> *truetrk_phi = 0;
    std::vector<float> *truetrk_pt = 0; 
    std::vector<float> *truetrk_charge = 0; 
    std::vector<float> *truetrk_vxz = 0; 
    std::vector<float> *truetrk_vxx = 0; 
    std::vector<float> *truetrk_vxy = 0; 

	std::vector<int> *truetrk_barcode = 0; 
	std::vector<int> *truetrk_pdg = 0; 
	std::vector<int> *truetrk_status = 0; 

	std::vector<const xAOD::IParticle*> myParticles;

	// track variables (reconstructed matched to truth)
	std::vector<float> *m_trks_eta_mc = 0; 
	std::vector<float> *m_trks_phi_mc = 0; 
	std::vector<float> *m_trks_pt_mc = 0; 
	std::vector<float> *m_trks_mcprob = 0; 
	std::vector<int> *m_trks_charge_mc = 0;
	std::vector<int> *m_trks_barcode = 0; 
	std::vector<int> *m_trks_pdg = 0; 
	std::vector<int> *m_trks_status = 0;

/*
	// Histograms
	TH2F *h_eta_SCTHits;
	TH2F *h_eta_SCTOutliers;
	TH2F *h_eta_SCTHoles;
	TH2F *h_eta_SCTDoubleHoles;
	TH2F *h_eta_SCTSharedHits;
	TH2F *h_eta_SCTDeadSensors;
	
	TH2F *h_eta_PixelHits;
	TH2F *h_eta_PixelOutliers;
	TH2F *h_eta_PixelHoles;
	TH2F *h_eta_PixelSharedHits;
	TH2F *h_eta_PixelSplitHits;
	TH2F *h_eta_PixelDeadSensors;

	TH2F *h_eta_TRTHits;
	TH2F *h_eta_TRTOutliers;
	TH2F *h_eta_TRTHoles;
	TH2F *h_eta_TRTSharedHits;

	TH2F *h_eta_SiHits;
	TH2F *h_eta_SiHole;

	TH2F *h_eta_IBLHits;
	TH2F *h_eta_BLHits;
	TH2F *h_eta_expIBLHits;
	TH2F *h_eta_expBLHits;*/

	std::vector<bool> *match_utps66_mu1 = 0;
	std::vector<bool> *match_utps64_mu1 = 0;
	std::vector<bool> *match_utps44_mu1 = 0;
	std::vector<bool> *match_utps44_tag_mu1 = 0;
	std::vector<bool> *match_utps66_tag_mu1 = 0;
	std::vector<bool> *match_utps10_10_mu1 = 0;
	std::vector<bool> *match_utps10_6_mu1 = 0;
	std::vector<bool> *match_utps11_6_mu1 = 0;
	std::vector<bool> *match_utps64_tag_mu1 = 0;
	std::vector<bool> *match_utps11_6_tag_mu1 = 0;

	std::vector<bool> *match_utps66_mu2 = 0;
	std::vector<bool> *match_utps64_mu2 = 0;
	std::vector<bool> *match_utps44_mu2 = 0;
	std::vector<bool> *match_utps44_tag_mu2 = 0;
	std::vector<bool> *match_utps66_tag_mu2 = 0;
	std::vector<bool> *match_utps10_10_mu2 = 0;
	std::vector<bool> *match_utps10_6_mu2 = 0;
	std::vector<bool> *match_utps11_6_mu2 = 0;
	std::vector<bool> *match_utps64_tag_mu2 = 0;
	std::vector<bool> *match_utps11_6_tag_mu2 = 0;

    std::vector<int> *mu1_loose = 0;
    std::vector<int> *mu1_medium = 0;
    std::vector<int> *mu1_tight = 0;

    std::vector<int> *mu2_loose = 0;
    std::vector<int> *mu2_medium = 0;
    std::vector<int> *mu2_tight = 0; 


	double getdR(double eta1,double phi1,double eta2,double phi2)
    {
        double dphi = phi1-phi2;
        if(dphi>TMath::Pi())dphi=dphi-2*TMath::Pi();
        else if(dphi<-TMath::Pi())dphi=dphi+2*TMath::Pi();
        double deta = eta1-eta2;
        return sqrt(deta*deta + dphi*dphi);
    }

    int getzposbin(double zpos)
    {
        if(zpos<-150) return 0;
        if(zpos>150) return 11;
        int bin=1;
        for(int i=0;i<10;i++) if(zpos > (i*30  -120) ) bin++; 
        return bin;
    }

};

#endif
