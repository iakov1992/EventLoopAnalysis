#include <AsgTools/MessageCheck.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <xAODEventInfo/EventInfo.h>
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TLorentzVector.h>
#include <vector>
#include <cmath>
#include <stdio.h>
// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
//Physics includes
#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"

#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "PATInterfaces/CorrectionCode.h" // to check the return correction code status of tools
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"

#include "AsgTools/ToolHandle.h"
#include "AsgTools/AnaToolHandle.h"
#include "AsgAnalysisInterfaces/IPileupReweightingTool.h"

#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"

//truth

#include "xAODTruth/xAODTruthHelpers.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"

using namespace std;


/*
List of triggers:

	Upsilon Triggers

		data15_13TeV, physics_Main stream:
			HLT_2mu10_bUpsimumu
			HLT_2mu4_bUpsimumu
			HLT_2mu6_bUpsimumu
			HLT_mu6_mu4_bUpsimumu
		data16_13TeV, physics_Main stream:
			HLT_2mu10_bUpsimumu
			HLT_2mu4_bUpsimumu_L1BPH-1M19-2MU4-BO_BPH-0DR34-2MU4
			HLT_2mu6_bUpsimumu
			HLT_mu10_mu6_bUpsimumu
			HLT_mu6_mu4_bUpsimumu
		data17_13TeV, physics_Main stream:
			HLT_2mu6_bUpsimumu_L1BPH-8M15-2MU6_BPH-0DR22-2MU6
		data17_13TeV, physics_BphysLS stream:
			HLT_2mu10_bUpsimumu
			HLT_2mu6_bUpsimumu
			HLT_mu10_mu6_bUpsimumu
			HLT_mu11_mu6_bUpsimumu
			HLT_mu6_mu4_bUpsimumu
			HLT_mu6_mu4_bUpsimumu_L1BPH-8M15-MU6MU4_BPH-0DR22-MU6MU4-BO
		data18_13TeV, physics_Main stream:
			HLT_2mu6_bUpsimumu_L1BPH-8M15-2MU6_BPH-0DR22-2MU6
		data18_13TeV, physics_BphysLS stream:
			HLT_2mu10_bUpsimumu
			HLT_mu11_mu6_bUpsimumu
			HLT_mu11_mu6_bUpsimumu_L1LFV-MU11
			HLT_mu6_mu4_bUpsimumu_L1BPH-8M15-MU6MU4_BPH-0DR22-MU6MU4-BO
	Zboson triggers:
		HLT_mu26_ivarmedium
		HLT_mu20_iloose_L1MU15
		HLT_mu50
		HLT_2mu10
		HLT_2mu14
*/

//ClassImp(MyxAODAnalysis)
//m_grl ("GoodRunsListSelectionTool/grl",this), 

MyxAODAnalysis :: MyxAODAnalysis (const std::string& name,
                                  ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm (name, pSvcLocator), m_jetCleaning_1("JetCleaningTool/JetCleaning_0", this),m_jetCleaning_2("JetCleaningTool/JetCleaning_1", this),m_JERTool("JERTool", this)
{

}

StatusCode MyxAODAnalysis :: initialize ()
{

    xAOD::TEvent *event = wk()->xaodEvent();

	// tools
	// GRL
	const char *GRLFilePath_1 = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data15_13TeV/20190708/data15_13TeV.periodAllYear_DetStatus-v105-pro22-13_Unknown_PHYS_StandardGRL_All_Good_25ns.xml";
	const char *GRLFilePath_2 = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data16_13TeV/20190708/data16_13TeV.periodAllYear_DetStatus-v105-pro22-13_Unknown_PHYS_StandardGRL_All_Good_25ns_WITH_IGNORES.xml";
	const char *GRLFilePath_3 = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data17_13TeV/20190708/data17_13TeV.periodAllYear_DetStatus-v105-pro22-13_Unknown_PHYS_StandardGRL_All_Good_25ns_Triggerno17e33prim.xml";
	const char *GRLFilePath_4 = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data18_13TeV/20190708/data18_13TeV.periodAllYear_DetStatus-v105-pro22-13_Unknown_PHYS_StandardGRL_All_Good_25ns_Triggerno17e33prim.xml";


	m_grl_1 = new GoodRunsListSelectionTool("GoodRunsListSelectionTool_1");
	const char *fullGRLFilePath_1 = gSystem->ExpandPathName(GRLFilePath_1);
	std::vector<std::string> vecStringGRL_1;
	vecStringGRL_1.push_back(fullGRLFilePath_1);
	ANA_CHECK(m_grl_1->setProperty("GoodRunsListVec", vecStringGRL_1));
    ANA_CHECK(m_grl_1->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    ANA_CHECK(m_grl_1->initialize());



    m_grl_2 = new GoodRunsListSelectionTool("GoodRunsListSelectionTool_2");
	const char *fullGRLFilePath_2 = gSystem->ExpandPathName(GRLFilePath_2);
	std::vector<std::string> vecStringGRL_2;
	vecStringGRL_2.push_back(fullGRLFilePath_2);
	ANA_CHECK(m_grl_2->setProperty("GoodRunsListVec", vecStringGRL_2));
    ANA_CHECK(m_grl_2->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    ANA_CHECK(m_grl_2->initialize());


    m_grl_3 = new GoodRunsListSelectionTool("GoodRunsListSelectionTool_3");
	const char *fullGRLFilePath_3 = gSystem->ExpandPathName(GRLFilePath_3);
	std::vector<std::string> vecStringGRL_3;
	vecStringGRL_3.push_back(fullGRLFilePath_3);
	ANA_CHECK(m_grl_3->setProperty("GoodRunsListVec", vecStringGRL_3));
    ANA_CHECK(m_grl_3->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    ANA_CHECK(m_grl_3->initialize());

    m_grl_4 = new GoodRunsListSelectionTool("GoodRunsListSelectionTool_4");
	const char *fullGRLFilePath_4 = gSystem->ExpandPathName(GRLFilePath_4);
	std::vector<std::string> vecStringGRL_4;
	vecStringGRL_4.push_back(fullGRLFilePath_4);
	ANA_CHECK(m_grl_4->setProperty("GoodRunsListVec", vecStringGRL_4));
    ANA_CHECK(m_grl_4->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    ANA_CHECK(m_grl_4->initialize());


	//initialize (and configure) the track selection tool
    m_trackSelection[0] = new InDet::InDetTrackSelectionTool("TrackSelection_0");
    m_trackSelection[0]->msg().setLevel(MSG::ERROR);
    ANA_CHECK(m_trackSelection[0]->setProperty("CutLevel","Loose"));
    ANA_CHECK(m_trackSelection[0]->setProperty("minPt",400));
    ANA_CHECK(m_trackSelection[0]->initialize());

    m_trackSelection[1] = new InDet::InDetTrackSelectionTool("TrackSelection_1");
    m_trackSelection[1]->msg().setLevel(MSG::ERROR);
    ANA_CHECK(m_trackSelection[1]->setProperty("CutLevel","TightPrimary"));
    ANA_CHECK(m_trackSelection[1]->setProperty("minPt",400));
    ANA_CHECK(m_trackSelection[1]->initialize());

    m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
    ANA_CHECK( m_trigConfigTool->initialize() );
    ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
    m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
    ANA_CHECK(m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ) ); // connect the TrigDecisionTool to the ConfigTool
    ANA_CHECK(m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision" ) );
    ANA_CHECK(m_trigDecisionTool->initialize());
    ToolHandle< Trig::TrigDecisionTool > m_trigDecision(m_trigDecisionTool);

	m_trigMatchingTool = new Trig::MatchingTool("TrigMatchingTool");
	//m_trigMatchingTool->setProperty("TrigDecisionTool", m_trigDecisionTool);
	ANA_CHECK(m_trigMatchingTool->initialize());


	//////////////
	// initialize and configure the jet cleaning tool
    ANA_CHECK (m_jetCleaning_1.setProperty("CutLevel", "LooseBad"));
    ANA_CHECK (m_jetCleaning_1.setProperty("DoUgly", false));
    ANA_CHECK (m_jetCleaning_1.initialize());

    ANA_CHECK (m_jetCleaning_2.setProperty("CutLevel", "TightBad"));
    ANA_CHECK (m_jetCleaning_2.setProperty("DoUgly", false));
    ANA_CHECK (m_jetCleaning_2.initialize());

    // instantiate and initialize the JER (using default configurations)
    ANA_CHECK(m_JERTool.initialize());

    // initialize and configure the jet calibration tool
	m_JetCalibration.setTypeAndName("JetCalibrationTool/MyJetCalibrationTool");
	const std::string name = "MyxAODAnalysis"; //string describing the current thread, for logging
    TString jetAlgo = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see above)
	//TString config = "JES_MC15cRecommendation_May2016.config";
	TString config = "JES_data2017_2016_2015_Recommendation_Feb2018_rel21.config";
    TString calibSeq = "JetArea_Residual_EtaJES_GSC_Insitu" ;
	TString calibArea = "00-04-81"; // Calibration Area tag (see below)
    bool isData = true; //bool describing if the events are data or from simulation
	if(!m_JetCalibration.isUserConfigured())
	{
        ANA_CHECK(m_JetCalibration.setProperty("JetCollection",jetAlgo.Data()));
        ANA_CHECK(m_JetCalibration.setProperty("ConfigFile",config.Data()) );
        ANA_CHECK(m_JetCalibration.setProperty("CalibSequence",calibSeq.Data()));
        ANA_CHECK(m_JetCalibration.setProperty("CalibArea",calibArea.Data()));
        ANA_CHECK(m_JetCalibration.setProperty("IsData",isData));
        ANA_CHECK(m_JetCalibration.retrieve());
    }
  
    // Configure the JVT tool.
    pjvtag = new JetVertexTaggerTool("jvtag");
    hjvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
    ANA_CHECK(pjvtag->setProperty("JVTFileName","JetMomentTools/JVTlikelihood_20140805.root"));
    ANA_CHECK(pjvtag->initialize());


	//ANA_CHECK(m_PRWTool.make("CP::PileupReweightingTool/prw"));//from egamma example
	m_PRWTool.setTypeAndName("CP::PileupReweightingTool/prw");
    //std::vector<std::string> m_ConfigFiles { "MyAnalysis/mc16a.user.iaizenbe.2015-2016.pp13.allUpsi.root" };
    std::vector<std::string> m_ConfigFiles { "MyAnalysis/mc16d.user.iaizenbe.2017.pp13.allUpsi.root" };
    ANA_CHECK(m_PRWTool.setProperty( "ConfigFiles", m_ConfigFiles));
    std::vector<std::string> lcalcFiles;

    //const char* LumiFilePath = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data15_13TeV/20190708/ilumicalc_histograms_None_276262-284484_OflLumi-13TeV-010.root";
	//const char* LumiFilePath = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data16_13TeV/20190708/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-010.root";
	const char* LumiFilePath = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data17_13TeV/20190708/ilumicalc_histograms_None_325713-340453_OflLumi-13TeV-010.root";
	//const char* LumiFilePath = "/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data18_13TeV/20190708/ilumicalc_histograms_None_348885-364292_OflLumi-13TeV-010.root";


    const char* fullLumiFilePath = gSystem->ExpandPathName (LumiFilePath);
    lcalcFiles.push_back(fullLumiFilePath);
    ANA_CHECK(m_PRWTool.setProperty("LumiCalcFiles", lcalcFiles));
    ANA_CHECK(m_PRWTool.setProperty("DataScaleFactorUP", 1.0/1.0));
    ANA_CHECK(m_PRWTool.setProperty("DataScaleFactor", 1.0/1.16));
    ANA_CHECK(m_PRWTool.setProperty("DataScaleFactorDOWN", 1.0/1.23));
    ANA_CHECK(m_PRWTool.retrieve()); 

    // Muon efficiency corrections
	m_effi_corr[0] = new CP::MuonEfficiencyScaleFactors("MuonEffTool_0");
    ANA_CHECK(m_effi_corr[0]->setProperty( "WorkingPoint", "Loose"));
    ANA_CHECK(m_effi_corr[0]->initialize());

    m_effi_corr[1] = new CP::MuonEfficiencyScaleFactors("MuonEffTool_1");
    ANA_CHECK(m_effi_corr[1]->setProperty( "WorkingPoint", "Medium"));
    ANA_CHECK(m_effi_corr[1]->initialize());

    m_effi_corr[2] = new CP::MuonEfficiencyScaleFactors("MuonEffTool_2");
    ANA_CHECK(m_effi_corr[2]->setProperty( "WorkingPoint", "Tight"));
    ANA_CHECK(m_effi_corr[2]->initialize());

    statup.insert(CP::SystematicVariation("MUON_EFF_STAT", 1));
    statdown.insert(CP::SystematicVariation("MUON_EFF_STAT", -1));
    sysup.insert(CP::SystematicVariation("MUON_EFF_SYS", 1));
    sysdown.insert(CP::SystematicVariation("MUON_EFF_SYS", -1));

    loptstatup.insert(CP::SystematicVariation("MUON_EFF_STAT_LOWPT", 1));
    loptstatdown.insert(CP::SystematicVariation("MUON_EFF_STAT_LOWPT", -1));
    loptsysup.insert(CP::SystematicVariation("MUON_EFF_SYS_LOWPT", 1));
    loptsysdown.insert(CP::SystematicVariation("MUON_EFF_SYS_LOWPT", -1));

    // initialize the muon calibration and smearing tool
    m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" );
    //m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::DEBUG );
    ANA_CHECK(m_muonCalibrationAndSmearingTool->initialize());


    m_muonSelection[0] = new CP::MuonSelectionTool("MuonSelectionTool_0");
	m_muonSelection[0]->msg().setLevel( MSG::ERROR ); 
	ANA_CHECK(m_muonSelection[0]->setProperty( "MaxEta", 2.5 )); 
	ANA_CHECK(m_muonSelection[0]->setProperty( "MuQuality", 2)); //Loose
	ANA_CHECK(m_muonSelection[0]->initialize());

    m_muonSelection[1] = new CP::MuonSelectionTool("MuonSelectionTool_1");
	m_muonSelection[1]->msg().setLevel( MSG::ERROR ); 
	ANA_CHECK(m_muonSelection[1]->setProperty( "MaxEta", 2.5 )); 
	ANA_CHECK(m_muonSelection[1]->setProperty( "MuQuality", 1)); //Medium
	ANA_CHECK(m_muonSelection[1]->initialize());

    m_muonSelection[2] = new CP::MuonSelectionTool("MuonSelectionTool_2");
	m_muonSelection[2]->msg().setLevel( MSG::ERROR ); 
	ANA_CHECK(m_muonSelection[2]->setProperty( "MaxEta", 2.5 )); 
	ANA_CHECK(m_muonSelection[2]->setProperty( "MuQuality", 0)); //Tight
	ANA_CHECK(m_muonSelection[2]->initialize());

    //initialize (and configure) the isolation selection tool
    m_IsoSelection[0] = new IsolationSelectionTool ("IsoSelection_0");
    m_IsoSelection[0]->msg().setLevel( MSG::ERROR );
    ANA_CHECK(m_IsoSelection[0]->setProperty( "MuonWP","PflowLoose_VarRad"));
    ANA_CHECK(m_IsoSelection[0]->initialize());

    m_IsoSelection[1] = new IsolationSelectionTool ("IsoSelection_1");
    m_IsoSelection[1]->msg().setLevel( MSG::ERROR );
    ANA_CHECK(m_IsoSelection[1]->setProperty( "MuonWP","PflowTight_VarRad"));
    ANA_CHECK(m_IsoSelection[1]->initialize());

    // Output Tree
	//TFile *outputFile = wk()->getOutputFile(outputName);
    pairs = new TTree ("pairs", "pairs");
    //pairs->SetDirectory (outputFile);

    pairs->Branch("RunNumber", &RunNumber);
    pairs->Branch("EventNumber", &EventNumber);
    pairs->Branch("lbn", &lbn);
    pairs->Branch("bcid", &bcid);
    pairs->Branch("averageIntPerXing", &averageIntPerXing);
    pairs->Branch("actualIntPerXing", &actualIntPerXing);
 	pairs->Branch("actualmufromtool", &actualmufromtool);
	pairs->Branch("averagemufromtool", &averagemufromtool);
	pairs->Branch("ScaledActualMuFromTool", &ScaledActualMuFromTool);
	pairs->Branch("ScaledAverageMuFromTool", &ScaledAverageMuFromTool);

    pairs->Branch("utps66", &utps66);
    pairs->Branch("utps64", &utps64);
    pairs->Branch("utps44", &utps44);
    pairs->Branch("utps66_tag", &utps66_tag);
    pairs->Branch("utps10_10", &utps10_10);
    pairs->Branch("utps10_6", &utps10_6);
    pairs->Branch("utps11_6", &utps11_6);
    pairs->Branch("utps64_tag", &utps64_tag);
    pairs->Branch("utps11_6_tag", &utps11_6_tag);
    pairs->Branch("zbt20", &zbt20);
    pairs->Branch("zbt26", &zbt26);
    pairs->Branch("zbt50", &zbt50);
    pairs->Branch("zbt10_10", &zbt10_10);
    pairs->Branch("zbt_14_14", &zbt_14_14);

	pairs->Branch("match_utps66_mu1", &match_utps66_mu1);
	pairs->Branch("match_utps64_mu1", &match_utps64_mu1);
	pairs->Branch("match_utps44_mu1", &match_utps44_mu1);
	pairs->Branch("match_utps44_tag_mu1", &match_utps44_tag_mu1);
	pairs->Branch("match_utps66_tag_mu1", &match_utps66_tag_mu1);
	pairs->Branch("match_utps10_10_mu1", &match_utps10_10_mu1);
	pairs->Branch("match_utps10_6_mu1", &match_utps10_6_mu1);
	pairs->Branch("match_utps11_6_mu1", &match_utps11_6_mu1);
	pairs->Branch("match_utps64_tag_mu1", &match_utps64_tag_mu1);
	pairs->Branch("match_utps11_6_tag_mu1", &match_utps11_6_tag_mu1);
	pairs->Branch("match_utps66_mu2", &match_utps66_mu2);
	pairs->Branch("match_utps64_mu2", &match_utps64_mu2);
	pairs->Branch("match_utps44_mu2", &match_utps44_mu2);
	pairs->Branch("match_utps44_tag_mu2", &match_utps44_tag_mu2);
	pairs->Branch("match_utps66_tag_mu2", &match_utps66_tag_mu2);
	pairs->Branch("match_utps10_10_mu2", &match_utps10_10_mu2);
	pairs->Branch("match_utps10_6_mu2", &match_utps10_6_mu2);
	pairs->Branch("match_utps11_6_mu2", &match_utps11_6_mu2);
	pairs->Branch("match_utps64_tag_mu2", &match_utps64_tag_mu2);
	pairs->Branch("match_utps11_6_tag_mu2", &match_utps11_6_tag_mu2);

	pairs->Branch("mu1_loose", &mu1_loose);
	pairs->Branch("mu1_medium", &mu1_medium);
	pairs->Branch("mu1_tight", &mu1_tight);
	pairs->Branch("mu2_loose", &mu2_loose);
	pairs->Branch("mu2_medium", &mu2_medium);
	pairs->Branch("mu2_tight", &mu2_tight);

    pairs->Branch("Noftracks", &Noftracks);
    pairs->Branch("Nofbkgs", &Nofbkgs);
    pairs->Branch("Nvtxs", &Nvtxs);
    pairs->Branch("primary_vx_z", &primary_vx_z);
    pairs->Branch("truth_primary_vx_z", &truth_primary_vx_z);
    pairs->Branch("closest_pup", &closest_pup);

    pairs->Branch("pairpt", &pairpt);
    pairs->Branch("pairy", &pairy);
    pairs->Branch("paireta", &paireta);
    pairs->Branch("pairm", &pairm);
    pairs->Branch("pairphi", &pairphi);

    pairs->Branch("recoeff1", &recoeff1);
    pairs->Branch("recotosub1", &recotosub1);
    pairs->Branch("recoerr1", &recoerr1);
    pairs->Branch("recopt1", &recopt1);
    pairs->Branch("recophi1", &recophi1);
    pairs->Branch("recoeta1", &recoeta1);
    pairs->Branch("recoch1", &recoch1);
    pairs->Branch("recod01", &recod01);
    pairs->Branch("recod0s1", &recod0s1);
    pairs->Branch("recotrig1", &recotrig1);
    pairs->Branch("recoutrig1", &recoutrig1);
    pairs->Branch("recoiso1", &recoiso1);

    pairs->Branch("recoeff2", &recoeff2);
    pairs->Branch("recotosub2", &recotosub2);
    pairs->Branch("recoerr2", &recoerr2);
    pairs->Branch("recopt2", &recopt2);
    pairs->Branch("recophi2", &recophi2);
    pairs->Branch("recoeta2", &recoeta2);
    pairs->Branch("recoch2", &recoch2);
    pairs->Branch("recod02", &recod02);
    pairs->Branch("recod0s2", &recod0s2);
    pairs->Branch("recotrig2", &recotrig2);
    pairs->Branch("recoutrig2", &recoutrig2);
    pairs->Branch("recoiso2", &recoiso2);

    pairs->Branch("vxs",&vxs);
    pairs->Branch("vxtype",&vxtype);
    pairs->Branch("vxsumpt",&vxsumpt);


    pairs->Branch("jet4_pt",&jet4_pt);
    pairs->Branch("jet4_uncalibpt",&jet4_uncalibpt);
    pairs->Branch("jet4_eta",&jet4_eta);
    pairs->Branch("jet4_phi",&jet4_phi);
	pairs->Branch("jet4_jvf",&jet4_jvf);
	pairs->Branch("jet_cleaning_1",&jet_cleaning_1);
	pairs->Branch("jet_cleaning_2",&jet_cleaning_2);

    pairs->Branch("trks_eta",&trks_eta);
    pairs->Branch("trks_sinthetaz0",&trks_sinthetaz0);
    pairs->Branch("trks_phi",&trks_phi);
    pairs->Branch("trks_pt",&trks_pt);
	pairs->Branch("trks_d0sig",&trks_d0sig);
    pairs->Branch("trks_true",&trks_true);

	// truth
	pairs->Branch("PileupWeight",&PileupWeight);
	pairs->Branch("Upspt",&Upspt);
	pairs->Branch("Upsy",&Upsy);
	pairs->Branch("Upsm",&Upsm);
	pairs->Branch("Upsphi",&Upsphi);
	pairs->Branch("UpsID",&UpsID);
	pairs->Branch("UpsParentID",&UpsParentID);
	pairs->Branch("truemupt1",&truemupt1);
	pairs->Branch("truemupt2",&truemupt2);
	pairs->Branch("truemueta1",&truemueta1);
	pairs->Branch("truemueta2",&truemueta2);
	pairs->Branch("truemuphi1",&truemuphi1);
	pairs->Branch("truemuphi2",&truemuphi2);
	pairs->Branch("truemuch1",&truemuch1);
	pairs->Branch("truemuch2",&truemuch2);
	pairs->Branch("truesib_recoq",&truesib_recoq);
	pairs->Branch("truesib_recoz0",&truesib_recoz0);
	pairs->Branch("truesib_eta",&truesib_eta);	
	pairs->Branch("truesib_phi",&truesib_phi);
	pairs->Branch("truesib_pt",&truesib_pt);
	pairs->Branch("truesib_ID",&truesib_ID);
	pairs->Branch("truetrk_recoq",&truetrk_recoq);
	pairs->Branch("truetrk_recoz0",&truetrk_recoz0);
	pairs->Branch("truetrk_eta",&truetrk_eta);
	pairs->Branch("truetrk_phi",&truetrk_phi);
	pairs->Branch("truetrk_pt",&truetrk_pt);
	pairs->Branch("truetrk_charge",&truetrk_charge);
	pairs->Branch("truetrk_vxz",&truetrk_vxz);
	pairs->Branch("truetrk_vxx",&truetrk_vxx);
	pairs->Branch("truetrk_vxy",&truetrk_vxy);

	pairs->Branch("truetrk_barcode",&truetrk_barcode);
	pairs->Branch("truetrk_pdg",&truetrk_pdg);
	pairs->Branch("truetrk_status",&truetrk_status);


	pairs->Branch("trks_eta_mc ",&m_trks_eta_mc);
	pairs->Branch("trks_phi_mc",&m_trks_phi_mc);
	pairs->Branch("trks_pt_mc",&m_trks_pt_mc);
	pairs->Branch("trks_mcprob",&m_trks_mcprob);
	pairs->Branch("trks_charge_mc",&m_trks_charge_mc);
	pairs->Branch("trks_barcode",&m_trks_barcode);
	pairs->Branch("trks_pdg",&m_trks_pdg);
	pairs->Branch("trks_status",&m_trks_status);

	wk()->addOutput(pairs);
	
/*
	h_eta_SCTHits = new TH2F("h_eta_SCTHits","h_eta_SCTHits; #eta; SCTHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_SCTOutliers = new TH2F("h_eta_SCTOutliers","h_eta_SCTOutliers; #eta; SCTOutliers",150,-3.,3.,50,-0.5,49.5);
	h_eta_SCTHoles = new TH2F("h_eta_SCTHoles","h_eta_SCTHoles; #eta; SCTHoles",150,-3.,3.,50,-0.5,49.5);
	h_eta_SCTDoubleHoles = new TH2F("h_eta_SCTDoubleHoles","h_eta_SCTDoubleHoles; #eta; SCTDoubleHoles",150,-3.,3.,50,-0.5,49.5);
	h_eta_SCTSharedHits = new TH2F("h_eta_SCTSharedHits","h_eta_SCTSharedHits; #eta; SCTSharedHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_SCTDeadSensors = new TH2F("h_eta_SCTDeadSensors","h_eta_SCTDeadSensors; #eta; SCTDeadSensors",150,-3.,3.,50,-0.5,49.5);
	
	h_eta_PixelHits = new TH2F("h_eta_PixelHits","h_eta_PixelHits; #eta; PixelHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_PixelOutliers = new TH2F("h_eta_PixelOutliers","h_eta_PixelOutliers; #eta; PixelOutliers",150,-3.,3.,50,-0.5,49.5);
	h_eta_PixelHoles = new TH2F("h_eta_PixelHoles","h_eta_PixelHoles; #eta; PixelHoles",150,-3.,3.,50,-0.5,49.5);
	h_eta_PixelSharedHits = new TH2F("h_eta_PixelSharedHits","h_eta_PixelSharedHits; #eta; SharedHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_PixelSplitHits = new TH2F("h_eta_PixelSplitHits","h_eta_PixelSplitHits; #eta; PixelSplitHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_PixelDeadSensors = new TH2F("h_eta_PixelDeadSensors","h_eta_PixelDeadSensors; #eta; PixelDeadSensors",150,-3.,3.,50,-0.5,49.5);

	h_eta_TRTHits = new TH2F("h_eta_TRTHits","h_eta_TRTHits; #eta; TRTHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_TRTOutliers = new TH2F("h_eta_TRTOutlier","h_eta_TRTOutlier; #eta; TRTOutliers",150,-3.,3.,50,-0.5,49.5);
	h_eta_TRTHoles = new TH2F("h_eta_TRTHoles","h_eta_TRTHoles; #eta; TRTHoles",150,-3.,3.,50,-0.5,49.5);
	h_eta_TRTSharedHits = new TH2F("h_eta_TRTSharedHits","h_eta_TRTSharedHits; #eta; TRTSharedHits",150,-3.,3.,50,-0.5,49.5);

	h_eta_SiHits = new TH2F("h_eta_SiHits","h_eta_SiHits; #eta; SiHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_SiHole = new TH2F("h_eta_SiHole","h_eta_SiHole; #eta; SiHoles",150,-3.,3.,50,-0.5,49.5);

	h_eta_IBLHits = new TH2F("h_eta_IBLHits","h_eta_IBLHits; #eta; IBLHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_BLHits = new TH2F("h_eta_BLHits","h_eta_BLHits; #eta; BLHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_expIBLHits = new TH2F("h_eta_expIBLHits","h_eta_expIBLHits; #eta; expIBLHits",150,-3.,3.,50,-0.5,49.5);
	h_eta_expBLHits = new TH2F("h_eta_expBLHits","h_eta_expBLHits; #eta; expBLHits",150,-3.,3.,50,-0.5,49.5);

	wk()->addOutput(h_eta_SCTHits);
	wk()->addOutput(h_eta_SCTOutliers);
	wk()->addOutput(h_eta_SCTHoles);
	wk()->addOutput(h_eta_SCTDoubleHoles);
	wk()->addOutput(h_eta_SCTSharedHits);
	wk()->addOutput(h_eta_SCTDeadSensors);
	wk()->addOutput(h_eta_PixelHits);
	wk()->addOutput(h_eta_PixelOutliers);
	wk()->addOutput(h_eta_PixelHoles);
	wk()->addOutput(h_eta_PixelSharedHits);
	wk()->addOutput(h_eta_PixelSplitHits);
	wk()->addOutput(h_eta_PixelDeadSensors);
	wk()->addOutput(h_eta_TRTHits);
	wk()->addOutput(h_eta_TRTOutliers);
	wk()->addOutput(h_eta_TRTHoles);
	wk()->addOutput(h_eta_TRTSharedHits);
	wk()->addOutput(h_eta_SiHits);
	wk()->addOutput(h_eta_SiHole);
	wk()->addOutput(h_eta_IBLHits);
	wk()->addOutput(h_eta_BLHits);
	wk()->addOutput(h_eta_expIBLHits);
	wk()->addOutput(h_eta_expBLHits);*/
	
	m_eventCounter = 0;
	m_passeventCounter = 0;
	m_numCleanEvents = 0;


    return StatusCode::SUCCESS;
}


StatusCode MyxAODAnalysis :: execute ()
{

    const xAOD::EventInfo *eventInfo = nullptr;
	xAOD::TEvent* event = wk()->xaodEvent();
  	ANA_CHECK (evtStore()->retrieve (eventInfo, "EventInfo"));
	// print every 1000 events, so we know where we are:
  	if( (m_eventCounter % 1000) ==0  ) Info("execute()", "Event number = %i", m_eventCounter);
  	m_eventCounter++;


	static SG::AuxElement::Decorator<unsigned int> dec_rnd("RandomRunNumber");

  	ANA_CHECK(event->retrieve( eventInfo, "EventInfo"));  
  	bool isMC=true;
  	if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION))
	{
		dec_rnd(*eventInfo) = eventInfo->runNumber();
		isMC=false;
	}

	m_PRWTool->apply(*eventInfo);
	if(isMC) PileupWeight = m_PRWTool->getCombinedWeight( *eventInfo );
    if(!eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) && !m_grl_1->passRunLB(*eventInfo) && !m_grl_2->passRunLB(*eventInfo)
                              && !m_grl_3->passRunLB(*eventInfo) && !m_grl_4->passRunLB(*eventInfo))
	{
    	return StatusCode::SUCCESS; // go to next event
  	}

	if((eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)) )
	{
  		return StatusCode::SUCCESS; // go to the next event
  	} // end if event flags check
  	m_numCleanEvents++;

	
	EventNumber = eventInfo->eventNumber();
	RunNumber = eventInfo->runNumber();
	lbn = eventInfo->lumiBlock();
	bcid = eventInfo->bcid();
	averageIntPerXing = eventInfo->averageInteractionsPerCrossing();
	actualIntPerXing =  eventInfo->actualInteractionsPerCrossing();


	utps64=-88;
  	utps44=-88;
  	utps66=-88;
	utps66_tag=-88;
	utps10_10=-88;
	utps10_6=-88;
	utps11_6=-88;
	utps64_tag=-88;
	utps11_6_tag=-88;

	zbt20=-88;
	zbt26=-88;
	zbt50=-88;
	zbt10_10=-88;
	zbt_14_14=-88;

	auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_mu.*");
	std::map<std::string,int> triggerCounts;
	for(auto &trig : chainGroup->getListOfTriggers()) 
	{
		auto cg = m_trigDecisionTool->getChainGroup(trig);
		std::string thisTrig = trig;
		if(thisTrig=="HLT_mu20_iloose_L1MU15") 
		{
			if(cg->isPassed())zbt20 = cg->getPrescale();
			else zbt20 =-98;
		}
		if(thisTrig=="HLT_mu26_ivarmedium") 
		{
			if(cg->isPassed())zbt26 = cg->getPrescale();
			else zbt26 =-98;
		}

		if(thisTrig=="HLT_mu50") 
		{
			if(cg->isPassed())zbt50 = cg->getPrescale();
			else zbt50 =-98;
		}

		if(thisTrig=="HLT_mu6_mu4_bUpsimumu") 
		{
			if(cg->isPassed())utps64 = cg->getPrescale();
			else utps64 =-98;
		}
		if(thisTrig=="HLT_mu6_mu4_bUpsimumu_L1BPH-8M15-MU6MU4_BPH-0DR22-MU6MU4-BO") 
    	{
			if(cg->isPassed())utps64_tag = cg->getPrescale();
			else utps64_tag =-98;
    	}
		if(thisTrig=="HLT_mu10_mu6_bUpsimumu") 
    	{
			if(cg->isPassed())utps10_6 = cg->getPrescale();
			else utps10_6 =-98;
    	}

		if(thisTrig=="HLT_mu11_mu6_bUpsimumu") 
    	{
			if(cg->isPassed())utps11_6 = cg->getPrescale();
			else utps11_6 =-98;
    	}

		if(thisTrig=="HLT_mu11_mu6_bUpsimumu_L1LFV-MU11") 
    	{
			if(cg->isPassed())utps11_6_tag = cg->getPrescale();
			else utps11_6_tag =-98;
    	}
  	} 


  	auto chainGroup2 = m_trigDecisionTool->getChainGroup("HLT_2mu.*");
  	std::map<std::string,int> triggerCounts2;
  	for(auto &trig2 : chainGroup2->getListOfTriggers()) 
	{
		auto cg = m_trigDecisionTool->getChainGroup(trig2);
		std::string thisTrig = trig2;

		if(thisTrig=="HLT_2mu10") 
	    {
			if(cg->isPassed())zbt10_10 = cg->getPrescale();
			else zbt10_10 =-98;
    	}

		if(thisTrig=="HLT_2mu14") 
    	{
			if(cg->isPassed())zbt_14_14 = cg->getPrescale();
			else zbt_14_14 =-98;
    	}

		if(thisTrig=="HLT_2mu4_bUpsimumu") 
    	{
			if(cg->isPassed())utps44 = cg->getPrescale();
			else utps44 =-98;
    	}

		if(thisTrig=="HLT_2mu4_bUpsimumu_L1BPH-1M19-2MU4-BO_BPH-0DR34-2MU4") 
    	{
			if(cg->isPassed())utps44_tag = cg->getPrescale();
			else utps44_tag =-98;
    	}

    	if(thisTrig=="HLT_2mu6_bUpsimumu")
    	{
			if(cg->isPassed())utps66 = cg->getPrescale();
			else utps66 =-98;
    	}

    	if(thisTrig=="HLT_2mu6_bUpsimumu_L1BPH-8M15-2MU6_BPH-0DR22-2MU6")
    	{
			if(cg->isPassed())utps66_tag = cg->getPrescale();
			else utps66_tag =-98;
    	}

    	if(thisTrig=="HLT_2mu10_bUpsimumu")
    	{
			if(cg->isPassed())utps10_10 = cg->getPrescale();
			else utps10_10 =-98;
    	}
  	}


	const xAOD::VertexContainer* pvxs = 0;
  	ANA_CHECK(event->retrieve( pvxs, "PrimaryVertices" ));
  	//loop over vertices
  	xAOD::VertexContainer::const_iterator pvx_itr = pvxs->begin();
  	xAOD::VertexContainer::const_iterator pvx_end = pvxs->end();

  	primary_vx_z=-9999;
  	truth_primary_vx_z=-9999;
  	Nvtxs=0;
  	for(const xAOD::Vertex* vx : *pvxs  ) 
  	{
  		if((vx)->vertexType() !=0) Nvtxs++;
    	else continue;

		vxs->push_back( vx->z());
		vxsumpt->push_back(sqrt(( vx)->auxdata< float >("sumPt2")));
    	if( (vx)->vertexType() == xAOD::VxType::PriVtx) vxtype->push_back(1);
    	else vxtype->push_back(0);

	}
  	//cout<<"first other count = "<<vxs->size()<<" vs "<<Nvtxs<<endl;
/*
  	if(primary_vx_z==-9999) 
	{
		cout<<"No primary vx reconstructed!"<<endl;
		return StatusCode::SUCCESS;
	} // go to the next event if there's no primary vertex*/



	//Look at the truth level stuff
	Upsm = -999;
	Upspt =-999;
	Upsy = -999;
	Upsphi = -999;
	UpsID = -999;
	UpsParentID = -999;

	truemupt1=-999;
	truemupt2=-999;
	truemueta1=-999;
	truemueta2=-999;
	truemuphi1=-999;
	truemuphi2=-999;
	truemuch1=-999;
	truemuch2=-999;


	//------------
	//TRACKS
	//------------
	const xAOD::TrackParticleContainer* tracks = 0;
	ANA_CHECK(event->retrieve( tracks, "InDetTrackParticles" ));
	xAOD::TrackParticleContainer::const_iterator track_itr;// = tracks->begin();
	xAOD::TrackParticleContainer::const_iterator track_end;// = tracks->end();

	int MC_state=0;
	int Ups_ID[3]={553,100553,200553};


	if(isMC)
	{
		const xAOD::TruthVertexContainer* true_vxs = 0;
    	ANA_CHECK(event->retrieve(true_vxs, "TruthVertices" ));

		xAOD::TruthVertexContainer::const_iterator tvItr = true_vxs->begin();
    	xAOD::TruthVertexContainer::const_iterator tvItrE = true_vxs->end();
    	const xAOD::TruthVertex* vertex = (*tvItr);
   		truth_primary_vx_z = vertex->z();

		const xAOD::TruthParticleContainer* trues = 0;
    	ANA_CHECK(event->retrieve(trues, "TruthParticles" ));

		xAOD::TruthParticleContainer::const_iterator tpItr = trues->begin();
		xAOD::TruthParticleContainer::const_iterator tpItrE = trues->end();
		xAOD::TruthParticleContainer::const_iterator tpJItr = trues->begin();
		xAOD::TruthParticleContainer::const_iterator tpJItrE = trues->end();

		for(; tpItr!=tpItrE; ++tpItr)
		{
			const xAOD::TruthParticle* particle = (*tpItr);
			if(particle->pdgId()!=Ups_ID[MC_state])continue;
			if(fabs(particle->rapidity())>2.5) continue;

			if (particle->nChildren()>0)
			{
				bool first=true;
				for (unsigned int i=0; i<particle->nChildren(); i++)
				{
					const xAOD::TruthParticle* Upschild = particle->child( i );
					if(Upschild->status()!=1) continue;//skip documentary particles
          			if(!Upschild->isMuon()) continue;//don't care about non-muon children

					if(first)
					{
						truemupt1 = Upschild->pt()/1000.;
						truemuphi1 = Upschild->phi();
						truemueta1 = Upschild->eta();
						truemuch1 = Upschild->charge();
						first=false;
					} // if first
					else
					{
						truemupt2 = Upschild->pt()/1000.;
						truemuphi2 = Upschild->phi();
						truemueta2 = Upschild->eta();
						truemuch2 = Upschild->charge();
					}
				} // for i
			} // if particle

			if(truemupt1<-900 || truemupt2<-900) continue;
			if(fabs(truth_primary_vx_z)>500) 
			{
				cout<<"Out of bounds vertex from truth upsilon!"<<endl;return EL::StatusCode::SUCCESS; 
			}

			Upsm = particle->m()/1000.;
			Upspt = particle->pt()/1000.;
			Upsy = particle->rapidity();
			Upsphi = particle->phi();//Upsphi = particle->nChildren();
			UpsID = particle->pdgId();//UpsID = particle->status();


			//look for Upsilon parent and siblings/cousins
			UpsParentID = -888.8;

			if(particle->nParents()>0)
			{
				for(unsigned int i=0; i<particle->nParents(); i++)
				{
					const xAOD::TruthParticle *parent = particle->parent(i);
					bool right_parent=false;
		    		for(int usi=MC_state+1;usi<3;usi++) 
					{
						if (parent->pdgId()==Ups_ID[usi]) right_parent=true;
					}
					if(!right_parent)continue;

					UpsParentID = parent->pdgId();	
				}

				if(UpsParentID>-888.7||UpsParentID<-888.9)
				{
					//if there's an excited upsilon parent for the 1S look for siblings
					
					for(; tpJItr!=tpJItrE; ++tpJItr)
					{
						const xAOD::TruthParticle* maybesibling = (*tpJItr);
						if(maybesibling->status()!=1 || maybesibling->isMuon() || maybesibling->pdgId()==22 ) continue;//if it's NOT final state or IS a muon/photon skip it
						
						if(maybesibling->nParents()>0)
						{
							bool is_real_sibling=false;
			    			for(unsigned int i=0; i < maybesibling->nParents(); i++) 
							{
								//const xAOD::TruthParticle *parent = maybesibling->parent(i);
				      			//if(parent->pdgId()!=100553&&parent->pdgId()!=200553)continue;
			    	  			bool aright_parent=false;
		    	  				for(int usi=MC_state+1;usi<3;usi++) 
								{	
									//cout<<"parent->pdgId() = "<<parent->pdgId()<<endl;
									//if (parent->pdgId() == Ups_ID[usi]) aright_parent=true;
								}
		      					if(!aright_parent)continue;
		      					is_real_sibling=true;
							}


							if( is_real_sibling)
							{
								truesib_ID->push_back(maybesibling->pdgId());
								truesib_pt->push_back(maybesibling->pt()/1000.);
								truesib_phi->push_back(maybesibling->phi());
								truesib_eta->push_back(maybesibling->eta());
							}//close if real sibling


						}//close if putative sibling has parents

					}//close jitr loop over truth



				}// if UpsParentID

			} // if particle

		} // for tpItr

		//now loop again for all tracks
    	tpItr = trues->begin();
    	tpItrE = trues->end();

		for(; tpItr!=tpItrE; ++tpItr)
		{
			const xAOD::TruthParticle* particle = (*tpItr);
			//now that I'm done looking for the upsilon let me get truth tracks
			if(particle->pt()<390) continue;
			//	if(fabs(particle->eta())>2.5) continue;
			if(particle->pdgId()==22)continue;
			if(particle->status()!=1 ) continue;
			if(particle->barcode()>=200000 ) continue;
			bool is_mu_from_ups = false;

			const xAOD::TruthParticle *tparent = particle->parent(0);
			if(tparent)
			{
				if(tparent->pdgId()==Ups_ID[MC_state] && particle->isMuon()) is_mu_from_ups=true;
			}

			if(!is_mu_from_ups)
			{
				double negflip =-1;//default to say it's pileup
	    		ElementLink< xAOD::TruthVertexContainer > truthverLink = particle->auxdata<ElementLink< xAOD::TruthVertexContainer > >("prodVtxLink");
				if(truthverLink.isValid())
				{
					if(truth_primary_vx_z==(*truthverLink)->z()) negflip=1;//if it's from the PV
					truetrk_vxz->push_back((*truthverLink)->z());
					truetrk_vxy->push_back((*truthverLink)->y());
					truetrk_vxx->push_back((*truthverLink)->x());
				}
				else 
				{
					negflip=0;
					cout<<"No vertex link!"<<endl;
					truetrk_vxz->push_back(-99999);
					truetrk_vxx->push_back(-99999);
					truetrk_vxy->push_back(-99999);
				}

				truetrk_pt->push_back(negflip*particle->pt()/1000.);
	    		truetrk_charge->push_back(particle->charge());
	    		truetrk_eta->push_back(particle->eta());
	    		truetrk_phi->push_back(particle->phi());

				truetrk_barcode->push_back(particle->barcode()); 
				truetrk_pdg->push_back(particle->pdgId()); 
				truetrk_status->push_back(particle->status()); 
			}//close chose a good truth particle - not muons from upsilons

		}//end (2nd) loop over truth particles

	}//close if(isMC)
	


	xAOD::TStore store;
	// get jet container of interest
  	const xAOD::JetContainer* jets = 0;
  	ANA_CHECK(event->retrieve(jets,"AntiKt4EMTopoJets"));
  	//  const xAOD::JetContainer* jets10 = 0;
  	//  ANA_CHECK(event->retrieve( jets10, "AntiKt10LCTopoJets" ));
  	//  cout<<"opened jet containers"<<endl;
  	//  int numGoodJets = 0;   
  	// loop over the jets in the container
  	float hi_jpt=10;
  	float hi_jeta=-9999;
  	float hi_jphi=-9999;

	xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	xAOD::JetContainer::const_iterator jet_end = jets->end();
	int NumGoodJets = 0;


	//for(auto *ijet : *jets)
  	for( ; jet_itr != jet_end; ++jet_itr ) 
	{
  		xAOD::Jet *this_calibjet = 0;
		// changed accept to keep
		//if(!m_jetCleaning->keep(*ijet)) continue;
    	if(!m_jetCleaning_1->keep(**jet_itr) && !m_jetCleaning_2->keep(**jet_itr)) continue; //only keep good clean jets
		if(m_jetCleaning_1->keep(**jet_itr))jet_cleaning_1->push_back(1);
		else jet_cleaning_1->push_back(0);

		if(m_jetCleaning_2->keep(**jet_itr))jet_cleaning_2->push_back(1);
		else jet_cleaning_2->push_back(0);
   		
		//m_JetCalibration->calibratedCopy(*ijet,this_calibjet);
    	m_JetCalibration->calibratedCopy(**jet_itr,this_calibjet); //make a calibrated copy
		
    	if(this_calibjet->pt() * 0.001 <18) 
		{
			delete this_calibjet;
			continue; 
		}//testf was run with 15
		
    //    if(getdR(recoeta1,recophi1,this_calibjet->eta(),this_calibjet->phi() )<0.4) continue;
    //    if(getdR(recoeta2,recophi2,this_calibjet->eta(),this_calibjet->phi() )<0.4) continue;

    	if(this_calibjet->pt() * 0.001 > hi_jpt) 
    	{
			hi_jpt=this_calibjet->pt() * 0.001;//the highest jet pT I've read so far
			hi_jeta=this_calibjet->eta();
			hi_jphi=this_calibjet->phi();
    	}
		
    	jet4_pt->push_back(this_calibjet->pt() * 0.001 );
    	jet4_uncalibpt->push_back((*jet_itr)->pt() * 0.001 );
    	jet4_eta->push_back(this_calibjet->eta() );
    	jet4_phi->push_back(this_calibjet->phi() );
    	jet4_jvf->push_back(hjvtagup->updateJvt(*this_calibjet)); 
		
		NumGoodJets++;
		delete this_calibjet;

  	} // end for loop over jets



  	//------------
	// MUONS
	//------------

	// get muon container of interest
	const xAOD::MuonContainer* muons = 0;
	ANA_CHECK(event->retrieve( muons, "Muons" ));

		// create a shallow copy of the muons container
	std::pair< xAOD::MuonContainer*,xAOD::ShallowAuxContainer* > muons_shallowCopy = xAOD::shallowCopyContainer( *muons );
	
	// iterate over our shallow copy
	xAOD::MuonContainer::iterator muon_itr = (muons_shallowCopy.first)->begin();
	xAOD::MuonContainer::iterator muon_end = (muons_shallowCopy.first)->end();
	
	int count=0;
	bool trigger_matched=false;
	const xAOD::TrackParticle* m1track =0;
	const xAOD::TrackParticle* m2track =0;
	double muon_etas[2]={-100,-100};
	double muon_pts[2]={-10000,-10000};
	double muon_phis[2]={-100,-100};
	
	for (;muon_itr != muon_end; ++muon_itr)
	{
		if(m_muonCalibrationAndSmearingTool->applyCorrection(**muon_itr) == CP::CorrectionCode::Error)
		{ 
			// apply correction and check return code
			// Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
			// If OutOfValidityRange is returned no modification is made and the original muon values are taken.
			Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
		}
    	if((*muon_itr)->pt() * 0.001 < 4) continue;
		if(!m_muonSelection[0]->accept(**muon_itr)) continue;

		if(m_muonSelection[0]->accept(**muon_itr)) mu1_loose->push_back(1);
		else mu1_loose->push_back(0);

		if(m_muonSelection[1]->accept(**muon_itr)) mu1_medium->push_back(1);
		else mu1_medium->push_back(0);

		if(m_muonSelection[2]->accept(**muon_itr)) mu1_tight->push_back(1);
		else mu1_tight->push_back(0);
	}
	




	//if(pairpt->size()>0 || UpsID>0) pairs->Fill();

	pairs->Fill();
	
	recotrig1->clear();
  	recotrig2->clear();
	recoutrig1->clear();
	recoutrig2->clear();
	recod01->clear();
	recod02->clear();
	recod0s1->clear();
	recod0s2->clear();
	recoiso1->clear();
	recoiso2->clear();
	recopt1->clear();
	recopt2->clear();
	recoeff1->clear();
	recoeff2->clear();
	recotosub1->clear();
	recotosub2->clear();
	recoerr1->clear();
	recoerr2->clear();
	recoeta1->clear();
	recoeta2->clear();
	recophi1->clear();
	recophi2->clear();
	recoch1->clear();
	recoch2->clear();
	pairpt->clear();
	pairy->clear();
	paireta->clear();
	pairm->clear();
	pairphi->clear();
	jet4_pt->clear();
	jet4_uncalibpt->clear();
	jet4_eta->clear();
	jet4_phi->clear();
	jet4_jvf->clear();
	jet_cleaning_1->clear();
	jet_cleaning_2->clear();
	vxs->clear();
	vxtype->clear();
	vxsumpt->clear();
	trks_eta->clear();
	trks_sinthetaz0->clear();
	trks_phi->clear();
	trks_pt->clear();
	trks_d0sig->clear();


	truesib_recoq->clear(); 
	truesib_recoz0->clear();
	truesib_eta->clear();
	truesib_phi->clear();
	truesib_pt->clear(); 
	truesib_ID->clear(); 
	truetrk_recoq->clear(); 
	truetrk_recoz0->clear();
	truetrk_eta->clear();
	truetrk_phi->clear();
	truetrk_pt->clear(); 
	truetrk_charge->clear(); 
	truetrk_vxz->clear(); 
	truetrk_vxx->clear(); 
	truetrk_vxy->clear(); 

	truetrk_barcode->clear(); 
	truetrk_pdg->clear(); 
	truetrk_status->clear(); 

	m_trks_eta_mc->clear();
	m_trks_phi_mc->clear();
	m_trks_pt_mc->clear();
	m_trks_mcprob->clear();
	m_trks_charge_mc->clear();
	m_trks_barcode->clear();
	m_trks_pdg->clear();
	m_trks_status->clear();


	match_utps66_mu1->clear();
	match_utps64_mu1->clear();
	match_utps44_mu1->clear();
	match_utps44_tag_mu1->clear();
	match_utps66_tag_mu1->clear();
	match_utps10_10_mu1->clear();
	match_utps10_6_mu1->clear();
	match_utps11_6_mu1->clear();
	match_utps64_tag_mu1->clear();
	match_utps11_6_tag_mu1->clear();
	match_utps66_mu2->clear();
	match_utps64_mu2->clear();
	match_utps44_mu2->clear();
	match_utps44_tag_mu2->clear();
	match_utps66_tag_mu2->clear();
	match_utps10_10_mu2->clear();
	match_utps10_6_mu2->clear();
	match_utps11_6_mu2->clear();
	match_utps64_tag_mu2->clear();
	match_utps11_6_tag_mu2->clear();

	mu1_loose->clear();
    mu1_medium->clear();
    mu1_tight->clear();

    mu2_loose->clear();
    mu2_medium->clear();
    mu2_tight->clear(); 

	return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: finalize ()
{

    ANA_CHECK_SET_TYPE (StatusCode); // set type of return code you are expecting (add to top of each function once) 
    xAOD::TEvent* event = wk()->xaodEvent();
    Info("finalize()", "Number of clean events = %i", m_numCleanEvents);
    cout<<"Number of pass mass  = "<<m_passeventCounter<<endl;

	if(m_grl_1)
	{
		delete m_grl_1;
		m_grl_1 = 0;
	}

    if(m_grl_2)
	{
		delete m_grl_2;
		m_grl_2 = 0;
	}

    if(m_grl_3)
	{
		delete m_grl_3;
		m_grl_3 = 0;
	}

    if(m_grl_4)
	{
		delete m_grl_4;
		m_grl_4 = 0;
	}

    for (size_t i = 0; i < 2; i++)
    {
        if(m_trackSelection[i])
        {
            delete m_trackSelection[i];
            m_trackSelection[i] = 0;
        }
    }
    

    for (size_t i = 0; i < 3; i++)
    {
        if(m_effi_corr[i])
	    {
	        delete m_effi_corr[i];
	        m_effi_corr[i] = 0;
	    }
    }

    for (size_t i = 0; i < 3; i++)
    {
        if(m_muonSelection[i])
	    {
		    delete m_muonSelection[i];
		    m_muonSelection[i] = 0;
	    }
    }

    if(m_muonCalibrationAndSmearingTool)
	{
  	    delete m_muonCalibrationAndSmearingTool;
   	    m_muonCalibrationAndSmearingTool = 0;
	}

    for (size_t i = 0; i < 2; i++)
    {
        if(m_IsoSelection[i])
        {
            delete m_IsoSelection[i];
            m_IsoSelection[i] = 0;
        }
    }


    if(m_trigDecisionTool)
	{
		delete m_trigDecisionTool;
		m_trigDecisionTool = 0;
	}

	if(m_trigConfigTool)
	{
		delete m_trigConfigTool;
		m_trigConfigTool = 0;
	}

	if(m_trigMatchingTool)
	{
		delete m_trigMatchingTool;
		m_trigMatchingTool = 0;
	}

    if(pjvtag)
	{
		delete pjvtag;
		pjvtag = 0;
	}

	return StatusCode::SUCCESS;
}
