# The name of the package:
atlas_subdir (MyAnalysis)

# Add the shared library:
atlas_add_library (MyAnalysisLib
  MyAnalysis/*.h Root/*.cxx
  PUBLIC_HEADERS MyAnalysis
  LINK_LIBRARIES AsgTools AnaAlgorithmLib xAODEventInfo AsgAnalysisInterfaces TrigDecisionToolLib TrigConfInterfaces AthContainers AthLinks GeoPrimitives EventLoopGrid EventLoop PathResolver GoodRunsListsLib xAODTracking xAODTruth TrigConfxAODLib SampleHandler InDetTrackSelectionToolLib xAODRootAccess xAODJet JetInterface xAODMuon MuonAnalysisInterfacesLib xAODTrigMuon xAODTrigger  MuonEfficiencyCorrectionsLib MuonMomentumCorrectionsLib MuonSelectorToolsLib TriggerMatchingToolLib EventLoopAlgs JetCalibToolsLib JetMomentToolsLib JetResolutionLib ElectronPhotonSelectorToolsLib IsolationSelectionLib JetSelectorToolsLib EventLoopGrid  xAODCaloEvent)

# The dependencies of the package:
atlas_depends_on_subdirs(
   PUBLIC
   Control/AthContainers
   Control/AthToolSupport/AsgTools
   Control/CxxUtils
   Control/xAODRootAccess
   Control/xAODRootAccessInterfaces
   DataQuality/GoodRunsLists/GoodRunsListSelectionTool
   Event/EventPrimitives
   Event/xAOD/xAODBase
   Event/xAOD/xAODCaloEvent
   Event/xAOD/xAODCore
   Event/xAOD/xAODEgamma
   Event/xAOD/xAODEventInfo
   Event/xAOD/xAODForward
   Event/xAOD/xAODHIEvent
   Event/xAOD/xAODJet
   Event/xAOD/xAODMissingET
   Event/xAOD/xAODMuon
   Event/xAOD/xAODRootAccess
   Event/xAOD/xAODTau
   Event/xAOD/xAODTracking
   Event/xAOD/xAODTrigger
   Event/xAOD/xAODTrigMinBias
   Event/xAOD/xAODTruth
   InnerDetector/InDetRecTools/InDetTrackSelectionTool
   PhysicsAnalysis/AnalysisCommon/AssociationUtils
   PhysicsAnalysis/AnalysisCommon/IsolationSelection
   PhysicsAnalysis/AnalysisCommon/PATInterfaces
   PhysicsAnalysis/AnalysisCommon/PileupReweighting
   PhysicsAnalysis/D3PDTools/EventLoop
   PhysicsAnalysis/D3PDTools/RootCoreUtils
   PhysicsAnalysis/D3PDTools/SampleHandler
   PhysicsAnalysis/ElectronPhotonID/ElectronEfficiencyCorrection
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonSelectorTools
   PhysicsAnalysis/ElectronPhotonID/ElectronPhotonShowerShapeFudgeTool
   PhysicsAnalysis/ElectronPhotonID/IsolationCorrections
   PhysicsAnalysis/ElectronPhotonID/PhotonEfficiencyCorrection
   PhysicsAnalysis/JetMissingEtID/JetSelectorTools
   PhysicsAnalysis/JetTagging/JetTagPerformanceCalibration/xAODBTaggingEfficiency
   PhysicsAnalysis/MuonID/MuonSelectorTools
   PhysicsAnalysis/MuonID/MuonsIDAnalysis/MuonEfficiencyCorrections
   PhysicsAnalysis/MuonID/MuonsIDAnalysis/MuonMomentumCorrections
   PhysicsAnalysis/TauID/TauAnalysisTools
   Reconstruction/Jet/JetCalibTools
   Reconstruction/Jet/JetJvtEfficiency
   Reconstruction/Jet/JetMomentTools
   Reconstruction/Jet/JetResolution
   Reconstruction/Jet/JetUncertainties
   Reconstruction/MET/METInterface
   Reconstruction/MET/METUtilities
   Tools/PathResolver
   Trigger/TrigAnalysis/TrigDecisionTool
   Trigger/TrigAnalysis/TriggerMatchingTool
   Trigger/TrigConfiguration/TrigConfxAOD )



if (XAOD_STANDALONE)
 # Add the dictionary (for AnalysisBase only):
 atlas_add_dictionary (MyAnalysisDict
  MyAnalysis/MyAnalysisDict.h
  MyAnalysis/selection.xml
  LINK_LIBRARIES MyAnalysisLib)
endif ()

if (NOT XAOD_STANDALONE)
  # Add a component library for AthAnalysis only:
  atlas_add_component (MyAnalysis
    src/components/*.cxx
    LINK_LIBRARIES MyAnalysisLib)
endif ()

# External(s) used by the package:
find_package( Boost COMPONENTS program_options )
find_package( ROOT COMPONENTS Core Hist RIO Gpad Tree MathCore Graf )
find_package( Eigen )


# Install files from the package:
atlas_install_data(share/*
	DESTINATION share)
atlas_install_xmls( data/*.xml )
atlas_install_joboptions( share/*_jobOptions.py )
atlas_install_scripts( share/*_eljob.py )
#atlas_add_executable(ATestSubmit
#   share/ATestSubmit.cxx
#   LINK_LIBRARIES MyAnalysisLib)
