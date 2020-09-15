#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/CondorDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopGrid/GridDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include "EventLoop/OutputStream.h"
#include "EventLoopAlgs/NTupleSvc.h"
#include <TSystem.h>
#include "SampleHandler/ScanDir.h"
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AnaAlgorithm/AnaAlgorithmConfig.h>

#include "MyAnalysis/MyxAODAnalysis.h"

void ATestSubmit(int vers,const std::string& submitDir)
{
  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // create a new sample handler to describe the data files we use
  SH::SampleHandler sh;

  // scan for datasets in the given directory
  // this works if you are on lxplus, otherwise you'd want to copy over files
  // to your local machine and use a local path.  if you do so, make sure
  // that you copy all subdirectories and point this to the directory
  // containing all the files, not the subdirectories.
	
	const int nFiles = 65;
	// 2015 RNs
	//													0				1				2				3				4
	int runNumber[nFiles] = {276262, 276329, 276336, 276416, 276511,
										 			 276689, 276778, 276790, 276952, 276954,
													 278880, 278912, 278968, 279169, 279259,
													 279279, 279284, 279345, 279515, 279598,
													 279685, 279813, 279867, 279928, 279932,
													 279984, 280231, 280273, 280319, 280368,
													 280423, 280464, 280500, 280520, 280614,
													 280673, 280753, 280853, 280862, 280950,
													 280977, 281070, 281074, 281075, 281317,
													 281385, 281411, 282625, 282631, 282712,
													 282784, 282992, 283074, 283155, 283270,
													 283429, 283608, 283780, 284006, 284154,
													 284213, 284285, 284420, 284427, 284484};

	char fileName[300];
	sprintf(fileName,"data15_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r9264_p3083_p3713/", runNumber[vers]);

	SH::scanRucio(sh, fileName);

  // set the name of the tree in our files
  // in the xAOD the TTree containing the EDM containers is "CollectionTree"
  sh.setMetaString ("nc_tree", "CollectionTree");
  sh.print();
  // this is the basic description of our job
  EL::Job job;
  job.sampleHandler (sh); // use SampleHandler in this job
  //job.options()->setDouble (EL::Job::optMaxEvents, 500); // for testing purposes, limit to run over the first 500 events only!
  // add our algorithm to the job
	EL::OutputStream output  ("myOutput");
  job.outputAdd (output);
  EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
  job.algsAdd (ntuple);
	job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);
	job.options()->setDouble(EL::Job::optMaxEvents, -1);
	job.options()->setDouble(EL::Job::optFilesPerWorker,1);

  EL::AnaAlgorithmConfig alg;
  alg.setType ("MyxAODAnalysis");
  alg.setName ("AnalysisAlg");
  job.algsAdd (alg);
	//alg->outputName = "myOutput";
  // make the driver we want to use:
  // this one works by running the algorithm directly:
	EL::PrunDriver driver;
	driver.options()->setString("nc_outputSampleName", "user.iaizenbe.data15_13TeV.%in:name[2]%.%in:name[6]%");
  // we can use other drivers to run things on the Grid, with PROOF, etc.
  // process the job using the driver
  //driver.submit (job, submitDir);

	driver.submitOnly(job, submitDir);
}
