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
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
void ATestSubmit (int vers, const std::string& submitDir)
{
  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // create a new sample handler to describe the data files we use
  SH::SampleHandler sh;


	const int nFiles = 589;
	// 2015 RNs 0-64
	// 2016 RNs 65-208
	// 2017 RNs 209-394
	// 2018 RNs 395-588
	// 													0				1				2				3				4
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
													 284213, 284285, 284420, 284427, 284484,

													 297730, 298595, 298633, 298690, 298771,
													 298773, 298862, 298967, 299055, 299144, 
													 299147, 299184, 299243, 299584, 300279,
													 300345, 300415, 300418, 300487, 300540,
													 300571, 300600, 300655, 300687, 300784, 
													 300800, 300863, 300908, 301912, 301918,
													 301932, 301973, 302053, 302137, 302265,
													 302269, 302300, 302347, 302380, 302391,
													 302393, 302737, 302831, 302872, 302919,
													 302925, 302956, 303007, 303079, 303201,
													 303208, 303264, 303266, 303291, 303304,
													 303338, 303421, 303499, 303560, 303638,
													 303832, 303846, 303892, 303943, 304006,
													 304128, 304178, 304198, 304211, 304308,
													 304337, 304409, 304431, 304494, 305380,
													 305543, 305571, 305618, 305671, 305674,
													 305723, 305727, 305735, 305777, 305811,
													 305920, 306269, 306278, 306310, 306384,
													 306419, 306442, 306448, 306451, 307126,
													 307195, 307259, 307306, 307354, 307358,
											 		 307394, 307454, 307539, 307569, 307601,
													 307619, 307656, 307710, 307716, 307732,
													 307861, 307935, 308047, 308084, 309375,
													 309390, 309440, 309516, 309674, 309759,
													 310015, 310247, 310249, 310341, 310370,
													 310405, 310468, 310473, 310634, 310691,
													 310738, 310809, 310863, 310872, 310969,
													 311071, 311170, 311244, 311287, 311321,
													 311365, 311402, 311473, 311481,

													 325713, 325790, 326439, 326446, 326468,
													 326551, 326657, 326695, 326834, 326870,
													 326923, 326945, 327057, 327103, 327265,
													 327342, 327490, 327582, 327636, 327662,
													 327745, 327761, 327764, 327860, 327862,
													 328042, 328099, 328221, 328263, 328333,
													 328374, 328393, 329716, 329778, 329780,
													 329829, 329835, 329869, 329964, 330025,
													 330074, 330079, 330101, 330160, 330166,
													 330203, 330294, 330470, 331033, 331082,
													 331085, 331129, 331215, 331239, 331697,
													 331710, 331742, 331772, 331804, 331825,
													 331860, 331875, 331905, 331951, 331975,
													 332303, 332304, 332720, 332896, 332915,
													 332953, 332955, 333181, 333192, 333367,
													 333380, 333426, 333469, 333487, 333519,
													 333650, 333707, 333778, 333828, 333853,
												   333904, 333979, 333994, 334264, 334317,
													 334350, 334384, 334413, 334443, 334455,
													 334487, 334564, 334580, 334588, 334637,
													 334678, 334737, 334779, 334842, 334849,
													 334878, 334890, 334907, 334960, 334993,
													 335016, 335022, 335056, 335082, 335083,
													 335131, 335170, 335177, 335222, 335282,
													 335290, 336506, 336548, 336567, 336630,
													 336678, 336719, 336782, 336832, 336852,
													 336915, 336927, 336944, 336998, 337005,
													 337052, 337107, 337156, 337176, 337215,
													 337263, 337335, 337371, 337404, 337451,
													 337491, 337542, 337662, 337705, 337833,
													 338183, 338220, 338259, 338263, 338349,
													 338377, 338480, 338498, 338608, 338675,
													 338712, 338767, 338834, 338846, 338897,
													 338933, 338967, 338987, 339037, 339070,
													 339205, 339346, 339387, 339396, 339435,
													 339500, 339535, 339562, 339590, 339758,
													 339849, 339957, 340030, 340072, 340368,
													 340453,

													 348885, 348894, 348895, 349011, 349014,
													 349033, 349051, 349111, 349114, 349169,
													 349268, 349309, 349327, 349335, 349451,
													 349481, 349498, 349526, 349533, 349534,
													 349582, 349592, 349637, 349646, 349693,
													 349841, 349842, 349944, 349977, 350013,
													 350067, 350121, 350144, 350160, 350184,
													 350220, 350310, 350361, 350440, 350479,
													 350531, 350682, 350749, 350751, 350803,
													 350842, 350848, 350880, 350923, 351062,
													 351160, 351223, 351296, 351325, 351359,
													 351364, 351455, 351550, 351628, 351636,
													 351671, 351698, 351832, 351894, 351969,
													 352056, 352107, 352274, 352340, 352394,
													 352436, 352448, 352494, 352514, 355261,
													 355273, 355529, 355544, 355563, 355599,
													 355650, 355651, 355754, 355848, 355861,
													 355877, 355995, 356077, 356095, 356124,
													 356177, 356205, 356250, 356259, 357193,
													 357283, 357293, 357355, 357409, 357451,
													 357500, 357539, 357620, 357679, 357713,
													 357750, 357772, 357821, 357887, 357962,
													 358031, 358096, 358115, 358175, 358215,
													 358233, 358300, 358325, 358333, 358395,
													 358516, 358541, 358577, 358615, 358656,
													 358985, 359010, 359058, 359124, 359170,
													 359171, 359191, 359279, 359286, 359310,
													 359355, 359398, 359441, 359472, 359541,
													 359586, 359593, 359623, 359677, 359678,
													 359717, 359735, 359766, 359823, 359872,
													 359918, 360026, 360063, 360129, 360161,
													 360209, 360244, 360293, 360309, 360348,
													 360373, 360402, 361738, 361795, 361862,
													 362204, 362297, 362345, 362354, 362388,
													 362445, 362552, 362619, 362661, 362776,
													 363033, 363096, 363129, 363198, 363262,
													 363400, 363664, 363710, 363738, 363830,
													 363910, 363947, 363979, 364030, 364076,
													 364098, 364160, 364214, 364292};




	char fileName[300];
	if(runNumber[vers] < 297730) sprintf(fileName,"data15_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r9264_p3083_p3713/", runNumber[vers]);
	if((runNumber[vers] >= 297730) && (runNumber[vers] < 325713)) sprintf(fileName,"data16_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r9264_p3083_p3713/", runNumber[vers]);
	if((runNumber[vers] >= 325713) && (runNumber[vers] < 348885)) 
	{
		//sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10202_p3399_p3713/", runNumber[vers]);
		//sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10250_p3399_p3713/", runNumber[vers]);
		//sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10203_p3399_p3713/", runNumber[vers]);
		//sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10258_p3399_p3713/", runNumber[vers]);
		//sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10259_p3399_p3713/", runNumber[vers]);
		//sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10260_p3399_p3713/", runNumber[vers]);
		sprintf(fileName,"data17_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.r10426_p3399_p3713/", runNumber[vers]);
	}
	if(runNumber[vers] >= 348885)
 	{
		//sprintf(fileName,"data18_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.f964_m2020_p3713/", runNumber[vers]);
		//sprintf(fileName,"data18_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.f1002_m2037_p3713/", runNumber[vers]);
		//sprintf(fileName,"data18_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.f937_m1972_p3713/", runNumber[vers]);
		//sprintf(fileName,"data18_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.f938_m1979_p3713/", runNumber[vers]);
		//sprintf(fileName,"data18_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.f947_m1993_p3713/", runNumber[vers]); // no samples ...
		sprintf(fileName,"data18_13TeV.00%i.physics_Main.deriv.DAOD_STDM7.f961_m2015_p3713/", runNumber[vers]);

	}
	SH::scanRucio (sh, fileName);
	sh.setMetaString ("nc_tree", "CollectionTree");

  // further sample handler configuration may go here

  // print out the samples we found
  sh.print ();

  // this is the basic description of our job
  EL::Job job;
  job.sampleHandler (sh); // use SampleHandler in this job
  //job.options()->setDouble (EL::Job::optMaxEvents, 500); // for testing purposes, limit to run over the first 500 events only!

  // add our algorithm to the job
  EL::AnaAlgorithmConfig alg;
  alg.setType ("MyxAODAnalysis");

  // set the name of the algorithm (this is the name use with
  // messages)
  alg.setName ("AnalysisAlg");

  // later on we'll add some configuration options for our algorithm that go here

  job.algsAdd (alg);
	job.options()->setDouble (EL::Job::optMaxEvents, 10000000);
	job.options()->setDouble (EL::Job::optCacheSize, 10*1024*1024);
	job.options()->setDouble (EL::Job::optCacheLearnEntries, 20);

  // make the driver we want to use:
  // this one works by running the algorithm directly:
  EL::PrunDriver driver;
	driver.options()->setString("nc_outputSampleName", "user.iaizenbe.test.%in:name[2]%.%in:name[6]%");
  // we can use other drivers to run things on the Grid, with PROOF, etc.

  // process the job using the driver
  driver.submitOnly(job, submitDir);
}
