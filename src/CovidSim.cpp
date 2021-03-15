/*
(c) 2004-20 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
*/

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "CovidSim.h"
#include "Rand.h"
#include "Error.h"
#include "Dist.h"
#include "Kernels.h"
#include "Bitmap.h"
#include "Model.h"
#include "Param.h"
#include "SetupModel.h"
#include "ModelMacros.h"
#include "InfStat.h"
#include "CalcInfSusc.h"
#include "Update.h"
#include "Sweep.h"
#include "Memory.h"
#include "CLI.h"

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

// Use the POSIX name for case-insensitive string comparison: strcasecmp.
#ifdef _WIN32
// Windows calls it _stricmp so make strcasecmp an alias.
#include <string.h>
#define strcasecmp _stricmp
#else
#include <strings.h>
#endif

void parse_bmp_option(std::string const&);
void parse_intervention_file_option(std::string const&);
void ReadParams(std::string const&, std::string const&, std::string const&);
void ReadInterventions(std::string const&);
int GetXMLNode(FILE*, const char*, const char*, char*, int);
void ReadAirTravel(std::string const&, std::string const&);
void InitModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void SeedInfection(double, int*, int, int); //adding run number as a parameter for event log: ggilani - 15/10/2014
int RunModel(int, std::string const&, std::string const&, std::string const&);

void SaveDistribs(std::string const&);
void SaveOriginDestMatrix(std::string const&); //added function to save origin destination matrix so it can be done separately to the main results: ggilani - 13/02/15
void SaveResults(std::string const&);
void SaveSummaryResults(std::string const&);
void SaveRandomSeeds(std::string const&); //added this function to save random seeds for each run: ggilani - 09/03/17
void SaveEvents(std::string const&); //added this function to save infection events from all realisations: ggilani - 15/10/14
void LoadSnapshot(std::string const&);
void SaveSnapshot(std::string const&);
void RecordInfTypes(void);

void RecordSample(double, int, std::string const&);
void CalibrationThresholdCheck(double, int);
void CalcLikelihood(int, std::string const&, std::string const&);
void CalcOriginDestMatrix_adunit(void); //added function to calculate origin destination matrix: ggilani 28/01/15

///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** /////
///// ***** ///// ***** ///// ***** ///// ***** ///// ***** GLOBAL VARIABLES (some structures in CovidSim.h file and some containers) - memory allocated later.
///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** /////

Param P;
Person* Hosts;
std::vector<PersonQuarantine> HostsQuarantine;
Household* Households;
PopVar State, StateT[MAX_NUM_THREADS];
Cell* Cells; // Cells[i] is the i'th cell
Cell ** CellLookup; // CellLookup[i] is a pointer to the i'th populated cell
Microcell* Mcells, ** McellLookup;
std::vector<uint16_t> mcell_country;
Place** Places;
AdminUnit AdUnits[MAX_ADUNITS];
//// Time Series defs:
//// TimeSeries is an array of type results, used to store (unsurprisingly) a time series of every quantity in results. Mostly used in RecordSample.
//// TSMeanNE and TSVarNE are the mean and variance of non-extinct time series. TSMeanE and TSVarE are the mean and variance of extinct time series. TSMean and TSVar are pointers that point to either extinct or non-extinct.
Results* TimeSeries, * TSMean, * TSVar, * TSMeanNE, * TSVarNE, * TSMeanE, * TSVarE; //// TimeSeries used in RecordSample, RecordInfTypes, SaveResults. TSMean and TSVar
Airport* Airports;
BitmapHeader* bmh;
//added declaration of pointer to events log: ggilani - 10/10/2014
Events* InfEventLog;
int nEvents;

double inftype[INFECT_TYPE_MASK], inftype_av[INFECT_TYPE_MASK], infcountry[MAX_COUNTRIES], infcountry_av[MAX_COUNTRIES], infcountry_num[MAX_COUNTRIES];
double indivR0[MAX_SEC_REC][MAX_GEN_REC], indivR0_av[MAX_SEC_REC][MAX_GEN_REC];
double inf_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], denom_household[MAX_HOUSEHOLD_SIZE + 1];
double inf_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], AgeDist[NUM_AGE_GROUPS], AgeDist2[NUM_AGE_GROUPS];
double case_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], case_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1];
double PropPlaces[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES];
double PropPlacesC[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES], AirTravelDist[MAX_DIST];
double PeakHeightSum, PeakHeightSS, PeakTimeSum, PeakTimeSS;

// These allow up to about 2 billion people per pixel, which should be ample.
int32_t *bmPopulation; // The population in each bitmap pixel. Special value -1 means "country boundary"
int32_t *bmInfected; // The number of infected people in each bitmap pixel.
int32_t *bmRecovered; // The number of recovered people in each bitmap pixel.
int32_t *bmTreated; // The number of treated people in each bitmap pixel.

int OutputTimeStepNumber; //// output timestep index. Global variable used in InitModel and RunModel
int DoInitUpdateProbs;
int InterruptRun = 0; // global variable set to zero at start of RunModel, and possibly modified in CalibrationThresholdCheck 
int PlaceDistDistrib[NUM_PLACE_TYPES][MAX_DIST], PlaceSizeDistrib[NUM_PLACE_TYPES][MAX_PLACE_SIZE];

/* int NumPC,NumPCD; */
const int MAXINTFILE = 10;
std::vector<std::string> InterventionFiles;

// default start value for icdf double arrays (was hardcoded as 100)
const double ICDF_START = 100.0;

void GetInverseCdf(ParamMap fallback, ParamMap params, const char* icdf_name, InverseCdf* inverseCdf, Param* P, double start_value = ICDF_START);

int main(int argc, char* argv[])
{
	Params::alloc_params(&P);
	///// Flags to ensure various parameters have been read; set to false as default.
	std::string pre_param_file, param_file, density_file, load_network_file, save_network_file, air_travel_file, school_file;
	std::string reg_demog_file, fit_file, data_file;
	std::string ad_unit_file, out_density_file, output_file_base;
	std::string snapshot_load_file, snapshot_save_file;

	int StopFit = 0;
	///// Flags to ensure various parameters have been read; set to false as default.
	int GotNR = 0;
	int	GotFI = 0;

	///// Read in command line arguments - lots of things, e.g. random number seeds; (pre)parameter files; binary files; population data; output directory? etc.
	P.FitIter = 0;
	auto parse_snapshot_save_option = [&snapshot_save_file](std::string const& input) {
		auto sep = input.find_first_of(',');
		if (sep == std::string::npos)
		{
			ERR_CRITICAL("Expected argument value to be in the format '<D>,<S>' where <D> is the "
							"timestep interval when a snapshot should be saved and <S> is the "
							"filename in which to save the snapshot");
		}
		parse_double(input.substr(0, sep), P.SnapshotSaveTime);
		parse_read_file(input.substr(sep + 1), snapshot_save_file);
	};

	double cl = (double) clock();

	// Default bitmap format is platform dependent.
#if defined(IMAGE_MAGICK) || defined(_WIN32)
	P.BitmapFormat = BitmapFormats::PNG;
#else
	P.BitmapFormat = BitmapFormats::BMP;
#endif

	// Set parameter defaults - read them in after
	P.PlaceCloseIndepThresh = P.MaxNumThreads = 0;
	P.CaseOrDeathThresholdBeforeAlert_CommandLine = 0;
	P.R0scale = 1.0;
	// added this so that kernel parameters are only changed if input from
	// the command line: ggilani - 15/10/2014
	P.KernelOffsetScale = P.KernelPowerScale = 1.0;
	P.DoLoadSnapshot = 0;

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
	//// **** PARSE COMMAND-LINE ARGS
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****


	CmdLineArgs args;
	args.add_string_option("A", parse_read_file, ad_unit_file, "Administrative Division");
	args.add_string_option("AP", parse_read_file, air_travel_file, "Air travel data file");
	args.add_custom_option("BM", parse_bmp_option, "Bitmap format to use [PNG,BMP]");
	args.add_integer_option("c", P.MaxNumThreads, "Number of threads to use");
	args.add_integer_option("C", P.PlaceCloseIndepThresh, "Sets the P.PlaceCloseIndepThresh parameter");

	/* Wes: Need to allow /CLPxx up to 99. I'll do this naively for now and prevent the help text from
			looking overly verybose. To satisfiy all behaviour, /CLP0: to /CLP9: should do the
			same as /CLP00: to /CLP09: - both valid, as are /CLP10: to /CLP99:
			I am not going to address here what happens if you specify both /CLP05: and /CLP5:
	*/

	for (int i = 0; i <= 99; i++)
	{
		std::string param = "CLP" + std::to_string(i);
		std::string description = "Overwrites #" + std::to_string(i) + " wildcard in parameter file";
		args.add_double_option(param, P.clP[i], description);
		if (i < 10)
		{
			param = "CLP0" + std::to_string(i);
			args.add_double_option(param, P.clP[i], description);
		}
	}

	args.add_string_option("d", parse_read_file, reg_demog_file, "Regional demography file");
	args.add_string_option("D", parse_read_file, density_file, "Population density file");
	args.add_string_option("DT", parse_read_file, data_file, "Likelihood data file");
	args.add_string_option("F", parse_string, fit_file, "Fitting file");
	args.add_integer_option("FI", GotFI, "Initial MCMC iteration");
	args.add_custom_option("I", parse_intervention_file_option, "Intervention file");
	// added Kernel Power and Offset scaling so that it can easily
	// be altered from the command line in order to vary the kernel
	// quickly: ggilani - 15/10/14
	args.add_double_option("KO", P.KernelOffsetScale, "Scales the P.KernelOffsetScale parameter");
	args.add_double_option("KP", P.KernelPowerScale, "Scales the P.KernelPowerScale parameter");
	args.add_string_option("L", parse_read_file, load_network_file, "Network file to load");
	args.add_string_option("LS", parse_read_file, snapshot_load_file, "Snapshot file to load");
	args.add_string_option("M", parse_write_dir, out_density_file, "Output density file");
	args.add_integer_option("NR", GotNR, "Number of realisations");
	args.add_string_option("O", parse_string, output_file_base, "Output file path prefix");
	args.add_string_option("P", parse_read_file, param_file, "Parameter file");
	args.add_string_option("PP", parse_read_file, pre_param_file, "Pre-Parameter file");
	args.add_double_option("R", P.R0scale, "R0 scaling");
	args.add_string_option("s", parse_read_file, school_file, "School file");
	args.add_string_option("S", parse_write_dir, save_network_file, "Network file to save");
	args.add_custom_option("SS", parse_snapshot_save_option, "Interval and file to save snapshots [double,string]");
	args.add_integer_option("T", P.CaseOrDeathThresholdBeforeAlert_CommandLine, "Sets the P.CaseOrDeathThresholdBeforeAlert parameter");
	args.parse(argc, argv, P);

    // Check if S and L options were both specified (can only be one)
	if (!save_network_file.empty() && !load_network_file.empty())
	{
		std::cerr << "Specifying both /L and /S is not allowed" << std::endl;
		args.print_detailed_help_and_exit();
	}

	// Check if P or O were not specified
	if (param_file.empty() || output_file_base.empty())
	{
		std::cerr << "Missing /P and /O arguments which are required" << std::endl;
		args.print_detailed_help_and_exit();
	}

	std::cerr << "Param=" << param_file << "\nOut=" << output_file_base << "\nDens=" << density_file << std::endl;
	Files::xfprintf_stderr("Bitmap Format = *.%s\n", P.BitmapFormat == BitmapFormats::PNG ? "png" : "bmp");
	Files::xfprintf_stderr("sizeof(int)=%i sizeof(long)=%i sizeof(float)=%i sizeof(double)=%i sizeof(unsigned short int)=%i sizeof(int *)=%i\n", (int)sizeof(int), (int)sizeof(long), (int)sizeof(float), (int)sizeof(double), (int)sizeof(unsigned short int), (int)sizeof(int*));

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
	//// **** SET UP OMP / THREADS
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

#ifdef _OPENMP
	P.NumThreads = omp_get_max_threads();
	if ((P.MaxNumThreads > 0) && (P.MaxNumThreads < P.NumThreads)) P.NumThreads = P.MaxNumThreads;
	if (P.NumThreads > MAX_NUM_THREADS)
	{
		Files::xfprintf_stderr("Assigned number of threads (%d) > MAX_NUM_THREADS (%d)\n", P.NumThreads, MAX_NUM_THREADS);
		P.NumThreads = MAX_NUM_THREADS;
	}
	Files::xfprintf_stderr("Using %d threads\n", P.NumThreads);
	omp_set_num_threads(P.NumThreads);

#pragma omp parallel default(shared)
	{
		Files::xfprintf_stderr("Thread %i initialised\n", omp_get_thread_num());
	}
	/* Files::xfprintf_stderr(int=%i\tfloat=%i\tdouble=%i\tint *=%i\n",(int) sizeof(int),(int) sizeof(float),(int) sizeof(double),(int) sizeof(int *));	*/
#else
	P.NumThreads = 1;
#endif
	if (pre_param_file.empty())
	{
		pre_param_file = std::string(".." DIRECTORY_SEPARATOR "Pre_") + param_file;
	}

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
	//// **** READ IN PARAMETERS, DATA ETC.
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

	P.NumRealisations = GotNR;
	ReadParams(param_file, pre_param_file, ad_unit_file);
	if (P.DoAirports)
	{
		if (air_travel_file.empty()) {
			ERR_CRITICAL("Parameter file indicated airports should be used but '/AP' file was not given");
		}
		ReadAirTravel(air_travel_file, output_file_base);
	}

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
	//// **** INITIALIZE
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

	///// initialize model (for all realisations).
	SetupModel(density_file, out_density_file, load_network_file, save_network_file, school_file, reg_demog_file, output_file_base);
	InitTransmissionCoeffs();
	for (int i = 0; i < MAX_ADUNITS; i++) AdUnits[i].NI = 0;
	for (auto const& int_file : InterventionFiles)
		ReadInterventions(int_file);

	Files::xfprintf_stderr("Model setup in %lf seconds\n", ((double) clock() - cl) / CLOCKS_PER_SEC);


	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
	//// **** RUN MODEL
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

	std::string output_file_base_f = output_file_base; // output_file_base_f remembers the original, as output_file_base changes with fitting.
	std::string output_file; // Historically, this was global, and was used for all save...(void) type functions.

	do
	{
		P.FitIter++;
		if (!fit_file.empty())
		{
			if (GotFI)
			{
				P.FitIter = GotFI;
				GotFI = 0;
				P.nextRunSeed1 = P.runSeed1;
				P.nextRunSeed2 = P.runSeed2;
			}
			StopFit = ReadFitIter(fit_file);
			if (!StopFit)
			{
				ReadParams(param_file, pre_param_file, ad_unit_file);
				if (!P.FixLocalBeta) InitTransmissionCoeffs();
				output_file_base = output_file_base_f + ".f" + std::to_string(P.FitIter);
			}
		}
		else
			StopFit = 1;
		if ((fit_file.empty()) || (!StopFit))
		{
			P.NRactE = P.NRactNE = 0;
			ResetTimeSeries();
			for (int Realisation = 0; (Realisation < P.NumRealisations) && (P.NRactNE < P.NumNonExtinctRealisations); Realisation++)
			{
				if (P.NumRealisations > 1)
				{
					output_file = output_file_base + "." + std::to_string(Realisation);
					Files::xfprintf_stderr("Realisation %i of %i  (time=%lf nr_ne=%i)\n", Realisation + 1, P.NumRealisations, ((double)(clock() - cl)) / CLOCKS_PER_SEC, P.NRactNE);
				}
				///// Set and save seeds
				if (((Realisation == 0) && (P.FitIter == 1)) || (P.ResetSeeds && P.KeepSameSeeds))
				{
					P.nextRunSeed1 = P.runSeed1;
					P.nextRunSeed2 = P.runSeed2;
				}
				if (P.ResetSeeds) {
					//save these seeds to file
					SaveRandomSeeds(output_file);
				}
				int32_t thisRunSeed1, thisRunSeed2;
				int ContCalib, ModelCalibLoop = 0;
				P.StopCalibration = P.ModelCalibIteration = ModelCalibLoop =  0;
				do
				{  // has been interrupted to reset holiday time. Note that this currently only happens in the first run, regardless of how many realisations are being run.
					if ((P.ModelCalibIteration % 14 == 0) && (ModelCalibLoop < 4))
					{
						thisRunSeed1 = P.nextRunSeed1;
						thisRunSeed2 = P.nextRunSeed2;
						setall(&P.nextRunSeed1, &P.nextRunSeed2);
						P.HolidaysStartDay_SimTime = 0; // needed for calibration to work for multiple realisations
						P.CaseOrDeathThresholdBeforeAlert = P.CaseOrDeathThresholdBeforeAlert_Fixed; // needed for calibration to work for multiple realisations
						if (!P.DoNoCalibration) P.SeedingScaling = 1.0; // needed for calibration to work for multiple realisations
						P.ModelCalibIteration = 0;  // needed for calibration to work for multiple realisations
						ModelCalibLoop++;
					}
					else
					{
						int32_t tmp1 = thisRunSeed1;
						int32_t tmp2 = thisRunSeed2;
						setall(&tmp1, &tmp2);  // reset random number seeds to generate same run again after calibration.
					}
					InitModel(Realisation);
					if (!snapshot_load_file.empty()) LoadSnapshot(snapshot_load_file);
					ContCalib = RunModel(Realisation, snapshot_save_file, snapshot_load_file, output_file_base);
				}
				while (ContCalib);
				if (!data_file.empty()) CalcLikelihood(Realisation, data_file, output_file_base);
				if (P.OutputNonSummaryResults)
				{
					if (((!TimeSeries[P.NumOutputTimeSteps - 1].extinct) || (!P.OutputOnlyNonExtinct)) && (P.OutputEveryRealisation))
					{
						SaveResults(output_file);
					}
				}
				if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun == 1))
				{
					SaveEvents(output_file);
				}
			}
			output_file = output_file_base + ".avNE";
			SaveSummaryResults(output_file);

			//Calculate origin destination matrix if needed
			if ((P.DoAdUnits) && (P.DoOriginDestinationMatrix))
			{
				CalcOriginDestMatrix_adunit();
				SaveOriginDestMatrix(output_file);
			}

			P.NRactual = P.NRactNE;
			TSMean = TSMeanNE; TSVar = TSVarNE;
			if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun == 0))
			{
				SaveEvents(output_file);
			}

			SaveSummaryResults(output_file);
			P.NRactual = P.NRactE;
			//TSMean = TSMeanE; TSVar = TSVarE;
			//Files::xsprintf(OutFile, "%s.avE", OutFileBase);
			//SaveSummaryResults();

			Bitmap_Finalise();

			Files::xfprintf_stderr("Extinction in %i out of %i runs\n", P.NRactE, P.NRactNE + P.NRactE);
			Files::xfprintf_stderr("Model ran in %lf seconds\n", ((double)clock() - cl) / CLOCKS_PER_SEC);
			Files::xfprintf_stderr("Model finished\n");
		}
	}
	while (!StopFit);
}

void parse_bmp_option(std::string const& input) {
	// make copy and convert input to lowercase
	std::string input_copy = input;
	std::transform(input_copy.begin(), input_copy.end(), input_copy.begin(), [](unsigned char c){ return std::tolower(c); });

	if (input_copy.compare("png") == 0) {
#if defined(IMAGE_MAGICK) || defined(_WIN32)
		P.BitmapFormat = BitmapFormats::PNG;
#else
		ERR_CRITICAL("PNG Bitmaps not supported - please build with Image Magic or WIN32 support");
#endif
	}
	else if (input_copy.compare("bmp") == 0) {
		P.BitmapFormat = BitmapFormats::BMP;
	}
	else {
		ERR_CRITICAL_FMT("Unrecognised bitmap format: %s", input_copy.c_str());
	}
}

void parse_intervention_file_option(std::string const& input)
{
	std::string output;
	parse_read_file(input, output);
	InterventionFiles.emplace_back(output);
}

void ReadParams(std::string const& ParamFile, std::string const& PreParamFile, std::string const& AdUnitFile)
{
	Files::xfprintf_stderr("Reading Params...\n");	fflush(stderr);

	double s, t, AgeSuscScale;
	int i, j, j1, j2, k, f, nc, na;

	char* CountryNameBuf = new char[128 * MAX_COUNTRIES];
	char* AdunitListNamesBuf = new char[128 * MAX_ADUNITS];
	char** CountryNames = new char*[MAX_COUNTRIES];
	char** AdunitListNames = new char*[MAX_ADUNITS];
	
	AgeSuscScale = 1.0;
	ParamMap params = Params::read_params_map(ParamFile.c_str());
	ParamMap pre_params = Params::read_params_map(PreParamFile.c_str());
	ParamMap adm_params = Params::read_params_map(AdUnitFile.c_str());
	
	if (P.FitIter == 0)
	{
		for (i = 0; i < MAX_COUNTRIES; i++) {
			CountryNames[i] = CountryNameBuf + INT64_C(128) * i;
			CountryNames[i][0] = 0;
		}
		for (i = 0; i < MAX_ADUNITS; i++) {
			AdunitListNames[i] = AdunitListNamesBuf + INT64_C(128) * i;
			AdunitListNames[i][0] = 0;
		}
		for (i = 0; i < 100; i++) P.clP_copies[i] = 0;
		P.LongitudeCutLine = Params::get_double(params, adm_params, "Longitude cut line", -360.0, &P);
		P.ModelTimeStep = Params::req_double(params, pre_params, "Update timestep", &P);
		P.OutputTimeStep = Params::req_double(params, pre_params, "Sampling timestep", &P);
		if (P.ModelTimeStep > P.OutputTimeStep) ERR_CRITICAL("Update step must be smaller than sampling step\n");
		t = ceil(P.OutputTimeStep / P.ModelTimeStep - 1e-6);
		P.NumModelTimeStepsPerOutputTimeStep = (int)t;
		P.ModelTimeStep = P.OutputTimeStep / t;
		P.TimeStepsPerDay = ceil(1.0 / P.ModelTimeStep - 1e-6);
		Files::xfprintf_stderr("Update step = %lf\nSampling step = %lf\nUpdates per sample=%i\nTimeStepsPerDay=%lf\n", P.ModelTimeStep, P.OutputTimeStep, P.NumModelTimeStepsPerOutputTimeStep, P.TimeStepsPerDay);
		P.SimulationDuration = Params::req_double(params, pre_params, "Sampling time", &P);
		P.NumOutputTimeSteps = 1 + (int)ceil(P.SimulationDuration / P.OutputTimeStep);
		P.PopSize = Params::req_int(pre_params, adm_params, "Population size", &P);
		if (P.NumRealisations == 0)
		{
			P.NumRealisations = Params::req_int(params, pre_params, "Number of realisations", &P);
			P.NumNonExtinctRealisations = Params::get_int(params, pre_params, "Number of non-extinct realisations", P.NumRealisations, &P);
		}
		else
		{
			P.NumNonExtinctRealisations = P.NumRealisations;
		}
		P.SmallEpidemicCases = Params::get_int(params, pre_params, "Maximum number of cases defining small outbreak", -1, &P);

		P.NumCells = -1;
		P.NMCL = Params::req_int(params, pre_params, "Number of micro-cells per spatial cell width", &P);
		//added parameter to reset seeds after every run
		P.ResetSeeds = Params::get_int(params, pre_params, "Reset seeds for every run", 0, &P);
		if (P.ResetSeeds != 0)
		{
			P.KeepSameSeeds = Params::get_int(params, pre_params, "Keep same seeds for every run", 0, &P); //added this to control which seeds are used: ggilani 27/11/19
		}
		P.ResetSeedsPostIntervention = Params::get_int(params, pre_params, "Reset seeds after intervention", 0, &P);
		if (P.ResetSeedsPostIntervention != 0)
		{
			P.TimeToResetSeeds = Params::get_int(params, pre_params, "Time to reset seeds after intervention", 1000000, &P);
		}
		P.DoHouseholds = Params::get_int(pre_params, adm_params, "Include households", 1, &P);

		P.KernelLookup.size_ = Params::get_int(pre_params, adm_params, "Kernel resolution", 4000000, &P);
		if (P.KernelLookup.size_ < 2000000)
		{
			ERR_CRITICAL_FMT("[Kernel resolution] needs to be at least 2000000 - not %d", P.KernelLookup.size_);
		}
		P.KernelLookup.expansion_factor_ = Params::get_int(pre_params, adm_params, "Kernel higher resolution factor", P.KernelLookup.size_ / 1600, &P);
		if (P.KernelLookup.expansion_factor_ < 1 || P.KernelLookup.expansion_factor_ >= P.KernelLookup.size_)
		{
			ERR_CRITICAL_FMT("[Kernel higher resolution factor] needs to be in range [1, P.NKR = %d) - not %d", P.KernelLookup.size_, P.KernelLookup.expansion_factor_);
		}
	}
	P.OutputAge = Params::get_int(params, pre_params, "OutputAge", 1, &P);
	P.OutputSeverity = Params::get_int(params, pre_params, "OutputSeverity", 1, &P);
	P.OutputSeverityAdminUnit= Params::get_int(params, pre_params, "OutputSeverityAdminUnit", 1, &P);
	P.OutputSeverityAge = Params::get_int(params, pre_params, "OutputSeverityAge", 1, &P);
	P.OutputAdUnitAge = Params::get_int(params, pre_params, "OutputAdUnitAge", 0, &P);
	P.OutputR0 = Params::get_int(params, pre_params, "OutputR0", 0, &P);
	P.OutputControls = Params::get_int(params, pre_params, "OutputControls", 0, &P);
	P.OutputCountry = Params::get_int(params, pre_params, "OutputCountry", 0, &P);
	P.OutputAdUnitVar = Params::get_int(params, pre_params, "OutputAdUnitVar", 0, &P);
	P.OutputHousehold = Params::get_int(params, pre_params, "OutputHousehold", 0, &P);
	P.OutputInfType = Params::get_int(params, pre_params, "OutputInfType", 0, &P);
	P.OutputNonSeverity = Params::get_int(params, pre_params, "OutputNonSeverity", 0, &P);
	P.OutputNonSummaryResults = Params::get_int(params, pre_params, "OutputNonSummaryResults", 0, &P);

	if (P.DoHouseholds != 0)
	{
		Params::req_double_vec(pre_params, adm_params, "Household size distribution", P.HouseholdSizeDistrib[0], MAX_HOUSEHOLD_SIZE, &P);
		P.HouseholdTrans = Params::req_double(params, pre_params, "Household attack rate", &P);
		P.HouseholdTransPow = Params::req_double(params, pre_params, "Household transmission denominator power", &P);
		P.DoCorrectAgeDist = Params::get_int(pre_params, adm_params, "Correct age distribution after household allocation to exactly match specified demography", 0, &P);
	}
	else
	{
		P.HouseholdTrans = 0.0;
		P.HouseholdTransPow = 1.0;
		P.HouseholdSizeDistrib[0][0] = 1.0;
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P.HouseholdSizeDistrib[0][i] = 0;
	}
	if (P.FitIter == 0)
	{
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P.HouseholdSizeDistrib[0][i] = P.HouseholdSizeDistrib[0][i] + P.HouseholdSizeDistrib[0][i - 1];
		P.HouseholdDenomLookup[0] = 1.0;
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
		    P.HouseholdDenomLookup[i] = 1 / pow(((double)(INT64_C(1) + i)), P.HouseholdTransPow);
		P.DoAdUnits = Params::get_int(pre_params, adm_params, "Include administrative units within countries", 1, &P);
		P.CountryDivisor = Params::get_int(pre_params, adm_params, "Divisor for countries", 1, &P);
		if (P.DoAdUnits != 0)
		{
		
			char** AdunitNames = (char**) Memory::xmalloc(3 * ADUNIT_LOOKUP_SIZE * sizeof(char*));
			char* AdunitNamesBuf = (char*) Memory::xmalloc(3 * ADUNIT_LOOKUP_SIZE * 360 * sizeof(char));

			for (i = 0; i < ADUNIT_LOOKUP_SIZE; i++)
			{
				P.AdunitLevel1Lookup[i] = -1;
				AdunitNames[3 * i] = AdunitNamesBuf + INT64_C(3) * i * 360;
				AdunitNames[3 * i + 1] = AdunitNamesBuf + INT64_C(3) * i * 360 + 60;
				AdunitNames[3 * i + 2] = AdunitNamesBuf + INT64_C(3) * i * 360 + 160;
			}
			P.AdunitLevel1Divisor = Params::get_int(pre_params, adm_params, "Divisor for level 1 administrative units", 1, &P);
			P.AdunitLevel1Mask = Params::get_int(pre_params, adm_params, "Mask for level 1 administrative units", 1000000000, &P);
			na = Params::req_string_vec(pre_params, adm_params, "Codes and country/province names for admin units", AdunitNames, 3 * ADUNIT_LOOKUP_SIZE, &P) / 3;
			nc = Params::get_int(pre_params, adm_params, "Number of countries to include", 0, &P);
			if ((na > 0) && (nc > 0))
			{
				P.DoAdunitBoundaries = (nc > 0);
				nc = abs(nc);
				Params::req_string_vec(pre_params, adm_params, "List of names of countries to include", CountryNames, nc, &P);
				P.NumAdunits = 0;
				for (i = 0; i < na; i++)
					for (j = 0; j < nc; j++)
						if ((AdunitNames[3 * i + 1][0]) && (!strcmp(AdunitNames[3 * i + 1], CountryNames[j])) && (atoi(AdunitNames[3 * i]) != 0))
						{
							AdUnits[P.NumAdunits].id = atoi(AdunitNames[3 * i]);
							P.AdunitLevel1Lookup[(AdUnits[P.NumAdunits].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor] = P.NumAdunits;
							if (strlen(AdunitNames[3 * i + 1]) < 100) strcpy(AdUnits[P.NumAdunits].cnt_name, AdunitNames[3 * i + 1]);
							if (strlen(AdunitNames[3 * i + 2]) < 200) strcpy(AdUnits[P.NumAdunits].ad_name, AdunitNames[3 * i + 2]);
							//						Files::xfprintf_stderr("%i %s %s ## ",AdUnits[P.NumAdunits].id,AdUnits[P.NumAdunits].cnt_name,AdUnits[P.NumAdunits].ad_name);
							P.NumAdunits++;
						}
			}
			else
			{
				P.NumAdunits = Params::get_int(pre_params, adm_params, "Number of level 1 administrative units to include", 0, &P);
				if (P.NumAdunits > 0)
				{
					P.DoAdunitBoundaries = 1;
					if (P.NumAdunits > MAX_ADUNITS) ERR_CRITICAL("MAX_ADUNITS too small.\n");
					Params::req_string_vec(pre_params, adm_params, "List of level 1 administrative units to include", AdunitListNames, P.NumAdunits, &P);
					na = P.NumAdunits;
					for (i = 0; i < P.NumAdunits; i++)
					{
						f = 0;
						if (na > 0)
						{
							for (j = 0; (j < na) && (!f); j++) f = (!strcmp(AdunitNames[3 * j + 2], AdunitListNames[i]));
							if (f) k = atoi(AdunitNames[3 * (j - 1)]);
						}
						if ((na == 0) || (!f)) k = atoi(AdunitListNames[i]);
						AdUnits[i].id = k;
						P.AdunitLevel1Lookup[(k % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor] = i;
						for (j = 0; j < na; j++)
							if (atoi(AdunitNames[3 * j]) == k)
							{
								if (strlen(AdunitNames[3 * j + 1]) < 100) strcpy(AdUnits[i].cnt_name, AdunitNames[3 * j + 1]);
								if (strlen(AdunitNames[3 * j + 2]) < 200) strcpy(AdUnits[i].ad_name, AdunitNames[3 * j + 2]);
								j = na;
							}
					}
				}
				else
					P.DoAdunitBoundaries = 0;
			}
			Memory::xfree(AdunitNames);
			Memory::xfree(AdunitNamesBuf);

			P.DoAdunitOutput = Params::get_int(params, pre_params, "Output incidence by administrative unit", 0, &P);
			P.DoAdunitBoundaryOutput = Params::get_int(pre_params, adm_params, "Draw administrative unit boundaries on maps", 0, &P);
			P.DoCorrectAdunitPop = Params::get_int(pre_params, adm_params, "Correct administrative unit populations", 0, &P);
			P.DoSpecifyPop = Params::get_int(pre_params, adm_params, "Fix population size at specified value", 0, &P);
			Files::xfprintf_stderr("Using %i administrative units\n", P.NumAdunits);
			P.AdunitBitmapDivisor = Params::get_int(pre_params, adm_params, "Divisor for administrative unit codes for boundary plotting on bitmaps", 1, &P);
			P.DoOutputPlaceDistForOneAdunit = Params::get_int(params, pre_params, "Only output household to place distance distribution for one administrative unit", 0, &P);
			if (P.DoOutputPlaceDistForOneAdunit != 0)
			{
				P.DoOutputPlaceDistForOneAdunit = Params::get_int(params, pre_params, "Administrative unit for which household to place distance distribution to be output", 0, &P);
			}
		}
		else

		{
			P.DoAdunitBoundaries = P.DoAdunitBoundaryOutput = P.DoAdunitOutput = P.DoCorrectAdunitPop = P.DoSpecifyPop = 0;
			P.AdunitLevel1Divisor = 1; P.AdunitLevel1Mask = 1000000000;
			P.AdunitBitmapDivisor = P.AdunitLevel1Divisor;
		}
	}

	Files::xfprintf_stderr("Reading Params 1\n"); fflush(stderr);


	P.DoAge = Params::get_int(pre_params, adm_params, "Include age", 1, &P);
	if (P.DoAge == 0)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++) {
			P.PropAgeGroup[0][i] = 1.0 / NUM_AGE_GROUPS;
			P.InitialImmunity[i] = 0;
			P.AgeInfectiousness[i] = P.AgeSusceptibility[i] = 1;
			P.RelativeSpatialContact[i] = P.RelativeTravelRate[i] = 1.0;
		}
	}
	else
	{

		P.DoPartialImmunity = Params::get_int(params, pre_params, "Initial immunity acts as partial immunity", 1, &P);
		if ((P.DoHouseholds != 0) && (P.DoPartialImmunity == 0))
		{
			P.DoWholeHouseholdImmunity = Params::get_int(params, pre_params, "Initial immunity applied to all household members", 0, &P);
		}
		else
			P.DoWholeHouseholdImmunity = 0;
		Params::get_double_vec(params, pre_params, "Initial immunity profile by age", P.InitialImmunity, NUM_AGE_GROUPS, 0, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "Relative spatial contact rates by age", P.RelativeSpatialContact, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, &P);
		if (Params::get_int(params, pre_params, "Apply spatial contact rates by age to susceptibles as well as infecteds", 0, &P) != 0)
			for (i = 0; i < NUM_AGE_GROUPS; i++)	P.RelativeSpatialContactSusc[i] = P.RelativeSpatialContact[i];
		else
			for (i = 0; i < NUM_AGE_GROUPS; i++)	P.RelativeSpatialContactSusc[i] = 1.0;
		Params::get_double_vec(params, pre_params, "Age-dependent infectiousness", P.AgeInfectiousness, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "Age-dependent susceptibility", P.AgeSusceptibility, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, &P);
		Params::req_double_vec(pre_params, adm_params, "Age distribution of population", P.PropAgeGroup[0], NUM_AGE_GROUPS, &P);
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			t += P.PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.PropAgeGroup[0][i] /= t;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			if (P.AgeSusceptibility[i] > t) t = P.AgeSusceptibility[i];  //peak susc has to be 1
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.AgeSusceptibility[i] /= t;
		AgeSuscScale = t;
		if (P.DoHouseholds) P.HouseholdTrans *= AgeSuscScale;
		Params::get_double_vec(pre_params, adm_params, "Relative travel rates by age", P.RelativeTravelRate, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, &P);

		//if (!GetInputParameter2(params, pre_params, "WAIFW matrix", "%lf", (void*)P.WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0))
		if (!Params::param_found(params, pre_params, "WAIFW matrix"))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					P.WAIFW_Matrix[i][j] = 1.0;
		}
		else
		{
		  Params::get_double_matrix(params, pre_params, "WAIFW matrix", P.WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 1.0, &P);

			/* WAIFW matrix needs to be scaled to have max value of 1.
			1st index of matrix specifies host being infected, second the infector.
			Overall age variation in infectiousness/contact rates/susceptibility should be factored
			out of WAIFW_matrix and put in Age dep infectiousness/susceptibility for efficiency. */
			t = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					if (P.WAIFW_Matrix[i][j] > t) t = P.WAIFW_Matrix[i][j];
			if (t > 0)
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P.WAIFW_Matrix[i][j] /= t;
			}
			else
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P.WAIFW_Matrix[i][j] = 1.0;
			}
		}
		if (!Params::param_found(params, pre_params, "WAIFW matrix spatial infections only"))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					P.WAIFW_Matrix_SpatialOnly[i][j] = 1.0;

			P.Got_WAIFW_Matrix_Spatial = 0; 
		}
		else
		{
		  Params::get_double_matrix(params, pre_params, "WAIFW matrix spatial infections only", P.WAIFW_Matrix_SpatialOnly, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0, &P);
			P.Got_WAIFW_Matrix_Spatial = 1;
			/* WAIFW matrix needs to be scaled to have max value of 1.
			1st index of matrix specifies host being infected, second the infector.
			Overall age variation in infectiousness/contact rates/susceptibility should be factored
			out of WAIFW_matrix and put in Age dep infectiousness/susceptibility for efficiency. */

			double Maximum = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					if (P.WAIFW_Matrix_SpatialOnly[i][j] > Maximum) Maximum = P.WAIFW_Matrix_SpatialOnly[i][j];
			if (Maximum > 0)
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P.WAIFW_Matrix_SpatialOnly[i][j] /= Maximum;
			}
			else
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P.WAIFW_Matrix_SpatialOnly[i][j] = 1.0;
			}

		}

		P.DoDeath = 0;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)	t += P.AgeInfectiousness[i] * P.PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)	P.AgeInfectiousness[i] /= t;
	}

	Files::xfprintf_stderr("Reading Params 2\n"); fflush(stderr);


	if (P.FitIter == 0)
	{
		P.DoSpatial = Params::get_int(pre_params, adm_params, "Include spatial transmission", 1, &P);
		P.MoveKernel.type_ = Params::req_int(pre_params, adm_params, "Kernel type", &P);
		P.MoveKernel.scale_ = Params::req_double(pre_params, adm_params, "Kernel scale", &P);
		if (P.KernelOffsetScale != 1)
		{
			P.MoveKernel.scale_ *= P.KernelOffsetScale;
		}
		P.MoveKernel.p3_ = Params::get_double(pre_params, adm_params, "Kernel 3rd param", 0, &P);
		P.MoveKernel.p4_ = Params::get_double(pre_params, adm_params, "Kernel 4th param", 0, &P);
		P.MoveKernel.shape_ = Params::get_double(pre_params, adm_params, "Kernel Shape", 1.0, &P);
		if (P.KernelPowerScale != 1)
		{
			P.MoveKernel.shape_ *= P.KernelPowerScale;
		}
		P.AirportKernel.type_ = Params::get_int(pre_params, adm_params, "Airport Kernel Type", P.MoveKernel.type_, &P);
		P.AirportKernel.scale_ = Params::get_double(pre_params, adm_params, "Airport Kernel Scale", P.MoveKernel.scale_, &P);
		P.AirportKernel.shape_ = Params::get_double(pre_params, adm_params, "Airport Kernel Shape", P.MoveKernel.shape_, &P);
		P.AirportKernel.p3_ = Params::get_double(pre_params, adm_params, "Airport Kernel 3rd param", P.MoveKernel.p3_, &P);
		P.AirportKernel.p4_ = Params::get_double(pre_params, adm_params, "Airport Kernel 4th param", P.MoveKernel.p4_, &P);

		P.DoPlaces = Params::get_int(pre_params, adm_params, "Include places", 1, &P);
		if (P.DoPlaces != 0)
		{
			P.PlaceTypeNum = Params::get_int(pre_params, adm_params, "Number of types of places", 0, &P);
			if (P.PlaceTypeNum == 0) P.DoPlaces = P.DoAirports = 0;
		}
		else
			P.PlaceTypeNum = P.DoAirports = 0;
	}

	Files::xfprintf_stderr("Reading Params 3\n"); fflush(stderr);


	if (P.DoPlaces != 0)
	{
		P.CareHomeResidentHouseholdScaling = Params::get_double(pre_params, adm_params, "Scaling of household contacts for care home residents", 1.0, &P);
		P.CareHomeResidentSpatialScaling = Params::get_double(pre_params, adm_params, "Scaling of spatial contacts for care home residents", 1.0, &P);
		P.CareHomeResidentPlaceScaling = Params::get_double(pre_params, adm_params, "Scaling of between group (home) contacts for care home residents", 1.0, &P);
		P.CareHomeWorkerGroupScaling = Params::get_double(pre_params, adm_params, "Scaling of within group (home) contacts for care home workers", 1.0, &P);
		P.CareHomeRelProbHosp = Params::get_double(pre_params, adm_params, "Relative probability that care home residents are hospitalised", 1.0, &P);

		if (P.FitIter == 0)
		{
			if (P.PlaceTypeNum > NUM_PLACE_TYPES) ERR_CRITICAL("Too many place types\n");
			P.CareHomePlaceType = Params::get_int(pre_params, adm_params, "Place type number for care homes", -1, &P);
			P.CareHomeAllowInitialInfections = Params::get_int(pre_params, adm_params, "Allow initial infections to be in care homes", 0, &P);
			P.CareHomeResidentMinimumAge = Params::get_int(pre_params, adm_params, "Minimum age of care home residents", 1000, &P);

			Params::req_int_vec(pre_params, adm_params, "Minimum age for age group 1 in place types", P.PlaceTypeAgeMin, P.PlaceTypeNum, &P);
			Params::req_int_vec(pre_params, adm_params, "Maximum age for age group 1 in place types", P.PlaceTypeAgeMax, P.PlaceTypeNum, &P);
			Params::req_double_vec(pre_params, adm_params, "Proportion of age group 1 in place types", P.PlaceTypePropAgeGroup, P.PlaceTypeNum, &P);

			Files::xfprintf_stderr("Reading Params 3.1\n"); fflush(stderr);


			if (!Params::param_found(pre_params, adm_params, "Proportion of age group 2 in place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P.PlaceTypePropAgeGroup2[i] = 0;
					P.PlaceTypeAgeMin2[i] = 0;
					P.PlaceTypeAgeMax2[i] = 1000;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Proportion of age group 2 in place types", P.PlaceTypePropAgeGroup2, P.PlaceTypeNum, &P);
			  Params::req_int_vec(pre_params, adm_params, "Minimum age for age group 2 in place types", P.PlaceTypeAgeMin2, P.PlaceTypeNum, &P);
				Params::req_int_vec(pre_params, adm_params, "Maximum age for age group 2 in place types", P.PlaceTypeAgeMax2, P.PlaceTypeNum, &P);
			}
			if (!Params::param_found(pre_params, adm_params, "Proportion of age group 3 in place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P.PlaceTypePropAgeGroup3[i] = 0;
					P.PlaceTypeAgeMin3[i] = 0;
					P.PlaceTypeAgeMax3[i] = 1000;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Proportion of age group 3 in place types", P.PlaceTypePropAgeGroup3, P.PlaceTypeNum, &P);
			  Params::req_int_vec(pre_params, adm_params, "Minimum age for age group 3 in place types", P.PlaceTypeAgeMin3, P.PlaceTypeNum, &P);
				Params::req_int_vec(pre_params, adm_params, "Maximum age for age group 3 in place types", P.PlaceTypeAgeMax3, P.PlaceTypeNum, &P);
			}
			if (!Params::param_found(pre_params, adm_params, "Kernel shape params for place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P.PlaceTypeKernelShape[i] = P.MoveKernel.shape_;
					P.PlaceTypeKernelScale[i] = P.MoveKernel.scale_;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Kernel shape params for place types", P.PlaceTypeKernelShape, P.PlaceTypeNum, &P);
				Params::req_double_vec(pre_params, adm_params, "Kernel scale params for place types", P.PlaceTypeKernelScale, P.PlaceTypeNum, &P);
			}
			if (!Params::param_found(pre_params, adm_params, "Kernel 3rd param for place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P.PlaceTypeKernelP3[i] = P.MoveKernel.p3_;
					P.PlaceTypeKernelP4[i] = P.MoveKernel.p4_;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Kernel 3rd param for place types", P.PlaceTypeKernelP3, P.PlaceTypeNum, &P);
				Params::req_double_vec(pre_params, adm_params, "Kernel 4th param for place types", P.PlaceTypeKernelP4, P.PlaceTypeNum, &P);
			}
			Params::get_int_vec(pre_params, adm_params, "Number of closest places people pick from (0=all) for place types", P.PlaceTypeNearestNeighb, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
			if (P.DoAdUnits != 0)
			{
				Params::get_double_vec(params, pre_params, "Degree to which crossing administrative unit boundaries to go to places is inhibited",
				                       P.InhibitInterAdunitPlaceAssignment, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
			}

			Files::xfprintf_stderr("Reading Params 3.2\n"); fflush(stderr);


			P.DoAirports = Params::get_int(params, pre_params, "Include air travel", 0, &P);
			if (P.DoAirports == 0)
			{
				// Airports disabled => all places are not to do with airports, and we
				// have no hotels.
				P.PlaceTypeNoAirNum = P.PlaceTypeNum;
				P.HotelPlaceType = P.PlaceTypeNum;
			}
			else
			{
				// When airports are activated we must have at least one airport place
				// // and a hotel type.
				P.PlaceTypeNoAirNum = Params::req_int(pre_params, adm_params, "Number of non-airport places", &P);
				P.HotelPlaceType = Params::req_int(pre_params, adm_params, "Hotel place type", &P);
				if (P.PlaceTypeNoAirNum >= P.PlaceTypeNum) {
					ERR_CRITICAL_FMT("[Number of non-airport places] parameter (%d) is greater than number of places (%d).\n", P.PlaceTypeNoAirNum, P.PlaceTypeNum);
				}
				if (P.HotelPlaceType < P.PlaceTypeNoAirNum || P.HotelPlaceType >= P.PlaceTypeNum) {
					ERR_CRITICAL_FMT("[Hotel place type] parameter (%d) not in the range [%d, %d)\n", P.HotelPlaceType, P.PlaceTypeNoAirNum, P.PlaceTypeNum);
				}

				P.AirportTrafficScale = Params::get_double(params, pre_params, "Scaling factor for input file to convert to daily traffic", 1.0, &P);
				P.HotelPropLocal = Params::get_double(params, pre_params, "Proportion of hotel attendees who are local", 0, &P);
				Params::get_double_vec(params, pre_params, "Distribution of duration of air journeys", P.JourneyDurationDistrib, MAX_TRAVEL_TIME, 0, MAX_TRAVEL_TIME, &P);
				if (!Params::param_found(params, pre_params, "Distribution of duration of air journeys"))
				{
					P.JourneyDurationDistrib[0] = 1;
				}
				Params::get_double_vec(params, pre_params, "Distribution of duration of local journeys", P.LocalJourneyDurationDistrib, MAX_TRAVEL_TIME, 0, MAX_TRAVEL_TIME, &P);
				if (!Params::param_found(params, pre_params, "Distribution of duration of local journeys"))
				{
					P.LocalJourneyDurationDistrib[0] = 1;
				}
				P.MeanJourneyTime = P.MeanLocalJourneyTime = 0;
				for (i = 0; i < MAX_TRAVEL_TIME; i++)
				{
					P.MeanJourneyTime += ((double)(i)) * P.JourneyDurationDistrib[i];
					P.MeanLocalJourneyTime += ((double)(i)) * P.LocalJourneyDurationDistrib[i];
				}
				Files::xfprintf_stderr("Mean duration of local journeys = %lf days\n", P.MeanLocalJourneyTime);
				for (i = 1; i < MAX_TRAVEL_TIME; i++)
				{
					P.JourneyDurationDistrib[i] += P.JourneyDurationDistrib[i - 1];
					P.LocalJourneyDurationDistrib[i] += P.LocalJourneyDurationDistrib[i - 1];
				}
				for (i = j1 = j2 = 0; i <= 1024; i++)
				{
					s = ((double)i) / 1024;
					while (P.JourneyDurationDistrib[j1] < s) j1++;
					P.InvJourneyDurationDistrib[i] = j1;
					while (P.LocalJourneyDurationDistrib[j2] < s) j2++;
					P.InvLocalJourneyDurationDistrib[i] = j2;
				}
			}

			Files::xfprintf_stderr("Reading Params 3.3\n"); fflush(stderr);


			Params::req_double_vec(pre_params, adm_params, "Mean size of place types", P.PlaceTypeMeanSize, P.PlaceTypeNum, &P);
			Params::req_double_vec(pre_params, adm_params, "Param 1 of place group size distribution", P.PlaceTypeGroupSizeParam1, P.PlaceTypeNum, &P);
			Params::get_double_vec(pre_params, adm_params, "Power of place size distribution", P.PlaceTypeSizePower, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);

			//added to enable lognormal distribution - ggilani 09/02/17
			Params::get_double_vec(pre_params, adm_params, "Standard deviation of place size distribution", P.PlaceTypeSizeSD, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
			Params::get_double_vec(pre_params, adm_params, "Offset of place size distribution", P.PlaceTypeSizeOffset, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
			Params::get_double_vec(pre_params, adm_params, "Maximum of place size distribution", P.PlaceTypeSizeMax, P.PlaceTypeNum, 1e20, NUM_PLACE_TYPES, &P);
			Params::get_double_vec(pre_params, adm_params, "Minimum of place size distribution", P.PlaceTypeSizeMin, P.PlaceTypeNum, 1.0, NUM_PLACE_TYPES, &P);
			Params::get_int_vec(pre_params, adm_params, "Kernel type for place types", P.PlaceTypeKernelType, P.PlaceTypeNum, P.MoveKernel.type_, NUM_PLACE_TYPES, &P);

			Files::xfprintf_stderr("Reading Params 3.4\n"); fflush(stderr);

			Params::get_double_vec(pre_params, adm_params, "Place overlap matrix", P.PlaceExclusivityMatrix, P.PlaceTypeNum * P.PlaceTypeNum, 0, P.PlaceTypeNum* P.PlaceTypeNum, &P);
			Files::xfprintf_stderr("Reading Params 3.5\n"); fflush(stderr);
			if (!Params::param_found(pre_params, adm_params, "Place overlap matrix"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)                                 // get_double_vec will set the zeroes if missing;
					P.PlaceExclusivityMatrix[i * (NUM_PLACE_TYPES + 1)] = 1;          // this line sets the diagonal to 1 (identity matrix)
			}
			Files::xfprintf_stderr("Reading Params 3.6\n"); fflush(stderr);

		}
		/* Note P.PlaceExclusivityMatrix not used at present - places assumed exclusive (each person belongs to 0 or 1 place) */

		Params::req_double_vec(params, pre_params, "Proportion of between group place links", P.PlaceTypePropBetweenGroupLinks, P.PlaceTypeNum, &P);
		Files::xfprintf_stderr("Reading Params 3.7\n"); fflush(stderr);
		Params::req_double_vec(params, pre_params, "Relative transmission rates for place types", P.PlaceTypeTrans, P.PlaceTypeNum, &P);
		Files::xfprintf_stderr("Reading Params 3.8\n"); fflush(stderr);
		for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceTypeTrans[i] *= AgeSuscScale;
	}
	Params::get_double_vec(params, pre_params, "Daily seasonality coefficients", P.Seasonality, DAYS_PER_YEAR, 1, DAYS_PER_YEAR, &P);
	if (!Params::param_found(params, pre_params, "Daily seasonality coefficients"))
	{
		P.DoSeasonality = 0;
	}
	else
	{
		P.DoSeasonality = 1;
		s = 0;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			s += P.Seasonality[i];
		s += 1e-20;
		s /= DAYS_PER_YEAR;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			P.Seasonality[i] /= s;
	}
	P.NumSeedLocations = Params::get_int(pre_params, adm_params, "Number of seed locations", 1, &P);
	if (P.NumSeedLocations > MAX_NUM_SEED_LOCATIONS)
	{
		Files::xfprintf_stderr("Too many seed locations\n");
		P.NumSeedLocations = MAX_NUM_SEED_LOCATIONS;
	}
	Files::xfprintf_stderr("Reading Params 3.9\n"); fflush(stderr);

	Params::req_int_vec(pre_params, adm_params, "Initial number of infecteds", P.NumInitialInfections, P.NumSeedLocations, &P);
	Files::xfprintf_stderr("Reading Params 3.10\n"); fflush(stderr);

	Params::get_double_matrix(pre_params, adm_params, "Location of initial infecteds", P.LocationInitialInfection, P.NumSeedLocations, 2, 0, &P);
	Files::xfprintf_stderr("Reading Params 3.11\n"); fflush(stderr);
	P.MinPopDensForInitialInfection = Params::get_int(pre_params, adm_params, "Minimum population in microcell of initial infection", 0, &P);
	P.MaxPopDensForInitialInfection = Params::get_int(pre_params, adm_params, "Maximum population in microcell of initial infection", 10000000, &P);
	P.MaxAgeForInitialInfection = Params::get_int(pre_params, adm_params, "Maximum age of initial infections", 1000, &P);
	P.DoRandomInitialInfectionLoc = Params::get_int(pre_params, adm_params, "Randomise initial infection location", 1, &P);
	P.DoAllInitialInfectioninSameLoc = Params::get_int(pre_params, adm_params, "All initial infections located in same microcell", 0, &P);

	Files::xfprintf_stderr("Reading Params 3 blah\n"); fflush(stderr);


	if (Params::param_found(pre_params, adm_params, "Day of year of start of seeding"))
	{
		P.InitialInfectionCalTime = Params::req_double(pre_params, adm_params, "Day of year of start of seeding", &P);
		if (Params::param_found(params, pre_params, "Scaling of infection seeding"))
		{
			P.SeedingScaling = Params::req_double(params, pre_params, "Scaling of infection seeding", &P);
		  P.DoNoCalibration = 1;
		}
		else
		{
			P.SeedingScaling = 1.0;
			P.DoNoCalibration = 0;
		}
	}
	else
	{
		P.SeedingScaling = 1.0;
		P.DoNoCalibration = 0;
		P.InitialInfectionCalTime = -1;
	}

	Files::xfprintf_stderr("Reading Params 4\n"); fflush(stderr);


	if (P.FitIter == 0)
	{
		if (P.DoAdUnits != 0)
		{
			if (!Params::param_found(pre_params, adm_params, "Administrative unit to seed initial infection into"))
			{
				Params::req_string_vec(pre_params, adm_params, "Administrative unit to seed initial infection into", AdunitListNames, P.NumSeedLocations, &P);
				for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
			}
			else
				for (i = 0; i < P.NumSeedLocations; i++)
				{
					f = 0;
					if (P.NumAdunits > 0)
					{
						for (j = 0; (j < P.NumAdunits) && (!f); j++) f = (!strcmp(AdUnits[j].ad_name, AdunitListNames[i]));
						if (f) k = AdUnits[j - 1].id;
					}
					if (!f) k = atoi(AdunitListNames[i]);
					P.InitialInfectionsAdminUnit[i] = k;
					P.InitialInfectionsAdminUnitId[i] = P.AdunitLevel1Lookup[(k % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor];
				}
			Params::get_double_vec(pre_params, adm_params, "Administrative unit seeding weights", P.InitialInfectionsAdminUnitWeight, P.NumSeedLocations, 1.0, P.NumSeedLocations, &P);
			s = 0;
			for (i = 0; i < P.NumSeedLocations; i++) s += P.InitialInfectionsAdminUnitWeight[i];
			for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnitWeight[i] /= s;
		}
		else
		{
			for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
		}
	}
	P.InfectionImportRate1 = Params::get_double(params, pre_params, "Initial rate of importation of infections", 0, &P);
	P.InfectionImportRate2 = Params::get_double(params, pre_params, "Changed rate of importation of infections", 0, &P);
	P.InfectionImportChangeTime = Params::get_double(params, pre_params, "Time when infection rate changes", 1e10, &P);
	P.DoImportsViaAirports = Params::get_int(params, pre_params, "Imports via air travel", 0, &P);
	P.DurImportTimeProfile = Params::get_int(params, pre_params, "Length of importation time profile provided", 0, &P);
	if (P.DurImportTimeProfile > 0)
	{
		if (P.DurImportTimeProfile >= MAX_DUR_IMPORT_PROFILE) ERR_CRITICAL("MAX_DUR_IMPORT_PROFILE too small\n");
		Params::req_double_vec(params, pre_params, "Daily importation time profile", P.ImportInfectionTimeProfile, P.DurImportTimeProfile, &P);
	}
	P.R0 = Params::req_double(params, pre_params, "Reproduction number", &P);
	if (Params::param_found(params, pre_params, "Beta for spatial transmission"))
	{
	  P.LocalBeta = Params::req_double(params, pre_params, "Beta for spatial transmission", &P);
	  P.FixLocalBeta = 1;
	}
	else
	{
		P.LocalBeta = -1.0;
		P.FixLocalBeta = 0;
	}
	P.InfectiousPeriod = Params::req_double(params, pre_params, "Infectious period", &P);
	P.SusceptibilitySD = Params::get_double(params, pre_params, "SD of individual variation in susceptibility", 0, &P);
	P.InfectiousnessSD = Params::get_double(params, pre_params, "SD of individual variation in infectiousness", 0, &P);
	if (Params::param_found(params, pre_params, "k of individual variation in infectiousness"))
    P.InfectiousnessSD = 1.0 / sqrt(Params::req_double(params, pre_params, "k of individual variation in infectiousness", &P));
	P.NoInfectiousnessSDinHH = Params::get_int(params, pre_params, "k does not apply in households", 0, &P);
	P.DoInfectiousnessProfile = Params::get_int(params, pre_params, "Model time varying infectiousness", 0, &P);
	P.R0DensityScalePower = Params::get_double(params, pre_params, "Power of scaling of spatial R0 with density", 0, &P);
	if (P.DoInfectiousnessProfile != 0)
	{
		Params::get_double_vec(params, pre_params, "Infectiousness profile", P.infectious_prof, INFPROF_RES, 1, INFPROF_RES, &P);
		k = (int)ceil(P.InfectiousPeriod / P.ModelTimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		s = 0;
		P.infectious_prof[INFPROF_RES] = 0;
		for (i = 0; i < MAX_INFECTIOUS_STEPS; i++)	P.infectiousness[i] = 0;
		for (i = 0; i < k; i++)
		{
			t = (((double)i) * P.ModelTimeStep / P.InfectiousPeriod * INFPROF_RES);
			j = (int)t;
			t -= (double)j;
			if (j < INFPROF_RES)
				s += (P.infectiousness[i] = P.infectious_prof[j] * (1 - t) + P.infectious_prof[j + 1] * t);
			else
				s += (P.infectiousness[i] = P.infectious_prof[INFPROF_RES]);
		}
		s /= ((double)k);
		for (i = 0; i <= k; i++) P.infectiousness[i] /= s;
		P.infectious_icdf.assign_exponent(-1.0);
	}
	else
	{
		if (Params::param_found(params, pre_params, "Infectious period inverse CDF"))
		{
			Params::req_double_vec(params, pre_params, "Infectious period inverse CDF", P.infectious_icdf.get_values(), CDF_RES + 1, &P);
			P.infectious_icdf.set_neg_log(ICDF_START);
		}
		k = (int)ceil(P.InfectiousPeriod * P.infectious_icdf[CDF_RES] / P.ModelTimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		for (i = 0; i < k; i++) P.infectiousness[i] = 1.0;
		P.infectiousness[k] = 0;
		P.infectious_icdf.assign_exponent();
	}

	Files::xfprintf_stderr("Reading Params 5\n"); fflush(stderr);


	P.DoLatent = Params::get_int(params, pre_params, "Include latent period", 0, &P);
	if (P.DoLatent != 0)
	{
		P.LatentPeriod = Params::req_double(params, pre_params, "Latent period", &P);
		GetInverseCdf(params, pre_params, "Latent period inverse CDF", &P.latent_icdf, &P, 1e10);
	}

	P.DoSymptoms = Params::get_int(params, pre_params, "Include symptoms", 0, &P);
	if (P.DoSymptoms == 0)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.ProportionSymptomatic[i] = 0;
		P.FalsePositiveRate = 0;
		P.SymptInfectiousness = P.AsymptInfectiousness = 1.0;
		P.LatentToSymptDelay = 0;
	}
	else
	{
		if (P.DoAge != 0)
			Params::req_double_vec(params, pre_params, "Proportion symptomatic by age group", P.ProportionSymptomatic, NUM_AGE_GROUPS, &P);
		else
		{
			P.ProportionSymptomatic[0] = Params::req_double(params, pre_params, "Proportion symptomatic", &P);
			for (i = 1; i < NUM_AGE_GROUPS; i++)
				P.ProportionSymptomatic[i] = P.ProportionSymptomatic[0];
		}
		P.LatentToSymptDelay = Params::req_double(params, pre_params, "Delay from end of latent period to start of symptoms", &P);
		P.SymptSpatialContactRate = Params::req_double(params, pre_params, "Relative rate of random contacts if symptomatic", &P);
		P.SymptInfectiousness = Params::get_double(params, pre_params, "Symptomatic infectiousness relative to asymptomatic", 1.0, &P);
		P.AsymptInfectiousness = Params::get_double(params, pre_params, "Asymptomatic infectiousness relative to symptomatic", 1.0, &P);
		P.DoRealSymptWithdrawal = Params::get_int(params, pre_params, "Model symptomatic withdrawal to home as true absenteeism", 0, &P);
		if (P.DoPlaces != 0)
		{
			Params::req_double_vec(params, pre_params, "Relative level of place attendance if symptomatic", P.SymptPlaceTypeContactRate, P.PlaceTypeNum, &P);
			if (P.DoRealSymptWithdrawal != 0)
			{
				for (j = 0; j < NUM_PLACE_TYPES; j++)
				{
					P.SymptPlaceTypeWithdrawalProp[j] = 1.0 - P.SymptPlaceTypeContactRate[j];
					P.SymptPlaceTypeContactRate[j] = 1.0;
				}
			}
			else
				for (j = 0; j < NUM_PLACE_TYPES; j++) P.SymptPlaceTypeWithdrawalProp[j] = 0.0;
		}
		P.CaseAbsentChildAgeCutoff = Params::get_int(params, pre_params, "Maximum age of child at home for whom one adult also stays at home", 0, &P);
		P.CaseAbsentChildPropAdultCarers = Params::get_double(params, pre_params, "Proportion of children at home for whom one adult also stays at home", 0, &P);
		P.PlaceCloseRoundHousehold = Params::get_int(params, pre_params, "Place close round household", 1, &P);
		P.AbsenteeismPlaceClosure = Params::get_int(params, pre_params, "Absenteeism place closure", 0, &P);
		if (P.AbsenteeismPlaceClosure != 0)
		{
			P.CaseAbsenteeismDelay = 0;  // Set to zero for tracking absenteeism
			P.MaxAbsentTime = Params::get_int(params, pre_params, "Max absent time", MAX_ABSENT_TIME, &P);
			if (P.MaxAbsentTime > MAX_ABSENT_TIME || P.MaxAbsentTime < 0)
			{
				ERR_CRITICAL_FMT("[Max absent time] out of range (%d), should be in range [0, %d]", P.MaxAbsentTime, MAX_ABSENT_TIME);
			}
		}
		else
		{
			P.CaseAbsenteeismDelay = Params::get_double(params, pre_params, "Delay in starting place absenteeism for cases who withdraw", 0, &P);
			P.MaxAbsentTime = 0; // Not used when !P.AbsenteeismPlaceClosure
		}
		P.CaseAbsenteeismDuration = Params::get_double(params, pre_params, "Duration of place absenteeism for cases who withdraw", 7, &P);

		P.FalsePositiveRate = Params::get_double(params, pre_params, "False positive rate", 0.0, &P);
		P.FalsePositivePerCapitaIncidence = Params::get_double(params, pre_params, "False positive per capita incidence", 0.0, &P);
		Params::get_double_vec(params, pre_params, "False positive relative incidence by age", P.FalsePositiveAgeRate, NUM_AGE_GROUPS, 1.0, NUM_AGE_GROUPS, &P);
	}

	Files::xfprintf_stderr("Reading Params 6\n"); fflush(stderr);


	P.SeroConvMaxSens = Params::get_double(params, pre_params, "Maximum sensitivity of serology assay", 1.0, &P);
	P.SeroConvP1 = Params::get_double(params, pre_params, "Seroconversion model parameter 1", 14.0, &P);
	P.SeroConvP2 = Params::get_double(params, pre_params, "Seroconversion model parameter 2", 3.0, &P);
	P.SeroConvSpec = Params::get_double(params, pre_params, "Specificity of serology assay", 1.0, &P);
	P.InfPrevSurveyScale = Params::get_double(params, pre_params, "Scaling of modelled infection prevalence to match surveys", 1.0, &P);

	P.DoSeverity = Params::get_int(params, pre_params, "Do Severity Analysis", 0, &P);
	if (P.DoSeverity != 0)
	{
		P.ScaleSymptProportions = Params::get_double(params, pre_params, "Factor to scale IFR", 1.0, &P);
		//// Means for icdf's.
		P.Mean_TimeToTest = Params::get_double(params, pre_params, "MeanTimeToTest", 0.0, &P);
		P.Mean_TimeToTestOffset = Params::get_double(params, pre_params, "MeanTimeToTestOffset", 1.0, &P);
		P.Mean_TimeToTestCriticalOffset = Params::get_double(params, pre_params, "MeanTimeToTestCriticalOffset", 1.0, &P);
		P.Mean_TimeToTestCritRecovOffset = Params::get_double(params, pre_params, "MeanTimeToTestCritRecovOffset", 1.0, &P);
		if (Params::get_int(params, pre_params, "Age dependent severity delays", 0, &P) == 0)
		{
			P.Mean_MildToRecovery[0] = Params::req_double(params, pre_params, "Mean_MildToRecovery", &P);
			P.Mean_ILIToRecovery[0] = Params::req_double(params, pre_params, "Mean_ILIToRecovery", &P);
			P.Mean_SARIToRecovery[0] = Params::req_double(params, pre_params, "Mean_SARIToRecovery", &P);
			P.Mean_CriticalToCritRecov[0] = Params::req_double(params, pre_params, "Mean_CriticalToCritRecov", &P);
			P.Mean_CritRecovToRecov[0] = Params::req_double(params, pre_params, "Mean_CritRecovToRecov", &P);
			P.Mean_ILIToSARI[0] = Params::req_double(params, pre_params, "Mean_ILIToSARI", &P);
			P.Mean_ILIToDeath[0] = Params::get_double(params, pre_params, "Mean_ILIToDeath", 7.0, &P);
			P.Mean_SARIToCritical[0] = Params::req_double(params, pre_params, "Mean_SARIToCritical", &P);
			P.Mean_SARIToDeath[0] = Params::req_double(params, pre_params, "Mean_SARIToDeath", &P);
			P.Mean_CriticalToDeath[0] = Params::req_double(params, pre_params, "Mean_CriticalToDeath", &P);
			for (int AgeGroup = 1; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
			{
				P.Mean_MildToRecovery		[AgeGroup] = P.Mean_MildToRecovery		[0];
				P.Mean_ILIToRecovery		[AgeGroup] = P.Mean_ILIToRecovery		[0];
				P.Mean_SARIToRecovery		[AgeGroup] = P.Mean_SARIToRecovery		[0];
				P.Mean_CriticalToCritRecov	[AgeGroup] = P.Mean_CriticalToCritRecov	[0];
				P.Mean_CritRecovToRecov		[AgeGroup] = P.Mean_CritRecovToRecov	[0];
				P.Mean_ILIToSARI			[AgeGroup] = P.Mean_ILIToSARI			[0];
				P.Mean_ILIToDeath			[AgeGroup] = P.Mean_ILIToDeath			[0];
				P.Mean_SARIToCritical		[AgeGroup] = P.Mean_SARIToCritical		[0];
				P.Mean_SARIToDeath			[AgeGroup] = P.Mean_SARIToDeath			[0];
				P.Mean_CriticalToDeath		[AgeGroup] = P.Mean_CriticalToDeath		[0];
			}
		}
		else
		{
			Params::req_double_vec(params, pre_params, "Mean_MildToRecovery", P.Mean_MildToRecovery, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_ILIToRecovery", P.Mean_ILIToRecovery, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_SARIToRecovery", P.Mean_SARIToRecovery, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_CriticalToCritRecov", P.Mean_CriticalToCritRecov, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_CritRecovToRecov", P.Mean_CritRecovToRecov, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_ILIToSARI", P.Mean_ILIToSARI, NUM_AGE_GROUPS, &P);
			Params::get_double_vec(params, pre_params, "Mean_ILIToDeath", P.Mean_ILIToDeath, NUM_AGE_GROUPS, 7.0, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_SARIToCritical", P.Mean_SARIToCritical, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_SARIToDeath", P.Mean_SARIToDeath, NUM_AGE_GROUPS, &P);
			Params::req_double_vec(params, pre_params, "Mean_CriticalToDeath", P.Mean_CriticalToDeath, NUM_AGE_GROUPS, &P);
		}

		//// Get InverseCDFs
		GetInverseCdf(params, pre_params, "MildToRecovery_icdf", &P.MildToRecovery_icdf, &P);
		GetInverseCdf(params, pre_params, "ILIToRecovery_icdf", &P.ILIToRecovery_icdf, &P);
		GetInverseCdf(params, pre_params, "ILIToDeath_icdf", &P.ILIToDeath_icdf, &P);
		GetInverseCdf(params, pre_params, "SARIToRecovery_icdf", &P.SARIToRecovery_icdf, &P);
		GetInverseCdf(params, pre_params, "CriticalToCritRecov_icdf", &P.CriticalToCritRecov_icdf, &P);
		GetInverseCdf(params, pre_params, "CritRecovToRecov_icdf", &P.CritRecovToRecov_icdf, &P);
		GetInverseCdf(params, pre_params, "ILIToSARI_icdf", &P.ILIToSARI_icdf, &P);
		GetInverseCdf(params, pre_params, "SARIToCritical_icdf", &P.SARIToCritical_icdf, &P);
		GetInverseCdf(params, pre_params, "SARIToDeath_icdf", &P.SARIToDeath_icdf, &P);
		GetInverseCdf(params, pre_params, "CriticalToDeath_icdf", &P.CriticalToDeath_icdf, &P);

		// If you decide to decompose Critical -> Death transition into Critical -> Stepdown and Stepdown -> Death, use the block below.
		P.IncludeStepDownToDeath = Params::get_int(params, pre_params, "IncludeStepDownToDeath", 0, &P);
		if (P.IncludeStepDownToDeath == 0) /// for backwards compatibility. If Stepdown to death not included (or if unspecified), set stepdown->death = stepdown->recovery.
		{
			for (int quantile = 0; quantile <= CDF_RES; quantile++)
				P.StepdownToDeath_icdf[quantile] = P.CritRecovToRecov_icdf[quantile];
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				P.Mean_StepdownToDeath[AgeGroup] = P.Mean_CritRecovToRecov[AgeGroup];
		}
		else
		{
			Params::req_double_vec(params, pre_params, "Mean_StepdownToDeath", P.Mean_StepdownToDeath, NUM_AGE_GROUPS, &P);
			GetInverseCdf(params, pre_params, "StepdownToDeath_icdf", &P.StepdownToDeath_icdf, &P);
		}

		Params::get_double_vec(params, pre_params, "Prop_Mild_ByAge", P.Prop_Mild_ByAge, NUM_AGE_GROUPS, 0.5, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "Prop_ILI_ByAge", P.Prop_ILI_ByAge, NUM_AGE_GROUPS, 0.3, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "Prop_SARI_ByAge", P.Prop_SARI_ByAge, NUM_AGE_GROUPS, 0.15, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "Prop_Critical_ByAge", P.Prop_Critical_ByAge, NUM_AGE_GROUPS, 0.05, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "CFR_SARI_ByAge", P.CFR_SARI_ByAge, NUM_AGE_GROUPS, 0.5, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "CFR_Critical_ByAge", P.CFR_Critical_ByAge, NUM_AGE_GROUPS, 0.5, NUM_AGE_GROUPS, &P);
		Params::get_double_vec(params, pre_params, "CFR_ILI_ByAge", P.CFR_ILI_ByAge, NUM_AGE_GROUPS, 0, NUM_AGE_GROUPS, &P);
			
		//Add param to allow severity to be uniformly scaled up or down.
		for (i = 0; i < NUM_AGE_GROUPS; i++)
		{
			P.Prop_SARI_ByAge[i] *= P.ScaleSymptProportions;
			P.Prop_Critical_ByAge[i] *= P.ScaleSymptProportions;
			P.Prop_ILI_ByAge[i] = 1.0 - P.Prop_Mild_ByAge[i] - P.Prop_SARI_ByAge[i] - P.Prop_Critical_ByAge[i];
		}
	}

	Files::xfprintf_stderr("Reading Params 7\n"); fflush(stderr);


	if (P.FitIter == 0)
	{
		if (Params::param_found(params, pre_params, "Bounding box for bitmap"))
		{
		  double* buf = new double[4];
			Params::req_double_vec(params, pre_params, "Bounding box for bitmap", buf, 4, &P);
			P.BoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(buf[0], buf[1]);
			P.BoundingBox.top_right() = CovidSim::Geometry::Vector2d(buf[2], buf[3]);
			delete[] buf;
		}
		else 
		{
			P.BoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(0.0, 0.0);
			P.BoundingBox.top_right()   = CovidSim::Geometry::Vector2d(1.0, 1.0);
		}
		if (Params::param_found(params, pre_params, "Spatial domain for simulation"))
		{
			double* buf = new double[4];
			Params::req_double_vec(params, pre_params, "Spatial domain for simulation", buf, 4, &P);
			P.SpatialBoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(buf[0], buf[1]);
			P.SpatialBoundingBox.top_right() = CovidSim::Geometry::Vector2d(buf[2], buf[3]);
			delete[] buf;
		}
		else
		{
			P.SpatialBoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(0.0, 0.0);
			P.SpatialBoundingBox.top_right()   = CovidSim::Geometry::Vector2d(1.0, 1.0);
		}
		P.in_cells_.width = Params::get_double(params, pre_params, "Grid size", 1.0 / 120.0, &P);
		P.DoUTM_coords = Params::get_int(params, pre_params, "Use long/lat coord system", 1, &P);
		P.BitmapScale = Params::get_double(params, pre_params, "Bitmap scale", 1, &P);
		P.BitmapAspectScale = Params::get_double(params, pre_params, "Bitmap y:x aspect scaling", 1, &P);
		P.BitmapMovieFrame = Params::get_int(params, pre_params, "Bitmap movie frame interval", 250, &P);
		P.OutputBitmap = Params::get_int(params, pre_params, "Output bitmap", 0, &P);
		P.OutputBitmapDetected = Params::get_int(params, pre_params, "Output bitmap detected", 0, &P);
		P.DoImmuneBitmap = Params::get_int(params, pre_params, "Output immunity on bitmap", 0, &P);
		P.DoInfectionTree = Params::get_int(params, pre_params, "Output infection tree", 0, &P);
		P.DoOneGen = Params::get_int(params, pre_params, "Do one generation", 0, &P);
		P.OutputEveryRealisation = Params::get_int(params, pre_params, "Output every realisation", 0, &P);
		P.MaxCorrSample = Params::get_int(params, pre_params, "Maximum number to sample for correlations", 1000000000, &P);
		P.DoSI = Params::get_int(params, pre_params, "Assume SI model", 0, &P);
		P.DoPeriodicBoundaries = Params::get_int(params, pre_params, "Assume periodic boundary conditions", 0, &P);
		P.OutputOnlyNonExtinct = Params::get_int(params, pre_params, "Only output non-extinct realisations", 0, &P);

		P.DoPerCapitaTriggers = Params::get_int(params, pre_params, "Use cases per thousand threshold for area controls", 0, &P);
		P.DoGlobalTriggers = Params::get_int(params, pre_params, "Use global triggers for interventions", 0, &P);
		P.DoAdminTriggers = Params::get_int(params, pre_params, "Use admin unit triggers for interventions", 0, &P);
		P.DoICUTriggers = Params::get_int(params, pre_params, "Use ICU case triggers for interventions", 0, &P);
		if (P.DoGlobalTriggers != 0)  P.DoAdminTriggers = 0;
		P.IncThreshPop = Params::get_int(params, pre_params, "Divisor for per-capita area threshold (default 1000)", 1000, &P);
		P.GlobalIncThreshPop = Params::get_int(params, pre_params, "Divisor for per-capita global threshold (default 1000)", 1000, &P);

		P.TriggersSamplingInterval = Params::get_int(params, pre_params, "Number of sampling intervals over which cumulative incidence measured for global trigger", 10000000, &P);
		P.PostAlertControlPropCasesId = Params::get_double(params, pre_params, "Proportion of cases detected for treatment", 1, &P);
		P.PreAlertControlPropCasesId = Params::get_double(params, pre_params, "Proportion of cases detected before outbreak alert", 1, &P);
		P.TriggerAlertOnDeaths = Params::get_int(params, pre_params, "Trigger alert on deaths", 0, &P);
	}

	Files::xfprintf_stderr("Reading Params 8\n"); fflush(stderr);


	if (P.TriggerAlertOnDeaths != 0)
	{
		P.CaseOrDeathThresholdBeforeAlert = Params::get_int(params, pre_params, "Number of deaths accummulated before alert", 0, &P);
	}
	else
	{
		P.CaseOrDeathThresholdBeforeAlert = Params::get_int(params, pre_params, "Number of detected cases needed before outbreak alert triggered", 0, &P);
	}

	if (P.CaseOrDeathThresholdBeforeAlert_CommandLine > 0) P.CaseOrDeathThresholdBeforeAlert = P.CaseOrDeathThresholdBeforeAlert_CommandLine;
	P.DoAlertTriggerAfterInterv = Params::get_int(params, pre_params, "Alert trigger starts after interventions", 0, &P);
	P.DateTriggerReached_CalTime = Params::get_double(params, pre_params, "Day of year trigger is reached", -1, &P);
	if (P.DoAlertTriggerAfterInterv != 0)
	{
		P.Interventions_StartDate_CalTime = Params::req_double(params, pre_params, "Day of year interventions start", &P);
		if (P.DateTriggerReached_CalTime <= P.Interventions_StartDate_CalTime)
			P.DoAlertTriggerAfterInterv = 0;
		else
		{
			P.AlertTriggerAfterIntervThreshold = P.CaseOrDeathThresholdBeforeAlert;
			P.CaseOrDeathThresholdBeforeAlert = 1000;
			Files::xfprintf_stderr("Threshold of %i deaths by day %lg\n", P.AlertTriggerAfterIntervThreshold, P.DateTriggerReached_CalTime);
		}
	}
	else
	{
		P.Interventions_StartDate_CalTime = P.DateTriggerReached_CalTime;
	}
	if (P.FitIter == 0)
	{
		P.CaseOrDeathThresholdBeforeAlert_Fixed = P.CaseOrDeathThresholdBeforeAlert;
	}

	P.WindowToEvaluateTriggerAlert = Params::get_int(params, pre_params, "Number of days to accummulate cases/deaths before alert", 1000, &P);

	P.DoPlaceGroupTreat = Params::get_int(params, pre_params, "Only treat mixing groups within places", 0, &P);

	P.TreatCellIncThresh = Params::get_double(params, pre_params, "Treatment trigger incidence per cell", INT32_MAX, &P);
	P.CaseIsolation_CellIncThresh = Params::get_double(params, pre_params, "Case isolation trigger incidence per cell", P.TreatCellIncThresh, &P);
	P.HHQuar_CellIncThresh = Params::get_double(params, pre_params, "Household quarantine trigger incidence per cell", P.TreatCellIncThresh, &P);

	P.TreatSuscDrop = Params::get_double(params, pre_params, "Relative susceptibility of treated individual", 1, &P);
	P.TreatInfDrop = Params::get_double(params, pre_params, "Relative infectiousness of treated individual", 1, &P);
	P.TreatDeathDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases resulting in death prevented by treatment", 0, &P);
	P.TreatSympDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases prevented by treatment", 0, &P);
	P.TreatDelayMean = Params::get_double(params, pre_params, "Delay to treat cell", 0, &P);
	P.TreatCaseCourseLength = Params::get_double(params, pre_params, "Duration of course of treatment", 5, &P);
	P.TreatProphCourseLength = Params::get_double(params, pre_params, "Duration of course of prophylaxis", 10, &P);
	P.TreatPropCases = Params::get_double(params, pre_params, "Proportion of detected cases treated", 1, &P);
	if (P.DoHouseholds != 0)
	{
		P.TreatPropCaseHouseholds = Params::get_double(params, pre_params, "Proportion of households of cases treated", 0, &P);
		P.TreatHouseholdsDuration = Params::get_double(params, pre_params, "Duration of household prophylaxis policy", USHRT_MAX / P.TimeStepsPerDay, &P);
	}
	// Check below - "Proportional treated" will always be ignored.
	//if (!GetInputParameter2(params, pre_params, "Proportion treated", "%lf", (void*) & (P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	//if (!GetInputParameter2(params, pre_params, "Proportion treated in radial prophylaxis", "%lf", (void*) & (P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;

	P.TreatPropRadial = Params::get_double(params, pre_params, "Proportion treated in radial prophylaxis", 1.0, &P);
	P.TreatRadius = Params::get_double(params, pre_params, "Treatment radius", 0, &P);
	P.TreatPlaceGeogDuration = Params::get_double(params, pre_params, "Duration of place/geographic prophylaxis policy", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.TreatTimeStartBase = Params::get_double(params, pre_params, "Treatment start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	if (P.DoPlaces != 0)
	{
		Params::get_double_vec(params, pre_params, "Proportion of places treated after case detected", P.TreatPlaceProbCaseId, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
		Params::get_double_vec(params, pre_params, "Proportion of people treated in targeted places", P.TreatPlaceTotalProp, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
	}
	P.TreatMaxCoursesBase = Params::get_double(params, pre_params, "Maximum number of doses available", 1e20, &P);
	P.TreatNewCoursesStartTime = Params::get_double(params, pre_params, "Start time of additional treatment production", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.TreatNewCoursesRate = Params::get_double(params, pre_params, "Rate of additional treatment production (courses per day)", 0, &P);
	P.TreatMaxCoursesPerCase = Params::get_int(params, pre_params, "Maximum number of people targeted with radial prophylaxis per case", INT32_MAX, &P);

	if (P.DoAdUnits != 0)
	{
		P.TreatByAdminUnit = Params::get_int(params, pre_params, "Treat administrative units rather than rings", 0, &P);
		P.TreatAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for treatment", 1, &P);
		if ((P.TreatAdminUnitDivisor == 0) || (P.TreatByAdminUnit == 0)) { P.TreatByAdminUnit = 0; P.TreatAdminUnitDivisor = 1; }
	}
	else
	{
		P.TreatAdminUnitDivisor = 1; P.TreatByAdminUnit = 0;
	}

	P.VaccCellIncThresh = Params::get_double(params, pre_params, "Vaccination trigger incidence per cell", 1000000000, &P);
	P.VaccSuscDrop = Params::get_double(params, pre_params, "Relative susceptibility of vaccinated individual", 1, &P);
	P.VaccSuscDrop2 = Params::get_double(params, pre_params, "Relative susceptibility of individual vaccinated after switch time", 1, &P);
	P.VaccTimeEfficacySwitch = Params::get_double(params, pre_params, "Switch time at which vaccine efficacy increases", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.VaccEfficacyDecay = Params::get_double(params, pre_params, "Decay rate of vaccine efficacy (per year)", 0, &P);
	P.VaccEfficacyDecay /= DAYS_PER_YEAR;
	P.VaccInfDrop = Params::get_double(params, pre_params, "Relative infectiousness of vaccinated individual", 1, &P);
	P.VaccMortDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases resulting in death prevented by vaccination", 0, &P);
	P.VaccSympDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases prevented by vaccination", 0, &P);
	P.VaccDelayMean = Params::get_double(params, pre_params, "Delay to vaccinate", 0, &P);

	P.VaccTimeToEfficacy = Params::get_double(params, pre_params, "Delay from vaccination to full protection", 0, &P);

	P.VaccCampaignInterval = Params::get_double(params, pre_params, "Years between rounds of vaccination", 1e10, &P);
	P.VaccDosePerDay = Params::get_int(params, pre_params, "Max vaccine doses per day", -1, &P);
	P.VaccCampaignInterval *= DAYS_PER_YEAR;
	P.VaccMaxRounds = Params::get_int(params, pre_params, "Maximum number of rounds of vaccination", 1, &P);
	if (P.DoHouseholds != 0)
	{
		P.VaccPropCaseHouseholds = Params::get_double(params, pre_params, "Proportion of households of cases vaccinated", 0, &P);
		P.VaccHouseholdsDuration = Params::get_double(params, pre_params, "Duration of household vaccination policy", USHRT_MAX / P.TimeStepsPerDay, &P);
	}

	P.VaccTimeStartBase = Params::get_double(params, pre_params, "Vaccination start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.VaccProp = Params::get_double(params, pre_params, "Proportion of population vaccinated", 0, &P);
	P.VaccCoverageIncreasePeriod = Params::get_double(params, pre_params, "Time taken to reach max vaccination coverage (in years)", 0, &P);
	P.VaccCoverageIncreasePeriod *= DAYS_PER_YEAR;
	P.VaccTimeStartGeo = Params::get_double(params, pre_params, "Time to start geographic vaccination", 1e10, &P);
	P.VaccRadius = Params::get_double(params, pre_params, "Vaccination radius", 0, &P);
	P.VaccMinRadius = Params::get_double(params, pre_params, "Minimum radius from case to vaccinate", 0, &P);
	P.VaccMaxCoursesBase = Params::get_double(params, pre_params, "Maximum number of vaccine courses available", 1e20, &P);
	P.VaccNewCoursesStartTime = Params::get_double(params, pre_params, "Start time of additional vaccine production", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.VaccNewCoursesEndTime = Params::get_double(params, pre_params, "End time of additional vaccine production", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.VaccNewCoursesRate = Params::get_double(params, pre_params, "Rate of additional vaccine production (courses per day)", 0, &P);
	P.DoMassVacc = Params::get_int(params, pre_params, "Apply mass rather than reactive vaccination", 0, &P);
	if (Params::param_found(params, pre_params, "Priority age range for mass vaccination")) {
	  Params::req_int_vec(params, pre_params, "Priority age range for mass vaccination", P.VaccPriorityGroupAge, 2, &P);
	}	else {
		P.VaccPriorityGroupAge[0] = 1; P.VaccPriorityGroupAge[1] = 0;
	}
	
	if (P.DoAdUnits != 0)
	{
		P.VaccByAdminUnit = Params::get_int(params, pre_params, "Vaccinate administrative units rather than rings", 0, &P);
		P.VaccAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for vaccination", 1, &P);
		if ((P.VaccAdminUnitDivisor == 0) || (P.VaccByAdminUnit == 0)) P.VaccAdminUnitDivisor = 1;
	}
	else
	{
		P.VaccAdminUnitDivisor = 1; P.VaccByAdminUnit = 0;
	}

	P.MoveRestrCellIncThresh = Params::get_int(params, pre_params, "Movement restrictions trigger incidence per cell", INT32_MAX, &P);
	P.MoveDelayMean = Params::get_double(params, pre_params, "Delay to start movement restrictions", 0, &P);
	P.MoveRestrDuration = Params::get_double(params, pre_params, "Duration of movement restrictions", 7, &P);
	P.MoveRestrEffect = Params::get_double(params, pre_params, "Residual movements after restrictions", 0, &P);
	P.MoveRestrRadius = Params::get_double(params, pre_params, "Minimum radius of movement restrictions", 0, &P);
	P.MoveRestrTimeStartBase = Params::get_double(params, pre_params, "Movement restrictions start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.DoBlanketMoveRestr = Params::get_int(params, pre_params, "Impose blanket movement restrictions", 0, &P);
	P.DoMoveRestrOnceOnly = Params::get_int(params, pre_params, "Movement restrictions only once", 0, &P);

	Files::xfprintf_stderr("Reading Params 9\n"); fflush(stderr);


	//if (P.DoMoveRestrOnceOnly) P.DoMoveRestrOnceOnly = 4; //// don't need this anymore with TreatStat option. Keep it as a boolean.
	if (P.DoAdUnits != 0)
	{
		P.MoveRestrByAdminUnit = Params::get_int(params, pre_params, "Movement restrictions in administrative units rather than rings", 0, &P);
		P.MoveRestrAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for movement restrictions", 1, &P);
		if ((P.MoveRestrAdminUnitDivisor == 0) || (P.MoveRestrByAdminUnit == 0)) P.MoveRestrAdminUnitDivisor = 1;
	}
	else
	{
		P.MoveRestrAdminUnitDivisor = 1; P.MoveRestrByAdminUnit = 0;
	}

	//Intervention delays and durations by admin unit: ggilani 16/03/20
	P.DoInterventionDelaysByAdUnit = Params::get_int(params, pre_params, "Include intervention delays by admin unit", 0, &P);
	if (P.DoInterventionDelaysByAdUnit)
	{
		//Set up arrays to temporarily store parameters per admin unit
		double AdunitDelayToSocialDistance	[MAX_ADUNITS];
		double AdunitDelayToHQuarantine		[MAX_ADUNITS];
		double AdunitDelayToCaseIsolation	[MAX_ADUNITS];
		double AdunitDelayToPlaceClose		[MAX_ADUNITS];
		double AdunitDurationSocialDistance	[MAX_ADUNITS];
		double AdunitDurationHQuarantine	[MAX_ADUNITS];
		double AdunitDurationCaseIsolation	[MAX_ADUNITS];
		double AdunitDurationPlaceClose		[MAX_ADUNITS];

		Params::get_double_vec(params, pre_params, "Delay to social distancing by admin unit", AdunitDelayToSocialDistance, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Delay to household quarantine by admin unit", AdunitDelayToHQuarantine, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Delay to case isolation by admin unit", AdunitDelayToCaseIsolation, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Delay to place closure by admin unit", AdunitDelayToPlaceClose, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Duration of social distancing by admin unit", AdunitDurationSocialDistance, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Duration of household quarantine by admin unit", AdunitDurationHQuarantine, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Duration of case isolation by admin unit", AdunitDurationCaseIsolation, P.NumAdunits, 0, P.NumAdunits, &P);
		Params::get_double_vec(params, pre_params, "Duration of place closure by admin unit", AdunitDurationPlaceClose, P.NumAdunits, 0, P.NumAdunits, &P);

		for (i = 0; i < P.NumAdunits; i++)
		{
			AdUnits[i].SocialDistanceDelay			= AdunitDelayToSocialDistance	[i];
			AdUnits[i].SocialDistanceDuration		= AdunitDurationSocialDistance	[i];
			AdUnits[i].HQuarantineDelay				= AdunitDelayToHQuarantine		[i];
			AdUnits[i].HQuarantineDuration			= AdunitDurationHQuarantine		[i];
			AdUnits[i].CaseIsolationDelay			= AdunitDelayToCaseIsolation	[i];
			AdUnits[i].CaseIsolationPolicyDuration	= AdunitDurationCaseIsolation	[i];
			AdUnits[i].PlaceCloseDelay				= AdunitDelayToPlaceClose		[i];
			AdUnits[i].PlaceCloseDuration			= AdunitDurationPlaceClose		[i];
		}
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** DIGITAL CONTACT TRACING
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	//New code for digital contact tracing - ggilani: 09/03/20
	P.DoDigitalContactTracing = Params::get_int(params, pre_params, "Include digital contact tracing", 0, &P);
	if (P.DoDigitalContactTracing != 0)
	{
		P.DigitalContactTracing_CellIncThresh = Params::get_double(params, pre_params, "Digital contact tracing trigger incidence per cell", 1000000000, &P);

		P.PropPopUsingDigitalContactTracing  = Params::get_double(params, pre_params, "Proportion of population or households covered by digital contact tracing", 1, &P);
		Params::get_double_vec(params, pre_params, "Proportion of smartphone users by age", P.ProportionSmartphoneUsersByAge, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, &P);
		if (P.DoPlaces != 0)
		{
			P.ClusterDigitalContactUsers = Params::get_int(params, pre_params, "Cluster digital app clusters by household", 0, &P); // by default, don't cluster by location
		}
		else
		{
			P.ClusterDigitalContactUsers = 0;
		}
		P.ProportionDigitalContactsIsolate = Params::get_double(params, pre_params, "Proportion of digital contacts who self-isolate", 0, &P);
		P.MaxDigitalContactsToTrace = Params::get_int(params, pre_params, "Maximum number of contacts to trace per index case", MAX_CONTACTS, &P);
		P.DigitalContactTracingDelay = Params::get_double(params, pre_params, "Delay between isolation of index case and contacts", P.ModelTimeStep, &P);
		//we really need one timestep between to make sure contact is not processed before index
		if (P.DigitalContactTracingDelay == 0) P.DigitalContactTracingDelay = P.ModelTimeStep;
		P.LengthDigitalContactIsolation = Params::get_double(params, pre_params, "Length of self-isolation for digital contacts", 0, &P);
		P.ScalingFactorSpatialDigitalContacts = Params::get_double(params, pre_params, "Spatial scaling factor - digital contact tracing", 1, &P);
		P.ScalingFactorPlaceDigitalContacts = Params::get_double(params, pre_params, "Place scaling factor - digital contact tracing", 1, &P);
		P.DigitalContactTracingTimeStartBase = Params::get_double(params, pre_params, "Digital contact tracing start time", USHRT_MAX / P.TimeStepsPerDay, &P);
		P.DigitalContactTracingPolicyDuration  = Params::get_double(params, pre_params, "Duration of digital contact tracing policy", 7, &P);
		P.OutputDigitalContactTracing = Params::get_int(params, pre_params, "Output digital contact tracing", 0, &P);
		P.OutputDigitalContactDist = Params::get_int(params, pre_params, "Output digital contact distribution", 0, &P);

		if (P.DoInterventionDelaysByAdUnit)
		{
			double AdunitDelayToDCT[MAX_ADUNITS];
			double AdunitDurationDCT[MAX_ADUNITS];

			Params::get_double_vec(params, pre_params, "Delay to digital contact tracing by admin unit", AdunitDelayToDCT, P.NumAdunits, 0, P.NumAdunits, &P);
			Params::get_double_vec(params, pre_params, "Duration of digital contact tracing by admin unit", AdunitDurationDCT, P.NumAdunits, 0, P.NumAdunits, &P);
			for (i = 0; i < P.NumAdunits; i++)
			{
				AdUnits[i].DCTDelay = AdunitDelayToDCT[i];
				AdUnits[i].DCTDuration = AdunitDurationDCT[i];
			}
		}
		P.DCTIsolateIndexCases = Params::get_int(params, pre_params, "Isolate index cases in digital contact tracing", 1, &P);
		P.DCTCaseIsolationEffectiveness = Params::get_double(params, pre_params, "Residual contacts after digital contact tracing isolation", P.CaseIsolationEffectiveness, &P);
		P.DCTCaseIsolationHouseEffectiveness = Params::get_double(params, pre_params, "Residual household contacts after digital contact tracing isolation", P.CaseIsolationHouseEffectiveness, &P);
		//initialise total number of users to 0
		P.NDigitalContactUsers = 0;
		P.NDigitalHouseholdUsers = 0;

		P.DelayFromIndexCaseDetectionToDCTIsolation = Params::get_double(params, pre_params, "Delay between symptom onset and isolation for index case", 0, &P);
		P.DoDCTTest = Params::get_int(params, pre_params, "Test index cases and contacts", 0, &P);
		P.DelayToTestIndexCase = Params::get_double(params, pre_params, "Delay to test index case", 1, &P);
		P.DelayToTestDCTContacts = Params::get_double(params, pre_params, "Delay to test DCT contacts", 7, &P);
		P.SpecificityDCT = Params::get_double(params, pre_params, "Testing specificity - DCT", 1, &P);
		P.SensitivityDCT = Params::get_double(params, pre_params, "Testing sensitivity - DCT", 1, &P);
		P.FindContactsOfDCTContacts = Params::get_int(params, pre_params, "Find contacts of digital contacts", 0, &P);
		P.RemoveContactsOfNegativeIndexCase = Params::get_int(params, pre_params, "Remove contacts of a negative index case", 0, &P);
	}
	else
	{
		//Set these to 1 so it doesn't interfere with code if we aren't using digital contact tracing.

		P.ScalingFactorSpatialDigitalContacts = 1;
		P.ScalingFactorPlaceDigitalContacts = 1;
	}

	Files::xfprintf_stderr("Reading Params 10\n"); fflush(stderr);


	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** PLACE CLOSURE
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****


	P.PlaceCloseCellIncThresh1 = Params::get_int(params, pre_params, "Trigger incidence per cell for place closure", 1000000000, &P);
	P.PlaceCloseCellIncThresh2 = Params::get_int(params, pre_params, "Trigger incidence per cell for second place closure", 1000000000, &P);
	if (P.PlaceCloseCellIncThresh1 < 0) P.PlaceCloseCellIncThresh1 = 1000000000;
	if (P.PlaceCloseCellIncThresh2 < 0) P.PlaceCloseCellIncThresh2 = 1000000000;
	P.PlaceCloseCellIncStopThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for end of place closure", 0, &P);
	P.PlaceCloseDelayMean = Params::get_double(params, pre_params, "Delay to start place closure", 0, &P);
	P.PlaceCloseDurationBase = Params::get_double(params, pre_params, "Duration of place closure", 7, &P);
	P.PlaceCloseDuration2 = Params::get_double(params, pre_params, "Duration of second place closure", 7, &P);
	if (P.DoPlaces != 0)
	{
		Params::get_double_vec(params, pre_params, "Proportion of places remaining open after closure by place type", P.PlaceCloseEffect, P.PlaceTypeNum, 1, NUM_PLACE_TYPES, &P);
		Params::get_double_vec(params, pre_params, "Proportional attendance after closure by place type", P.PlaceClosePropAttending, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
	}
	if (P.DoHouseholds != 0)
		P.PlaceCloseHouseholdRelContact = Params::get_double(params, pre_params, "Relative household contact rate after closure", 1, &P);
	P.PlaceCloseSpatialRelContact = Params::get_double(params, pre_params, "Relative spatial contact rate after closure", 1, &P);

	P.DoHolidays = Params::get_int(pre_params, adm_params, "Include holidays", 0, &P);
	if (P.DoHolidays != 0)
	{
		Params::get_double_vec(pre_params, adm_params, "Proportion of places remaining open during holidays by place type", P.HolidayEffect, P.PlaceTypeNum, 1, NUM_PLACE_TYPES, &P);
		P.NumHolidays = Params::get_int(pre_params, adm_params, "Number of holidays", 0, &P);
		if (P.NumHolidays > DAYS_PER_YEAR) P.NumHolidays = DAYS_PER_YEAR;
		if (P.NumHolidays > 0)
		{
			Params::req_double_vec(pre_params, adm_params, "Holiday start times", P.HolidayStartTime, P.NumHolidays, &P);
			Params::req_double_vec(pre_params, adm_params, "Holiday durations", P.HolidayDuration, P.NumHolidays, &P);
		}
	}
	else
		P.NumHolidays = 0;
	P.PlaceCloseRadius = Params::get_double(params, pre_params, "Minimum radius for place closure", 0, &P);
	P.PlaceCloseTimeStartBase = Params::get_double(params, pre_params, "Place closure start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.PlaceCloseTimeStartBase2 = Params::get_double(params, pre_params, "Place closure second start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.DoPlaceCloseOnceOnly = Params::get_int(params, pre_params, "Places close only once", 0, &P);
	//if (P.DoPlaceCloseOnceOnly) P.DoPlaceCloseOnceOnly = 4; //// don't need this anymore with TreatStat option. Keep it as a boolean.
	P.PlaceCloseIncTrig1 = Params::get_int(params, pre_params, "Place closure incidence threshold", 1, &P);
	P.PlaceCloseIncTrig2 = Params::get_int(params, pre_params, "Place closure second incidence threshold", P.PlaceCloseIncTrig1, &P);
	P.PlaceCloseFracIncTrig = Params::get_double(params, pre_params, "Place closure fractional incidence threshold", 0, &P);
	if ((P.DoAdUnits != 0) && (P.DoPlaces != 0))
	{
		P.PlaceCloseByAdminUnit = Params::get_int(params, pre_params, "Place closure in administrative units rather than rings", 0, &P);
		P.PlaceCloseAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for place closure", 1, &P);
		Params::get_int_vec(params, pre_params, "Place types to close for admin unit closure (0/1 array)", P.PlaceCloseAdunitPlaceTypes, P.PlaceTypeNum, 0, P.PlaceTypeNum, &P);
		P.PlaceCloseCasePropThresh = Params::get_double(params, pre_params, "Cumulative proportion of place members needing to become sick for admin unit closure", 2, &P);
		P.PlaceCloseAdunitPropThresh = Params::get_double(params, pre_params, "Proportion of places in admin unit needing to pass threshold for place closure", 2, &P);
		if ((P.PlaceCloseAdminUnitDivisor < 1) || (P.PlaceCloseByAdminUnit == 0)) P.PlaceCloseAdminUnitDivisor = 1;
	}
	else
	{
		P.PlaceCloseAdminUnitDivisor = 1; P.PlaceCloseByAdminUnit = 0;
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** SOCIAL DISTANCING
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	P.SocDistCellIncThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for social distancing", 1000000000, &P);
	P.SocDistCellIncStopThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for end of social distancing", 0, &P);
	P.SocDistDuration = Params::get_double(params, pre_params, "Duration of social distancing", 7, &P);
	P.SocDistDuration2 = Params::get_double(params, pre_params, "Duration of social distancing after change", 7, &P);
	if (P.DoPlaces != 0)
	{
	  Params::get_double_vec(params, pre_params, "Relative place contact rate given social distancing by place type", P.SocDistPlaceEffect, P.PlaceTypeNum, 1, NUM_PLACE_TYPES, &P);
		Params::get_double_vec(params, pre_params, "Relative place contact rate given enhanced social distancing by place type", P.EnhancedSocDistPlaceEffect, P.PlaceTypeNum, 1, NUM_PLACE_TYPES, &P);
		if (Params::param_found(params, pre_params, "Relative place contact rate given social distancing by place type after change"))
		{
		  Params::req_double_vec(params, pre_params, "Relative place contact rate given social distancing by place type after change", P.SocDistPlaceEffect2, P.PlaceTypeNum, &P);
		} else {
		  for (i = 0; i < NUM_PLACE_TYPES; i++) P.SocDistPlaceEffect2[i] = P.SocDistPlaceEffect[i];
		}

		if (Params::param_found(params, pre_params, "Relative place contact rate given enhanced social distancing by place type after change"))
		{
			Params::req_double_vec(params, pre_params, "Relative place contact rate given enhanced social distancing by place type after change", P.EnhancedSocDistPlaceEffect2, P.PlaceTypeNum, &P);
		}
		else
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.EnhancedSocDistPlaceEffect2[i] = P.EnhancedSocDistPlaceEffect[i];
		}
	}
	if (P.DoHouseholds != 0)
	{
		P.SocDistHouseholdEffect = Params::get_double(params, pre_params, "Relative household contact rate given social distancing", 1, &P);
		P.EnhancedSocDistHouseholdEffect = Params::get_double(params, pre_params, "Relative household contact rate given enhanced social distancing", 1, &P);
		P.SocDistHouseholdEffect2 = Params::get_double(params, pre_params, "Relative household contact rate given social distancing  after change", P.SocDistHouseholdEffect, &P);
		P.EnhancedSocDistHouseholdEffect2 = Params::get_double(params, pre_params, "Relative household contact rate given enhanced social distancing after change", P.EnhancedSocDistHouseholdEffect, &P);
		P.EnhancedSocDistClusterByHousehold = Params::get_int(params, pre_params, "Cluster compliance with enhanced social distancing by household", 0, &P);
	}
	else
		P.EnhancedSocDistClusterByHousehold = 0;
	P.SocDistSpatialEffect = Params::get_double(params, pre_params, "Relative spatial contact rate given social distancing", 1, &P);
	P.SocDistSpatialEffect2 = Params::get_double(params, pre_params, "Relative spatial contact rate given social distancing after change", P.SocDistSpatialEffect, &P);
	P.SocDistRadius = Params::get_double(params, pre_params, "Minimum radius for social distancing", 0, &P);
	P.SocDistTimeStartBase = Params::get_double(params, pre_params, "Social distancing start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.SocDistChangeDelay = Params::get_double(params, pre_params, "Delay for change in effectiveness of social distancing", USHRT_MAX / P.TimeStepsPerDay, &P);
	if (Params::param_found(params, pre_params, "Proportion compliant with enhanced social distancing by age group")) {
		Params::req_double_vec(params, pre_params, "Proportion compliant with enhanced social distancing by age group", P.EnhancedSocDistProportionCompliant, NUM_AGE_GROUPS, &P);

	} else {
		t = Params::get_double(params, pre_params, "Proportion compliant with enhanced social distancing", 0, &P);
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.EnhancedSocDistProportionCompliant[i] = t;
	}

	P.EnhancedSocDistSpatialEffect = Params::get_double(params, pre_params, "Relative spatial contact rate given enhanced social distancing", 1, &P);
	P.EnhancedSocDistSpatialEffect2 = Params::get_double(params, pre_params, "Relative spatial contact rate given enhanced social distancing after change", P.EnhancedSocDistSpatialEffect, &P);

	P.DoSocDistOnceOnly = Params::get_int(params, pre_params, "Social distancing only once", 0, &P);
	//if (P.DoSocDistOnceOnly) P.DoSocDistOnceOnly = 4;  //// don't need this anymore with TreatStat option. Keep it as a boolean.

	P.AirportCloseEffectiveness = Params::get_double(params, pre_params, "Airport closure effectiveness", 0, &P);
	P.AirportCloseEffectiveness = 1.0 - P.AirportCloseEffectiveness;
	P.AirportCloseTimeStartBase = Params::get_double(params, pre_params, "Airport closure start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.AirportCloseDuration = Params::get_double(params, pre_params, "Airport closure duration", USHRT_MAX / P.TimeStepsPerDay, &P);

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** HOUSEHOLD QUARANTINE
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	if (P.DoHouseholds != 0)
	{
		P.DoHQretrigger = Params::get_int(params, pre_params, "Retrigger household quarantine with each new case in quarantine window", 0, &P);
		P.HQuarantineTimeStartBase = Params::get_double(params, pre_params, "Household quarantine start time", USHRT_MAX / P.TimeStepsPerDay, &P);
		P.HQuarantineDelay = Params::get_double(params, pre_params, "Delay to start household quarantine", 0, &P);
		P.HQuarantineHouseDuration = Params::get_double(params, pre_params, "Length of time households are quarantined", 0, &P);
		P.HQuarantinePolicyDuration = Params::get_double(params, pre_params, "Duration of household quarantine policy", USHRT_MAX / P.TimeStepsPerDay, &P);
		P.HQuarantineHouseEffect = Params::get_double(params, pre_params, "Relative household contact rate after quarantine", 1, &P);
		if (P.DoPlaces != 0)
		{
			Params::get_double_vec(params, pre_params, "Residual place contacts after household quarantine by place type", P.HQuarantinePlaceEffect, P.PlaceTypeNum, 1,NUM_PLACE_TYPES, &P);
		}
		P.HQuarantineSpatialEffect = Params::get_double(params, pre_params, "Residual spatial contacts after household quarantine", 1, &P);
		P.HQuarantinePropHouseCompliant = Params::get_double(params, pre_params, "Household level compliance with quarantine", 1, &P);
		P.HQuarantinePropIndivCompliant = Params::get_double(params, pre_params, "Individual level compliance with quarantine", 1, &P);
	}
	else
		P.HQuarantineTimeStartBase = 1e10;
	P.CaseIsolationTimeStartBase = Params::get_double(params, pre_params, "Case isolation start time", USHRT_MAX / P.TimeStepsPerDay, &P);
	P.CaseIsolationProp = Params::get_double(params, pre_params, "Proportion of detected cases isolated", 0, &P);
	P.CaseIsolationDelay = Params::get_double(params, pre_params, "Delay to start case isolation", 0, &P);
	P.CaseIsolationDuration = Params::get_double(params, pre_params, "Duration of case isolation", 0, &P);
	P.CaseIsolationPolicyDuration = Params::get_double(params, pre_params, "Duration of case isolation policy", 1e10, &P);
	P.CaseIsolationEffectiveness = Params::get_double(params, pre_params, "Residual contacts after case isolation", 1, &P);
	if (P.DoHouseholds != 0)
	{
		P.CaseIsolationHouseEffectiveness = Params::get_double(params, pre_params, "Residual household contacts after case isolation", P.CaseIsolationEffectiveness, &P);
	}

	Files::xfprintf_stderr("Reading Params 11\n"); fflush(stderr);


	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** VARIABLE EFFICACIES OVER TIME
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	P.VaryEfficaciesOverTime  = Params::get_int(params, pre_params, "Vary efficacies over time", 0, &P);
	//// **** number of change times
	if (P.VaryEfficaciesOverTime == 0)
	{
		P.Num_SD_ChangeTimes = 1;
		P.Num_CI_ChangeTimes = 1;
		P.Num_HQ_ChangeTimes = 1;
		P.Num_PC_ChangeTimes = 1;
		P.Num_DCT_ChangeTimes = 1;
	}
	else
	{
		P.Num_SD_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of social distancing", 1, &P);
		P.Num_CI_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of case isolation", 1, &P);
		P.Num_HQ_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of household quarantine", 1, &P);
		P.Num_PC_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of place closure", 1, &P);
		P.Num_DCT_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of digital contact tracing", 1, &P);
	}

	//// **** change times:
	//// By default, initialize first change time to zero and all subsequent change times to occur after simulation time, i.e. single value e.g. of efficacy for social distancing.
	P.SD_ChangeTimes			[0] = 0;
	P.CI_ChangeTimes			[0] = 0;
	P.HQ_ChangeTimes			[0] = 0;
	P.PC_ChangeTimes			[0] = 0;
	P.DCT_ChangeTimes			[0] = 0;
	for (int ChangeTime = 1; ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES; ChangeTime++)
	{
		P.SD_ChangeTimes			[ChangeTime] = 1e10;
		P.CI_ChangeTimes			[ChangeTime] = 1e10;
		P.HQ_ChangeTimes			[ChangeTime] = 1e10;
		P.PC_ChangeTimes			[ChangeTime] = 1e10;
		P.DCT_ChangeTimes			[ChangeTime] = 1e10;
		P.CFR_ChangeTimes_CalTime	[ChangeTime] = INT32_MAX; // Out of bounds for int
	}
	//// Get real values from (pre)param file

	Params::get_double_vec(params, pre_params, "Change times for levels of social distancing", P.SD_ChangeTimes, P.Num_SD_ChangeTimes, 0, P.Num_SD_ChangeTimes, &P);
	Params::get_double_vec(params, pre_params, "Change times for levels of case isolation", P.CI_ChangeTimes, P.Num_CI_ChangeTimes, 0, P.Num_CI_ChangeTimes, &P);
	Params::get_double_vec(params, pre_params, "Change times for levels of household quarantine", P.HQ_ChangeTimes, P.Num_HQ_ChangeTimes, 0, P.Num_HQ_ChangeTimes, &P);
	Params::get_double_vec(params, pre_params, "Change times for levels of place closure", P.PC_ChangeTimes, P.Num_PC_ChangeTimes, 0, P.Num_PC_ChangeTimes, &P);
	Params::get_double_vec(params, pre_params, "Change times for levels of digital contact tracing", P.DCT_ChangeTimes, P.Num_DCT_ChangeTimes, 0, P.Num_DCT_ChangeTimes, &P);

	// initialize to zero (regardless of whether doing places or households).
	for (int ChangeTime = 0; ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES; ChangeTime++)
	{
		//// **** "efficacies"
		//// spatial
		P.SD_SpatialEffects_OverTime				[ChangeTime] = 0;
		P.Enhanced_SD_SpatialEffects_OverTime		[ChangeTime] = 0;
		P.CI_SpatialAndPlaceEffects_OverTime		[ChangeTime] = 0;
		P.HQ_SpatialEffects_OverTime				[ChangeTime] = 0;
		P.PC_SpatialEffects_OverTime				[ChangeTime] = 0;
		P.DCT_SpatialAndPlaceEffects_OverTime		[ChangeTime] = 0;

		//// Household
		P.SD_HouseholdEffects_OverTime			[ChangeTime] = 0;
		P.Enhanced_SD_HouseholdEffects_OverTime	[ChangeTime] = 0;
		P.CI_HouseholdEffects_OverTime			[ChangeTime] = 0;
		P.HQ_HouseholdEffects_OverTime			[ChangeTime] = 0;
		P.PC_HouseholdEffects_OverTime			[ChangeTime] = 0;
		P.DCT_HouseholdEffects_OverTime			[ChangeTime] = 0;

		//// place
		for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
		{
			P.SD_PlaceEffects_OverTime			[ChangeTime][PlaceType] = 0;
			P.Enhanced_SD_PlaceEffects_OverTime	[ChangeTime][PlaceType] = 0;
			P.HQ_PlaceEffects_OverTime			[ChangeTime][PlaceType] = 0;
			P.PC_PlaceEffects_OverTime			[ChangeTime][PlaceType] = 0;
		}
		P.PC_Durs_OverTime[ChangeTime] = 0;

		//// **** compliance
		P.CI_Prop_OverTime					[ChangeTime] = 0;
		P.HQ_Individual_PropComply_OverTime	[ChangeTime] = 0;
		P.HQ_Household_PropComply_OverTime	[ChangeTime] = 0;
		P.DCT_Prop_OverTime					[ChangeTime] = 0;
	}


	//// **** "efficacies": by default, initialize to values read in previously.
	///// spatial contact rates rates over time (and place too for CI and DCT)
	bool v = (P.VaryEfficaciesOverTime == 0);
	Params::get_double_vec_ff(v, params, pre_params, "Relative spatial contact rates over time given social distancing", P.SD_SpatialEffects_OverTime, P.Num_SD_ChangeTimes, P.SocDistSpatialEffect, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Relative spatial contact rates over time given enhanced social distancing", P.Enhanced_SD_SpatialEffects_OverTime, P.Num_SD_ChangeTimes, P.EnhancedSocDistSpatialEffect, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Residual contacts after case isolation over time", P.CI_SpatialAndPlaceEffects_OverTime, P.Num_CI_ChangeTimes, P.CaseIsolationEffectiveness, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Residual spatial contacts over time after household quarantine", P.HQ_SpatialEffects_OverTime, P.Num_HQ_ChangeTimes, P.HQuarantineSpatialEffect, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Relative spatial contact rates over time after place closure", P.PC_SpatialEffects_OverTime, P.Num_PC_ChangeTimes, P.PlaceCloseSpatialRelContact, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Residual contacts after digital contact tracing isolation over time", P.DCT_SpatialAndPlaceEffects_OverTime, P.Num_DCT_ChangeTimes, P.DCTCaseIsolationEffectiveness, &P);

	///// Household contact rates over time
	if (P.DoHouseholds != 0)
	{
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time given social distancing", P.SD_HouseholdEffects_OverTime, P.Num_SD_ChangeTimes, P.SocDistHouseholdEffect, &P);
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time given enhanced social distancing", P.Enhanced_SD_HouseholdEffects_OverTime, P.Num_SD_ChangeTimes, P.EnhancedSocDistHouseholdEffect, &P);
		Params::get_double_vec_ff(v, params, pre_params, "Residual household contacts after case isolation over time", P.CI_HouseholdEffects_OverTime, P.Num_CI_ChangeTimes, P.CaseIsolationHouseEffectiveness, &P);
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time after quarantine", P.HQ_HouseholdEffects_OverTime, P.Num_HQ_ChangeTimes, P.HQuarantineHouseEffect, &P);
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time after place closure", P.PC_HouseholdEffects_OverTime, P.Num_PC_ChangeTimes, P.PlaceCloseHouseholdRelContact, &P);
		Params::get_double_vec_ff(v, params, pre_params, "Residual household contacts after digital contact tracing isolation over time", P.DCT_HouseholdEffects_OverTime, P.Num_DCT_ChangeTimes, P.DCTCaseIsolationHouseEffectiveness, &P);
	}

	///// place contact rates over time
	if (P.DoPlaces != 0)
	{ 
		//// soc dist
		Params::get_double_matrix(params, pre_params, "Relative place contact rates over time given social distancing by place type", P.SD_PlaceEffects_OverTime, P.Num_SD_ChangeTimes, P.PlaceTypeNum, 0, &P);
		if (v || !Params::param_found(params, pre_params, "Relative place contact rates over time given social distancing by place type"))
			for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) //// by default populate to values of P.SocDistPlaceEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.SocDistPlaceEffect[PlaceType];

		//// enhanced soc dist
		Params::get_double_matrix(params, pre_params, "Relative place contact rates over time given enhanced social distancing by place type", P.Enhanced_SD_PlaceEffects_OverTime, P.Num_SD_ChangeTimes, P.PlaceTypeNum, 0, &P);
		if (v || !Params::param_found(params, pre_params, "Relative place contact rates over time given enhanced social distancing by place type"))
			for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) //// by default populate to values of P.EnhancedSocDistPlaceEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.Enhanced_SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.EnhancedSocDistPlaceEffect[PlaceType];

		//// household quarantine
		Params::get_double_matrix(params, pre_params, "Residual place contacts over time after household quarantine by place type", P.HQ_PlaceEffects_OverTime, P.Num_HQ_ChangeTimes, P.PlaceTypeNum, 0, &P);
		if (v || !Params::param_found(params, pre_params, "Residual place contacts over time after household quarantine by place type"))
			for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) //// by default populate to values of P.HQuarantinePlaceEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.HQ_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.HQuarantinePlaceEffect[PlaceType];

		//// place closure
		Params::get_double_matrix(params, pre_params, "Proportion of places remaining open after closure by place type over time", P.PC_PlaceEffects_OverTime, P.Num_PC_ChangeTimes, P.PlaceTypeNum, 0, &P);
		if (v || !Params::param_found(params, pre_params, "Proportion of places remaining open after closure by place type over time"))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) //// by default populate to values of P.PlaceCloseEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.PC_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.PlaceCloseEffect[PlaceType];

		Params::get_double_matrix(params, pre_params, "Proportional attendance after closure by place type over time", P.PC_PropAttending_OverTime, P.Num_PC_ChangeTimes, P.PlaceTypeNum, 0, &P);
		if (v || !Params::param_found(params, pre_params, "Proportional attendance after closure by place type over time"))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) //// by default populate to values of P.PlaceClosePropAttending
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.PC_PropAttending_OverTime[ChangeTime][PlaceType] = P.PlaceClosePropAttending[PlaceType];
	}


	//// ****  compliance
	Params::get_double_vec_ff(v, params, pre_params, "Proportion of detected cases isolated over time", P.CI_Prop_OverTime, P.Num_CI_ChangeTimes, P.CaseIsolationProp, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Individual level compliance with quarantine over time", P.HQ_Individual_PropComply_OverTime, P.Num_HQ_ChangeTimes, P.HQuarantinePropIndivCompliant, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Household level compliance with quarantine over time", P.HQ_Household_PropComply_OverTime, P.Num_HQ_ChangeTimes, P.HQuarantinePropHouseCompliant, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Proportion of digital contacts who self-isolate over time", P.DCT_Prop_OverTime, P.Num_DCT_ChangeTimes, P.ProportionDigitalContactsIsolate, &P);
	Params::get_int_vec_ff(v, params, pre_params, "Maximum number of contacts to trace per index case over time", P.DCT_MaxToTrace_OverTime, P.Num_DCT_ChangeTimes, P.MaxDigitalContactsToTrace, &P);
	if (P.DoPlaces != 0)
	{
		//// ****  thresholds
		//// place closure (global threshold)
		Params::get_int_vec_ff(v, params, pre_params, "Place closure incidence threshold over time", P.PC_IncThresh_OverTime, P.Num_PC_ChangeTimes, P.PlaceCloseIncTrig1, &P);
		Params::get_double_vec_ff(v, params, pre_params, "Place closure fractional incidence threshold over time", P.PC_FracIncThresh_OverTime, P.Num_PC_ChangeTimes, P.PlaceCloseFracIncTrig, &P);
		Params::get_int_vec_ff(v, params, pre_params, "Trigger incidence per cell for place closure over time", P.PC_CellIncThresh_OverTime, P.Num_PC_ChangeTimes, P.PlaceCloseCellIncThresh1, &P);
		for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++)
		  if (P.PC_CellIncThresh_OverTime[ChangeTime]<0) P.PC_CellIncThresh_OverTime[ChangeTime] = 1000000000; // allows -1 to be used as a proxy for no cell-based triggering
	}
	//// household quarantine
	Params::get_double_vec_ff(v, params, pre_params, "Household quarantine trigger incidence per cell over time", P.HQ_CellIncThresh_OverTime, P.Num_HQ_ChangeTimes, P.HHQuar_CellIncThresh, &P);
	Params::get_double_vec_ff(v, params, pre_params, "Case isolation trigger incidence per cell over time", P.CI_CellIncThresh_OverTime, P.Num_CI_ChangeTimes, P.CaseIsolation_CellIncThresh, &P);
	Params::get_int_vec_ff(v, params, pre_params, "Trigger incidence per cell for social distancing over time", P.SD_CellIncThresh_OverTime, P.Num_SD_ChangeTimes, P.SocDistCellIncThresh, &P);

	//// **** Durations (later add Case isolation and Household quarantine)
	// place closure
	Params::get_double_vec_ff(v, params, pre_params, "Duration of place closure over time", P.PC_Durs_OverTime, P.Num_PC_ChangeTimes, P.PlaceCloseDurationBase, &P);

	//// Guards: make unused change values in array equal to final used value
	if (P.VaryEfficaciesOverTime)
	{
		//// soc dist
		for (int SD_ChangeTime = P.Num_SD_ChangeTimes; SD_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; SD_ChangeTime++)
		{
			//// non-enhanced
			P.SD_SpatialEffects_OverTime	[SD_ChangeTime] = P.SD_SpatialEffects_OverTime		[P.Num_SD_ChangeTimes - 1];
			P.SD_HouseholdEffects_OverTime	[SD_ChangeTime] = P.SD_HouseholdEffects_OverTime	[P.Num_SD_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
				P.SD_PlaceEffects_OverTime[SD_ChangeTime][PlaceType] = P.SD_PlaceEffects_OverTime[P.Num_SD_ChangeTimes - 1][PlaceType];
			//// enhanced
			P.Enhanced_SD_SpatialEffects_OverTime	[SD_ChangeTime] = P.Enhanced_SD_SpatialEffects_OverTime		[P.Num_SD_ChangeTimes - 1];
			P.Enhanced_SD_HouseholdEffects_OverTime	[SD_ChangeTime] = P.Enhanced_SD_HouseholdEffects_OverTime	[P.Num_SD_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
				P.Enhanced_SD_PlaceEffects_OverTime[SD_ChangeTime][PlaceType] = P.Enhanced_SD_PlaceEffects_OverTime[P.Num_SD_ChangeTimes - 1][PlaceType];

			P.SD_CellIncThresh_OverTime				[SD_ChangeTime] = P.SD_CellIncThresh_OverTime				[P.Num_SD_ChangeTimes - 1];
		}

		//// case isolation
		for (int CI_ChangeTime = P.Num_CI_ChangeTimes; CI_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; CI_ChangeTime++)
		{
			P.CI_SpatialAndPlaceEffects_OverTime[CI_ChangeTime] = P.CI_SpatialAndPlaceEffects_OverTime	[P.Num_CI_ChangeTimes - 1];
			P.CI_HouseholdEffects_OverTime		[CI_ChangeTime] = P.CI_HouseholdEffects_OverTime		[P.Num_CI_ChangeTimes - 1];
			P.CI_Prop_OverTime					[CI_ChangeTime] = P.CI_Prop_OverTime					[P.Num_CI_ChangeTimes - 1];
			P.CI_CellIncThresh_OverTime			[CI_ChangeTime] = P.CI_CellIncThresh_OverTime			[P.Num_CI_ChangeTimes - 1];
		}

		//// household quarantine
		for (int HQ_ChangeTime = P.Num_HQ_ChangeTimes; HQ_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; HQ_ChangeTime++)
		{
			P.HQ_SpatialEffects_OverTime	[HQ_ChangeTime] = P.HQ_SpatialEffects_OverTime	[P.Num_HQ_ChangeTimes - 1];
			P.HQ_HouseholdEffects_OverTime	[HQ_ChangeTime] = P.HQ_HouseholdEffects_OverTime[P.Num_HQ_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
				P.HQ_PlaceEffects_OverTime[HQ_ChangeTime][PlaceType] = P.HQ_PlaceEffects_OverTime[P.Num_HQ_ChangeTimes - 1][PlaceType];

			P.HQ_Individual_PropComply_OverTime	[HQ_ChangeTime] = P.HQ_Individual_PropComply_OverTime	[P.Num_HQ_ChangeTimes - 1];
			P.HQ_Household_PropComply_OverTime	[HQ_ChangeTime] = P.HQ_Household_PropComply_OverTime	[P.Num_HQ_ChangeTimes - 1];

			P.HQ_CellIncThresh_OverTime			[HQ_ChangeTime] = P.HQ_CellIncThresh_OverTime			[P.Num_HQ_ChangeTimes - 1];
		}

		//// place closure
		for (int PC_ChangeTime = P.Num_PC_ChangeTimes; PC_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; PC_ChangeTime++)
		{
			P.PC_SpatialEffects_OverTime	[PC_ChangeTime] = P.PC_SpatialEffects_OverTime	[P.Num_PC_ChangeTimes - 1];
			P.PC_HouseholdEffects_OverTime	[PC_ChangeTime] = P.PC_HouseholdEffects_OverTime[P.Num_PC_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
			{
				P.PC_PlaceEffects_OverTime[PC_ChangeTime][PlaceType] = P.PC_PlaceEffects_OverTime[P.Num_PC_ChangeTimes - 1][PlaceType];
				P.PC_PropAttending_OverTime[PC_ChangeTime][PlaceType] = P.PC_PropAttending_OverTime[P.Num_PC_ChangeTimes - 1][PlaceType];
			}

			P.PC_IncThresh_OverTime			[PC_ChangeTime]	= P.PC_IncThresh_OverTime		[P.Num_PC_ChangeTimes - 1];
			P.PC_FracIncThresh_OverTime		[PC_ChangeTime]	= P.PC_FracIncThresh_OverTime	[P.Num_PC_ChangeTimes - 1];
			P.PC_CellIncThresh_OverTime		[PC_ChangeTime]	= P.PC_CellIncThresh_OverTime	[P.Num_PC_ChangeTimes - 1];
		}

		//// digital contact tracing
		for (int DCT_ChangeTime = P.Num_DCT_ChangeTimes; DCT_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; DCT_ChangeTime++)
		{
			P.DCT_SpatialAndPlaceEffects_OverTime	[DCT_ChangeTime] = P.DCT_SpatialAndPlaceEffects_OverTime[P.Num_DCT_ChangeTimes - 1];
			P.DCT_HouseholdEffects_OverTime			[DCT_ChangeTime] = P.DCT_HouseholdEffects_OverTime		[P.Num_DCT_ChangeTimes - 1];
			P.DCT_Prop_OverTime						[DCT_ChangeTime] = P.DCT_Prop_OverTime					[P.Num_DCT_ChangeTimes - 1];
			P.DCT_MaxToTrace_OverTime				[DCT_ChangeTime] = P.DCT_MaxToTrace_OverTime			[P.Num_DCT_ChangeTimes - 1];
		}
	}

	Files::xfprintf_stderr("Reading Params 12\n"); fflush(stderr);


	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** CFR SCALINGS OVER TIME (logic to set this up is the same as for VARIABLE EFFICACIES OVER TIME)
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	P.Num_CFR_ChangeTimes = Params::get_int(params, pre_params, "Num_CFR_ChangeTimes", 1, &P);

	//// By default, initialize first change time to zero and all subsequent change times to occur after simulation time, i.e. single value e.g. of Critical CFR.
	P.CFR_ChangeTimes_CalTime[0] = 0;

  for (int ChangeTime = 1; ChangeTime < MAX_NUM_CFR_CHANGE_TIMES; ChangeTime++) P.CFR_ChangeTimes_CalTime[ChangeTime] = INT32_MAX;
	Params::get_int_vec(params, pre_params, "CFR_ChangeTimes_CalTime", P.CFR_ChangeTimes_CalTime, P.Num_CFR_ChangeTimes, INT32_MAX, P.Num_CFR_ChangeTimes, &P);

	// Get various CFR scalings. 
	Params::get_double_vec(params, pre_params, "CFR_TimeScaling_Critical", P.CFR_TimeScaling_Critical, P.Num_CFR_ChangeTimes, 1, P.Num_CFR_ChangeTimes, &P);
	Params::get_double_vec(params, pre_params, "CFR_TimeScaling_SARI", P.CFR_TimeScaling_SARI, P.Num_CFR_ChangeTimes, 1, P.Num_CFR_ChangeTimes, &P);
	Params::get_double_vec(params, pre_params, "CFR_TimeScaling_ILI", P.CFR_TimeScaling_ILI, P.Num_CFR_ChangeTimes, 1, P.Num_CFR_ChangeTimes, &P);

	Files::xfprintf_stderr("Reading Params 13\n"); fflush(stderr);


	if(P.FitIter == 0)
	{
		if (P.DoPlaces != 0)
		{
			P.KeyWorkerPopNum = Params::get_int(params, pre_params, "Number of key workers randomly distributed in the population", 0, &P);
			Params::get_int_vec(params, pre_params, "Number of key workers in different places by place type", P.KeyWorkerPlaceNum, P.PlaceTypeNum, 0, NUM_PLACE_TYPES, &P);
			Params::get_double_vec(params, pre_params, "Proportion of staff who are key workers per chosen place by place type", P.KeyWorkerPropInKeyPlaces, P.PlaceTypeNum, 1, NUM_PLACE_TYPES, &P);
			P.KeyWorkerProphCellIncThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for key worker prophylaxis", 1000000000, &P);
			P.KeyWorkerProphTimeStartBase = Params::get_double(params, pre_params, "Key worker prophylaxis start time", USHRT_MAX / P.TimeStepsPerDay, &P);
			P.KeyWorkerProphDuration = Params::get_double(params, pre_params, "Duration of key worker prophylaxis", 0, &P);
			P.KeyWorkerProphRenewalDuration = Params::get_double(params, pre_params, "Time interval from start of key worker prophylaxis before policy restarted", P.KeyWorkerProphDuration, &P);
			if (P.DoHouseholds != 0)
			{
				P.KeyWorkerHouseProp = Params::get_double(params, pre_params, "Proportion of key workers whose households are also treated as key workers", 0, &P);
			}
			P.KeyWorkerProphRadius = Params::get_double(params, pre_params, "Minimum radius for key worker prophylaxis", 0, &P);
		}
		else
		{
			P.KeyWorkerPopNum = 0;
			P.KeyWorkerProphTimeStartBase = 1e10;
		}

		//Added this to parameter list so that recording infection events (and the number to record) can easily be turned off and on: ggilani - 10/10/2014
		P.DoRecordInfEvents = Params::get_int(params, pre_params, "Record infection events", 0, &P);
		if (P.DoRecordInfEvents != 0)
		{
			P.MaxInfEvents = Params::get_int(params, pre_params, "Max number of infection events to record", 1000, &P);
			P.RecordInfEventsPerRun = Params::get_int(params, pre_params, "Record infection events per run", 0, &P);
		}
		else
		{
			P.MaxInfEvents = 0;
		}
		//Include a limit to the number of infections to simulate, if this happens before time runs out
		P.LimitNumInfections = Params::get_int(params, pre_params, "Limit number of infections", 0, &P);
		if (P.LimitNumInfections != 0)
		{
			P.MaxNumInfections = Params::get_int(params, pre_params, "Max number of infections", 60000, &P);
		}
		//Add origin-destination matrix parameter
		P.DoOriginDestinationMatrix = Params::get_int(params, pre_params, "Output origin destination matrix", 0, &P);

		P.MeanChildAgeGap = Params::get_int(params, pre_params, adm_params, "Mean child age gap", 2, &P);
		P.MinAdultAge = Params::get_int(params, pre_params, adm_params, "Min adult age", 19, &P);
		P.MaxMFPartnerAgeGap = Params::get_int(params, pre_params, adm_params, "Max MF partner age gap", 5, &P);
		P.MaxFMPartnerAgeGap = Params::get_int(params, pre_params, adm_params, "Max FM partner age gap", 5, &P);
		P.MinParentAgeGap = Params::get_int(params, pre_params, adm_params, "Min parent age gap", 19, &P);
		P.MaxParentAgeGap = Params::get_int(params, pre_params, adm_params, "Max parent age gap", 44, &P);
		P.MaxChildAge = Params::get_int(params, pre_params, adm_params, "Max child age", 20, &P);
		P.OneChildTwoPersProb = Params::get_double(params, pre_params, adm_params, "One Child Two Pers Prob", 0.08, &P);
		P.TwoChildThreePersProb = Params::get_double(params, pre_params, adm_params, "Two Child Three Pers Prob", 0.11, &P);
		P.OnePersHouseProbOld = Params::get_double(params, pre_params, adm_params, "One Pers House Prob Old", 0.5, &P);
		P.TwoPersHouseProbOld = Params::get_double(params, pre_params, adm_params, "Two Pers House Prob Old", 0.5, &P);
		P.OnePersHouseProbYoung = Params::get_double(params, pre_params, adm_params, "One Pers House Prob Young", 0.23, &P);
		P.TwoPersHouseProbYoung = Params::get_double(params, pre_params, adm_params, "Two Pers House Prob Young", 0.23, &P);
		P.OneChildProbYoungestChildUnderFive = Params::get_double(params, pre_params, adm_params, "One Child Prob Youngest Child Under Five", 0.5, &P);
		P.TwoChildrenProbYoungestUnderFive = Params::get_double(params, pre_params, adm_params, "Two Children Prob Youngest Under Five", 0.0, &P);
		P.TwoChildrenProbYoungestUnderFive = Params::get_double(params, pre_params,adm_params, "Prob Youngest Child Under Five", 0, &P);
		P.ZeroChildThreePersProb = Params::get_double(params, pre_params, adm_params, "Zero Child Three Pers Prob",0.25, &P);
		P.OneChildFourPersProb = Params::get_double(params, pre_params, adm_params, "One Child Four Pers Prob", 0.2, &P);
		P.YoungAndSingleSlope = Params::get_double(params, pre_params, adm_params, "Young And Single Slope", 0.7, &P);
		P.YoungAndSingle = Params::get_int(params, pre_params, adm_params, "Young And Single", 36, &P);
		P.NoChildPersAge = Params::get_int(params, pre_params, adm_params, "No Child Pers Age", 44, &P);
		P.OldPersAge = Params::get_int(params, pre_params, adm_params, "Old Pers Age", 60, &P);
		P.ThreeChildFivePersProb = Params::get_double(params, pre_params, adm_params, "Three Child Five Pers Prob", 0.5, &P);
		P.OlderGenGap = Params::get_int(params, pre_params, adm_params, "Older Gen Gap", 19, &P);
	}

	Files::xfprintf_stderr("Reading Params 13\n"); fflush(stderr);


	if (P.DoOneGen != 0) P.DoOneGen = 1;
	P.ColourPeriod = 2000;
	P.MoveRestrRadius2 = P.MoveRestrRadius * P.MoveRestrRadius;
	P.SocDistRadius2 = P.SocDistRadius * P.SocDistRadius;
	P.VaccRadius2 = P.VaccRadius * P.VaccRadius;
	P.VaccMinRadius2 = P.VaccMinRadius * P.VaccMinRadius;
	P.TreatRadius2 = P.TreatRadius * P.TreatRadius;
	P.PlaceCloseRadius2 = P.PlaceCloseRadius * P.PlaceCloseRadius;
	P.KeyWorkerProphRadius2 = P.KeyWorkerProphRadius * P.KeyWorkerProphRadius;
	if (P.TreatRadius2 == 0) P.TreatRadius2 = -1;
	if (P.VaccRadius2 == 0) P.VaccRadius2 = -1;
	if (P.PlaceCloseRadius2 == 0) P.PlaceCloseRadius2 = -1;
	if (P.MoveRestrRadius2 == 0) P.MoveRestrRadius2 = -1;
	if (P.SocDistRadius2 == 0) P.SocDistRadius2 = -1;
	if (P.KeyWorkerProphRadius2 == 0) P.KeyWorkerProphRadius2 = -1;
/*	if (P.TreatCellIncThresh < 1) P.TreatCellIncThresh = 1;
	if (P.CaseIsolation_CellIncThresh < 1) P.CaseIsolation_CellIncThresh = 1;
	if (P.DigitalContactTracing_CellIncThresh < 1) P.DigitalContactTracing_CellIncThresh = 1;
	if (P.HHQuar_CellIncThresh < 1) P.HHQuar_CellIncThresh = 1;
	if (P.MoveRestrCellIncThresh < 1) P.MoveRestrCellIncThresh = 1;
	if (P.PlaceCloseCellIncThresh < 1) P.PlaceCloseCellIncThresh = 1;
	if (P.KeyWorkerProphCellIncThresh < 1) P.KeyWorkerProphCellIncThresh = 1;
*/

	//// Make unsigned short versions of various intervention variables. And scaled them by number of timesteps per day
	P.usHQuarantineHouseDuration = ((unsigned short int) (P.HQuarantineHouseDuration * P.TimeStepsPerDay));
	P.usVaccTimeToEfficacy = ((unsigned short int) (P.VaccTimeToEfficacy * P.TimeStepsPerDay));
	P.usVaccTimeEfficacySwitch = ((unsigned short int) (P.VaccTimeEfficacySwitch * P.TimeStepsPerDay));
	P.usCaseIsolationDelay = ((unsigned short int) (P.CaseIsolationDelay * P.TimeStepsPerDay));
	P.usCaseIsolationDuration = ((unsigned short int) (P.CaseIsolationDuration * P.TimeStepsPerDay));
	P.usCaseAbsenteeismDuration = ((unsigned short int) (P.CaseAbsenteeismDuration * P.TimeStepsPerDay));
	P.usCaseAbsenteeismDelay = ((unsigned short int) (P.CaseAbsenteeismDelay * P.TimeStepsPerDay));
	if (P.DoUTM_coords)
	{
		for (i = 0; i <= 1000; i++)
		{
			asin2sqx[i] = asin(sqrt(i / 1000.0));
			asin2sqx[i] = asin2sqx[i] * asin2sqx[i];
		}
		for (i = 0; i <= DEGREES_PER_TURN; i++)
		{
			t = PI * i / 180;
			sinx[i] = sin(t);
			cosx[i] = cos(t);
		}
	}
	if (P.R0scale != 1.0)
	{
		P.HouseholdTrans *= P.R0scale;
		if (P.FixLocalBeta) P.LocalBeta *= P.R0scale;
		P.R0 *= P.R0scale;
		for (int place_type = 0; place_type < P.PlaceTypeNum; place_type++) {
			P.PlaceTypeTrans[place_type] *= P.R0scale;
		}
		Files::xfprintf_stderr("Rescaled transmission coefficients by factor of %lg\n", P.R0scale);
	}
	Files::xfprintf_stderr("Parameters read\n");	fflush(stderr);

	delete [] CountryNameBuf;
	delete [] AdunitListNamesBuf;
	delete [] CountryNames;
	delete [] AdunitListNames;
}

void ReadInterventions(std::string const& IntFile)
{
	FILE* dat;
	double r, s, startt, stopt;
	int j, k, au, ni, f, nsr;
	char buf[65536], txt[65536];
	Intervention CurInterv;

	Files::xfprintf_stderr("Reading intervention file.\n");
	dat = Files::xfopen(IntFile.c_str(), "rb");

	Files::xfscanf(dat, 0, "%*[^<]");
	Files::xfscanf(dat, 1, "<%[^>]", txt);

	if (strcmp(txt, "\?xml version=\"1.0\" encoding=\"ISO-8859-1\"\?") != 0) ERR_CRITICAL("Intervention file not XML.\n");

	Files::xfscanf(dat, 1, "%*[^<]<%[^>]", txt);
	if (strcmp(txt, "InterventionSettings") != 0) ERR_CRITICAL("Intervention has no top level.\n");
	ni = 0;
	while (!feof(dat))
	{
		Files::xfscanf(dat, 1, "%*[^<]<%[^>]", txt);
		if (strcmp(txt, "intervention") == 0)
		{
			ni++;
			Files::xfscanf(dat, 1, "%*[^<]<%[^>]", txt);
			if (strcmp(txt, "parameters") != 0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if (!GetXMLNode(dat, "Type", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if (strcmp(txt, "Treatment") == 0)
				CurInterv.InterventionType = 0;
			else if (strcmp(txt, "Vaccination") == 0)
				CurInterv.InterventionType = 1;
			else if (strcmp(txt, "ITN") == 0)
				CurInterv.InterventionType = 2;
			else if (strcmp(txt, "IRS") == 0)
				CurInterv.InterventionType = 3;
			else if (strcmp(txt, "GM") == 0)
				CurInterv.InterventionType = 4;
			else if (strcmp(txt, "MSAT") == 0)
				CurInterv.InterventionType = 5;
			else
				Files::xsscanf(txt, 1, "%i", &CurInterv.InterventionType);
			if (!GetXMLNode(dat, "AUThresh", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%i", &CurInterv.DoAUThresh);
			if (!GetXMLNode(dat, "StartTime", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.StartTime);
			startt = CurInterv.StartTime;
			if (!GetXMLNode(dat, "StopTime", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
		  Files::xsscanf(txt, 1, "%lf", &CurInterv.StopTime);
			stopt = CurInterv.StopTime;
			if (!GetXMLNode(dat, "MinDuration", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.MinDuration);
			CurInterv.MinDuration *= DAYS_PER_YEAR;
			if (!GetXMLNode(dat, "RepeatInterval", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.RepeatInterval);
			CurInterv.RepeatInterval *= DAYS_PER_YEAR;
			if (!GetXMLNode(dat, "MaxPrevAtStart", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.StartThresholdHigh);
			if (!GetXMLNode(dat, "MinPrevAtStart", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.StartThresholdLow);
			if (!GetXMLNode(dat, "MaxPrevAtStop", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.StopThreshold);
			if (GetXMLNode(dat, "NoStartAfterMinDur", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%i", &CurInterv.NoStartAfterMin);
			else
				CurInterv.NoStartAfterMin = 0;
			if (!GetXMLNode(dat, "Level", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%lf", &CurInterv.Level);
			if (GetXMLNode(dat, "LevelCellVar", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%lf", &CurInterv.LevelCellVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelAUVar", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%lf", &CurInterv.LevelAUVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelCountryVar", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%lf", &CurInterv.LevelCountryVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelClustering", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%lf", &CurInterv.LevelClustering);
			else
				CurInterv.LevelClustering = 0;
			if (GetXMLNode(dat, "ControlParam", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%lf", &CurInterv.ControlParam);
			else
				CurInterv.ControlParam = 0;
			if (GetXMLNode(dat, "TimeOffset", "parameters", txt, 1))
				Files::xsscanf(txt, 1, "%lf", &CurInterv.TimeOffset);
			else
				CurInterv.TimeOffset = 0;

			if (!GetXMLNode(dat, "MaxRounds", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%u", &CurInterv.MaxRounds);
			if (!GetXMLNode(dat, "MaxResource", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xsscanf(txt, 1, "%u", &CurInterv.MaxResource);
			if (GetXMLNode(dat, "NumSequentialReplicas", "parameters", txt, 1))
			Files::xsscanf(txt, 1, "%i", &nsr);
			else
				nsr = 0;
			do {
			    Files::xfscanf(dat, 1, "%*[^<]<%[^>]", txt);
      } while ((strcmp(txt, "/intervention") != 0) && (strcmp(txt, "/parameters") != 0) && (!feof(dat)));
			if (strcmp(txt, "/parameters") != 0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			Files::xfscanf(dat, 1, "%*[^<]<%[^>]", txt);
			if ((strcmp(txt, "adunits") != 0) && (strcmp(txt, "countries") != 0)) ERR_CRITICAL("Incomplete adunits/countries specification in intervention file\n");
			if (strcmp(txt, "adunits") == 0)
			{
				while (GetXMLNode(dat, "A", "adunits", buf, 0))
				{
					Files::xsscanf(buf, 1, "%s", txt);
					j = atoi(txt);
					if (j == 0)
					{
						f = 1; au = -1;
						do
						{
							au++; f = strcmp(txt, AdUnits[au].ad_name);
						} while ((f) && (au < P.NumAdunits));
						if (!f)
						{
							r = fabs(CurInterv.Level) + (2.0 * ranf() - 1) * CurInterv.LevelAUVar;
							if ((CurInterv.Level < 1) && (r > 1))
								r = 1;
							else if (r < 0)
								r = 0;
							for (k = 0; k <= nsr; k++)
							{
								AdUnits[au].InterventionList[AdUnits[au].NI] = CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level = r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime = startt + ((double)k) * (stopt - startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime = stopt + ((double)k) * (stopt - startt);
								AdUnits[au].NI++;
							}
						}
					}
					else
					{
						k = (j % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor;
						au = P.AdunitLevel1Lookup[k];
						if ((au >= 0) && (AdUnits[au].id / P.AdunitLevel1Divisor == j / P.AdunitLevel1Divisor))
						{
							r = CurInterv.Level + (2.0 * ranf() - 1) * CurInterv.LevelAUVar;
							if ((CurInterv.Level < 1) && (r > 1))
								r = 1;
							else if (r < 0)
								r = 0;
							for (k = 0; k <= nsr; k++)
							{
								AdUnits[au].InterventionList[AdUnits[au].NI] = CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level = r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime = startt + ((double)k) * (stopt - startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime = stopt + ((double)k) * (stopt - startt);
								AdUnits[au].NI++;
							}
						}
					}
				}
			}
			else
			{
				while (GetXMLNode(dat, "C", "countries", buf, 0))
				{
					s = (2.0 * ranf() - 1) * CurInterv.LevelCountryVar;
					Files::xsscanf(buf, 1, "%s", txt);
					j = atoi(txt);
					for (au = 0; au < P.NumAdunits; au++)
						if (((j == 0) && (strcmp(txt, AdUnits[au].cnt_name) == 0)) || ((j > 0) && (j == AdUnits[au].cnt_id)))
						{
							r = CurInterv.Level + (2.0 * ranf() - 1) * CurInterv.LevelAUVar + s;
							if ((CurInterv.Level < 1) && (r > 1))
								r = 1;
							else if (r < 0)
								r = 0;
							for (k = 0; k <= nsr; k++)
							{
								AdUnits[au].InterventionList[AdUnits[au].NI] = CurInterv;
								AdUnits[au].InterventionList[AdUnits[au].NI].Level = r;
								AdUnits[au].InterventionList[AdUnits[au].NI].StartTime = startt + ((double)k) * (stopt - startt);
								AdUnits[au].InterventionList[AdUnits[au].NI].StopTime = stopt + ((double)k) * (stopt - startt);
								AdUnits[au].NI++;
							}
						}
				}
			}
			Files::xfscanf(dat, 1, "%*[^<]<%[^>]", txt);
			if (strcmp(txt, "/intervention") != 0) ERR_CRITICAL("Incorrect intervention specification in intervention file\n");
		}
	}
	if (strcmp(txt, "/InterventionSettings") != 0) ERR_CRITICAL("Intervention has no top level closure.\n");
	Files::xfprintf_stderr("%i interventions read\n", ni);
	Files::xfclose(dat);
}
int GetXMLNode(FILE* dat, const char* NodeName, const char* ParentName, char* Value, int ResetFilePos)
{
	// ResetFilePos=1 leaves dat cursor in same position as when function was called. 0 leaves it at end of NodeName closure
	// GetXMLNode returns 1 if NodeName found, 0 otherwise. If NodeName not found, ParentName closure must be

	char buf[65536], CloseNode[2048], CloseParent[2048];
	int CurPos, ret;

	Files::xsprintf(CloseParent, "/%s", ParentName);
	CurPos = ftell(dat);
	do
	{
		Files::xfscanf(dat, 1, "%*[^<]<%[^>]", buf);
	} while ((strcmp(buf, CloseParent) != 0) && (strcmp(buf, NodeName) != 0) && (!feof(dat)));
	if (strcmp(buf, CloseParent) == 0)
		ret = 0;
	else
	{
		if (strcmp(buf, NodeName) != 0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		Files::xfscanf(dat, 1, ">%[^<]", buf);
		if (strlen(buf) < 2048) strcpy(Value, buf);
		//		Files::xfprintf_stderr("# %s=%s\n",NodeName,Value);
		Files::xfscanf(dat, 1, "<%[^>]", buf);
		Files::xsprintf(CloseNode, "/%s", NodeName);
		if (strcmp(buf, CloseNode) != 0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		ret = 1;
	}
	if (ResetFilePos) fseek(dat, CurPos, 0);
	return ret;
}

void ReadAirTravel(std::string const& air_travel_file, std::string const& output_file_base)
{
	int i, j, k, l;
	float sc, t, t2;
	float* buf;
	double traf;
	std::string outname;
	FILE* dat;

	Files::xfprintf_stderr("Reading airport data...\nAirports with no connections = ");
	dat = Files::xfopen(air_travel_file.c_str(), "rb");
	Files::xfscanf(dat, 2, "%i %i", &P.Nairports, &P.Air_popscale);
	sc = (float)((double)P.PopSize / (double)P.Air_popscale);
	if (P.Nairports > MAX_AIRPORTS) ERR_CRITICAL("Too many airports\n");
	if (P.Nairports < 2) ERR_CRITICAL("Too few airports\n");
	buf = (float*)Memory::xcalloc(_I64(P.Nairports) + 1, sizeof(float));
	Airports = (Airport*)Memory::xcalloc(P.Nairports, sizeof(Airport));
	for (i = 0; i < P.Nairports; i++)
	{
		Files::xfscanf(dat, 3, "%f %f %lf", &(Airports[i].loc.x), &(Airports[i].loc.y), &traf);
		traf *= (P.AirportTrafficScale * sc);
		if (!P.SpatialBoundingBox.inside(CovidSim::Geometry::Vector2d(Airports[i].loc)))
		{
			Airports[i].loc.x = Airports[i].loc.y = -1;
			Airports[i].total_traffic = 0;
		}
		else
		{
			Airports[i].loc.x -= (float)P.SpatialBoundingBox.bottom_left().x;
			Airports[i].loc.y -= (float)P.SpatialBoundingBox.bottom_left().y;
			Airports[i].total_traffic = (float)traf;
		}
		t = 0;
		for (j = k = 0; j < P.Nairports; j++)
		{
			Files::xfscanf(dat, 1, "%f", buf + j);
			if (buf[j] > 0) { k++; t += buf[j]; }
		}
		Airports[i].num_connected = k;
		if (Airports[i].num_connected > 0)
		{
			Airports[i].prop_traffic = (float*)Memory::xcalloc(Airports[i].num_connected, sizeof(float));
			Airports[i].conn_airports = (unsigned short int*)Memory::xcalloc(Airports[i].num_connected, sizeof(unsigned short int));
			for (j = k = 0; j < P.Nairports; j++)
				if (buf[j] > 0)
				{
					Airports[i].conn_airports[k] = j;
					Airports[i].prop_traffic[k] = buf[j] / t;
					k++;
				}
		}
		else
		{
			if (Airports[i].total_traffic > 0)
				Files::xfprintf_stderr("#%i# ", i);
			else
				Files::xfprintf_stderr("%i ", i);
		}
	}
	Files::xfclose(dat);
	Memory::xfree(buf);
	Files::xfprintf_stderr("\nAirport data read OK.\n");
	for (i = 0; i < P.Nairports; i++)
	{
		/*		Files::xfprintf_stderr("(%f %i|",Airports[i].total_traffic,Airports[i].num_connected);
		*/		t = 0; k = 0;
	for (j = Airports[i].num_connected - 1; j >= 0; j--)
	{
		if ((Airports[i].prop_traffic[j] > 0) && (Airports[Airports[i].conn_airports[j]].total_traffic == 0))
		{
			t += Airports[i].prop_traffic[j];
			Airports[i].num_connected--;
			if (j < Airports[i].num_connected)
			{
				Airports[i].prop_traffic[j] = Airports[i].prop_traffic[Airports[i].num_connected];
				Airports[i].conn_airports[j] = Airports[i].conn_airports[Airports[i].num_connected];
			}
			Airports[i].prop_traffic[Airports[i].num_connected] = 0;
			Airports[i].conn_airports[Airports[i].num_connected] = 0;
		}
		else if (Airports[i].prop_traffic[j] > 0)
			k = 1;
	}
	/*		Files::xfprintf_stderr("%f %i ",t,k);
	*/		t = 1.0f - t;
	if (k)
	{
		Airports[i].total_traffic *= t;
		t2 = 0;
		for (j = 0; j < Airports[i].num_connected; j++)
		{
			Airports[i].prop_traffic[j] = t2 + Airports[i].prop_traffic[j];
			t2 = Airports[i].prop_traffic[j];
		}
		for (j = 0; j < Airports[i].num_connected; j++)
			Airports[i].prop_traffic[j] /= t2;
		/*			if((Airports[i].num_connected>0)&&(Airports[i].prop_traffic[Airports[i].num_connected-1]!=1))
						Files::xfprintf_stderr("<%f> ",Airports[i].prop_traffic[Airports[i].num_connected-1]);
		*/
	}
	else
	{
		Airports[i].total_traffic = 0; Airports[i].num_connected = 0;
	}
	if (Airports[i].num_connected > 0)
	{
		for (j = k = 0; k < 128; k++)
		{
			t = (float)((double)k / 128);
			while (Airports[i].prop_traffic[j] < t) j++;
			Airports[i].Inv_prop_traffic[k] = j;
		}
		Airports[i].Inv_prop_traffic[128] = Airports[i].num_connected - 1;
	}
	/*		Files::xfprintf_stderr("%f) ",Airports[i].total_traffic);
	*/
	}
	Files::xfprintf_stderr("Airport data clipped OK.\n");
	for (i = 0; i < MAX_DIST; i++) AirTravelDist[i] = 0;
	for (i = 0; i < P.Nairports; i++)
		if (Airports[i].total_traffic > 0)
		{
			for (j = 0; j < Airports[i].num_connected; j++)
			{
				k = (int)Airports[i].conn_airports[j];
				traf = floor(sqrt(dist2_raw(Airports[i].loc.x, Airports[i].loc.y, Airports[k].loc.x, Airports[k].loc.y)) / OUTPUT_DIST_SCALE);
				l = (int)traf;
				//Files::xfprintf_stderr("%(%i) ",l);
				if (l < MAX_DIST)
					AirTravelDist[l] += (double) Airports[i].total_traffic * Airports[i].prop_traffic[j];
			}
		}
	outname = output_file_base + ".airdist.xls";
	dat = Files::xfopen(outname.c_str(), "wb");
	Files::xfprintf(dat, "dist\tfreq\n");
	for (i = 0; i < MAX_DIST; i++)
		Files::xfprintf(dat, "%i\t%.10f\n", i, AirTravelDist[i]);
	Files::xfclose(dat);
}

void InitModel(int run) // passing run number so we can save run number in the infection event log: ggilani - 15/10/2014
{
	int nim;
	int NumSeedingInfections_byLocation[MAX_NUM_SEED_LOCATIONS];

	if (P.OutputBitmap)
	{
#ifdef _WIN32
		//if (P.OutputBitmap == 1)
		//{
		//	char buf[200];
		//	Files::xsprintf(buf, "%s.ge" DIRECTORY_SEPARATOR "%s.avi", OutFile, OutFile);
		//	avi = CreateAvi(buf, P.BitmapMovieFrame, NULL);
		//}
#endif
		for (unsigned p = 0; p < bmh->imagesize; p++)
		{
			bmInfected[p] = bmRecovered[p] = bmTreated[p] = 0;
		}
	}
	OutputTimeStepNumber = 0;
	State.S = P.PopSize;
	State.L = State.I = State.R = State.D = 0;
	State.cumI = State.cumR = State.cumC = State.cumFC = State.cumCT = State.cumCC = State.cumTC = State.cumD = State.cumDC = State.trigDetectedCases = State.DCT = State.cumDCT
		= State.cumTG = State.cumSI = State.nTG = State.cumHQ = State.cumAC = State.cumAH = State.cumAA = State.cumACS = State.cumAPC = State.cumAPA = State.cumAPCS = 0;
	State.cumT = State.cumUT = State.cumTP = State.cumV = State.sumRad2 = State.maxRad2 = State.cumV_daily = State.cumVG = 0; //added State.cumVG
	State.mvacc_cum = 0;
	if (P.DoSeverity)
	{
		State.Mild		= State.ILI			= State.SARI	= State.Critical	= State.CritRecov		= 0;
		State.cumMild	= State.cumILI		= State.cumSARI = State.cumCritical = State.cumCritRecov	= 0;
		State.cumDeath_ILI = State.cumDeath_SARI = State.cumDeath_Critical = 0;

		if (P.DoAdUnits)
			for (int AdminUnit = 0; AdminUnit <= P.NumAdunits; AdminUnit++)
			{
				State.Mild_adunit[AdminUnit] = State.ILI_adunit[AdminUnit] =
				State.SARI_adunit[AdminUnit] = State.Critical_adunit[AdminUnit] = State.CritRecov_adunit[AdminUnit] =
				State.cumMild_adunit[AdminUnit] = State.cumILI_adunit[AdminUnit] =
				State.cumSARI_adunit[AdminUnit] = State.cumCritical_adunit[AdminUnit] = State.cumCritRecov_adunit[AdminUnit] =
				State.cumDeath_ILI_adunit[AdminUnit] = State.cumDeath_SARI_adunit[AdminUnit] = State.cumDeath_Critical_adunit[AdminUnit] =
				State.cumD_adunit[AdminUnit] = 0;
			}
		for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
		{
			State.Mild_age[AgeGroup] = State.ILI_age[AgeGroup] =
				State.SARI_age[AgeGroup] = State.Critical_age[AgeGroup] = State.CritRecov_age[AgeGroup] =
				State.cumMild_age[AgeGroup] = State.cumILI_age[AgeGroup] =
				State.cumSARI_age[AgeGroup] = State.cumCritical_age[AgeGroup] = State.cumCritRecov_age[AgeGroup] =
				State.cumDeath_ILI_age[AgeGroup] = State.cumDeath_SARI_age[AgeGroup] = State.cumDeath_Critical_age[AgeGroup] = 0;
		}
	}
	if (P.DoAdUnits && P.OutputAdUnitAge)
		for (int Adunit = 0; Adunit < P.NumAdunits; Adunit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
			{
				State.prevInf_age_adunit[AgeGroup][Adunit] = 0;
				State.cumInf_age_adunit	[AgeGroup][Adunit] = 0;
			}

	for (int i = 0; i < NUM_AGE_GROUPS; i++) State.cumCa[i] = State.cumIa[i] = State.cumDa[i] = 0;
	for (int i = 0; i < 2; i++) State.cumC_keyworker[i] = State.cumI_keyworker[i] = State.cumT_keyworker[i] = 0;
	for (int i = 0; i < NUM_PLACE_TYPES; i++) State.NumPlacesClosed[i] = 0;
	for (int i = 0; i < INFECT_TYPE_MASK; i++) State.cumItype[i] = 0;
	//initialise cumulative case counts per country to zero: ggilani 12/11/14
	for (int i = 0; i < MAX_COUNTRIES; i++) State.cumC_country[i] = 0;
	if (P.DoAdUnits)
		for (int i = 0; i <= P.NumAdunits; i++)
		{
			State.cumI_adunit[i] = State.cumC_adunit[i] = State.cumD_adunit[i] = State.cumT_adunit[i] = State.cumH_adunit[i] =
				State.cumDC_adunit[i] = State.cumCT_adunit[i] = State.cumCC_adunit[i] = State.trigDC_adunit[i] = State.DCT_adunit[i] = State.cumDCT_adunit[i] = 0; //added hospitalisation, added detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
			AdUnits[i].place_close_trig = 0;
			AdUnits[i].CaseIsolationTimeStart = AdUnits[i].HQuarantineTimeStart = AdUnits[i].DigitalContactTracingTimeStart = AdUnits[i].SocialDistanceTimeStart = AdUnits[i].PlaceCloseTimeStart = 1e10;
			AdUnits[i].ndct = 0; //noone being digitally contact traced at beginning of run
		}

	//update state variables for storing contact distribution
	for (int i = 0; i < MAX_CONTACTS+1; i++) State.contact_dist[i] = 0;

	for (int j = 0; j < MAX_NUM_THREADS; j++)
	{
		StateT[j].L = StateT[j].I = StateT[j].R = StateT[j].D = 0;
		StateT[j].cumI = StateT[j].cumR = StateT[j].cumC = StateT[j].cumFC = StateT[j].cumCT = StateT[j].cumCC = StateT[j].DCT = StateT[j].cumDCT
			= StateT[j].cumTG = StateT[j].cumSI = StateT[j].nTG = StateT[j].cumTC = StateT[j].cumD = StateT[j].cumDC = StateT[j].cumHQ = StateT[j].cumAC
			= StateT[j].cumACS = StateT[j].cumAH = StateT[j].cumAA = StateT[j].cumAPC = StateT[j].cumAPA = StateT[j].cumAPCS = 0;
		StateT[j].cumT = StateT[j].cumUT = StateT[j].cumTP = StateT[j].cumV = StateT[j].sumRad2 = StateT[j].maxRad2 = StateT[j].cumV_daily =  0;
		for (int i = 0; i < NUM_AGE_GROUPS; i++) StateT[j].cumCa[i] = StateT[j].cumIa[i] = StateT[j].cumDa[i] = 0;
		for (int i = 0; i < 2; i++) StateT[j].cumC_keyworker[i] = StateT[j].cumI_keyworker[i] = StateT[j].cumT_keyworker[i] = 0;
		for (int i = 0; i < NUM_PLACE_TYPES; i++) StateT[j].NumPlacesClosed[i] = 0;
		for (int i = 0; i < INFECT_TYPE_MASK; i++) StateT[j].cumItype[i] = 0;
		//initialise cumulative case counts per country per thread to zero: ggilani 12/11/14
		for (int i = 0; i < MAX_COUNTRIES; i++) StateT[j].cumC_country[i] = 0;
		if (P.DoAdUnits)
			for (int i = 0; i <= P.NumAdunits; i++)
				StateT[j].cumI_adunit[i] = StateT[j].cumC_adunit[i] = StateT[j].cumD_adunit[i] = StateT[j].cumT_adunit[i] = StateT[j].cumH_adunit[i] = StateT[j].cumDC_adunit[i] =
				StateT[j].cumCT_adunit[i] = StateT[j].cumCC_adunit[i] = StateT[j].nct_queue[i] = StateT[j].cumDCT_adunit[i] = StateT[j].DCT_adunit[i] = StateT[j].ndct_queue[i] = 0; //added hospitalisation, detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17

		if (P.DoSeverity)
		{
			StateT[j].Mild		= StateT[j].ILI		= StateT[j].SARI	= StateT[j].Critical	= StateT[j].CritRecov		= 0;
			StateT[j].cumMild	= StateT[j].cumILI	= StateT[j].cumSARI = StateT[j].cumCritical = StateT[j].cumCritRecov	= 0;
			StateT[j].cumDeath_ILI = StateT[j].cumDeath_SARI = StateT[j].cumDeath_Critical = 0;

			for (int AdminUnit = 0; AdminUnit <= P.NumAdunits; AdminUnit++)
			{
				StateT[j].Mild_adunit[AdminUnit] = StateT[j].ILI_adunit[AdminUnit] =
				StateT[j].SARI_adunit[AdminUnit] = StateT[j].Critical_adunit[AdminUnit] = StateT[j].CritRecov_adunit[AdminUnit] =
				StateT[j].cumMild_adunit[AdminUnit] = StateT[j].cumILI_adunit[AdminUnit] =
				StateT[j].cumSARI_adunit[AdminUnit] = StateT[j].cumCritical_adunit[AdminUnit] = StateT[j].cumCritRecov_adunit[AdminUnit] =
				StateT[j].cumDeath_ILI_adunit[AdminUnit] = StateT[j].cumDeath_SARI_adunit[AdminUnit] = StateT[j].cumDeath_Critical_adunit[AdminUnit] =
				StateT[j].cumD_adunit[AdminUnit] =  0;
			}
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
			{
				StateT[j].Mild_age[AgeGroup] = StateT[j].ILI_age[AgeGroup] =
					StateT[j].SARI_age[AgeGroup] = StateT[j].Critical_age[AgeGroup] = StateT[j].CritRecov_age[AgeGroup] =
					StateT[j].cumMild_age[AgeGroup] = StateT[j].cumILI_age[AgeGroup] =
					StateT[j].cumSARI_age[AgeGroup] = StateT[j].cumCritical_age[AgeGroup] = StateT[j].cumCritRecov_age[AgeGroup] =
					StateT[j].cumDeath_ILI_age[AgeGroup] = StateT[j].cumDeath_SARI_age[AgeGroup] = StateT[j].cumDeath_Critical_age[AgeGroup] = 0;
			}

		}
		//resetting thread specific parameters for storing contact distribution
		for (int i = 0; i < MAX_CONTACTS+1; i++) StateT[j].contact_dist[i] = 0;

	}
	nim = 0;

	if (P.DoAdUnits && P.OutputAdUnitAge)
		for (int Thread = 0; Thread < P.NumThreads; Thread++)
			for (int Adunit = 0; Adunit < P.NumAdunits; Adunit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				{
					StateT[Thread].prevInf_age_adunit[AgeGroup][Adunit] = 0;
					StateT[Thread].cumInf_age_adunit [AgeGroup][Adunit] = 0;
				}

	std::fill(HostsQuarantine.begin(), HostsQuarantine.end(), PersonQuarantine());
#pragma omp parallel for schedule(static,1) default(none) \
		shared(P, Hosts)
	for (int tn = 0; tn < P.NumThreads; tn++)
		for (int k = tn; k < P.PopSize; k+= P.NumThreads)
		{
			Hosts[k].absent_start_time = USHRT_MAX - 1;
			Hosts[k].absent_stop_time = 0;
			if (P.DoAirports) Hosts[k].PlaceLinks[P.HotelPlaceType] = -1;
			Hosts[k].vacc_start_time = Hosts[k].treat_start_time = Hosts[k].isolation_start_time = Hosts[k].absent_start_time = Hosts[k].dct_start_time = Hosts[k].dct_trigger_time = USHRT_MAX - 1;
			Hosts[k].treat_stop_time = Hosts[k].absent_stop_time = Hosts[k].dct_end_time = 0;
			Hosts[k].to_die = 0;
			Hosts[k].Travelling = 0;
			Hosts[k].detected = 0; //set detected to zero initially: ggilani - 19/02/15
			Hosts[k].detected_time = 0;
			Hosts[k].digitalContactTraced = 0;
			Hosts[k].set_susceptible();
			Hosts[k].num_treats = 0;
			Hosts[k].latent_time = Hosts[k].recovery_or_death_time = 0; //also set hospitalisation time to zero: ggilani 28/10/2014
			Hosts[k].infector = -1;
			Hosts[k].infect_type = 0;
			Hosts[k].index_case_dct = 0;
			Hosts[k].ProbAbsent =(float) ranf_mt(tn);
			Hosts[k].ProbCare = (float) ranf_mt(tn);
			Hosts[k].susc = (float)((P.DoPartialImmunity) ? (1.0 - P.InitialImmunity[HOST_AGE_GROUP(k)]) : 1.0);
			if(P.SusceptibilitySD > 0) Hosts[k].susc *= (float) gen_gamma_mt(1 / (P.SusceptibilitySD * P.SusceptibilitySD), 1 / (P.SusceptibilitySD * P.SusceptibilitySD), tn);
			if (P.DoSeverity)
			{
				Hosts[k].SARI_time = USHRT_MAX - 1; //// think better to set to initialize to maximum possible value, but keep this way for now.
				Hosts[k].Critical_time = USHRT_MAX - 1;
				Hosts[k].Stepdown_time = USHRT_MAX - 1;
				Hosts[k].Severity_Current = Severity::Asymptomatic;
				Hosts[k].Severity_Final = Severity::Asymptomatic;
				Hosts[k].set_susceptible();
			}
		}

#pragma omp parallel for reduction(+:nim) schedule(static,1) default(none) \
		shared(P, Cells, Hosts, Households)
	for (int tn = 0; tn < P.NumThreads; tn++)
	{
		for (int i = tn; i < P.NumCells; i+=P.NumThreads)
		{
			if ((Cells[i].tot_treat != 0) || (Cells[i].tot_vacc != 0) || (Cells[i].S != Cells[i].n) || (Cells[i].D > 0) || (Cells[i].R > 0))
			{
				for (int j = 0; j < Cells[i].n; j++)
				{
					int k = Cells[i].members[j];
					Cells[i].susceptible[j] = k; //added this in here instead
					Hosts[k].listpos = j;
				}
				Cells[i].S = Cells[i].n;
				Cells[i].L = Cells[i].I = Cells[i].R = Cells[i].cumTC = Cells[i].D = 0;
				Cells[i].infected = Cells[i].latent = Cells[i].susceptible + Cells[i].S;
				Cells[i].tot_treat = Cells[i].tot_vacc = 0;
				for (int l = 0; l < MAX_INTERVENTION_TYPES; l++) Cells[i].CurInterv[l] = -1;

				// Next loop needs to count down for DoImmune host list reordering to work
				if(!P.DoPartialImmunity)
					for (int j = Cells[i].n - 1; j >= 0; j--)
					{
						int k = Cells[i].members[j];
						if (P.DoWholeHouseholdImmunity)
						{
	// note that this breaks determinism of runs if executed due to reordering of Cell members list each realisation
							if (P.InitialImmunity[0] != 0)
							{
								if (Households[Hosts[k].hh].FirstPerson == k)
								{
									if ((P.InitialImmunity[0] == 1) || (ranf_mt(tn) < P.InitialImmunity[0]))
									{
										nim += Households[Hosts[k].hh].nh;
										for (int m = Households[Hosts[k].hh].nh - 1; m >= 0; m--)
											DoImmune(k + m);
									}
								}
							}
						}
						else
						{
							int m = HOST_AGE_GROUP(k);
							if ((P.InitialImmunity[m] == 1) || ((P.InitialImmunity[m] > 0) && (ranf_mt(tn) < P.InitialImmunity[m])))
							{
								DoImmune(k); nim += 1;
							}
						}
					}
			}
		}
	}

#pragma omp parallel for schedule(static,500) default(none) \
		shared(P, Mcells, McellLookup)
	for (int l = 0; l < P.NumPopulatedMicrocells; l++)
	{
		int i = (int)(McellLookup[l] - Mcells);
		Mcells[i].vacc_start_time = Mcells[i].treat_start_time = USHRT_MAX - 1;
		Mcells[i].treat_end_time = 0;
		Mcells[i].treat_trig = Mcells[i].vacc_trig = 0;
		Mcells[i].vacc = Mcells[i].treat = Mcells[i].placeclose = Mcells[i].socdist = Mcells[i].moverest = TreatStat::Untreated;
		Mcells[i].place_trig = Mcells[i].move_trig = Mcells[i].socdist_trig = Mcells[i].keyworkerproph_trig = Mcells[i].keyworkerproph = 0;
		Mcells[i].move_start_time = USHRT_MAX - 1;
		Mcells[i].place_end_time = Mcells[i].move_end_time =
			Mcells[i].socdist_end_time = Mcells[i].keyworkerproph_end_time = 0;
	}
	if (P.DoPlaces)
#pragma omp parallel for schedule(static,1) default(none) \
			shared(P, Places)
		for (int m = 0; m < P.PlaceTypeNum; m++)
		{
			for(int l = 0; l < P.Nplace[m]; l++)
			{
				Places[m][l].close_start_time = USHRT_MAX - 1;
				Places[m][l].treat = Places[m][l].control_trig = 0;
				Places[m][l].treat_end_time = Places[m][l].close_end_time = 0;
				Places[m][l].ProbClose = (float) ranf_mt(m);
				if (P.AbsenteeismPlaceClosure)
				{
					Places[m][l].AbsentLastUpdateTime = 0;
					for (int i2 = 0; i2 < P.MaxAbsentTime; i2++) Places[m][l].Absent[i2] = 0;
				}
			}
		}

	//// **** //// **** //// **** Initialize Current effects
	//// **** soc dist
	P.SocDistDurationCurrent			= P.SocDistDuration;
	P.SocDistSpatialEffectCurrent		= P.SD_SpatialEffects_OverTime	[0];				//// spatial
	P.SocDistHouseholdEffectCurrent		= P.SD_HouseholdEffects_OverTime[0];				//// household
	for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
		P.SocDistPlaceEffectCurrent[PlaceType] = P.SD_PlaceEffects_OverTime[0][PlaceType];	//// place
	P.SocDistCellIncThresh				= P.SD_CellIncThresh_OverTime	[0];				//// cell incidence threshold

	//// **** enhanced soc dist
	P.EnhancedSocDistSpatialEffectCurrent		= P.Enhanced_SD_SpatialEffects_OverTime		[0];	//// spatial
	P.EnhancedSocDistHouseholdEffectCurrent		= P.Enhanced_SD_HouseholdEffects_OverTime	[0];	//// household
	for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
		P.EnhancedSocDistPlaceEffectCurrent[PlaceType] = P.Enhanced_SD_PlaceEffects_OverTime[0][PlaceType];	//// place

	//// **** case isolation
	P.CaseIsolationEffectiveness		= P.CI_SpatialAndPlaceEffects_OverTime	[0];	//// spatial / place
	P.CaseIsolationHouseEffectiveness	= P.CI_HouseholdEffects_OverTime		[0];	//// household
	P.CaseIsolationProp					= P.CI_Prop_OverTime					[0];	//// compliance
	P.CaseIsolation_CellIncThresh		= P.CI_CellIncThresh_OverTime			[0];	//// cell incidence threshold


	//// **** household quarantine
	P.HQuarantineSpatialEffect	= P.HQ_SpatialEffects_OverTime	[0];	//// spatial
	P.HQuarantineHouseEffect	= P.HQ_HouseholdEffects_OverTime[0];	//// household
	for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
		P.HQuarantinePlaceEffect[PlaceType] = P.HQ_PlaceEffects_OverTime	[0][PlaceType];	//// place
	P.HQuarantinePropIndivCompliant = P.HQ_Individual_PropComply_OverTime	[0]; //// individual compliance
	P.HQuarantinePropHouseCompliant = P.HQ_Household_PropComply_OverTime	[0]; //// household compliance
	P.HHQuar_CellIncThresh			= P.HQ_CellIncThresh_OverTime			[0]; //// cell incidence threshold


	//// **** place closure
	P.PlaceCloseSpatialRelContact	= P.PC_SpatialEffects_OverTime	[0];			//// spatial
	P.PlaceCloseHouseholdRelContact = P.PC_HouseholdEffects_OverTime[0];			//// household
	for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
	{
		P.PlaceCloseEffect[PlaceType] = P.PC_PlaceEffects_OverTime[0][PlaceType];	//// place
		P.PlaceClosePropAttending[PlaceType] = P.PC_PropAttending_OverTime[0][PlaceType];
	}
	P.PlaceCloseIncTrig1			= P.PC_IncThresh_OverTime		[0];			//// global incidence threshold
	P.PlaceCloseFracIncTrig			= P.PC_FracIncThresh_OverTime	[0];			//// fractional incidence threshold
	P.PlaceCloseCellIncThresh1		= P.PC_CellIncThresh_OverTime	[0];			//// cell incidence threshold
	P.PlaceCloseDurationBase = P.PC_Durs_OverTime[0]; //// duration of place closure


	//// **** digital contact tracing
	P.DCTCaseIsolationEffectiveness			= P.DCT_SpatialAndPlaceEffects_OverTime	[0];	//// spatial / place
	P.DCTCaseIsolationHouseEffectiveness	= P.DCT_HouseholdEffects_OverTime		[0];	//// household
	P.ProportionDigitalContactsIsolate		= P.DCT_Prop_OverTime					[0];	//// compliance
	P.MaxDigitalContactsToTrace				= P.DCT_MaxToTrace_OverTime				[0];

	// Initialize CFR scalings
	P.CFR_Critical_Scale_Current	= P.CFR_TimeScaling_Critical[0];
	P.CFR_SARI_Scale_Current		= P.CFR_TimeScaling_SARI	[0];
	P.CFR_ILI_Scale_Current			= P.CFR_TimeScaling_ILI		[0];


	for (int i = 0; i < MAX_NUM_THREADS; i++)
	{
		for (int j = 0; j < MAX_NUM_THREADS; j++)	StateT[i].n_queue[j] = 0;
		for (int j = 0; j < P.PlaceTypeNum; j++)	StateT[i].np_queue[j] = 0;
		StateT[i].host_closure_queue_size = 0;
	}
	if (DoInitUpdateProbs)
	{
		UpdateProbs(0);
		DoInitUpdateProbs = 0;
	}
	//initialise event log to zero at the beginning of every run: ggilani - 10/10/2014. UPDATE: 15/10/14 - we are now going to store all events from all realisations in one file
	if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun))
	{
		nEvents = 0;
		for (int i = 0; i < P.MaxInfEvents; i++)
		{
			InfEventLog[i].t = InfEventLog[i].infectee_x = InfEventLog[i].infectee_y = InfEventLog[i].t_infector = 0.0;
			InfEventLog[i].infectee_ind = InfEventLog[i].infector_ind = 0;
			InfEventLog[i].infectee_adunit = InfEventLog[i].listpos = InfEventLog[i].infectee_cell = InfEventLog[i].infector_cell = InfEventLog[i].thread = 0;
		}
	}

	for (int i = 0; i < P.NumSeedLocations; i++) NumSeedingInfections_byLocation[i] = (int) (((double) P.NumInitialInfections[i]) * P.InitialInfectionsAdminUnitWeight[i]* P.SeedingScaling +0.5);
	SeedInfection(0, NumSeedingInfections_byLocation, 0, run);
	P.ControlPropCasesId = P.PreAlertControlPropCasesId;
	P.TreatTimeStart = 1e10;

	P.VaccTimeStart = 1e10;
	P.MoveRestrTimeStart = 1e10;
	P.PlaceCloseTimeStart = 1e10;
	P.PlaceCloseTimeStart2 = 2e10;
	P.SocDistTimeStart = 1e10;
	P.AirportCloseTimeStart = 1e10;
	//P.DigitalContactTracingTimeStart = 1e10;
	P.HQuarantineTimeStart = 1e10;
	P.KeyWorkerProphTimeStart = 1e10;
	P.TreatMaxCourses = P.TreatMaxCoursesBase;
	P.VaccMaxCourses = P.VaccMaxCoursesBase;
	P.PlaceCloseDuration = P.PlaceCloseDurationBase; //// duration of place closure
	P.PlaceCloseIncTrig = P.PlaceCloseIncTrig1;
	P.PlaceCloseTimeStartPrevious = 1e10;
	P.PlaceCloseCellIncThresh = P.PlaceCloseCellIncThresh1;
	P.ResetSeedsFlag = 0; //added this to allow resetting seeds part way through run: ggilani 27/11/2019
	if (!P.StopCalibration) P.DateTriggerReached_SimTime = 0;
	if (P.InitialInfectionCalTime > 0)
	{
		P.HolidaysStartDay_SimTime = -P.InitialInfectionCalTime;
		P.DateTriggerReached_SimTime = P.Epidemic_StartDate_CalTime = P.Interventions_StartDate_CalTime - P.InitialInfectionCalTime;
	}
	Files::xfprintf_stderr("Finished InitModel.\n");
}

void SeedInfection(double t, int* NumSeedingInfections_byLocation, int AlreadyInitialized, int run) //adding run number to pass it to event log
{
	/* *NumSeedingInfections_byLocation is an array of the number of seeding infections by location. During runtime, usually just a single int (given by a poisson distribution)*/
	/*AlreadyInitialized set to 0 when initializing model, otherwise set to 1 during runtime. i.e. when initial seeds are being set, not when imported cases are being set*/

	int mcellnum /*microcell number*/;
	int Person = 0, Mcell_x = 0, Mcell_y = 0;
	int NumMCellSeedingChoices = 0; // Although inconsistently, zero coded as choice successful and no more choices needed. if NumMCellSeedingChoices = 0, RunningTotalGuesses = 0. 
	int RunningTotalGuesses = 0 /*range = {0, 1000}*/;
	int NumSeedLocations = ((AlreadyInitialized == 0) ? P.NumSeedLocations : 1); /*number of seed locations?*/;

	for (int SeedLocation = 0; SeedLocation < NumSeedLocations; SeedLocation++)
	{
		if ((!P.DoRandomInitialInfectionLoc) || ((P.DoAllInitialInfectioninSameLoc) && (AlreadyInitialized))) //// either non-random locations, doing all initial infections in same location, and not initializing.
		{
			Mcell_x = (int)(P.LocationInitialInfection[SeedLocation][0] / P.in_microcells_.width);
			Mcell_y = (int)(P.LocationInitialInfection[SeedLocation][1] / P.in_microcells_.height);
			mcellnum = Mcell_x * P.total_microcells_high_ + Mcell_y;
			NumMCellSeedingChoices = 0;
			for (int Infection = 0; (Infection < NumSeedingInfections_byLocation[SeedLocation]) && (NumMCellSeedingChoices < 10000); Infection++)
			{
				int Person = Mcells[mcellnum].members[(int)(ranf() * ((double)Mcells[mcellnum].n))]; //// randomly choose member of microcell mcellnum. Name this member l
				if (Hosts[Person].is_susceptible()) //// If Host l is uninfected.
				{
					if ((CalcPersonSusc(Person, 0, 0) > 0) && (Hosts[Person].age <= P.MaxAgeForInitialInfection) &&
						(P.CareHomeAllowInitialInfections || P.CareHomePlaceType < 0 || Hosts[Person].PlaceLinks[P.CareHomePlaceType] < 0))
					{
						//only reset the initial location if AlreadyInitialized == 0, i.e. when initial seeds are being set, not when imported cases are being set
						if (AlreadyInitialized == 0)
						{
							P.LocationInitialInfection[SeedLocation][0] = Households[Hosts[Person].hh].loc.x;
							P.LocationInitialInfection[SeedLocation][1] = Households[Hosts[Person].hh].loc.y;
						}
						Hosts[Person].infector = -2;
						Hosts[Person].infect_type = INFECT_TYPE_MASK - 1;
						DoInfect(Person, t, 0, run); ///// guessing this updates a number of things about person l at time t in thread 0 for this run.
						NumMCellSeedingChoices = 0;
					}
				}
				else { Infection--; NumMCellSeedingChoices++; } //// k-- means if person l chosen is already infected, go again. The NumMCellSeedingChoices < 10000 is a guard against a) too many infections; b) an infinite loop if no more uninfected people left.
			}
		}
		else if (P.DoAllInitialInfectioninSameLoc) // difference between this block and block below is that here you chosse a single microcell to distribute all seeding infections for this location. Block below chooses different microcells. 
		{
			RunningTotalGuesses = 0; // initialize to zero.
			do // while RunningTotalGuesses between 1 and 999, i.e. while still need to choose person to seed infection into, but also while loop hasn't made more than 1000 attempts.
			{
				NumMCellSeedingChoices = 0; // initialize to zero.

				// choose Person (in order to choose microcell, really just sampling microcells weighted by population. Peron within microcell will be re-chosen later)
				do
				{
					Person = (int)(ranf() * ((double)P.PopSize));
					mcellnum = Hosts[Person].mcell;

				} while ((Mcells[mcellnum].n < NumSeedingInfections_byLocation[SeedLocation]) // choose person again if their microcell has fewer people than number seeded in this location
					|| (Mcells[mcellnum].n > P.MaxPopDensForInitialInfection) // choose person again if their microcell has too many people
					|| (Mcells[mcellnum].n < P.MinPopDensForInitialInfection) // choose person again if their microcell doesn't have enough people
					|| ((P.InitialInfectionsAdminUnit[SeedLocation] > 0) && ((AdUnits[Mcells[mcellnum].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[SeedLocation] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));

				/// having chosen microcell, choose person to potentially infect
				for (int Infection = 0; (Infection < NumSeedingInfections_byLocation[SeedLocation]) && (NumMCellSeedingChoices < 10000); Infection++)
				{
					// choose Peron within microcell
					Person = Mcells[mcellnum].members[(int)(ranf() * ((double)Mcells[mcellnum].n))];
					if (Hosts[Person].is_susceptible())
					{
						if ((CalcPersonSusc(Person, 0, 0) > 0) && // if person has non-zero susceptibility
							(Hosts[Person].age <= P.MaxAgeForInitialInfection) && // and they're not too young
							(P.CareHomeAllowInitialInfections || P.CareHomePlaceType < 0 || Hosts[Person].PlaceLinks[P.CareHomePlaceType] < 0))
						{
							P.LocationInitialInfection[SeedLocation][0] = Households[Hosts[Person].hh].loc.x;
							P.LocationInitialInfection[SeedLocation][1] = Households[Hosts[Person].hh].loc.y;
							Hosts[Person].infector = -2; Hosts[Person].infect_type = INFECT_TYPE_MASK - 1;
							DoInfect(Person, t, 0, run);
							NumMCellSeedingChoices = 0; // can move on from do-while loop
						}
					}
					else
					{
						/// Choose person again. Increment number of choices 
						Infection--;
						NumMCellSeedingChoices++;
					}
				}
				if (NumMCellSeedingChoices)
					RunningTotalGuesses++;
				else
					RunningTotalGuesses = 0;
			} while ((RunningTotalGuesses > 0) && (RunningTotalGuesses < 1000));
		}
		else // difference between this block and block above is that above you chosse a single microcell to distribute all seeding infections for this location. Block here chooses different microcells. 
		{
			NumMCellSeedingChoices = 0;
			for (int Infection = 0; (Infection < NumSeedingInfections_byLocation[SeedLocation]) && (NumMCellSeedingChoices < 10000); Infection++)
			{
				// choose Person (in order to choose microcell, really just sampling microcells weighted by population. Peron within microcell will be re-chosen later)
				do
				{
					Person = (int)(ranf() * ((double)P.PopSize));
					mcellnum = Hosts[Person].mcell;

				} while ((Mcells[mcellnum].n == 0) || (Mcells[mcellnum].n > P.MaxPopDensForInitialInfection) // choose again if microcell is unpopulated, not populated enough or too popuated
					|| (Mcells[mcellnum].n < P.MinPopDensForInitialInfection)
					|| ((P.InitialInfectionsAdminUnit[SeedLocation] > 0) && ((AdUnits[Mcells[mcellnum].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[SeedLocation] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));

				/// having chosen microcell, choose person to potentially infect
				Person = Mcells[mcellnum].members[(int)(ranf() * ((double)Mcells[mcellnum].n))];
				if (Hosts[Person].is_susceptible())
				{
					if ((CalcPersonSusc(Person, 0, 0) > 0) && // if person has non-zero susceptibility
						(Hosts[Person].age <= P.MaxAgeForInitialInfection) && // and they're not too young
						(P.CareHomeAllowInitialInfections || P.CareHomePlaceType < 0 || Hosts[Person].PlaceLinks[P.CareHomePlaceType] < 0))
					{
						P.LocationInitialInfection[SeedLocation][0] = Households[Hosts[Person].hh].loc.x;
						P.LocationInitialInfection[SeedLocation][1] = Households[Hosts[Person].hh].loc.y;
						Hosts[Person].infector = -2; Hosts[Person].infect_type = INFECT_TYPE_MASK - 1;
						DoInfect(Person, t, 0, run);
						NumMCellSeedingChoices = 0;
					}
					else
					{
						/// Choose person again. Increment number of choices
						Infection--;
						NumMCellSeedingChoices++;
					}
				}
				else
				{
					/// Choose person again. Increment number of choices
					Infection--;
					NumMCellSeedingChoices++;
				}
			}
		}
	}
	if (NumMCellSeedingChoices > 0) Files::xfprintf_stderr("### Seeding error ###\n");
}

int RunModel(int run, std::string const& snapshot_save_file, std::string const& snapshot_load_file, std::string const& output_file_base)
{
	//// **** Structure of function is as follows. For each timestep: 
		// i) Seed Infections with SeedInfection function
		// ii) Have infected people infect other people via InfectSweep function
		// iii) Move people along their disease progression via IncubRecoverySweep function
		// iv) Digitially contact trace them with DigitalContactTracingSweep function
		// v) Treat infected people with TreatSweep function
		// Calculate/Record Model output if timestep is an output timestep. 

	int KeepRunning = 1, IsEpidemicStillGoing = 0, NumSeedingInfections, NumSeedingInfections_byLocation[MAX_NUM_SEED_LOCATIONS] /*Denotes either Num imported Infections given rate ir, or number false positive "infections"*/;
	double InfectionImportRate; // infection import rate?;
	double CurrSimTime, ProportionSusceptible = 1, PreviousProportionSusceptible = 1, t2;
	unsigned short int CurrTimeStep; //// Timestep in simulation time.
	int continueEvents = 1;

	InterruptRun = 0; // global variable set to zero at start of RunModel, and possibly modified in CalibrationThresholdCheck
	if (snapshot_load_file.empty())
	{
		CurrSimTime = 0;
		P.ts_age = 0;
	}
	else
	{
		P.ts_age = (int)(P.SnapshotLoadTime * P.TimeStepsPerDay);
		CurrSimTime = ((double)P.ts_age) * P.ModelTimeStep;
	}

	for (OutputTimeStepNumber = 1; ((OutputTimeStepNumber < P.NumOutputTimeSteps) && (!InterruptRun)); OutputTimeStepNumber++) // OutputTimeStepNumber starts from 1 here as is zero in InitModel
	{
		RecordSample				(CurrSimTime, OutputTimeStepNumber - 1, output_file_base);
		CalibrationThresholdCheck	(CurrSimTime, OutputTimeStepNumber - 1);
		UpdateCFRs					(CurrSimTime - P.Epidemic_StartDate_CalTime); 

		/// print various quantities to console
		Files::xfprintf_stderr("\r    t=%lg   %i    %i|%i    %i     %i [=%i]  %i (%lg %lg %lg)   %lg    ", CurrSimTime,
			State.S, State.L, State.I, State.R, State.D, State.S + State.L + State.I + State.R + State.D, State.cumD, State.cumT, State.cumV, State.cumVG, sqrt(State.maxRad2) / 1000); //added State.cumVG

		if (!InterruptRun)
		{
			//Only run to a certain number of infections: ggilani 28/10/14
			if (P.LimitNumInfections) continueEvents = (State.cumI < P.MaxNumInfections);

			for (int ModelTimeStep = 0; ((ModelTimeStep < P.NumModelTimeStepsPerOutputTimeStep) && (!InterruptRun) && (continueEvents)); ModelTimeStep++) // local (int)ModelTimeStep not used, but CurrSimTime is updated by P.ModelTimeStep.
			{
				CurrTimeStep = (unsigned short int) (P.TimeStepsPerDay * CurrSimTime);

				//if we are to reset random numbers after an intervention event, specific time
				if (P.ResetSeedsPostIntervention)
					if ((P.ResetSeedsFlag == 0) && (CurrTimeStep >= (P.TimeToResetSeeds * P.TimeStepsPerDay)))
					{
						setall(&P.nextRunSeed1, &P.nextRunSeed2);
						P.ResetSeedsFlag = 1;
					}

				if (KeepRunning)
				{
					if (P.DoAirports) TravelDepartSweep(CurrSimTime);

					// calculate importation rate (if appropriate)
					if (P.DurImportTimeProfile > 0)
					{
						if ((int)CurrSimTime < P.DurImportTimeProfile)
							InfectionImportRate = P.ImportInfectionTimeProfile[(int)CurrSimTime] * ((CurrSimTime > P.InfectionImportChangeTime) ? (P.InfectionImportRate2 / P.InfectionImportRate1) : 1.0);
						else
							InfectionImportRate = 0;
					}
					else	InfectionImportRate = (CurrSimTime > P.InfectionImportChangeTime) ? P.InfectionImportRate2 : P.InfectionImportRate1;

					// Calculated number of seeding infections (and seed infections).
					if (InfectionImportRate > 0) //// if infection import rate > 0, seed some infections
					{
						for (int SeedLoc = NumSeedingInfections = 0; SeedLoc < P.NumSeedLocations; SeedLoc++)
						{
							// sample number imported infections in this location from from Poisson distribution.
							NumSeedingInfections_byLocation[SeedLoc] = (int)ignpoi(P.ModelTimeStep * InfectionImportRate * P.InitialInfectionsAdminUnitWeight[SeedLoc] * P.SeedingScaling); 
							// Add to total
							NumSeedingInfections += NumSeedingInfections_byLocation[SeedLoc];
						}
						// ** // ** SeedInfection
						if (NumSeedingInfections > 0)	SeedInfection(CurrSimTime, NumSeedingInfections_byLocation, 1, run);
					}

					if (P.FalsePositivePerCapitaIncidence > 0)
					{
						NumSeedingInfections = (int)ignpoi(P.ModelTimeStep * P.FalsePositivePerCapitaIncidence * ((double)P.PopSize));

						if (NumSeedingInfections > 0)
						{
							int Person = 0; 
							for (int SeedLoc = 0; SeedLoc < NumSeedingInfections; SeedLoc++)
							{
								do
								{
									Person = (int)(((double)P.PopSize) * ranf()); //// choose person lPerson randomly from entire population. (but change person if they're dead or not a false positive)
								} while (Hosts[Person].is_dead() || (ranf() > P.FalsePositiveAgeRate[HOST_AGE_GROUP(Person)]));
								DoFalseCase(Person, CurrSimTime, CurrTimeStep, 0);
							}
						}
					}

					// ** // ** Infected Sweep: loops over all infectious people and decides which susceptible people to infect (at household, place and spatial level), and adds them to queue. Then changes each person's various characteristics using DoInfect function.  adding run number as a parameter to infect sweep so we can track run number: ggilani - 15/10/14
					InfectSweep(CurrSimTime, run);  //
					// ** // ** IncubRecoverySweep: loops over all infecteds (either latent or infectious). If CurrSimTime is the right time, latent people moved to being infected, and infectious people moved to being clinical cases. Possibly also adds them to recoveries or deaths. Add them to hospitalisation & hospitalisation discharge queues.
					if (!P.DoSI) IncubRecoverySweep(CurrSimTime);
					// ** // ** DigitalContactTracingSweep: If doing new contact tracing, update numbers of people under contact tracing after each time step
					if (P.DoDigitalContactTracing)
						DigitalContactTracingSweep(CurrSimTime);

					IsEpidemicStillGoing = ((P.DoDeath) || (State.L + State.I > 0) || (InfectionImportRate > 0) || (P.FalsePositivePerCapitaIncidence > 0));

					// ** // ** TreatSweep loops over microcells to decide which cells are treated (either with treatment, vaccine, social distancing, movement restrictions etc.). Calls DoVacc, DoPlaceClose, DoProphNoDelay etc. to change (threaded) State variables
					if (!TreatSweep(CurrSimTime)) // TreatSweep will return zero if no treatments are used at CurrSimTime
						if ((!IsEpidemicStillGoing) && (State.L + State.I == 0) && (P.FalsePositivePerCapitaIncidence == 0)) // i.e. if no more infections and no false positives
							if ((InfectionImportRate == 0) && (((int)CurrSimTime) > P.DurImportTimeProfile)) KeepRunning = 0;

					if (P.DoAirports) TravelReturnSweep(CurrSimTime);
					UpdateHostClosure();
				}
				CurrSimTime += P.ModelTimeStep;

				if (P.DoDeath) P.ts_age++;
				// save snapshot (possibly)
				if (!snapshot_save_file.empty() && (CurrSimTime <= P.SnapshotSaveTime) && (CurrSimTime + P.ModelTimeStep > P.SnapshotSaveTime)) SaveSnapshot(snapshot_save_file);
				// Add to Maximum number of treatment courses
				if (CurrSimTime > P.TreatNewCoursesStartTime) P.TreatMaxCourses += P.ModelTimeStep * P.TreatNewCoursesRate;
				// Add to Maximum number of vaccine courses
				if ((CurrSimTime > P.VaccNewCoursesStartTime) && (CurrSimTime < P.VaccNewCoursesEndTime)) P.VaccMaxCourses += P.ModelTimeStep * P.VaccNewCoursesRate;

				ProportionSusceptible = ((double)(State.S)) / ((double)P.PopSize);
				if ((PreviousProportionSusceptible - ProportionSusceptible) > 0.2)
				{
					PreviousProportionSusceptible = ProportionSusceptible;
					UpdateProbs(0);
					DoInitUpdateProbs = 1;
				}
			}
		}
	}
	if (!InterruptRun) RecordSample(CurrSimTime, P.NumOutputTimeSteps - 1, output_file_base);
	Files::xfprintf_stderr("\nEnd of run\n");
	t2 = CurrSimTime + P.SimulationDuration;
	while (KeepRunning)
	{
		KeepRunning = TreatSweep(t2);
		t2 += P.OutputTimeStep;
	}
	//	Files::xfprintf_stderr(,"End RunModel\n");
	if (P.DoAirports)
	{
		t2 = CurrSimTime;
		for (t2 = CurrSimTime; t2 <= CurrSimTime + MAX_TRAVEL_TIME; t2 += P.ModelTimeStep)
			TravelReturnSweep(t2);
	}

	if(!InterruptRun) RecordInfTypes();
	return (InterruptRun);
}

void SaveDistribs(std::string const& output_file_base)
{
	int i, j, k;
	FILE* dat;
	std::string outname;
	double s;

	if (P.DoPlaces)
	{
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != P.HotelPlaceType)
			{
				for (i = 0; i < P.Nplace[j]; i++)
					Places[j][i].n = 0;
				for (i = 0; i < P.PopSize; i++)
				{
					if (Hosts[i].PlaceLinks[j] >= P.Nplace[j])
						Files::xfprintf_stderr("*%i %i: %i %i", i, j, Hosts[i].PlaceLinks[j], P.Nplace[j]);
					else if (Hosts[i].PlaceLinks[j] >= 0)
						Places[j][Hosts[i].PlaceLinks[j]].n++;
				}
			}
		for (j = 0; j < P.PlaceTypeNum; j++)
			for (i = 0; i < MAX_DIST; i++)
				PlaceDistDistrib[j][i] = 0;
		for (i = 0; i < P.PopSize; i++)
			for (j = 0; j < P.PlaceTypeNum; j++)
				if ((j != P.HotelPlaceType) && (Hosts[i].PlaceLinks[j] >= 0))
				{
					if (Hosts[i].PlaceLinks[j] >= P.Nplace[j])
						Files::xfprintf_stderr("*%i %i: %i ", i, j, Hosts[i].PlaceLinks[j]);
					else if ((!P.DoOutputPlaceDistForOneAdunit) ||
						((AdUnits[Mcells[Hosts[i].mcell].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor == (P.OutputPlaceDistAdunit % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor))
					{
						k = Hosts[i].PlaceLinks[j];
						s = sqrt(dist2_raw(Households[Hosts[i].hh].loc.x, Households[Hosts[i].hh].loc.y, Places[j][k].loc.x, Places[j][k].loc.y)) / OUTPUT_DIST_SCALE;
						k = (int)s;
						if (k < MAX_DIST) PlaceDistDistrib[j][k]++;
					}
				}
		for (j = 0; j < P.PlaceTypeNum; j++)
			for (i = 0; i < MAX_PLACE_SIZE; i++)
				PlaceSizeDistrib[j][i] = 0;
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != P.HotelPlaceType)
				for (i = 0; i < P.Nplace[j]; i++)
					if (Places[j][i].n < MAX_PLACE_SIZE)
						PlaceSizeDistrib[j][Places[j][i].n]++;
		outname = output_file_base + ".placedist.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "dist");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != P.HotelPlaceType)
				Files::xfprintf(dat, "\tfreq_p%i", j);
		Files::xfprintf(dat, "\n");
		for (i = 0; i < MAX_DIST; i++)
		{
			Files::xfprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != P.HotelPlaceType)
					Files::xfprintf(dat, "\t%i", PlaceDistDistrib[j][i]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
		outname = output_file_base + ".placesize.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "size");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != P.HotelPlaceType)
				Files::xfprintf(dat, "\tfreq_p%i", j);
		Files::xfprintf(dat, "\n");
		for (i = 0; i < MAX_PLACE_SIZE; i++)
		{
			Files::xfprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != P.HotelPlaceType)
					Files::xfprintf(dat, "\t%i", PlaceSizeDistrib[j][i]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}
}
void SaveOriginDestMatrix(std::string const& output_file_base)
{
	/** function: SaveOriginDestMatrix(std::string const&)
	 *
	 * purpose: to save the calculated origin destination matrix to file
	 * parameters: name of output file
	 * returns: none
	 *
	 * author: ggilani, 13/02/15
	 */
	int i, j;
	FILE* dat;

	std::string outname = output_file_base + ".origdestmat.xls";
	dat = Files::xfopen(outname.c_str(), "wb");
	Files::xfprintf(dat, "0,");
	for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "%i,", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
	Files::xfprintf(dat, "\n");
	for (i = 0; i < P.NumAdunits; i++)
	{
		Files::xfprintf(dat, "%i,", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
		for (j = 0; j < P.NumAdunits; j++)
		{
			Files::xfprintf(dat, "%.10f,", AdUnits[i].origin_dest[j]);
		}
		Files::xfprintf(dat, "\n");
	}
	Files::xfclose(dat);
}

void SaveResults(std::string const& output_file_base)
{
	int i, j;
	FILE* dat;
	std::string outname;

	if (P.OutputNonSeverity)
	{
		outname = output_file_base + ".xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t\tS\tL\tI\tR\tD\tincI\tincR\tincFC\tincC\tincDC\tincTC\tincCT\tincCC\tcumT\tcumTP\tcumV\tcumVG\tExtinct\trmsRad\tmaxRad\n");//\t\t%.10f\t%.10f\t%.10f\n",P.R0household,P.R0places,P.R0spatial);
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10ft%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].L, TimeSeries[i].I,
				TimeSeries[i].R, TimeSeries[i].D, TimeSeries[i].incI,
				TimeSeries[i].incR, TimeSeries[i].incFC, TimeSeries[i].incC, TimeSeries[i].incDC, TimeSeries[i].incTC, TimeSeries[i].incCT, TimeSeries[i].incCC,
				TimeSeries[i].cumT, TimeSeries[i].cumTP, TimeSeries[i].cumV, TimeSeries[i].cumVG, TimeSeries[i].extinct, TimeSeries[i].rmsRad, TimeSeries[i].maxRad);
		}
		Files::xfclose(dat);
	}

	if ((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		outname = output_file_base + ".adunit.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tI_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tC_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tDC_%s", AdUnits[i].ad_name);

		Files::xfprintf(dat, "\n");
		for (i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f", TimeSeries[i].t);
			for (j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incI_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incC_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incDC_adunit[j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}

	if ((P.DoDigitalContactTracing) && (P.DoAdUnits) && (P.OutputDigitalContactTracing))
	{
		outname = output_file_base + ".digitalcontacttracing.xls"; //modifying to csv file
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++)
		{
			Files::xfprintf(dat, "\tincDCT_%s", AdUnits[i].ad_name);
		}
		for (i = 0; i < P.NumAdunits; i++)
		{
			Files::xfprintf(dat, "\tDCT_%s", AdUnits[i].ad_name);
		}
		Files::xfprintf(dat, "\n");
		//print actual output
		for(i=0; i<P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10lf", TimeSeries[i].t);
			for (j = 0; j < P.NumAdunits; j++)
			{
				Files::xfprintf(dat, "\t%.10lf", TimeSeries[i].incDCT_adunit[j]);
			}
			for (j = 0; j < P.NumAdunits; j++)
			{
				Files::xfprintf(dat, "\t%.10lf", TimeSeries[i].DCT_adunit[j]);
			}
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);

	}

	if (P.DoDigitalContactTracing && P.OutputDigitalContactDist)
	{
		outname = output_file_base + ".digitalcontactdist.xls"; //modifying to csv file
		dat = Files::xfopen(outname.c_str(), "wb");
		//print headers
		Files::xfprintf(dat, "nContacts\tFrequency\n");
		for (i = 0; i < (MAX_CONTACTS + 1); i++)
		{
			Files::xfprintf(dat, "%i\t%i\n", i, State.contact_dist[i]);
		}
		Files::xfclose(dat);
	}

	if(P.KeyWorkerProphTimeStartBase < P.SimulationDuration)
	{
		outname = output_file_base + ".keyworker.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tI%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tC%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tT%i", i);
		Files::xfprintf(dat, "\t%i\t%i\n", P.KeyWorkerNum, P.KeyWorkerIncHouseNum);
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f", TimeSeries[i].t);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incI_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incC_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumT_keyworker[j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}

	if(P.DoInfectionTree)
	{
		outname = output_file_base + "%s.tree.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		for(i = 0; i < P.PopSize; i++)
			if(Hosts[i].infect_type % INFECT_TYPE_MASK > 0)
				Files::xfprintf(dat, "%i\t%i\t%i\t%i\n", i, Hosts[i].infector, Hosts[i].infect_type % INFECT_TYPE_MASK, (int)HOST_AGE_YEAR(i));
		Files::xfclose(dat);
	}
#if defined(_WIN32) || defined(IMAGE_MAGICK)
	static int dm[13] ={0,31,28,31,30,31,30,31,31,30,31,30,31};
	int d, m, y, dml, f;
#ifdef _WIN32
	//if(P.OutputBitmap == 1) CloseAvi(avi);
	//if((TimeSeries[P.NumOutputTimeSteps - 1].extinct) && (P.OutputOnlyNonExtinct))
	//	{
	//	outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base + ".avi";
	//	DeleteFile(outname);
	//	}
#endif
	if(P.OutputBitmap >= 1 && P.BitmapFormat == BitmapFormats::PNG)
	{
		// Generate Google Earth .kml file
		outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base + ".ge.kml"; // outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base ".kml";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.2\">\n<Document>\n");
		Files::xfprintf(dat, "<name>%s</name>\n", output_file_base.c_str());
		y = 2009;
		m = 1;
		d = 1;
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "<GroundOverlay>\n<name>Snapshot %i</name>\n", i + 1);
			Files::xfprintf(dat, "<TimeSpan>\n<begin>%i-%02i-%02iT00:00:00Z</begin>\n", y, m, d);
			d += (int)P.OutputTimeStep; // OutputTimeStep has to be an integer here.
			do
			{
				f = 1;
				dml = dm[m];
				if((m == 2) && (y % 4 == 0)) dml = 29;
				if(d > dml)
				{
					m++;
					if(m > 12)
					{
						m -= 12;
						y++;
					}
					d -= dml;
					f = 0;
				}
			} while(!f);
			Files::xfprintf(dat, "<end>%i-%02i-%02iT00:00:00Z</end>\n</TimeSpan>\n", y, m, d);
			outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base + "." + std::to_string(i + 1) + ".png";
			Files::xfprintf(dat, "<Icon>\n<href>%s</href>\n</Icon>\n", outname.c_str());
			Files::xfprintf(dat, "<LatLonBox>\n<north>%.10f</north>\n<south>%.10f</south>\n<east>%.10f</east>\n<west>%.10f</west>\n</LatLonBox>\n",
				P.SpatialBoundingBox.top_right().y, P.SpatialBoundingBox.bottom_left().y, P.SpatialBoundingBox.top_right().x, P.SpatialBoundingBox.bottom_left().x);
			Files::xfprintf(dat, "</GroundOverlay>\n");
		}
		Files::xfprintf(dat, "</Document>\n</kml>\n");
		Files::xfclose(dat);
	}
#endif


	if((P.DoSeverity)&&(P.OutputSeverity))
	{
		outname = output_file_base + ".severity.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t\tRt\tTG\tSI\tS\tI\tR\tincI\tMild\tILI\tSARI\tCritical\tCritRecov\tincMild\tincILI\tincSARI\tincCritical\tincCritRecov\tincDeath\tincDeath_ILI\tincDeath_SARI\tincDeath_Critical\tcumMild\tcumILI\tcumSARI\tcumCritical\tcumCritRecov\tcumDeath\tcumDeath_ILI\tcumDeath_SARI\tcumDeath_Critical\n");//\t\t%.10f\t%.10f\t%.10f\n",P.R0household,P.R0places,P.R0spatial);
		for (i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				TimeSeries[i].t, TimeSeries[i].Rdenom, TimeSeries[i].meanTG, TimeSeries[i].meanSI, TimeSeries[i].S, TimeSeries[i].I, TimeSeries[i].R, TimeSeries[i].incI,
				TimeSeries[i].Mild		, TimeSeries[i].ILI		, TimeSeries[i].SARI	, TimeSeries[i].Critical	, TimeSeries[i].CritRecov	,
				TimeSeries[i].incMild	, TimeSeries[i].incILI	, TimeSeries[i].incSARI	, TimeSeries[i].incCritical	, TimeSeries[i].incCritRecov,
				TimeSeries[i].incD,	TimeSeries[i].incDeath_ILI, TimeSeries[i].incDeath_SARI, TimeSeries[i].incDeath_Critical,
				TimeSeries[i].cumMild	, TimeSeries[i].cumILI	, TimeSeries[i].cumSARI	, TimeSeries[i].cumCritical	, TimeSeries[i].cumCritRecov, TimeSeries[i].D	,
				TimeSeries[i].cumDeath_ILI, TimeSeries[i].cumDeath_SARI, TimeSeries[i].cumDeath_Critical);
		}
		Files::xfclose(dat);

		if((P.DoAdUnits) && (P.OutputSeverityAdminUnit))
		{
			//// output severity results by admin unit
			outname = output_file_base + ".severity.adunit.xls";
			dat = Files::xfopen(outname.c_str(), "wb");
			Files::xfprintf(dat, "t");

			/////// ****** /////// ****** /////// ****** COLNAMES
			const std::string colnames[] = {
				// prevalence
				"Mild_", "ILI_", "SARI_", "Critical_", "CritRecov_",
				// incidence
				"incI_", "incMild_", "incILI_", "incSARI_", "incCritical_", "incCritRecov_",
				"incDeath_adu", "incDeath_ILI_adu", "incDeath_SARI_adu", "incDeath_Critical_adu"
				// cumulative incidence
				"cumMild_", "cumILI_", "cumSARI_", "cumCritical_", "cumCritRecov_", "cumDeaths_",
				"cumDeath_ILI_", "cumDeath_SARI_", "cumDeath_Critical_"
			};
			for (auto colname : colnames)
			{
				for (i = 0; i < P.NumAdunits; i++)
				{
					Files::xfprintf(dat, "\t%s%s", colname.c_str(), AdUnits[i].ad_name);
				}
			}

			Files::xfprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table.
			for(i = 0; i < P.NumOutputTimeSteps; i++)
			{
				Files::xfprintf(dat, "%.10f", TimeSeries[i].t);

				//// prevalence
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].Mild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].Critical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].CritRecov_adunit[j]);

				//// incidence
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].incDeath_Critical_adunit[j]);

				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", TimeSeries[i].cumDeath_Critical_adunit[j]);

				if(i != P.NumOutputTimeSteps - 1) Files::xfprintf(dat, "\n");
			}
			Files::xfclose(dat);
		}
	}

	if (P.DoAdUnits && P.OutputAdUnitAge)
	{
		//// output infections by age and admin unit
		outname = output_file_base + ".age.adunit.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");

		// colnames
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				Files::xfprintf(dat, "\tincInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// incidence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				Files::xfprintf(dat, "\tprevInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// prevalence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				Files::xfprintf(dat, "\tcumInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// cumulative incidence
		Files::xfprintf(dat, "\n");

		// Populate
		for (int Time = 0; Time < P.NumOutputTimeSteps; Time++)
		{
			Files::xfprintf(dat, "%.10f", TSMean[Time].t);
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					Files::xfprintf(dat, "\t%.10f", TimeSeries[Time].incInf_age_adunit[AgeGroup][AdUnit]);	// incidence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					Files::xfprintf(dat, "\t%.10f", TimeSeries[Time].prevInf_age_adunit[AgeGroup][AdUnit]);	// prevalence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					Files::xfprintf(dat, "\t%.10f", TimeSeries[Time].cumInf_age_adunit[AgeGroup][AdUnit]);	// cumulative incidence
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}
}

void SaveSummaryResults(std::string const& output_file_base) //// calculates and saves summary results (called for average of extinct and non-extinct realisation time series - look in main)
{
	int i, j;
	double c, t;
	FILE* dat;
	std::string outname;

	c = 1 / ((double)(_I64(P.NRactE) + P.NRactNE));

	if (P.OutputNonSeverity)
	{
		outname = output_file_base + ".xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		//// set colnames
		Files::xfprintf(dat, "t\tS\tL\tI\tR\tD\tincI\tincR\tincD\tincC\tincDC\tincTC\tcumT\tcumTmax\tcumTP\tcumV\tcumVmax\tExtinct\trmsRad\tmaxRad\tvS\tvI\tvR\tvD\tvincI\tvincR\tvincFC\tvincC\tvincDC\tvincTC\tvrmsRad\tvmaxRad\t\t%i\t%i\t%.10f\t%.10f\t%.10f\t\t%.10f\t%.10f\t%.10f\t%.10f\n",
			P.NRactNE, P.NRactE, P.R0household, P.R0places, P.R0spatial, c * PeakHeightSum, c * PeakHeightSS - c * c * PeakHeightSum * PeakHeightSum, c * PeakTimeSum, c * PeakTimeSS - c * c * PeakTimeSum * PeakTimeSum);
		c = 1 / ((double)P.NRactual);

		//// populate table
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].L, c * TSMean[i].I, c * TSMean[i].R,
				c * TSMean[i].D, c * TSMean[i].incI, c * TSMean[i].incR, c * TSMean[i].incFC, c * TSMean[i].incC, c * TSMean[i].incDC, c * TSMean[i].incTC,
				c * TSMean[i].cumT, TSMean[i].cumTmax, c * TSMean[i].cumTP, c * TSMean[i].cumV, TSMean[i].cumVmax, c * TSMean[i].extinct, c * TSMean[i].rmsRad, c * TSMean[i].maxRad);
			Files::xfprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				c * TSVar[i].S		- c * c * TSMean[i].S		* TSMean[i].S,
				c * TSVar[i].I		- c * c * TSMean[i].I		* TSMean[i].I,
				c * TSVar[i].R		- c * c * TSMean[i].R		* TSMean[i].R,
				c * TSVar[i].D		- c * c * TSMean[i].D		* TSMean[i].D,
				c * TSVar[i].incI	- c * c * TSMean[i].incI	* TSMean[i].incI,
				c * TSVar[i].incR	- c * c * TSMean[i].incR	* TSMean[i].incR,
				c * TSVar[i].incD	- c * c * TSMean[i].incD	* TSMean[i].incD,
				c * TSVar[i].incC	- c * c * TSMean[i].incC	* TSMean[i].incC,
				c * TSVar[i].incDC	- c * c * TSMean[i].incDC	* TSMean[i].incDC, //added detected cases
				c * TSVar[i].incTC	- c * c * TSMean[i].incTC	* TSMean[i].incTC,
				c * TSVar[i].rmsRad - c * c * TSMean[i].rmsRad	* TSMean[i].rmsRad,
				c * TSVar[i].maxRad - c * c * TSMean[i].maxRad	* TSMean[i].maxRad);
		}
		Files::xfclose(dat);
	}

	if (P.OutputControls)
	{
		outname = output_file_base + ".controls.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t\tS\tincC\tincTC\tincFC\tcumT\tcumUT\tcumTP\tcumV\tincHQ\tincAC\tincAH\tincAA\tincACS\tincAPC\tincAPA\tincAPCS\tpropSocDist");
		for(j = 0; j < NUM_PLACE_TYPES; j++) Files::xfprintf(dat, "\tprClosed_%i", j);
		Files::xfprintf(dat, "t\tvS\tvincC\tvincTC\tvincFC\tvcumT\tvcumUT\tvcumTP\tvcumV");
		for(j = 0; j < NUM_PLACE_TYPES; j++) Files::xfprintf(dat, "\tvprClosed_%i", j);
		Files::xfprintf(dat, "\n");
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].incC, c * TSMean[i].incTC, c * TSMean[i].incFC,
				c * TSMean[i].cumT, c * TSMean[i].cumUT, c * TSMean[i].cumTP, c * TSMean[i].cumV, c * TSMean[i].incHQ,
				c * TSMean[i].incAC, c * TSMean[i].incAH, c * TSMean[i].incAA, c * TSMean[i].incACS,
				c * TSMean[i].incAPC, c * TSMean[i].incAPA, c * TSMean[i].incAPCS,c*TSMean[i].PropSocDist);
			for(j = 0; j < NUM_PLACE_TYPES; j++) Files::xfprintf(dat, "\t%lf", c * TSMean[i].PropPlacesClosed[j]);
			Files::xfprintf(dat, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				c * TSVar[i].S - c * c * TSMean[i].S * TSMean[i].S,
				c * TSVar[i].incC - c * c * TSMean[i].incC * TSMean[i].incC,
				c * TSVar[i].incTC - c * c * TSMean[i].incTC * TSMean[i].incTC,
				c * TSVar[i].incFC - c * c * TSMean[i].incFC * TSMean[i].incFC,
				c * TSVar[i].cumT - c * c * TSMean[i].cumT * TSMean[i].cumT,
				c * TSVar[i].cumUT - c * c * TSMean[i].cumUT * TSMean[i].cumUT,
				c * TSVar[i].cumTP - c * c * TSMean[i].cumTP * TSMean[i].cumTP,
				c * TSVar[i].cumV - c * c * TSMean[i].cumV * TSMean[i].cumV);
			for(j = 0; j < NUM_PLACE_TYPES; j++) Files::xfprintf(dat, "\t%lf", TSVar[i].PropPlacesClosed[j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);

	}

	if (P.OutputAge)
	{
		outname = output_file_base + ".age.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for(i = 0; i < NUM_AGE_GROUPS; i++)
			Files::xfprintf(dat, "\tI%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for(i = 0; i < NUM_AGE_GROUPS; i++)
			Files::xfprintf(dat, "\tC%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for(i = 0; i < NUM_AGE_GROUPS; i++)
			Files::xfprintf(dat, "\tD%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		Files::xfprintf(dat, "\n");
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f", c * TSMean[i].t);
			for(j = 0; j < NUM_AGE_GROUPS; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incIa[j]);
			for(j = 0; j < NUM_AGE_GROUPS; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incCa[j]);
			for(j = 0; j < NUM_AGE_GROUPS; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDa[j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfprintf(dat, "dist");
		for(j = 0; j < NUM_AGE_GROUPS; j++)
			Files::xfprintf(dat, "\t%.10f", AgeDist[j]);
		Files::xfprintf(dat, "\n");
		Files::xfclose(dat);
	}

	if((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		outname = output_file_base + ".adunit.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for(i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tI_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tC_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tDC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for(i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tT_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\t%i", AdUnits[i].n);
		for(i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\t%.10f", P.PopByAdunit[i][1]);
		Files::xfprintf(dat, "\n");
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f", c * TSMean[i].t);
			for(j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incI_adunit[j]);
			for(j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incC_adunit[j]);
			for(j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
			for(j = 0; j < P.NumAdunits; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumT_adunit[j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);

		if (P.OutputAdUnitVar)
		{
			outname = output_file_base + ".adunitVar.xls";
			dat = Files::xfopen(outname.c_str(), "wb");
			Files::xfprintf(dat, "t");
			for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tI_%s", AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tC_%s", AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tDC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
			for (i = 0; i < P.NumAdunits; i++) Files::xfprintf(dat, "\tT_%s", AdUnits[i].ad_name);
			Files::xfprintf(dat, "\n");
			for (i = 0; i < P.NumOutputTimeSteps; i++)
			{
				Files::xfprintf(dat, "%.10f", c * TSMean[i].t);
				for (j = 0; j < P.NumAdunits; j++)
					Files::xfprintf(dat, "\t%.10f", c * TSVar[i].incI_adunit[j] - c * c * TSMean[i].incI_adunit[j] * TSMean[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)
					Files::xfprintf(dat, "\t%.10f", c * TSVar[i].incC_adunit[j] - c * c * TSMean[i].incC_adunit[j] * TSMean[i].incC_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)
					Files::xfprintf(dat, "\t%.10f", c * TSVar[i].incDC_adunit[j] - c * c * TSMean[i].incDC_adunit[j] * TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
				for (j = 0; j < P.NumAdunits; j++)
					Files::xfprintf(dat, "\t%.10f", c * TSVar[i].cumT_adunit[j] - c * c * TSMean[i].cumT_adunit[j] * TSMean[i].cumT_adunit[j]);
				Files::xfprintf(dat, "\n");
			}
			Files::xfclose(dat);
		}
	}

	if ((P.DoDigitalContactTracing) && (P.DoAdUnits) && (P.OutputDigitalContactTracing))
	{
		outname = output_file_base + ".digitalcontacttracing.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++)
		{
			Files::xfprintf(dat, "\tincDCT_%s", AdUnits[i].ad_name); // //printing headers for inc per admin unit
		}
		for (i = 0; i < P.NumAdunits; i++)
		{
			Files::xfprintf(dat, "\tDCT_%s", AdUnits[i].ad_name); // //printing headers for prevalence of digital contact tracing per admin unit
		}
		Files::xfprintf(dat, "\n");
		//print actual output
		for (i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10lf", c* TSMean[i].t);
			for (j = 0; j < P.NumAdunits; j++)
			{
				Files::xfprintf(dat, "\t%.10lf", c * TSMean[i].incDCT_adunit[j]);
			}
			for (j = 0; j < P.NumAdunits; j++)
			{
				Files::xfprintf(dat, "\t%.10lf", c * TSMean[i].DCT_adunit[j]);
			}
			Files::xfprintf(dat, "\n");
		}

		Files::xfclose(dat);

	}

	if(P.KeyWorkerProphTimeStartBase < P.SimulationDuration)
	{
		outname = output_file_base + ".keyworker.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tI%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tC%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tT%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tvI%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tvC%i", i);
		for(i = 0; i < 2; i++) Files::xfprintf(dat, "\tvT%i", i);
		Files::xfprintf(dat, "\t%i\t%i\n", P.KeyWorkerNum, P.KeyWorkerIncHouseNum);
		for(i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%.10f", c * TSMean[i].t);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incI_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incC_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumT_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSVar[i].incI_keyworker[j] - c * c * TSMean[i].incI_keyworker[j] * TSMean[i].incI_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSVar[i].incC_keyworker[j] - c * c * TSMean[i].incC_keyworker[j] * TSMean[i].incC_keyworker[j]);
			for(j = 0; j < 2; j++)
				Files::xfprintf(dat, "\t%.10f", c * TSVar[i].cumT_keyworker[j] - c * c * TSMean[i].cumT_keyworker[j] * TSMean[i].cumT_keyworker[j]);
			Files::xfprintf(dat,  "\n");
		}
		Files::xfclose(dat);
	}

	if (P.OutputInfType)
	{
		outname = output_file_base + ".inftype.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t\tR\tTG\tSI");
		for (j = 0; j < INFECT_TYPE_MASK; j++) Files::xfprintf(dat, "\tRtype_%i", j);
		for (j = 0; j < INFECT_TYPE_MASK; j++) Files::xfprintf(dat, "\tincItype_%i", j);
		for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\tRage_%i", j);
		Files::xfprintf(dat, "\n");
		for (i = 0; i < P.NumOutputTimeSteps; i++)
		{
			Files::xfprintf(dat, "%lf\t%lf\t%lf\t%lf", c * TSMean[i].t, c * TSMean[i].Rdenom, c* TSMean[i].meanTG, c* TSMean[i].meanSI);
			for (j = 0; j < INFECT_TYPE_MASK; j++) Files::xfprintf(dat, "\t%lf", c * TSMean[i].Rtype[j]);
			for (j = 0; j < INFECT_TYPE_MASK; j++) Files::xfprintf(dat, "\t%lf", c * TSMean[i].incItype[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%lf", c * TSMean[i].Rage[j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}

	if (P.OutputR0)
	{
		outname = output_file_base + ".R0.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		for (i = 0; i < MAX_SEC_REC; i++)
		{
			Files::xfprintf(dat, "%i", i);
			for (j = 0; j < MAX_GEN_REC; j++)
				Files::xfprintf(dat, "\t%.10f", c * indivR0_av[i][j]);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}

	if (P.OutputHousehold)
	{
		outname = output_file_base + ".household.xls";
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			t = 0;
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				t += inf_household_av[i][j];
			inf_household_av[i][0] = denom_household[i] / c - t;
		}
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			t = 0;
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				t += case_household_av[i][j];
			case_household_av[i][0] = denom_household[i] / c - t;
		}
		dat = Files::xfopen(outname.c_str(), "wb");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			Files::xfprintf(dat, "\t%i", i);
		Files::xfprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			Files::xfprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				Files::xfprintf(dat, "\t%.10f", inf_household_av[j][i] * c);
			Files::xfprintf(dat, "\n");
		}
		Files::xfprintf(dat, "\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			Files::xfprintf(dat, "\t%i", i);
		Files::xfprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			Files::xfprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				Files::xfprintf(dat, "\t%.10f", case_household_av[j][i] * c);
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}

	if (P.OutputCountry)
	{
		outname = output_file_base + ".country.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		for (i = 0; i < MAX_COUNTRIES; i++)
			Files::xfprintf(dat, "%i\t%.10f\t%.10f\n", i, infcountry_av[i] * c, infcountry_num[i] * c);
		Files::xfclose(dat);
	}

	if ((P.DoSeverity)&&(P.OutputSeverity))
	{
		//// output separate severity file (can integrate with main if need be)
		outname = output_file_base + ".severity.xls";

		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t\tPropSocDist\tRt\tTG\tSI\tS\tI\tR\tincI\tincC\tMild\tILI\tSARI\tCritical\tCritRecov\tSARIP\tCriticalP\tCritRecovP\tprevQuarNotInfected\tprevQuarNotSymptomatic\tincMild\tincILI\tincSARI\tincCritical\tincCritRecov\tincSARIP\tincCriticalP\tincCritRecovP\tincDeath\tincDeath_ILI\tincDeath_SARI\tincDeath_Critical\tcumMild\tcumILI\tcumSARI\tcumCritical\tcumCritRecov\tcumDeath\tcumDeath_ILI\tcumDeath_SARI\tcumDeath_Critical\t");
		Files::xfprintf(dat, "PropSocDist_v\tRt_v\tTG_v\tSI_v\tS_v\tI_v\tR_v\tincI_v\tincC_v\tMild_v\tILI_v\tSARI_v\tCritical_v\tCritRecov_v\tincMild_v\tincILI_v\tincSARI_v\tincCritical_v\tincCritRecov_v\tincDeath_v\tincDeath_ILI_v\tincDeath_SARI_v\tincDeath_Critical_v\tcumMild_v\tcumILI_v\tcumSARI_v\tcumCritical_v\tcumCritRecov_v\tcumDeath_v\tcumDeath_ILI_v\tcumDeath_SARI_v\tcumDeath_Critical_v\n");
		double SARI, Critical, CritRecov, incSARI, incCritical, incCritRecov, sc1, sc2,sc3,sc4; //this stuff corrects bed prevalence for exponentially distributed time to test results in hospital
		sc1 = (P.Mean_TimeToTest > 0) ? exp(-1.0 / P.Mean_TimeToTest) : 0.0;
		sc2 = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestOffset / P.Mean_TimeToTest) : 0.0;
		sc3 = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCriticalOffset / P.Mean_TimeToTest) : 0.0;
		sc4 = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCritRecovOffset / P.Mean_TimeToTest) : 0.0;
		incSARI = incCritical = incCritRecov = 0;
		for (i = 0; i < P.NumOutputTimeSteps; i++)
		{
			if (i > 0)
			{
				SARI = (TSMean[i].SARI - TSMean[i - 1].SARI) * sc2 + SARI * sc1;
				Critical = (TSMean[i].Critical - TSMean[i - 1].Critical) * sc3 + Critical * sc1;
				CritRecov = (TSMean[i].CritRecov - TSMean[i - 1].CritRecov) * sc4 + CritRecov * sc1;
				incSARI = TSMean[i].incSARI * (1.0 - sc2) + incSARI * sc1;
				incCritical = TSMean[i].incCritical * (1.0 - sc3) + incCritical * sc1;
				incCritRecov = TSMean[i].incCritRecov * (1.0 - sc4) + incCritRecov * sc1;
			}
			else
			{
				SARI = TSMean[i].SARI * sc2;
				Critical = TSMean[i].Critical * sc3;
				CritRecov = TSMean[i].CritRecov * sc4;
			}

			Files::xfprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.17f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t",
				c* TSMean[i].t, c* TSMean[i].PropSocDist, c* TSMean[i].Rdenom, c* TSMean[i].meanTG, c* TSMean[i].meanSI, c* TSMean[i].S, c* TSMean[i].I, c* TSMean[i].R, c* TSMean[i].incI, c* TSMean[i].incC,
				c* TSMean[i].Mild, c* TSMean[i].ILI, c* TSMean[i].SARI,c* TSMean[i].Critical, c* TSMean[i].CritRecov,c* (TSMean[i].SARI - SARI), c* (TSMean[i].Critical - Critical), c* (TSMean[i].CritRecov - CritRecov),
				c * TSMean[i].prevQuarNotInfected, c * TSMean[i].prevQuarNotSymptomatic,
				c * TSMean[i].incMild, c * TSMean[i].incILI, c * TSMean[i].incSARI, c * TSMean[i].incCritical, c * TSMean[i].incCritRecov, c * incSARI, c * incCritical, c * incCritRecov, c * TSMean[i].incD,
				c * TSMean[i].incDeath_ILI, c * TSMean[i].incDeath_SARI, c * TSMean[i].incDeath_Critical,
				c * TSMean[i].cumMild, c * TSMean[i].cumILI, c * TSMean[i].cumSARI, c * TSMean[i].cumCritical, c * TSMean[i].cumCritRecov, c*TSMean[i].D,
				c * TSMean[i].cumDeath_ILI, c * TSMean[i].cumDeath_SARI, c * TSMean[i].cumDeath_Critical);
			Files::xfprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				c* TSVar[i].PropSocDist- c * TSMean[i].PropSocDist* c * TSMean[i].PropSocDist,
				c* TSVar[i].Rdenom - c * TSMean[i].Rdenom * c * TSMean[i].Rdenom,
				c* TSVar[i].meanTG - c * TSMean[i].meanTG * c * TSMean[i].meanTG,
				c* TSVar[i].meanSI - c * TSMean[i].meanSI * c * TSMean[i].meanSI,
				c* TSVar[i].S- c * TSMean[i].S* c * TSMean[i].S,
				c* TSVar[i].I- c * TSMean[i].I* c * TSMean[i].I,
				c* TSVar[i].R- c * TSMean[i].R* c * TSMean[i].R,
				c* TSVar[i].incI- c * TSMean[i].incI* c * TSMean[i].incI,
				c* TSVar[i].incC- c * TSMean[i].incC* c * TSMean[i].incC,
				c* TSVar[i].Mild- c * TSMean[i].Mild* c * TSMean[i].Mild,
				c* TSVar[i].ILI- c * TSMean[i].incILI* c * TSMean[i].incILI,
				c* TSVar[i].SARI- c * TSMean[i].SARI* c * TSMean[i].SARI,
				c* TSVar[i].Critical- c * TSMean[i].Critical* c * TSMean[i].Critical,
				c* TSVar[i].CritRecov- c * TSMean[i].CritRecov* c * TSMean[i].CritRecov,
				c* TSVar[i].incMild- c * TSMean[i].incMild* c * TSMean[i].incMild,
				c* TSVar[i].incILI- c * TSMean[i].incILI* c * TSMean[i].incILI,
				c* TSVar[i].incSARI- c * TSMean[i].incSARI* c * TSMean[i].incSARI,
				c* TSVar[i].incCritical- c * TSMean[i].incCritical* c * TSMean[i].incCritical,
				c* TSVar[i].incCritRecov- c * TSMean[i].incCritRecov* c * TSMean[i].incCritRecov,
				c* TSVar[i].incD- c * TSMean[i].incD* c * TSMean[i].incD,
				c* TSVar[i].incDeath_ILI- c * TSMean[i].incDeath_ILI* c * TSMean[i].incDeath_ILI,
				c* TSVar[i].incDeath_SARI- c * TSMean[i].incDeath_SARI* c * TSMean[i].incDeath_SARI,
				c* TSVar[i].incDeath_Critical- c * TSMean[i].incDeath_Critical* c * TSMean[i].incDeath_Critical,
				c* TSVar[i].cumMild- c * TSMean[i].cumMild* c * TSMean[i].cumMild,
				c* TSVar[i].cumILI- c * TSMean[i].cumILI* c * TSMean[i].cumILI,
				c* TSVar[i].cumSARI- c * TSMean[i].cumSARI* c * TSMean[i].cumSARI,
				c* TSVar[i].cumCritical- c * TSMean[i].cumCritical* c * TSMean[i].cumCritical,
				c* TSVar[i].cumCritRecov- c * TSMean[i].cumCritRecov* c * TSMean[i].cumCritRecov,
				c* TSVar[i].D- c * TSMean[i].D* c * TSMean[i].D,
				c* TSVar[i].cumDeath_ILI- c * TSMean[i].cumDeath_ILI* c * TSMean[i].cumDeath_ILI,
				c* TSVar[i].cumDeath_SARI- c * TSMean[i].cumDeath_SARI* c * TSMean[i].cumDeath_SARI,
				c* TSVar[i].cumDeath_Critical- c * TSMean[i].cumDeath_Critical* c * TSMean[i].cumDeath_Critical);
		}
		Files::xfclose(dat);

		const std::string colnames[] = {
			// prevalence
			"Mild", "ILI", "SARI", "Critical", "CritRecov", "SARIP", "CriticalP", "CritRecovP",
			// incidence
			"incI", "incMild", "incILI", "incSARI", "incCritical", "incCritRecov", "incSARIP",
			"incCriticalP", "incCritRecovP", "incDeath", "incDeath_ILI", "incDeath_SARI", "incDeath__Critical",
			// cumulative incidence
			"cumMild", "cumILI", "cumSARI", "cumCritical", "cumCritRecov", "cumDeaths",
			"cumDeaths_ILI", "cumDeaths_SARI", "cumDeaths_Critical"
		};

		if (P.OutputSeverityAge)
		{
			double* SARI_a, * Critical_a, * CritRecov_a, * incSARI_a, * incCritical_a, * incCritRecov_a, sc1a, sc2a, sc3a, sc4a; //this stuff corrects bed prevalence for exponentially distributed time to test results in hospital

			SARI_a = (double*)Memory::xcalloc(NUM_AGE_GROUPS, sizeof(double));
			Critical_a = (double*)Memory::xcalloc(NUM_AGE_GROUPS, sizeof(double));
			CritRecov_a = (double*)Memory::xcalloc(NUM_AGE_GROUPS, sizeof(double));
			incSARI_a = (double*)Memory::xcalloc(NUM_AGE_GROUPS, sizeof(double));
			incCritical_a = (double*)Memory::xcalloc(NUM_AGE_GROUPS, sizeof(double));
			incCritRecov_a = (double*)Memory::xcalloc(NUM_AGE_GROUPS, sizeof(double));
			sc1a = (P.Mean_TimeToTest > 0) ? exp(-1.0 / P.Mean_TimeToTest) : 0.0;
			sc2a = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestOffset / P.Mean_TimeToTest) : 0.0;
			sc3a = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCriticalOffset / P.Mean_TimeToTest) : 0.0;
			sc4a = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCritRecovOffset / P.Mean_TimeToTest) : 0.0;
			for (i = 0; i < NUM_AGE_GROUPS; i++) incSARI_a[i] = incCritical_a[i] = incCritRecov_a[i] = 0;
			//// output severity results by age group
			outname = output_file_base + ".severity.age.xls";
			dat = Files::xfopen(outname.c_str(), "wb");
			Files::xfprintf(dat, "t");

			for (auto colname : colnames) {
				for (i = 0; i < NUM_AGE_GROUPS; i++)
				{
					Files::xfprintf(dat, "\t%s_%i", colname.c_str(), i);
				}
			}

			Files::xfprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table.
			for (i = 0; i < P.NumOutputTimeSteps; i++)
			{
				for (j = 0; j < NUM_AGE_GROUPS; j++)
				{
					if (i > 0)
					{
						SARI_a[j] = (TSMean[i].SARI_age[j] - TSMean[i - 1].SARI_age[j]) * sc2a + SARI_a[j] * sc1a;
						Critical_a[j] = (TSMean[i].Critical_age[j] - TSMean[i - 1].Critical_age[j]) * sc3a + Critical_a[j] * sc1a;
						CritRecov_a[j] = (TSMean[i].CritRecov_age[j] - TSMean[i - 1].CritRecov_age[j]) * sc4a + CritRecov_a[j] * sc1a;
						incSARI_a[j] = TSMean[i].incSARI_age[j] * (1.0 - sc2a) + incSARI_a[j] * sc1a;
						incCritical_a[j] = TSMean[i].incCritical_age[j] * (1.0 - sc3a) + incCritical_a[j] * sc1a;
						incCritRecov_a[j] = TSMean[i].incCritRecov_age[j] * (1.0 - sc4a) + incCritRecov_a[j] * sc1a;
					}
					else
					{
						SARI_a[j] = TSMean[i].SARI_age[j] * sc2a;
						Critical_a[j] = TSMean[i].Critical_age[j] * sc3a;
						CritRecov_a[j] = TSMean[i].CritRecov_age[j] * sc4a;
					}
				}
				Files::xfprintf(dat, "%.10f", c * TSMean[i].t);
				//// prevalance
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].Mild_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].ILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].SARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].Critical_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].CritRecov_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].SARI_age[j] - SARI_a[j]));
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].Critical_age[j] - Critical_a[j]));
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].CritRecov_age[j] - CritRecov_a[j]));

				//// incidence
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incIa[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incMild_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incSARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incCritical_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incCritRecov_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * incSARI_a[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * incCritical_a[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * incCritRecov_a[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDa[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDeath_ILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDeath_SARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDeath_Critical_age[j]);

				//// cumulative incidence
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumMild_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumSARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumCritical_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumCritRecov_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].cumDeath_ILI_age[j] + TSMean[i].cumDeath_SARI_age[j] + TSMean[i].cumDeath_Critical_age[j]));
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_ILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_SARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_Critical_age[j]);
				Files::xfprintf(dat, "\n");
			}
			Files::xfclose(dat);
			Memory::xfree(SARI_a); Memory::xfree(Critical_a); Memory::xfree(CritRecov_a);
			Memory::xfree(incSARI_a); Memory::xfree(incCritical_a); Memory::xfree(incCritRecov_a);
		}
		if ((P.DoAdUnits) && (P.OutputSeverityAdminUnit))
		{
			double* SARI_a, * Critical_a, * CritRecov_a, * incSARI_a, * incCritical_a, * incCritRecov_a, sc1a, sc2a,sc3a,sc4a; //this stuff corrects bed prevalence for exponentially distributed time to test results in hospital

			SARI_a = (double*)Memory::xcalloc(MAX_ADUNITS, sizeof(double));
			Critical_a = (double*)Memory::xcalloc(MAX_ADUNITS, sizeof(double));
			CritRecov_a = (double*)Memory::xcalloc(MAX_ADUNITS, sizeof(double));
			incSARI_a = (double*)Memory::xcalloc(MAX_ADUNITS, sizeof(double));
			incCritical_a = (double*)Memory::xcalloc(MAX_ADUNITS, sizeof(double));
			incCritRecov_a = (double*)Memory::xcalloc(MAX_ADUNITS, sizeof(double));
			sc1a = (P.Mean_TimeToTest > 0) ? exp(-1.0 / P.Mean_TimeToTest) : 0.0;
			sc2a = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestOffset / P.Mean_TimeToTest) : 0.0;
			sc3a = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCriticalOffset / P.Mean_TimeToTest) : 0.0;
			sc4a = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCritRecovOffset / P.Mean_TimeToTest) : 0.0;
			for (i = 0; i < P.NumAdunits; i++) incSARI_a[i] = incCritical_a[i] = incCritRecov_a[i] = 0;
			//// output severity results by admin unit
			outname = output_file_base + ".severity.adunit.xls";
			dat = Files::xfopen(outname.c_str(), "wb");
			Files::xfprintf(dat, "t");
			for (auto colname : colnames) {
				for (i = 0; i < P.NumAdunits; i++)
				{
					Files::xfprintf(dat, "\t%s_%s", colname.c_str(), AdUnits[i].ad_name);
				}
			}
			Files::xfprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table.
			for (i = 0; i < P.NumOutputTimeSteps; i++)
			{
				for (j = 0; j < P.NumAdunits; j++)
				{
					if (i > 0)
					{
						SARI_a[j] = (TSMean[i].SARI_adunit[j] - TSMean[i - 1].SARI_adunit[j]) * sc2a + SARI_a[j] * sc1a;
						Critical_a[j] = (TSMean[i].Critical_adunit[j] - TSMean[i - 1].Critical_adunit[j]) * sc3a + Critical_a[j] * sc1a;
						CritRecov_a[j] = (TSMean[i].CritRecov_adunit[j] - TSMean[i - 1].CritRecov_adunit[j]) * sc4a + CritRecov_a[j] * sc1a;
						incSARI_a[j] = TSMean[i].incSARI_adunit[j] * (1.0 - sc2a) + incSARI_a[j] * sc1a;
						incCritical_a[j] = TSMean[i].incCritical_adunit[j] * (1.0 - sc3a) + incCritical_a[j] * sc1a;
						incCritRecov_a[j] = TSMean[i].incCritRecov_adunit[j] * (1.0 - sc4a) + incCritRecov_a[j] * sc1a;
					}
					else
					{
						SARI_a[j] = TSMean[i].SARI_adunit[j] * sc2a;
						Critical_a[j] = TSMean[i].Critical_adunit[j] * sc3a;
						CritRecov_a[j] = TSMean[i].CritRecov_adunit[j] * sc4a;
					}
				}
				Files::xfprintf(dat, "%.10f", c*TSMean[i].t);
				//// prevalance
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].Mild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].Critical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].CritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].SARI_adunit[j] - SARI_a[j]));
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].Critical_adunit[j] - Critical_a[j]));
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * (TSMean[i].CritRecov_adunit[j] - CritRecov_a[j]));

				//// incidence
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * incSARI_a[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * incCritical_a[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * incCritRecov_a[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].incDeath_Critical_adunit[j]);

				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		Files::xfprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_Critical_adunit[j]);

				Files::xfprintf(dat, "\n");
			}
			Files::xfclose(dat);
			Memory::xfree(SARI_a); Memory::xfree(Critical_a); Memory::xfree(CritRecov_a);
			Memory::xfree(incSARI_a); Memory::xfree(incCritical_a); Memory::xfree(incCritRecov_a);
		}
	}

	if (P.DoAdUnits && P.OutputAdUnitAge)
	{
		//// output infections by age and admin unit
		outname = output_file_base + ".age.adunit.xls";
		dat = Files::xfopen(outname.c_str(), "wb");
		Files::xfprintf(dat, "t");

		// colnames
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				Files::xfprintf(dat, "\tincInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// incidence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				Files::xfprintf(dat, "\tprevInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// prevalence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				Files::xfprintf(dat, "\tcumInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// cumulative incidence
		Files::xfprintf(dat, "\n");

		// Populate
		for (int Time = 0; Time < P.NumOutputTimeSteps; Time++)
		{
			Files::xfprintf(dat, "%.10f", c * TSMean[Time].t);
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					Files::xfprintf(dat, "\t%.10f", c * TSMean[Time].incInf_age_adunit[AgeGroup][AdUnit]);	// incidence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					Files::xfprintf(dat, "\t%.10f", c * TSMean[Time].prevInf_age_adunit[AgeGroup][AdUnit]);	// prevalence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					Files::xfprintf(dat, "\t%.10f", c * TSMean[Time].cumInf_age_adunit[AgeGroup][AdUnit]);	// cumulative incidence
			Files::xfprintf(dat, "\n");
		}
		Files::xfclose(dat);
	}

}

void SaveRandomSeeds(std::string const& output_file_base)
{
	/* function: SaveRandomSeeds(std::string const&)
	 *
	 * Purpose: outputs the random seeds used for each run to a file
	 * Parameter: name of output file
	 * Returns: none
	 *
	 * Author: ggilani, 09/03/17
	 */
	std::string outname = output_file_base + ".seeds.xls";
	FILE* dat = Files::xfopen(outname.c_str(), "wb");
	Files::xfprintf(dat, "%i\t%i\n", P.nextRunSeed1, P.nextRunSeed2);
	Files::xfclose(dat);
}

void SaveEvents(std::string const& output_file_base)
{
	/* function: SaveEvents(std::string const&)
	 *
	 * Purpose: outputs event log to a csv file if required
	 * Parameters: name of output file
	 * Returns: none
	 *
	 * Author: ggilani, 15/10/2014
	 */
	int i;
	
	std::string outname = output_file_base + ".infevents.xls";
	FILE* dat = Files::xfopen(outname.c_str(), "wb");
	Files::xfprintf(dat, "type,t,thread,ind_infectee,cell_infectee,listpos_infectee,adunit_infectee,x_infectee,y_infectee,t_infector,ind_infector,cell_infector\n");
	for (i = 0; i < nEvents; i++)
	{
		Files::xfprintf(dat, "%i\t%.10f\t%i\t%i\t%i\t%i\t%i\t%.10f\t%.10f\t%.10f\t%i\t%i\n",
			InfEventLog[i].type, InfEventLog[i].t, InfEventLog[i].thread, InfEventLog[i].infectee_ind, InfEventLog[i].infectee_cell, InfEventLog[i].listpos, InfEventLog[i].infectee_adunit, InfEventLog[i].infectee_x, InfEventLog[i].infectee_y, InfEventLog[i].t_infector, InfEventLog[i].infector_ind, InfEventLog[i].infector_cell);
	}
	Files::xfclose(dat);
}

void LoadSnapshot(std::string const& snapshot_load_file)
{
	int i, j, * CellMemberArray, * CellSuscMemberArray;
	int32_t l;
	long long CM_offset, CSM_offset;
	double t;
	int** Array_InvCDF;
	float* Array_tot_prob, ** Array_cum_trans, ** Array_max_trans;

	FILE* dat = Files::xfopen(snapshot_load_file.c_str(), "rb");
	Files::xfprintf_stderr("Loading snapshot.");
	Array_InvCDF = (int**)Memory::xcalloc(P.NumPopulatedCells, sizeof(int*));
	Array_max_trans = (float**)Memory::xcalloc(P.NumPopulatedCells, sizeof(float*));
	Array_cum_trans = (float**)Memory::xcalloc(P.NumPopulatedCells, sizeof(float*));
	Array_tot_prob = (float*)Memory::xcalloc(P.NumPopulatedCells, sizeof(float));
	for (i = 0; i < P.NumPopulatedCells; i++)
	{
		Array_InvCDF[i] = Cells[i].InvCDF;
		Array_max_trans[i] = Cells[i].max_trans;
		Array_cum_trans[i] = Cells[i].cum_trans;
		Array_tot_prob[i] = Cells[i].tot_prob;
	}

	Files::fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.PopSize) ERR_CRITICAL_FMT("Incorrect N (%i %i) in snapshot file.\n", P.PopSize, i);
	Files::fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.NumHouseholds) ERR_CRITICAL("Incorrect NH in snapshot file.\n");
	Files::fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.NumCells) ERR_CRITICAL_FMT("## %i neq %i\nIncorrect NC in snapshot file.", i, P.NumCells);
	Files::fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.NumPopulatedCells) ERR_CRITICAL("Incorrect NCP in snapshot file.\n");
	Files::fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.ncw) ERR_CRITICAL("Incorrect ncw in snapshot file.\n");
	Files::fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.nch) ERR_CRITICAL("Incorrect nch in snapshot file.\n");
	Files::fread_big((void*)& l, sizeof(int32_t), 1, dat); if (l != P.setupSeed1) ERR_CRITICAL("Incorrect setupSeed1 in snapshot file.\n");
	Files::fread_big((void*)& l, sizeof(int32_t), 1, dat); if (l != P.setupSeed2) ERR_CRITICAL("Incorrect setupSeed2 in snapshot file.\n");
	Files::fread_big((void*)& t, sizeof(double), 1, dat); if (t != P.ModelTimeStep) ERR_CRITICAL("Incorrect ModelTimeStep in snapshot file.\n");
	Files::fread_big((void*) & (P.SnapshotLoadTime), sizeof(double), 1, dat);
	P.NumOutputTimeSteps = 1 + (int)ceil((P.SimulationDuration - P.SnapshotLoadTime) / P.OutputTimeStep);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)& CellMemberArray, sizeof(int*), 1, dat);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)& CellSuscMemberArray, sizeof(int*), 1, dat);
	Files::xfprintf_stderr(".");
	CM_offset = State.CellMemberArray - CellMemberArray;
	CSM_offset = State.CellSuscMemberArray - CellSuscMemberArray;

	Files::fread_big((void*)Hosts, sizeof(Person), (size_t)P.PopSize, dat);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)Households, sizeof(Household), (size_t)P.NumHouseholds, dat);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)Cells, sizeof(Cell), (size_t)P.NumCells, dat);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)Mcells, sizeof(Microcell), (size_t)P.NumMicrocells, dat);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	Files::xfprintf_stderr(".");
	Files::fread_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	Files::xfprintf_stderr(".");
	for (i = 0; i < P.NumCells; i++)
	{
		if (Cells[i].n > 0)
		{
			Cells[i].members += CM_offset;
			Cells[i].susceptible += CSM_offset;
			Cells[i].latent += CSM_offset;
			Cells[i].infected += CSM_offset;
		}
		for (j = 0; j < MAX_INTERVENTION_TYPES; j++) Cells[i].CurInterv[j] = -1; // turn interventions off in loaded image
	}
	for (i = 0; i < P.NumMicrocells; i++)
		if (Mcells[i].n > 0)
			Mcells[i].members += CM_offset;

	for (i = 0; i < P.NumPopulatedCells; i++)
	{
		Cells[i].InvCDF = Array_InvCDF[i];
		Cells[i].max_trans = Array_max_trans[i];
		Cells[i].cum_trans = Array_cum_trans[i];
		Cells[i].tot_prob = Array_tot_prob[i];
	}
	Memory::xfree(Array_tot_prob);
	Memory::xfree(Array_cum_trans);
	Memory::xfree(Array_max_trans);
	Memory::xfree(Array_InvCDF);
	Files::xfprintf_stderr("\n");
	Files::xfclose(dat);
}

void SaveSnapshot(std::string const& snapshot_save_file)
{
	int i = 1;

	FILE* dat = Files::xfopen(snapshot_save_file.c_str(), "wb");

	Files::fwrite_big((void*) & (P.PopSize), sizeof(int), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.NumHouseholds), sizeof(int), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.NumCells), sizeof(int), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.NumPopulatedCells), sizeof(int), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.ncw), sizeof(int), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.nch), sizeof(int), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.setupSeed1), sizeof(int32_t), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.setupSeed2), sizeof(int32_t), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.ModelTimeStep), sizeof(double), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (P.SnapshotSaveTime), sizeof(double), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (State.CellMemberArray), sizeof(int*), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*) & (State.CellSuscMemberArray), sizeof(int*), 1, dat);
	Files::xfprintf_stderr("## %i\n", i++);

	Files::fwrite_big((void*)Hosts, sizeof(Person), (size_t)P.PopSize, dat);

	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*)Households, sizeof(Household), (size_t)P.NumHouseholds, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*)Cells, sizeof(Cell), (size_t)P.NumCells, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*)Mcells, sizeof(Microcell), (size_t)P.NumMicrocells, dat);
	Files::xfprintf_stderr("## %i\n", i++);

	Files::fwrite_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	Files::xfprintf_stderr("## %i\n", i++);
	Files::fwrite_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	Files::xfprintf_stderr("## %i\n", i++);

	Files::xfclose(dat);
}

void UpdateProbs(int DoPlace)
{
	if (!DoPlace)
	{
#pragma omp parallel for schedule(static,500) default(none) \
			shared(P, CellLookup)
		for (int j = 0; j < P.NumPopulatedCells; j++)
		{
			CellLookup[j]->tot_prob = 0;
			CellLookup[j]->S0 = CellLookup[j]->S + CellLookup[j]->L + CellLookup[j]->I;
			if (P.DoDeath)
			{
				CellLookup[j]->S0 += CellLookup[j]->n / 5;
				if ((CellLookup[j]->n < 100) || (CellLookup[j]->S0 > CellLookup[j]->n)) CellLookup[j]->S0 = CellLookup[j]->n;
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(static,500) default(none) \
			shared(P, CellLookup)
		for (int j = 0; j < P.NumPopulatedCells; j++)
		{
			CellLookup[j]->S0 = CellLookup[j]->S;
			CellLookup[j]->tot_prob = 0;
		}
	}
#pragma omp parallel for schedule(static,500) default(none) \
		shared(P, CellLookup)
	for (int j = 0; j < P.NumPopulatedCells; j++)
	{
		int m, k;
		float t;
		CellLookup[j]->cum_trans[0] = ((float)(CellLookup[0]->S0)) * CellLookup[j]->max_trans[0];
		t = ((float)CellLookup[0]->n) * CellLookup[j]->max_trans[0];
		for (m = 1; m < P.NumPopulatedCells; m++)
		{
				CellLookup[j]->cum_trans[m] = CellLookup[j]->cum_trans[m - 1] + ((float)(CellLookup[m]->S0)) * CellLookup[j]->max_trans[m];
				t += ((float)CellLookup[m]->n) * CellLookup[j]->max_trans[m];
		}
		CellLookup[j]->tot_prob = CellLookup[j]->cum_trans[P.NumPopulatedCells - 1];
		for (m = 0; m < P.NumPopulatedCells; m++)
			CellLookup[j]->cum_trans[m] /= CellLookup[j]->tot_prob;
		CellLookup[j]->tot_prob /= t;
		for (k = m = 0; k <= 1024; k++)
		{
			while (CellLookup[j]->cum_trans[m] * 1024 < ((float)k)) m++;
			CellLookup[j]->InvCDF[k] = m;
		}
	}
}

int ChooseTriggerVariableAndValue(int AdUnit)
{
	int VariableAndValue = 0;
	if (P.DoGlobalTriggers)
	{
		if (P.DoPerCapitaTriggers)
			VariableAndValue = (int)floor(((double)State.trigDetectedCases) * P.GlobalIncThreshPop / ((double)P.PopSize));
		else
			VariableAndValue = State.trigDetectedCases;
	}
	else if (P.DoAdminTriggers) VariableAndValue = State.trigDC_adunit[AdUnit];
	else VariableAndValue = INT32_MAX; //// i.e. if not doing triggering (at either admin or global level) then set value to be arbitrarily large so that it will surpass any trigger threshold. Probably other ways around this if anybody wants to correct?

	return VariableAndValue;
}
double ChooseThreshold(int AdUnit, double WhichThreshold) //// point is that this threshold needs to be generalised, so this is likely insufficient.
{
	double Threshold = 0;
	if (P.DoGlobalTriggers) Threshold = WhichThreshold;
	else if (P.DoAdminTriggers)
	{
		if (P.DoPerCapitaTriggers)
			Threshold = (int)ceil(((double)(AdUnits[AdUnit].n * WhichThreshold)) / P.IncThreshPop);
		else
			Threshold = WhichThreshold;
	}
	return Threshold;
}

void UpdateCFRs(double t_CalTime)
{
	for (int ChangeTime = 0; ChangeTime < P.Num_CFR_ChangeTimes; ChangeTime++)
		if (t_CalTime == P.CFR_ChangeTimes_CalTime[ChangeTime])
		{
			P.CFR_Critical_Scale_Current	= P.CFR_TimeScaling_Critical[ChangeTime];
			P.CFR_SARI_Scale_Current		= P.CFR_TimeScaling_SARI	[ChangeTime];
			P.CFR_ILI_Scale_Current			= P.CFR_TimeScaling_ILI		[ChangeTime];
		}
}


void UpdateCurrentInterventionParams(double t_CalTime)
{
	//// **** social distancing
	for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++)
		if (t_CalTime == P.SD_ChangeTimes[ChangeTime])
		{
			//// **** non-enhanced
			P.SocDistHouseholdEffectCurrent = P.SD_HouseholdEffects_OverTime[ChangeTime];	//// household
			P.SocDistSpatialEffectCurrent	= P.SD_SpatialEffects_OverTime	[ChangeTime];	//// spatial
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
				P.SocDistPlaceEffectCurrent[PlaceType] = P.SD_PlaceEffects_OverTime[ChangeTime][PlaceType]; ///// place

			//// **** enhanced
			P.EnhancedSocDistHouseholdEffectCurrent = P.Enhanced_SD_HouseholdEffects_OverTime	[ChangeTime];	//// household
			P.EnhancedSocDistSpatialEffectCurrent	= P.Enhanced_SD_SpatialEffects_OverTime		[ChangeTime];	//// spatial
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
				P.EnhancedSocDistPlaceEffectCurrent[PlaceType] = P.Enhanced_SD_PlaceEffects_OverTime[ChangeTime][PlaceType]; ///// place

			P.SocDistCellIncThresh = P.SD_CellIncThresh_OverTime[ChangeTime];				//// cell incidence threshold
		}

	//// **** case isolation
	for (int ChangeTime = 0; ChangeTime < P.Num_CI_ChangeTimes; ChangeTime++)
		if (t_CalTime == P.CI_ChangeTimes[ChangeTime])
		{
			P.CaseIsolationEffectiveness		= P.CI_SpatialAndPlaceEffects_OverTime	[ChangeTime]; //// spatial / place
			P.CaseIsolationHouseEffectiveness	= P.CI_HouseholdEffects_OverTime		[ChangeTime]; //// household

			P.CaseIsolationProp					= P.CI_Prop_OverTime					[ChangeTime]; //// compliance
			P.CaseIsolation_CellIncThresh		= P.CI_CellIncThresh_OverTime			[ChangeTime]; //// cell incidence threshold
		}

	////// **** household quarantine
	if (P.DoHouseholds)
		for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++)
			if (t_CalTime == P.HQ_ChangeTimes[ChangeTime])
			{
				P.HQuarantineSpatialEffect	= P.HQ_SpatialEffects_OverTime				[ChangeTime];				//// spatial
				P.HQuarantineHouseEffect 	= P.HQ_HouseholdEffects_OverTime			[ChangeTime];				//// household
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.HQuarantinePlaceEffect[PlaceType] = P.HQ_PlaceEffects_OverTime	[ChangeTime][PlaceType];	//// place

				P.HQuarantinePropIndivCompliant = P.HQ_Individual_PropComply_OverTime	[ChangeTime]; //// individual compliance
				P.HQuarantinePropHouseCompliant = P.HQ_Household_PropComply_OverTime	[ChangeTime]; //// household compliance

				P.HHQuar_CellIncThresh			= P.HQ_CellIncThresh_OverTime			[ChangeTime]; //// cell incidence threshold
			}

	//// **** place closure
	if (P.DoPlaces)
	{
		for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++)
			if (t_CalTime == P.PC_ChangeTimes[ChangeTime])
			{
				//// First open all the places - keep commented out in case becomes necessary but avoid if possible to avoid runtime costs.
//				unsigned short int ts = (unsigned short int) (P.TimeStepsPerDay * t);
//				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
//#pragma omp parallel for schedule(static,1)
//					for (int ThreadNum = 0; ThreadNum < P.NumThreads; ThreadNum++)
//						for (int PlaceNum = ThreadNum; PlaceNum < P.Nplace[PlaceType]; PlaceNum += P.NumThreads)
//							DoPlaceOpen(PlaceType, PlaceNum, ts, ThreadNum);

				P.PlaceCloseSpatialRelContact	= P.PC_SpatialEffects_OverTime	[ChangeTime];					//// spatial
				P.PlaceCloseHouseholdRelContact = P.PC_HouseholdEffects_OverTime[ChangeTime];					//// household
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
				{
					P.PlaceCloseEffect[PlaceType] = P.PC_PlaceEffects_OverTime[ChangeTime][PlaceType];			//// place
					P.PlaceClosePropAttending[PlaceType] = P.PC_PropAttending_OverTime[ChangeTime][PlaceType];	//// place
				}

				P.PlaceCloseIncTrig				= P.PC_IncThresh_OverTime		[ChangeTime];					//// global incidence threshold
				P.PlaceCloseFracIncTrig			= P.PC_FracIncThresh_OverTime	[ChangeTime];					//// fractional incidence threshold
				P.PlaceCloseCellIncThresh		= P.PC_CellIncThresh_OverTime	[ChangeTime];					//// cell incidence threshold
				P.PlaceCloseDuration			= P.PC_Durs_OverTime			[ChangeTime];					//// duration of place closure

				//// reset place close time start - has been set to 9e9 in event of no triggers. m
				if(P.PlaceCloseTimeStart < 1e10) P.PlaceCloseTimeStart = t_CalTime;

				// ensure that new duration doesn't go over next change time. Judgement call here - talk to Neil if this is what he wants.
				if ((ChangeTime < P.Num_PC_ChangeTimes - 1) && (P.PlaceCloseTimeStart + P.PlaceCloseDuration >= P.PC_ChangeTimes[ChangeTime + 1]))
					P.PlaceCloseDuration = P.PC_ChangeTimes[ChangeTime + 1] - P.PC_ChangeTimes[ChangeTime]; // -1;
				//Files::xfprintf_stderr("\nt=%lf, n=%i (%i)  PlaceCloseDuration = %lf  (%lf) \n", t, ChangeTime, P.Num_PC_ChangeTimes, P.PlaceCloseDuration, P.PC_ChangeTimes[ChangeTime+1]);
			}
	}

	//// **** digital contact tracing
	for (int ChangeTime = 0; ChangeTime < P.Num_DCT_ChangeTimes; ChangeTime++)
		if (t_CalTime == P.DCT_ChangeTimes[ChangeTime])
		{
			P.DCTCaseIsolationEffectiveness			= P.DCT_SpatialAndPlaceEffects_OverTime	[ChangeTime];	//// spatial / place
			P.DCTCaseIsolationHouseEffectiveness	= P.DCT_HouseholdEffects_OverTime		[ChangeTime];	//// household
			P.ProportionDigitalContactsIsolate		= P.DCT_Prop_OverTime					[ChangeTime];	//// compliance
			P.MaxDigitalContactsToTrace				= P.DCT_MaxToTrace_OverTime				[ChangeTime];
		}
}

void RecordAdminAgeBreakdowns(int t_int)
{
	//// **** Infections by age and admin unit (can parallelize later)
	for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
		{
			// Record incidence. Need new total minus old total (same as minus old total plus new total).
			// First subtract old total before reset, collated and updated.
			TimeSeries[t_int].incInf_age_adunit[AgeGroup][AdUnit] = (double)(-State.cumInf_age_adunit[AgeGroup][AdUnit]);

			// reset totals
			State.prevInf_age_adunit[AgeGroup][AdUnit] = 0;
			State.cumInf_age_adunit [AgeGroup][AdUnit] = 0;

			// sum totals
			for (int Thread = 0; Thread < P.NumThreads; Thread++)
			{
				State.prevInf_age_adunit[AgeGroup][AdUnit] += StateT[Thread].prevInf_age_adunit[AgeGroup][AdUnit];
				State.cumInf_age_adunit [AgeGroup][AdUnit] += StateT[Thread].cumInf_age_adunit [AgeGroup][AdUnit];
			}

			// record in time series.
			TimeSeries[t_int].prevInf_age_adunit[AgeGroup][AdUnit] = State.prevInf_age_adunit[AgeGroup][AdUnit];
			TimeSeries[t_int].cumInf_age_adunit [AgeGroup][AdUnit] = State.cumInf_age_adunit [AgeGroup][AdUnit];

			// Record incidence. Need new total minus old total. Add new total
			TimeSeries[t_int].incInf_age_adunit[AgeGroup][AdUnit] += (double)(State.cumInf_age_adunit[AgeGroup][AdUnit]);
		}
}

void RecordQuarNotInfected(int n, unsigned short int ts)
{
	int QuarNotInfected = 0, QuarNotSymptomatic = 0;
#pragma omp parallel for schedule(static,1) reduction(+:QuarNotInfected, QuarNotSymptomatic)
	for (int thread_no = 0; thread_no < P.NumThreads; thread_no++)
		for (int Person = thread_no; Person < P.PopSize; Person += P.NumThreads)
			if (HOST_QUARANTINED(Person))
			{
				if (Hosts[Person].is_susceptible() || Hosts[Person].is_recovered()) QuarNotInfected++;
				if (Hosts[Person].is_never_symptomatic()) QuarNotSymptomatic++;
			}

	TimeSeries[n].prevQuarNotInfected		= (double) QuarNotInfected;
	TimeSeries[n].prevQuarNotSymptomatic	= (double) QuarNotSymptomatic;
}

void RecordSample(double t, int n, std::string const& output_file_base)
{
	int j, k, S = 0, L = 0, I = 0, R = 0, D = 0, N, cumC = 0, cumTC = 0, cumI = 0, cumR, cumD = 0, cumDC = 0, cumFC = 0, cumTG = 0, cumSI = 0, nTG = 0;
	int cumCT = 0; //added cumulative number of contact traced: ggilani 15/06/17
	int cumCC = 0; //added cumulative number of cases who are contacts: ggilani 28/05/2019
	int cumDCT = 0; //added cumulative number of cases who are digitally contact traced: ggilani 11/03/20
	int cumHQ = 0, cumAC = 0, cumAH = 0, cumAA = 0, cumACS = 0, cumAPC = 0, cumAPA = 0, cumAPCS = 0, numPC, trigDetectedCases;
	int cumC_country[MAX_COUNTRIES]; //add cumulative cases per country
	unsigned short int ts = (unsigned short int) (P.TimeStepsPerDay * t);

	//// Severity quantities
	int Mild = 0, ILI = 0, SARI = 0, Critical = 0, CritRecov = 0, cumMild = 0, cumILI = 0, cumSARI = 0, cumCritical = 0;
	int cumCritRecov = 0, cumDeath_ILI = 0, cumDeath_SARI = 0, cumDeath_Critical = 0;

	//// initialize to zero

	for (int i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] = 0;
	
#pragma omp parallel for schedule(static,10000) reduction(+:S,L,I,R,D,cumTC) default(none) \
		shared(P, CellLookup)
	for (int i = 0; i < P.NumPopulatedCells; i++)
	{
		Cell* ct = CellLookup[i];
		S += (int)ct->S;
		L += (int)ct->L;
		I += (int)ct->I;
		R += (int)ct->R;
		D += (int)ct->D;
		cumTC += (int)ct->cumTC;
	}
	cumR = R;
	cumD = D;
	//cumD = 0;
	N = S + L + I + R + D;
	if (N != P.PopSize) Files::xfprintf_stderr("## %i #\n", P.PopSize - N);
	State.sumRad2 = 0;
	for (j = 0; j < P.NumThreads; j++)
	{
		cumI += StateT[j].cumI;
		cumC += StateT[j].cumC;
		cumDC += StateT[j].cumDC;
		cumTG += StateT[j].cumTG;
		cumSI += StateT[j].cumSI;
		nTG += StateT[j].nTG;
		StateT[j].cumTG = StateT[j].cumSI = StateT[j].nTG = 0;
		cumFC += StateT[j].cumFC;
		cumCT += StateT[j].cumCT; //added contact tracing
		cumCC += StateT[j].cumCC; //added cases who are contacts
		cumDCT += StateT[j].cumDCT; //added cases who are digitally contact traced
		State.sumRad2 += StateT[j].sumRad2;
		State.sumRad2 += StateT[j].sumRad2;
		cumHQ += StateT[j].cumHQ;
		cumAC += StateT[j].cumAC;
		cumAA += StateT[j].cumAA;
		cumAPC += StateT[j].cumAPC;
		cumAPA += StateT[j].cumAPA;
		cumAPCS += StateT[j].cumAPCS;
		cumAH += StateT[j].cumAH;
		cumACS += StateT[j].cumACS;
		//cumD += StateT[j].cumD;

		if (P.DoSeverity)
		{
			///// severity states by thread
			Mild += StateT[j].Mild;
			ILI += StateT[j].ILI;
			SARI += StateT[j].SARI;
			Critical += StateT[j].Critical;
			CritRecov += StateT[j].CritRecov;

			///// cumulative severity states by thread
			cumMild += StateT[j].cumMild;
			cumILI += StateT[j].cumILI;
			cumSARI += StateT[j].cumSARI;
			cumCritical += StateT[j].cumCritical;
			cumCritRecov += StateT[j].cumCritRecov;
			cumDeath_ILI += StateT[j].cumDeath_ILI;
			cumDeath_SARI += StateT[j].cumDeath_SARI;
			cumDeath_Critical += StateT[j].cumDeath_Critical;
		}

		//add up cumulative country counts: ggilani - 12/11/14
		for (int i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] += StateT[j].cumC_country[i];
		if (State.maxRad2 < StateT[j].maxRad2) State.maxRad2 = StateT[j].maxRad2;
	}
	for (j = 0; j < P.NumThreads; j++)
		StateT[j].maxRad2 = State.maxRad2;
	TimeSeries[n].t = t;
	TimeSeries[n].S = (double)S;
	TimeSeries[n].L = (double)L;
	TimeSeries[n].I = (double)I;
	TimeSeries[n].R = (double)R;
	TimeSeries[n].D = (double)D;
	TimeSeries[n].incI = (double)(_I64(cumI) - State.cumI);
	TimeSeries[n].incC = (double)(_I64(cumC) - State.cumC);
	TimeSeries[n].incFC = (double)(_I64(cumFC) - State.cumFC);
	TimeSeries[n].incCT = (double)(_I64(cumCT) - State.cumCT); // added contact tracing
	TimeSeries[n].incCC = (double)(_I64(cumCC) - State.cumCC); // added cases who are contacts
	TimeSeries[n].incDCT = (double)(_I64(cumDCT) - State.cumDCT); //added cases who are digitally contact traced
	TimeSeries[n].incDC = (double)(_I64(cumDC) - State.cumDC); //added incidence of detected cases
	TimeSeries[n].incTC = (double)(_I64(cumTC) - State.cumTC);
	TimeSeries[n].incR = (double)(_I64(cumR) - State.cumR);
	TimeSeries[n].incD = (double)(_I64(cumD) - State.cumD);
	TimeSeries[n].incHQ = (double)(_I64(cumHQ) - State.cumHQ);
	TimeSeries[n].incAC = (double)(_I64(cumAC) - State.cumAC);
	TimeSeries[n].incAH = (double)(_I64(cumAH) - State.cumAH);
	TimeSeries[n].incAA = (double)(_I64(cumAA) - State.cumAA);
	TimeSeries[n].incACS = (double)(_I64(cumACS) - State.cumACS);
	TimeSeries[n].incAPC = (double)(_I64(cumAPC) - State.cumAPC);
	TimeSeries[n].incAPA = (double)(_I64(cumAPA) - State.cumAPA);
	TimeSeries[n].incAPCS = (double)(_I64(cumAPCS) - State.cumAPCS);
	TimeSeries[n].cumT = State.cumT;
	TimeSeries[n].cumUT = State.cumUT;
	TimeSeries[n].cumTP = State.cumTP;
	TimeSeries[n].cumV = State.cumV;
	TimeSeries[n].cumVG = State.cumVG; //added VG;
	TimeSeries[n].cumDC = cumDC;
	TimeSeries[n].meanTG = TimeSeries[n].meanSI = 0;
	if (nTG > 0) // record mean generation time and serial interval in timestep
	{
		TimeSeries[n].meanTG = P.ModelTimeStep * ((double)cumTG) / ((double)nTG);
		TimeSeries[n].meanSI = P.ModelTimeStep * ((double)cumSI) / ((double)nTG);
	}
	//Files::xfprintf_stderr("\ncumD=%i last_cumD=%i incD=%lg\n ", cumD, State.cumD, TimeSeries[n].incD);
	//incidence per country
	for (int i = 0; i < MAX_COUNTRIES; i++) TimeSeries[n].incC_country[i] = (double)(_I64(cumC_country[i]) - State.cumC_country[i]);
	if (P.DoICUTriggers)
	{
		trigDetectedCases = cumCritical;
		if (n >= P.TriggersSamplingInterval) trigDetectedCases -= (int)TimeSeries[n - P.TriggersSamplingInterval].cumCritical;
	}
	else
	{
		trigDetectedCases = cumDC;
		if (n >= P.TriggersSamplingInterval) trigDetectedCases -= (int)TimeSeries[n - P.TriggersSamplingInterval].cumDC;
	}
	State.trigDetectedCases = trigDetectedCases;

	//// update State with new totals from threads.
	State.S = S;
	State.L = L;
	State.I = I;
	State.R = R;
	State.D = D;
	State.cumI = cumI;
	State.cumDC = cumDC;
	State.cumTC = cumTC;
	State.cumFC = cumFC;
	State.cumCT = cumCT; //added cumulative contact tracing
	State.cumCC = cumCC; //added cumulative cases who are contacts
	State.cumDCT = cumDCT; //added cumulative cases who are digitally contact traced
	State.cumC = cumC;
	State.cumR = cumR;
	State.cumD = cumD;
	State.cumHQ = cumHQ;
	State.cumAC = cumAC;
	State.cumAH = cumAH;
	State.cumAA = cumAA;
	State.cumACS = cumACS;
	State.cumAPC = cumAPC;
	State.cumAPA = cumAPA;
	State.cumAPCS = cumAPCS;

	//// **** Infections by age and admin unit (can parallelize later)
	if (P.DoAdUnits && P.OutputAdUnitAge)
		RecordAdminAgeBreakdowns(n);

	RecordQuarNotInfected(n, ts);

	if (P.DoSeverity)
	{

		//// Record incidence. (Must be done with old State totals)
		TimeSeries[n].incMild = (double)(_I64(cumMild) - State.cumMild);
		TimeSeries[n].incILI = (double)(_I64(cumILI) - State.cumILI);
		TimeSeries[n].incSARI = (double)(_I64(cumSARI) - State.cumSARI);
		TimeSeries[n].incCritical = (double)(_I64(cumCritical) - State.cumCritical);
		TimeSeries[n].incCritRecov = (double)(_I64(cumCritRecov) - State.cumCritRecov);
		TimeSeries[n].incDeath_ILI = (double)(_I64(cumDeath_ILI) - State.cumDeath_ILI);
		TimeSeries[n].incDeath_SARI = (double)(_I64(cumDeath_SARI) - State.cumDeath_SARI);
		TimeSeries[n].incDeath_Critical = (double)(_I64(cumDeath_Critical) - State.cumDeath_Critical);

		/////// update state with totals
		State.Mild = Mild;
		State.ILI = ILI;
		State.SARI = SARI;
		State.Critical = Critical;
		State.CritRecov = CritRecov;
		State.cumMild = cumMild;
		State.cumILI = cumILI;
		State.cumSARI = cumSARI;
		State.cumCritical = cumCritical;
		State.cumCritRecov = cumCritRecov;
		State.cumDeath_ILI = cumDeath_ILI;
		State.cumDeath_SARI = cumDeath_SARI;
		State.cumDeath_Critical = cumDeath_Critical;

		//// Record new totals for time series. (Must be done with old State totals)
		TimeSeries[n].Mild = Mild;
		TimeSeries[n].ILI = ILI;
		TimeSeries[n].SARI = SARI;
		TimeSeries[n].Critical = Critical;
		TimeSeries[n].CritRecov = CritRecov;
		TimeSeries[n].cumMild = cumMild;
		TimeSeries[n].cumILI = cumILI;
		TimeSeries[n].cumSARI = cumSARI;
		TimeSeries[n].cumCritical = cumCritical;
		TimeSeries[n].cumCritRecov = cumCritRecov;
		TimeSeries[n].cumDeath_ILI = cumDeath_ILI;
		TimeSeries[n].cumDeath_SARI = cumDeath_SARI;
		TimeSeries[n].cumDeath_Critical = cumDeath_Critical;

		for (int i = 0; i < NUM_AGE_GROUPS; i++)
		{
			//// Record incidence. Need new total minus old total (same as minus old total plus new total).
			//// First subtract old total while unchanged.
			TimeSeries[n].incMild_age[i] = (double)(-State.cumMild_age[i]);
			TimeSeries[n].incILI_age[i] = (double)(-State.cumILI_age[i]);
			TimeSeries[n].incSARI_age[i] = (double)(-State.cumSARI_age[i]);
			TimeSeries[n].incCritical_age[i] = (double)(-State.cumCritical_age[i]);
			TimeSeries[n].incCritRecov_age[i] = (double)(-State.cumCritRecov_age[i]);
			TimeSeries[n].incDeath_ILI_age[i] = (double)(-State.cumDeath_ILI_age[i]);
			TimeSeries[n].incDeath_SARI_age[i] = (double)(-State.cumDeath_SARI_age[i]);
			TimeSeries[n].incDeath_Critical_age[i] = (double)(-State.cumDeath_Critical_age[i]);

			State.Mild_age[i] = 0;
			State.ILI_age[i] = 0;
			State.SARI_age[i] = 0;
			State.Critical_age[i] = 0;
			State.CritRecov_age[i] = 0;
			State.cumMild_age[i] = 0;
			State.cumILI_age[i] = 0;
			State.cumSARI_age[i] = 0;
			State.cumCritical_age[i] = 0;
			State.cumCritRecov_age[i] = 0;
			State.cumDeath_ILI_age[i] = 0;
			State.cumDeath_SARI_age[i] = 0;
			State.cumDeath_Critical_age[i] = 0;

			for (j = 0; j < P.NumThreads; j++)
			{
				//// collate from threads
				State.Mild_age[i] += StateT[j].Mild_age[i];
				State.ILI_age[i] += StateT[j].ILI_age[i];
				State.SARI_age[i] += StateT[j].SARI_age[i];
				State.Critical_age[i] += StateT[j].Critical_age[i];
				State.CritRecov_age[i] += StateT[j].CritRecov_age[i];
				State.cumMild_age[i] += StateT[j].cumMild_age[i];
				State.cumILI_age[i] += StateT[j].cumILI_age[i];
				State.cumSARI_age[i] += StateT[j].cumSARI_age[i];
				State.cumCritical_age[i] += StateT[j].cumCritical_age[i];
				State.cumCritRecov_age[i] += StateT[j].cumCritRecov_age[i];
				State.cumDeath_ILI_age[i] += StateT[j].cumDeath_ILI_age[i];
				State.cumDeath_SARI_age[i] += StateT[j].cumDeath_SARI_age[i];
				State.cumDeath_Critical_age[i] += StateT[j].cumDeath_Critical_age[i];
			}

			//// Record incidence. Need new total minus old total. Add new total
			TimeSeries[n].incMild_age[i] += (double)(State.cumMild_age[i]);
			TimeSeries[n].incILI_age[i] += (double)(State.cumILI_age[i]);
			TimeSeries[n].incSARI_age[i] += (double)(State.cumSARI_age[i]);
			TimeSeries[n].incCritical_age[i] += (double)(State.cumCritical_age[i]);
			TimeSeries[n].incCritRecov_age[i] += (double)(State.cumCritRecov_age[i]);
			TimeSeries[n].incDeath_ILI_age[i] += (double)(State.cumDeath_ILI_age[i]);
			TimeSeries[n].incDeath_SARI_age[i] += (double)(State.cumDeath_SARI_age[i]);
			TimeSeries[n].incDeath_Critical_age[i] += (double)(State.cumDeath_Critical_age[i]);

			//// Record new totals for time series. (Must be done with old State totals)
			TimeSeries[n].Mild_age[i] = State.Mild_age[i];
			TimeSeries[n].ILI_age[i] = State.ILI_age[i];
			TimeSeries[n].SARI_age[i] = State.SARI_age[i];
			TimeSeries[n].Critical_age[i] = State.Critical_age[i];
			TimeSeries[n].CritRecov_age[i] = State.CritRecov_age[i];
			TimeSeries[n].cumMild_age[i] = State.cumMild_age[i];
			TimeSeries[n].cumILI_age[i] = State.cumILI_age[i];
			TimeSeries[n].cumSARI_age[i] = State.cumSARI_age[i];
			TimeSeries[n].cumCritical_age[i] = State.cumCritical_age[i];
			TimeSeries[n].cumCritRecov_age[i] = State.cumCritRecov_age[i];
			TimeSeries[n].cumDeath_ILI_age[i] = State.cumDeath_ILI_age[i];
			TimeSeries[n].cumDeath_SARI_age[i] = State.cumDeath_SARI_age[i];
			TimeSeries[n].cumDeath_Critical_age[i] = State.cumDeath_Critical_age[i];
		}
		if (P.DoAdUnits)
			for (int i = 0; i <= P.NumAdunits; i++)
			{
				//// Record incidence. Need new total minus old total (same as minus old total plus new total).
				//// First subtract old total while unchanged.
				TimeSeries[n].incMild_adunit[i] = (double)(-State.cumMild_adunit[i]);
				TimeSeries[n].incILI_adunit[i] = (double)(-State.cumILI_adunit[i]);
				TimeSeries[n].incSARI_adunit[i] = (double)(-State.cumSARI_adunit[i]);
				TimeSeries[n].incCritical_adunit[i] = (double)(-State.cumCritical_adunit[i]);
				TimeSeries[n].incCritRecov_adunit[i] = (double)(-State.cumCritRecov_adunit[i]);
				TimeSeries[n].incD_adunit[i] = (double)(-State.cumD_adunit[i]);
				TimeSeries[n].incDeath_ILI_adunit[i] = (double)(-State.cumDeath_ILI_adunit[i]);
				TimeSeries[n].incDeath_SARI_adunit[i] = (double)(-State.cumDeath_SARI_adunit[i]);
				TimeSeries[n].incDeath_Critical_adunit[i] = (double)(-State.cumDeath_Critical_adunit[i]);

				//// reset State (not StateT) to zero. Don't need to do this with non-admin unit as local variables Mild, cumSARI etc. initialized to zero at beginning of function. Check with Gemma
				State.Mild_adunit[i] = 0;
				State.ILI_adunit[i] = 0;
				State.SARI_adunit[i] = 0;
				State.Critical_adunit[i] = 0;
				State.CritRecov_adunit[i] = 0;
				State.cumMild_adunit[i] = 0;
				State.cumILI_adunit[i] = 0;
				State.cumSARI_adunit[i] = 0;
				State.cumCritical_adunit[i] = 0;
				State.cumCritRecov_adunit[i] = 0;
				State.cumD_adunit[i] = 0;
				State.cumDeath_ILI_adunit[i] = 0;
				State.cumDeath_SARI_adunit[i] = 0;
				State.cumDeath_Critical_adunit[i] = 0;

				for (j = 0; j < P.NumThreads; j++)
				{
					//// collate from threads
					State.Mild_adunit[i] += StateT[j].Mild_adunit[i];
					State.ILI_adunit[i] += StateT[j].ILI_adunit[i];
					State.SARI_adunit[i] += StateT[j].SARI_adunit[i];
					State.Critical_adunit[i] += StateT[j].Critical_adunit[i];
					State.CritRecov_adunit[i] += StateT[j].CritRecov_adunit[i];
					State.cumMild_adunit[i] += StateT[j].cumMild_adunit[i];
					State.cumILI_adunit[i] += StateT[j].cumILI_adunit[i];
					State.cumSARI_adunit[i] += StateT[j].cumSARI_adunit[i];
					State.cumCritical_adunit[i] += StateT[j].cumCritical_adunit[i];
					State.cumCritRecov_adunit[i] += StateT[j].cumCritRecov_adunit[i];
					State.cumD_adunit[i] += StateT[j].cumD_adunit[i];
					State.cumDeath_ILI_adunit[i] += StateT[j].cumDeath_ILI_adunit[i];
					State.cumDeath_SARI_adunit[i] += StateT[j].cumDeath_SARI_adunit[i];
					State.cumDeath_Critical_adunit[i] += StateT[j].cumDeath_Critical_adunit[i];
				}

				//// Record incidence. Need new total minus old total. Add new total
				TimeSeries[n].incMild_adunit[i] += (double)(State.cumMild_adunit[i]);
				TimeSeries[n].incILI_adunit[i] += (double)(State.cumILI_adunit[i]);
				TimeSeries[n].incSARI_adunit[i] += (double)(State.cumSARI_adunit[i]);
				TimeSeries[n].incCritical_adunit[i] += (double)(State.cumCritical_adunit[i]);
				TimeSeries[n].incCritRecov_adunit[i] += (double)(State.cumCritRecov_adunit[i]);
				TimeSeries[n].incD_adunit[i] += (double)(State.cumD_adunit[i]);
				TimeSeries[n].incDeath_ILI_adunit[i] += (double)(State.cumDeath_ILI_adunit[i]);
				TimeSeries[n].incDeath_SARI_adunit[i] += (double)(State.cumDeath_SARI_adunit[i]);
				TimeSeries[n].incDeath_Critical_adunit[i] += (double)(State.cumDeath_Critical_adunit[i]);

				//// Record new totals for time series. (Must be done with old State totals)
				TimeSeries[n].Mild_adunit[i] = State.Mild_adunit[i];
				TimeSeries[n].ILI_adunit[i] = State.ILI_adunit[i];
				TimeSeries[n].SARI_adunit[i] = State.SARI_adunit[i];
				TimeSeries[n].Critical_adunit[i] = State.Critical_adunit[i];
				TimeSeries[n].CritRecov_adunit[i] = State.CritRecov_adunit[i];
				TimeSeries[n].cumMild_adunit[i] = State.cumMild_adunit[i];
				TimeSeries[n].cumILI_adunit[i] = State.cumILI_adunit[i];
				TimeSeries[n].cumSARI_adunit[i] = State.cumSARI_adunit[i];
				TimeSeries[n].cumCritical_adunit[i] = State.cumCritical_adunit[i];
				TimeSeries[n].cumCritRecov_adunit[i] = State.cumCritRecov_adunit[i];
				TimeSeries[n].cumD_adunit[i] = State.cumD_adunit[i];
				TimeSeries[n].cumDeath_ILI_adunit[i] = State.cumDeath_ILI_adunit[i];
				TimeSeries[n].cumDeath_SARI_adunit[i] = State.cumDeath_SARI_adunit[i];
				TimeSeries[n].cumDeath_Critical_adunit[i] = State.cumDeath_Critical_adunit[i];
			}
	}

	//update cumulative cases per country
	for (int i = 0; i < MAX_COUNTRIES; i++) State.cumC_country[i] = cumC_country[i];
	//update overall state variable for cumulative cases per adunit

	TimeSeries[n].rmsRad = (State.cumI > 0) ? sqrt(State.sumRad2 / ((double)State.cumI)) : 0;
	TimeSeries[n].maxRad = sqrt(State.maxRad2);
	TimeSeries[n].extinct = ((((P.SmallEpidemicCases >= 0) && (State.R <= P.SmallEpidemicCases)) || (P.SmallEpidemicCases < 0)) && (State.I + State.L == 0)) ? 1 : 0;
	for (int i = 0; i < NUM_AGE_GROUPS; i++)
	{
		TimeSeries[n].incCa[i] = TimeSeries[n].incIa[i] = TimeSeries[n].incDa[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incCa[i] += (double)StateT[j].cumCa[i];
			TimeSeries[n].incIa[i] += (double)StateT[j].cumIa[i];
			TimeSeries[n].incDa[i] += (double)StateT[j].cumDa[i];
		}
	}

	for (int i = 0; i < 2; i++)
	{
		TimeSeries[n].incC_keyworker[i] = TimeSeries[n].incI_keyworker[i] = TimeSeries[n].cumT_keyworker[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incC_keyworker[i] += (double)StateT[j].cumC_keyworker[i];
			TimeSeries[n].incI_keyworker[i] += (double)StateT[j].cumI_keyworker[i];
			TimeSeries[n].cumT_keyworker[i] += (double)StateT[j].cumT_keyworker[i];
			StateT[j].cumC_keyworker[i] = StateT[j].cumI_keyworker[i] = 0;
		}
	}

	for (int i = 0; i < INFECT_TYPE_MASK; i++)
	{
		TimeSeries[n].incItype[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incItype[i] += (double)StateT[j].cumItype[i];
			StateT[j].cumItype[i] = 0;
		}
	}
	if (P.DoAdUnits)
		for (int i = 0; i <= P.NumAdunits; i++)
		{
			TimeSeries[n].incI_adunit[i] = TimeSeries[n].incC_adunit[i] = TimeSeries[n].cumT_adunit[i] = TimeSeries[n].incH_adunit[i] = TimeSeries[n].incDC_adunit[i] = TimeSeries[n].incCT_adunit[i] = TimeSeries[n].incDCT_adunit[i] = 0; //added detected cases: ggilani 03/02/15
			for (j = 0; j < P.NumThreads; j++)
			{
				TimeSeries[n].incI_adunit[i] += (double)StateT[j].cumI_adunit[i];
				TimeSeries[n].incC_adunit[i] += (double)StateT[j].cumC_adunit[i];
				TimeSeries[n].incDC_adunit[i] += (double)StateT[j].cumDC_adunit[i]; //added detected cases: ggilani 03/02/15
				TimeSeries[n].incH_adunit[i] += (double)StateT[j].cumH_adunit[i]; //added hospitalisation
				TimeSeries[n].incCT_adunit[i] += (double)StateT[j].cumCT_adunit[i]; //added contact tracing: ggilani 15/06/17
				TimeSeries[n].incCC_adunit[i] += (double)StateT[j].cumCC_adunit[i]; //added cases who are contacts: ggilani 28/05/2019
				TimeSeries[n].incDCT_adunit[i] += (double)StateT[j].cumDCT_adunit[i]; //added cases who are digitally contact traced: ggilani 11/03/20
				TimeSeries[n].cumT_adunit[i] += (double)StateT[j].cumT_adunit[i];
				State.cumC_adunit[i] += StateT[j].cumC_adunit[i];
				State.cumDC_adunit[i] += StateT[j].cumDC_adunit[i];
				StateT[j].cumI_adunit[i] = StateT[j].cumC_adunit[i] = StateT[j].cumH_adunit[i] = StateT[j].cumDC_adunit[i] = StateT[j].cumCT_adunit[i] = StateT[j].cumCC_adunit[i] = StateT[j].cumDCT_adunit[i] = 0; //added hospitalisation, detected cases, contact tracing: ggilani 03/02/15, 15/06/17
			}
			if (P.DoICUTriggers)
			{
				State.trigDC_adunit[i] += (int)TimeSeries[n].incCritical_adunit[i];
				if (n >= P.TriggersSamplingInterval) State.trigDC_adunit[i] -= (int)TimeSeries[n - P.TriggersSamplingInterval].incCritical_adunit[i];
			}
			else
			{
				State.trigDC_adunit[i] += (int)TimeSeries[n].incDC_adunit[i];
				if (n >= P.TriggersSamplingInterval) State.trigDC_adunit[i] -= (int)TimeSeries[n - P.TriggersSamplingInterval].incDC_adunit[i];
			}
		}
	if (P.DoDigitalContactTracing)
		for (int i = 0; i < P.NumAdunits; i++)
			TimeSeries[n].DCT_adunit[i] = (double)AdUnits[i].ndct; //added total numbers of contacts currently isolated due to digital contact tracing: ggilani 11/03/20
	if (P.DoPlaces)
		for (int i = 0; i < NUM_PLACE_TYPES; i++)
		{
			numPC = 0;
			for (j = 0; j < P.Nplace[i]; j++)
				if (PLACE_CLOSED(i, j)) numPC++;
			State.NumPlacesClosed[i] = numPC;
			TimeSeries[n].PropPlacesClosed[i] = ((double)numPC) / ((double)P.Nplace[i]);
		}
	for (int i = k = 0; i < P.NumMicrocells; i++) if (Mcells[i].socdist == TreatStat::Treated) k++;
	TimeSeries[n].PropSocDist = ((double)k) / ((double)P.NumMicrocells);

	//update contact number distribution in State
	for (int i = 0; i < (MAX_CONTACTS + 1); i++)
	{
		for (j = 0; j < P.NumThreads; j++)
		{
			State.contact_dist[i] += StateT[j].contact_dist[i];
			StateT[j].contact_dist[i] = 0;
		}
	}
	if (P.OutputBitmap >= 1)
	{
		TSMean = TSMeanNE; TSVar = TSVarNE;
		CaptureBitmap();
		OutputBitmap(0, output_file_base);
	}
}

void CalibrationThresholdCheck(double t,int n)
{
	int k;
	int trigAlert, trigAlertCases;
	/* Never used
	unsigned short int ts = (unsigned short int) (P.TimeStepsPerDay * t);
	*/

	trigAlertCases = State.cumDC;
	if (n >= P.WindowToEvaluateTriggerAlert) 
		trigAlertCases -= (int)TimeSeries[n - P.WindowToEvaluateTriggerAlert].cumDC;

	if (P.TriggerAlertOnDeaths)  //// if using deaths as trigger (as opposed to detected cases)
	{
		trigAlert = (int)TimeSeries[n].D;
		if (n >= P.WindowToEvaluateTriggerAlert) trigAlert -= (int) TimeSeries[n - P.WindowToEvaluateTriggerAlert].D;
	}
	else
		trigAlert = trigAlertCases;

	double RatioPredictedObserved, DesiredAccuracy; // calibration variables.

	if (((P.DoNoCalibration) && (t >= P.Epidemic_StartDate_CalTime)) ||
		((!P.DoNoCalibration) && (((!P.DoAlertTriggerAfterInterv) && (trigAlert >= P.CaseOrDeathThresholdBeforeAlert)) ||
		((P.DoAlertTriggerAfterInterv) && (((trigAlertCases >= P.CaseOrDeathThresholdBeforeAlert) && (P.ModelCalibIteration < 2))
			|| ((t >= P.Epidemic_StartDate_CalTime) && ((P.ModelCalibIteration >= 2)||(P.InitialInfectionCalTime > 0) )))))))
	{
		if ((!P.DoNoCalibration) && (!P.StopCalibration) && (!InterruptRun))
		{
			if ((P.DateTriggerReached_SimTime == 0) && (P.InitialInfectionCalTime <=0))
			{
				P.Epidemic_StartDate_CalTime = P.DateTriggerReached_SimTime = t; // initialize Epidemic_StartDate_CalTime & DateTriggerReached_SimTime to now (i.e. simulation time that trigger was reached)
				if (P.DateTriggerReached_CalTime >= 0)
				{
					P.HolidaysStartDay_SimTime = P.DateTriggerReached_SimTime - P.Interventions_StartDate_CalTime; /// initialize holiday offset to time difference between DateTriggerReached_SimTime and day of year interventions start.
//					Files::xfprintf_stderr("@@## trigAlertCases=%i P.HolidaysStartDay_SimTime=%lg \n",trigAlertCases, P.HolidaysStartDay_SimTime);
				}
			}
			if ((P.DateTriggerReached_CalTime >= 0) && (!P.DoAlertTriggerAfterInterv) && (P.InitialInfectionCalTime <= 0))
			{
				P.StopCalibration = 1;
				InterruptRun = 1;
			}
			if ((P.DoAlertTriggerAfterInterv) && (t == P.DateTriggerReached_SimTime + P.DateTriggerReached_CalTime - P.Interventions_StartDate_CalTime))
			{
				if ((trigAlert > 0) && (P.ModelCalibIteration < 20))
				{
					RatioPredictedObserved = ((double)trigAlert)/((double)P.AlertTriggerAfterIntervThreshold);
					DesiredAccuracy = 1.1 / sqrt((double)P.AlertTriggerAfterIntervThreshold);
					if (DesiredAccuracy < 0.05) DesiredAccuracy = 0.05;
					Files::xfprintf_stderr("\n** %i %lf %lf | %lg / %lg \t", P.ModelCalibIteration, t, P.DateTriggerReached_SimTime + P.DateTriggerReached_CalTime - P.Interventions_StartDate_CalTime, P.HolidaysStartDay_SimTime, RatioPredictedObserved);
					Files::xfprintf_stderr("| %i %i %i %i -> ", trigAlert, trigAlertCases, P.AlertTriggerAfterIntervThreshold, P.CaseOrDeathThresholdBeforeAlert);

					if (P.InitialInfectionCalTime > 0)
					{
						if ((P.ModelCalibIteration >= 2) &&
							((((RatioPredictedObserved - 1.0) <= DesiredAccuracy) && (RatioPredictedObserved >= 1)) ||
							(((1.0 - RatioPredictedObserved) <= DesiredAccuracy) && (RatioPredictedObserved < 1))))
							P.StopCalibration = 1;
						else if (P.ModelCalibIteration == 1)
							P.SeedingScaling /= pow(RatioPredictedObserved, 0.7);
						else if (P.ModelCalibIteration == 2)
							P.SeedingScaling /= pow(RatioPredictedObserved, 0.6);
						else if (P.ModelCalibIteration > 2)
							P.SeedingScaling /= pow(RatioPredictedObserved, 0.3 + 0.2 * ranf()); // include random number to prevent infinite loops
					}
					else
					{
						if ((P.ModelCalibIteration >= 2) &&
							((((RatioPredictedObserved - 1.0) <= DesiredAccuracy) && (RatioPredictedObserved >= 1)) ||
							(((1.0 - RatioPredictedObserved) <= DesiredAccuracy) && (RatioPredictedObserved < 1))))
							P.StopCalibration = 1;
						else if (P.ModelCalibIteration == 0)
						{
							k = (int)(((double)P.CaseOrDeathThresholdBeforeAlert) / RatioPredictedObserved);
							if (k > 0) P.CaseOrDeathThresholdBeforeAlert = k;
						}
						else if ((P.ModelCalibIteration >= 2) && ((P.ModelCalibIteration) % 3 < 2))
						{
							if (RatioPredictedObserved > 1)
							{
								P.Epidemic_StartDate_CalTime--;
								P.HolidaysStartDay_SimTime--;
							}
							else if (RatioPredictedObserved < 1)
							{
								P.Epidemic_StartDate_CalTime++;
								P.HolidaysStartDay_SimTime++;
							}
						}
						else if ((P.ModelCalibIteration >= 2) && ((P.ModelCalibIteration) % 3 == 2))
						{
							P.SeedingScaling /= pow(RatioPredictedObserved, 0.2 + 0.3 * ranf()); // include random number to prevent loops
						}
					}
					P.ModelCalibIteration++;
					Files::xfprintf_stderr("%i : %lg\n", P.CaseOrDeathThresholdBeforeAlert, P.SeedingScaling);

					if(P.StopCalibration)
						Files::xfprintf_stderr("Calibration ended.\n");
					else
						InterruptRun = 1;
				}
				else
					P.StopCalibration = 1;
			}
		}
		P.ControlPropCasesId = P.PostAlertControlPropCasesId;

		if (P.VaryEfficaciesOverTime)
			UpdateCurrentInterventionParams(t - P.Epidemic_StartDate_CalTime); // t - P.Epidemic_StartDate_CalTime converts simulation time (t) into calendar time. 

// changed to a define for speed (though always likely inlined anyway) and to avoid clang compiler warnings re double alignment
#define DO_OR_DONT_AMEND_START_TIME(X,Y) if(X>=1e10) X=Y;

		//// Set Case isolation start time (by admin unit)
		for (int i = 0; i < P.NumAdunits; i++)
			if (ChooseTriggerVariableAndValue(i) > ChooseThreshold(i, P.CaseIsolation_CellIncThresh)) //// a little wasteful if doing Global trigs as function called more times than necessary, but worth it for much simpler code. Also this function is small portion of runtime.
				DO_OR_DONT_AMEND_START_TIME(AdUnits[i].CaseIsolationTimeStart, t + ((P.DoInterventionDelaysByAdUnit)?AdUnits[i].CaseIsolationDelay: P.CaseIsolationTimeStartBase))
		//// Set Household Quarantine start time (by admin unit)
		for (int i = 0; i < P.NumAdunits; i++)
			if (ChooseTriggerVariableAndValue(i) > ChooseThreshold(i, P.HHQuar_CellIncThresh)) //// a little wasteful if doing Global trigs as function called more times than necessary, but worth it for much simpler code. Also this function is small portion of runtime.
					DO_OR_DONT_AMEND_START_TIME(AdUnits[i].HQuarantineTimeStart, t + ((P.DoInterventionDelaysByAdUnit)?AdUnits[i].HQuarantineDelay: P.HQuarantineTimeStartBase));

		//// Set DigitalContactTracingTimeStart
		if (P.DoDigitalContactTracing)
			for (int i = 0; i < P.NumAdunits; i++)
				if (ChooseTriggerVariableAndValue(i) > ChooseThreshold(i, P.DigitalContactTracing_CellIncThresh)) //// a little wasteful if doing Global trigs as function called more times than necessary, but worth it for much simpler code. Also this function is small portion of runtime.
					DO_OR_DONT_AMEND_START_TIME(AdUnits[i].DigitalContactTracingTimeStart, t + ((P.DoInterventionDelaysByAdUnit)?AdUnits[i].DCTDelay: P.DigitalContactTracingTimeStartBase));

		if (P.DoGlobalTriggers)
		{
			int TriggerValue = ChooseTriggerVariableAndValue(0);
			if (TriggerValue >= ChooseThreshold(0, P.TreatCellIncThresh))
				DO_OR_DONT_AMEND_START_TIME((P.TreatTimeStart), t + P.TreatTimeStartBase);
			if (TriggerValue >= P.VaccCellIncThresh) DO_OR_DONT_AMEND_START_TIME(P.VaccTimeStart, t + P.VaccTimeStartBase);
			if (TriggerValue >= P.SocDistCellIncThresh)
			{
				DO_OR_DONT_AMEND_START_TIME(P.SocDistTimeStart, t + P.SocDistTimeStartBase);
				//added this for admin unit based intervention delays based on a global trigger: ggilani 17/03/20
				if (P.DoInterventionDelaysByAdUnit)
					for (int i = 0; i < P.NumAdunits; i++)
						DO_OR_DONT_AMEND_START_TIME(AdUnits[i].SocialDistanceTimeStart, t + AdUnits[i].SocialDistanceDelay);
			}
			if (TriggerValue >= P.PlaceCloseCellIncThresh)
			{
				DO_OR_DONT_AMEND_START_TIME(P.PlaceCloseTimeStart, t + P.PlaceCloseTimeStartBase);
				if (P.DoInterventionDelaysByAdUnit)
					for (int i = 0; i < P.NumAdunits; i++)
						DO_OR_DONT_AMEND_START_TIME(AdUnits[i].PlaceCloseTimeStart, t + AdUnits[i].PlaceCloseDelay);
			}
			if (TriggerValue >= P.MoveRestrCellIncThresh)
				DO_OR_DONT_AMEND_START_TIME(P.MoveRestrTimeStart, t + P.MoveRestrTimeStartBase);
			if (TriggerValue >= P.KeyWorkerProphCellIncThresh)
				DO_OR_DONT_AMEND_START_TIME(P.KeyWorkerProphTimeStart, t + P.KeyWorkerProphTimeStartBase);
		}
		else
		{
		  DO_OR_DONT_AMEND_START_TIME(P.TreatTimeStart, t + P.TreatTimeStartBase);
			DO_OR_DONT_AMEND_START_TIME(P.VaccTimeStart	, t + P.VaccTimeStartBase);
			DO_OR_DONT_AMEND_START_TIME(P.SocDistTimeStart, t + P.SocDistTimeStartBase);
			DO_OR_DONT_AMEND_START_TIME(P.PlaceCloseTimeStart, t + P.PlaceCloseTimeStartBase);
			DO_OR_DONT_AMEND_START_TIME(P.MoveRestrTimeStart, t + P.MoveRestrTimeStartBase);
			DO_OR_DONT_AMEND_START_TIME(P.KeyWorkerProphTimeStart, t + P.KeyWorkerProphTimeStartBase);
		}
		DO_OR_DONT_AMEND_START_TIME(P.AirportCloseTimeStart, t + P.AirportCloseTimeStartBase);
	}
	if ((P.PlaceCloseIndepThresh > 0) && (((double)State.cumDC) >= P.PlaceCloseIndepThresh))
		DO_OR_DONT_AMEND_START_TIME(P.PlaceCloseTimeStart, t + P.PlaceCloseTimeStartBase);

	if (t > P.SocDistTimeStart + P.SocDistChangeDelay)
	{
		P.SocDistDurationCurrent = P.SocDistDuration2;
		P.SocDistHouseholdEffectCurrent = P.SocDistHouseholdEffect2;
		P.SocDistSpatialEffectCurrent = P.SocDistSpatialEffect2;
		P.EnhancedSocDistHouseholdEffectCurrent = P.EnhancedSocDistHouseholdEffect2;
		P.EnhancedSocDistSpatialEffectCurrent = P.EnhancedSocDistSpatialEffect2;
		for (int i = 0; i < P.PlaceTypeNum; i++)
		{
			P.SocDistPlaceEffectCurrent[i] = P.SocDistPlaceEffect2[i];
			P.EnhancedSocDistPlaceEffectCurrent[i] = P.EnhancedSocDistPlaceEffect2[i];
		}
	}
	//fix to switch off first place closure after P.PlaceCloseDuration has elapsed, if there are no school or cell-based triggers set
	if (t == P.PlaceCloseTimeStart + P.PlaceCloseDuration)
	{
		P.PlaceCloseTimeStartPrevious = P.PlaceCloseTimeStart;
		if ((P.PlaceCloseIncTrig == 0) && (P.PlaceCloseFracIncTrig == 0) && (P.PlaceCloseCellIncThresh == 0)) P.PlaceCloseTimeStart = 9e9;
	}

	if (!P.VaryEfficaciesOverTime)
	{
		if ((P.PlaceCloseTimeStart2 > P.PlaceCloseTimeStartPrevious) && //// if second place closure start time after previous start time AND
			(t >= P.PlaceCloseTimeStartPrevious + P.PlaceCloseDuration) &&	//// if now after previous place closure period has finished AND
			(t >= P.PlaceCloseTimeStartPrevious + P.PlaceCloseTimeStartBase2 - P.PlaceCloseTimeStartBase))	//// if now after previous start time + plus difference between 1st and 2nd base start times
		{
			Files::xfprintf_stderr("\nSecond place closure period (t=%lg)\n", t);
			P.PlaceCloseTimeStartPrevious = P.PlaceCloseTimeStart2 = P.PlaceCloseTimeStart = t;
			P.PlaceCloseDuration = P.PlaceCloseDuration2;
			P.PlaceCloseIncTrig = P.PlaceCloseIncTrig2;
			P.PlaceCloseCellIncThresh = P.PlaceCloseCellIncThresh2;
		}
	}
}

void CalcLikelihood(int run, std::string const& DataFile, std::string const& OutFileBase)
{
	FILE* dat;

	static int DataAlreadyRead = 0, ncols, nrows, * ColTypes; // static i.e. won't be reset upon subsequent calls to CalcLikelihood
	static double** Data, NegBinK, sumL;

	if (!DataAlreadyRead)
	{
		char FieldName[1024];
		dat = Files::xfopen(DataFile.c_str(), "r");
		// Extract numbers of rows and columns, and overdispersion parameter of negative bionomial distribution, from Data file
		Files::xfscanf(dat, 3, "%i %i %lg", &nrows, &ncols, &NegBinK);

		// allocate memory
		ColTypes = (int*) Memory::xcalloc(ncols, sizeof(int));
		Data = (double**) Memory::xcalloc(nrows, sizeof(double *));
		for (int i = 0; i < nrows; i++)
			Data[i] = (double*) Memory::xcalloc(ncols, sizeof(double));

		// cycle through columns assigning an int label to each data/column type in data file. Essentially renaming column names to integers. 
		for (int i = 0; i < ncols; i++)
		{
			ColTypes[i] = -100;
			Files::xfscanf(dat, 1, "%s", FieldName);
			if (!strcmp(FieldName, "day"))
			{
				ColTypes[i] = -1;
				if (i != 0) ERR_CRITICAL("'day' must be first column in data file\n");
			}
			if (!strcmp(FieldName, "all_deaths"))
				ColTypes[i] = 0;
			else if (!strcmp(FieldName, "hospital_deaths"))
				ColTypes[i] = 1;
			else if (!strcmp(FieldName, "care_home_deaths"))
				ColTypes[i] = 2;
			else if (!strcmp(FieldName, "seroprevalence_numerator"))
				ColTypes[i] = 3;
			else if (!strcmp(FieldName, "seroprevalence_denominator"))
			{
				ColTypes[i] = 4;
				if (ColTypes[i - 1] != 3) ERR_CRITICAL("Seroprevalence denominator must be next column after numerator in data file\n");
			}
			else if (!strcmp(FieldName, "infection_prevalence_numerator"))
				ColTypes[i] = 5;
			else if (!strcmp(FieldName, "infection_prevalence_denominator"))
			{
				ColTypes[i] = 6;
				if (ColTypes[i - 1] != 5) ERR_CRITICAL("Infection prevalence denominator must be next column after numerator in data file\n");
			}
		}

		// extract data into Data array.
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				Files::xfscanf(dat, 1, "%lg", &(Data[i][j]));
		Files::xfclose(dat);
		DataAlreadyRead = 1;
	}

	// calculate likelihood function
	double c, LL = 0.0;
	double kp = (P.clP[99] > 0) ? P.clP[99] : NegBinK; // clP[99] reserved for fitting overdispersion. If not positive-definite assign to NegBinK extracted above. 
	c = 1.0; // 1 / ((double)(P.NRactE + P.NRactNE));
	int offset = (P.Interventions_StartDate_CalTime > 0) ? ((int)(P.Interventions_StartDate_CalTime - P.DateTriggerReached_SimTime)) : 0;

	for (int col = 1; col < ncols; col++) /// cycle through columns (different sources of data contributing to likelihood), and add to log likelihood (LL) accordingly. 
	{
		if ((ColTypes[col] >= 0) && (ColTypes[col] <= 2)) // i.e. "all deaths", "hospital deaths", "care home deaths"
		{
			double ModelValueSum = 0.0;
			for (int row = 0; row < nrows; row++)
			{
				int day = (int)Data[row][0]; // day is day of year - directly indexes TimeSeries[]
				if ((Data[row][col] >= -1) && (day < P.NumOutputTimeSteps)) // data is not NA (-ve) and within time range of model run
				{
					double ModelValue;
					if (ColTypes[col] == 0)
						ModelValue = c * TimeSeries[day - offset].incD; // all deaths by date of death
					else if (ColTypes[col] == 1)
						ModelValue = c * (TimeSeries[day - offset].incDeath_Critical + TimeSeries[day - offset].incDeath_SARI); // hospital deaths (SARI and Critical) by date of death
					else if (ColTypes[col] == 2)
						ModelValue = c * TimeSeries[day - offset].incDeath_ILI; // care home deaths (ILI) by date of death
					ModelValueSum += ModelValue;
					if (Data[row][col] >= 0)
					{
						if ((row > 0) && (Data[row - 1][col] == -1)) // cumulative column: -1 means sum column up to first >=0 value
						{
							ModelValue = ModelValueSum;
							ModelValueSum = 0.0;  // reset cumulative sum
						}
						if (NegBinK >= 10000)
							//prob model and data from same underlying poisson
							LL += lgamma(2 * (Data[row][col] + ModelValue) + 1) - lgamma(Data[row][col] + ModelValue + 1) - lgamma(Data[row][col] + 1) - lgamma(ModelValue + 1) - (3 * (Data[row][col] + ModelValue) + 1) * log(2);
						else
						{
							//neg bin LL (NegBinK=1 implies no over-dispersion. >1 implies more)
							double knb = 1.0 + ModelValue / kp;
							double pnb = kp / (1.0 + kp);
							LL += lgamma(Data[row][col] + knb) - lgamma(Data[row][col] + 1) - lgamma(knb) + knb * log(1.0 - pnb) + Data[row][col] * log(pnb);
						}
					}
				}
			}
		}
		else if (ColTypes[col] == 3) // seroprevalence by date of sample
		{
			for (int row = 0; row < nrows; row++)
			{
				int day = (int)Data[row][0]; // day is day of year - directly indexes TimeSeries[]
				if ((Data[row][col] >= 0) && (day < P.NumOutputTimeSteps)) // data is not NA (-ve) and within time range of model run
				{
					double m = Data[row][col]; // numerator
					double N = Data[row][col + 1]; // denominator
					double ModelValue = 0.0;
					for (int k = offset; k < day; k++) // loop over all days of infection up to day of sample
					{
						double prob_seroconvert = P.SeroConvMaxSens * (1.0 - 0.5 * ((exp(-((double)(_I64(day) - k)) * P.SeroConvP1) + 1.0) * exp(-((double)(_I64(day) - k)) * P.SeroConvP2))); // add P1 to P2 to prevent degeneracy
						ModelValue += c * TimeSeries[k - offset].incI * prob_seroconvert;
					}
					ModelValue += c * TimeSeries[day - offset].S * (1.0 - P.SeroConvSpec);
					ModelValue /= ((double)P.PopSize);
					LL += m * log((ModelValue + 1e-20) / (m / N + 1e-20)) + (N - m) * log((1.0 - ModelValue + 1e-20) / (1.0 - m / N + 1e-20));  // subtract saturated likelihood
				}
			}
		}
		else if (ColTypes[col] == 5) // infection prevalence by date of sample
		{
			for (int row = 0; row < nrows; row++)
			{
				int day = (int)Data[row][0]; // day is day of year - directly indexes TimeSeries[]
				if ((Data[row][col] >= 0) && (day < P.NumOutputTimeSteps)) // data is not NA (-ve) and within time range of model run
				{
					double m = Data[row][col]; // numerator
					double N = Data[row][col + 1]; // denominator
					double ModelValue = P.InfPrevSurveyScale * c * TimeSeries[day - offset].I / ((double)P.PopSize);
					LL += m * log(ModelValue + 1e-20) + (N - m) * log(1.0 - ModelValue);
				}
			}
		}
	}
	Files::xfprintf_stderr("Log-likelihood = %lg\n", LL);
	if (run == 0)
		sumL = LL;
	else
	{
		double maxLL = LL;
		if (sumL > maxLL) maxLL = sumL;
		sumL = maxLL + log(exp(sumL - maxLL) + exp(LL - maxLL));
	}

	if (run + 1 == P.NumRealisations) // at final realisation, output log-likelihood
	{
		LL = sumL - log((double)P.NumRealisations);
		std::string TmpFile = OutFileBase + ".ll.tmp";
		std::string OutFile = OutFileBase + ".ll.txt";
		dat = Files::xfopen(TmpFile.c_str(), "w");
		Files::xfprintf(dat, "%i\t%.8lg\n", P.FitIter, LL);
		Files::xfclose(dat);
		Files::xrename(TmpFile.c_str(), OutFile.c_str()); // rename only when file is complete and closed
	}
}

void RecordInfTypes(void)
{
	int i, j, k, l, lc, lc2, b, c, n, i2;
	unsigned int nf;
	double* res, * res_av, * res_var, t, s = 0;

	for (n = 0; n < P.NumOutputTimeSteps; n++)
	{
		for (i = 0; i < INFECT_TYPE_MASK; i++) TimeSeries[n].Rtype[i] = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++) TimeSeries[n].Rage[i] = 0;
		TimeSeries[n].Rdenom = 0;
	}
	for (i = 0; i < INFECT_TYPE_MASK; i++) inftype[i] = 0;
	for (i = 0; i < MAX_COUNTRIES; i++) infcountry[i] = 0;
	for (i = 0; i < MAX_SEC_REC; i++)
		for (j = 0; j < MAX_GEN_REC; j++)
			indivR0[i][j] = 0;
	for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		for (j = 0; j <= MAX_HOUSEHOLD_SIZE; j++)
			inf_household[i][j] = case_household[i][j] = 0;
	for (b = 0; b < P.NumCells; b++)
		if ((Cells[b].S != Cells[b].n) || (Cells[b].R > 0))
			for (c = 0; c < Cells[b].n; c++)
				Hosts[Cells[b].members[c]].listpos = 0;
	//	for(b=0;b<P.NumCells;b++)
	//		if((Cells[b].S!=Cells[b].n)||(Cells[b].R>0))
	{
		j = k = l = lc = lc2 = 0; t = 1e10;
		//			for(c=0;c<Cells[b].n;c++)
		for (i = 0; i < P.PopSize; i++)
		{
			//				i=Cells[b].members[c];
			if (j == 0) j = k = Households[Hosts[i].hh].nh;
			if (!Hosts[i].is_susceptible() && !Hosts[i].is_immune_at_start())
			{
				if (Hosts[i].latent_time * P.ModelTimeStep <= P.SimulationDuration)
					TimeSeries[(int)(Hosts[i].latent_time * P.ModelTimeStep / P.OutputTimeStep)].Rdenom++;
				infcountry[mcell_country[Hosts[i].mcell]]++;
				if (Hosts[i].is_susceptible_or_infected())
				{
					l = -1;
				}
				else if (l >= 0)
				{
					l++;
				}
				if ((l >= 0) && (Hosts[i].is_recovered_symp() || Hosts[i].is_dead_was_symp()))
				{
					lc2++;
					if (Hosts[i].latent_time * P.ModelTimeStep <= t) // This convoluted logic is to pick up households where the index is symptomatic
					{
						lc = 1; t = Hosts[i].latent_time * P.ModelTimeStep;
					}
				}
				else if ((l > 0) && (Hosts[i].latent_time * P.ModelTimeStep < t))
				{
					lc = 0; t = Hosts[i].latent_time * P.ModelTimeStep;
				}
				i2 = Hosts[i].infector;
				if (i2 >= 0)
				{
					Hosts[i2].listpos++;
					if (Hosts[i2].latent_time * P.ModelTimeStep <= P.SimulationDuration)
					{
						TimeSeries[(int)(Hosts[i2].latent_time * P.ModelTimeStep / P.OutputTimeStep)].Rtype[Hosts[i].infect_type % INFECT_TYPE_MASK]++;
						TimeSeries[(int)(Hosts[i2].latent_time * P.ModelTimeStep / P.OutputTimeStep)].Rage[HOST_AGE_GROUP(i)]++;
					}
				}
			}
			inftype[Hosts[i].infect_type % INFECT_TYPE_MASK]++;
			j--;
			if (j == 0)
			{
				if (l < 0) l = 0;
				inf_household[k][l]++;
				case_household[k][lc2]++; //now recording total symptomatic cases, rather than infections conditional on symptomatic index
				l = lc = lc2 = 0; t = 1e10;
			}
		}
	}
	for (b = 0; b < P.NumCells; b++)
		if ((Cells[b].S != Cells[b].n) || (Cells[b].R > 0))
			for (c = 0; c < Cells[b].n; c++)
			{
				i = Cells[b].members[c];
				if (Hosts[i].is_recovered() || Hosts[i].is_dead())
				{
					l = Hosts[i].infect_type / INFECT_TYPE_MASK;
					if ((l < MAX_GEN_REC) && (Hosts[i].listpos < MAX_SEC_REC)) indivR0[Hosts[i].listpos][l]++;
				}
			}
	/* 	if(!TimeSeries[P.NumOutputTimeSteps-1].extinct) */
	{
		for (i = 0; i < INFECT_TYPE_MASK; i++) inftype_av[i] += inftype[i];
		for (i = 0; i < MAX_COUNTRIES; i++)
		{
			infcountry_av[i] += infcountry[i];
			if (infcountry[i] > 0) infcountry_num[i]++;
		}
		for (i = 0; i < MAX_SEC_REC; i++)
			for (j = 0; j < MAX_GEN_REC; j++)
				indivR0_av[i][j] += indivR0[i][j];
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
			for (j = 0; j <= MAX_HOUSEHOLD_SIZE; j++)
			{
				inf_household_av[i][j] += inf_household[i][j];
				case_household_av[i][j] += case_household[i][j];
			}
	}
	k = (P.Interventions_StartDate_CalTime > 0) ? ((int)(P.Interventions_StartDate_CalTime - P.DateTriggerReached_SimTime)) : 0;
	for (n = 0; n < P.NumOutputTimeSteps; n++)
	{
		TimeSeries[n].t += k;
		s = 0;
		if (TimeSeries[n].Rdenom == 0) TimeSeries[n].Rdenom = 1e-10;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			TimeSeries[n].Rage[i] /= TimeSeries[n].Rdenom;
		for (i = 0; i < INFECT_TYPE_MASK; i++)
			s += (TimeSeries[n].Rtype[i] /= TimeSeries[n].Rdenom);
		TimeSeries[n].Rdenom = s;
	}
	nf = sizeof(Results) / sizeof(double);
	if (!P.DoAdUnits) nf -= MAX_ADUNITS; // TODO: This still processes most of the AdUnit arrays; just not the last one

	if (TimeSeries[P.NumOutputTimeSteps - 1].extinct)
	{
		TSMean = TSMeanE; TSVar = TSVarE; P.NRactE++;
	}
	else
	{
		TSMean = TSMeanNE; TSVar = TSVarNE; P.NRactNE++;
	}
	lc = -k;

	// This calculates sum and sum of squares of entire TimeSeries array
	for (n = 0; n < P.NumOutputTimeSteps; n++)
	{
		if ((n + lc >= 0) && (n + lc < P.NumOutputTimeSteps))
		{
			if (s < TimeSeries[n + lc].incC) { s = TimeSeries[n + lc].incC; t = P.OutputTimeStep * ((double)(_I64(n) + lc)); }
			res = (double*)&TimeSeries[n + lc];
			res_av = (double*)&TSMean[n];
			res_var = (double*)&TSVar[n];
			for (std::size_t i3 = ResultsDoubleOffsetStart /* skip over initial fields */; i3 < nf; i3++)
			{
				res_av[i3] += res[i3];
				res_var[i3] += res[i3] * res[i3];
			}
			if (P.DoAdUnits && P.OutputAdUnitAge)
				for (std::size_t age = 0; age < NUM_AGE_GROUPS; ++age)
					for (std::size_t adunit = 0; adunit < (size_t) P.NumAdunits; ++adunit)
					{
						TSMean[n].prevInf_age_adunit[age][adunit] = TimeSeries[n + lc].prevInf_age_adunit[age][adunit];
						TSMean[n].cumInf_age_adunit [age][adunit] = TimeSeries[n + lc].cumInf_age_adunit [age][adunit];
						TSMean[n].incInf_age_adunit [age][adunit] = TimeSeries[n + lc].incInf_age_adunit [age][adunit];
					}

			if (TSMean[n].cumTmax < TimeSeries[n + lc].cumT) TSMean[n].cumTmax = TimeSeries[n + lc].cumT;
			if (TSMean[n].cumVmax < TimeSeries[n + lc].cumV) TSMean[n].cumVmax = TimeSeries[n + lc].cumV;
		}
		TSMean[n].t += ((double) n )* P.OutputTimeStep;
	}
	PeakHeightSum += s;
	PeakHeightSS += s * s;
	PeakTimeSum += t;
	PeakTimeSS += t * t;
}


void CalcOriginDestMatrix_adunit()
{
	/** function: CalcOriginDestMatrix_adunit()
	 *
	 * purpose: to output the origin destination matrix between admin units
	 *
	 * parameters: none
	 *
	 * returns: none
	 *
	 * author: ggilani, date: 28/01/15
	 */

#pragma omp parallel for schedule(static) default(none) \
		shared(P, Cells, CellLookup, Mcells, StateT)
	for (int tn = 0; tn < P.NumThreads; tn++)
	{
		for (int i = tn; i < P.NumPopulatedCells; i += P.NumThreads)
		{
			//reset pop density matrix to zero
			double pop_dens_from[MAX_ADUNITS] = {};

			//find index of cell from which flow travels
			ptrdiff_t cl_from = CellLookup[i] - Cells;
			ptrdiff_t cl_from_mcl = (cl_from / P.nch) * P.NMCL * P.total_microcells_high_ + (cl_from % P.nch) * P.NMCL;

			//loop over microcells in these cells to find populations in each admin unit and so flows
			for (int k = 0; k < P.NMCL; k++)
			{
				for (int l = 0; l < P.NMCL; l++)
				{
					//get index of microcell
					ptrdiff_t mcl_from = cl_from_mcl + l + _I64(k) * P.total_microcells_high_;
					if (Mcells[mcl_from].n > 0)
					{
						//get proportion of each population of cell that exists in each admin unit
						pop_dens_from[Mcells[mcl_from].adunit] += (((double)Mcells[mcl_from].n) / ((double)Cells[cl_from].n));
					}
				}
			}

			for (int j = i; j < P.NumPopulatedCells; j++)
			{
				//reset pop density matrix to zero
				double pop_dens_to[MAX_ADUNITS] = {};

				//find index of cell which flow travels to
				ptrdiff_t cl_to = CellLookup[j] - Cells;
				ptrdiff_t cl_to_mcl = (cl_to / P.nch) * P.NMCL * P.total_microcells_high_ + (cl_to % P.nch) * P.NMCL;
				//calculate distance and kernel between the cells
				//total_flow=Cells[cl_from].max_trans[j]*Cells[cl_from].n*Cells[cl_to].n;
				double total_flow;
				if (j == 0)
				{
					total_flow = (double)Cells[cl_from].cum_trans[j] * Cells[cl_from].n;
				}
				else
				{
					total_flow = ((double)Cells[cl_from].cum_trans[j] - Cells[cl_from].cum_trans[j - 1]) * Cells[cl_from].n;
				}

				//loop over microcells within destination cell
				for (int m = 0; m < P.NMCL; m++)
				{
					for (int p = 0; p < P.NMCL; p++)
					{
						//get index of microcell
						ptrdiff_t mcl_to = cl_to_mcl + p + _I64(m) * P.total_microcells_high_;
						if (Mcells[mcl_to].n > 0)
						{
							//get proportion of each population of cell that exists in each admin unit
							pop_dens_to[Mcells[mcl_to].adunit] += (((double)Mcells[mcl_to].n) / ((double)Cells[cl_to].n));
						}
					}
				}

				for (int m = 0; m < P.NumAdunits; m++)
				{
					for (int p = 0; p < P.NumAdunits; p++)
					{
						if (m != p)
						{
							double flow = total_flow * pop_dens_from[m] * pop_dens_to[p]; //updated to remove reference to cross-border flows: ggilani 26/03/20
							StateT[tn].origin_dest[m][p] += flow;
							StateT[tn].origin_dest[p][m] += flow;
						}
					}
				}
			}
		}
	}

	//Sum up flow between adunits across threads
	for (int i = 0; i < P.NumAdunits; i++)
	{
		for (int j = 0; j < P.NumAdunits; j++)
		{
			for (int k = 0; k < P.NumThreads; k++)
			{
				AdUnits[i].origin_dest[j] += StateT[k].origin_dest[i][j];
			}
		}
	}
}

void GetInverseCdf(ParamMap fallback, ParamMap params, const char* icdf_name, InverseCdf* inverseCdf, Param* P, double start_value)
{
	Params::get_double_vec(fallback, params, icdf_name, inverseCdf->get_values(), CDF_RES + 1, 0, CDF_RES + 1, P);
  if (!Params::param_found(fallback, params, icdf_name)) {
		inverseCdf->set_neg_log(start_value);
	}
	inverseCdf->assign_exponent();
}
