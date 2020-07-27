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
#include "BinIO.h"
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

int GetInputParameter(FILE*, FILE*, const char*, const char*, void*, int, int, int);
int GetInputParameter2(FILE*, FILE*, const char*, const char*, void*, int, int, int);
int GetInputParameter2all(FILE*, FILE*, FILE*, const char*, const char*, void*, int, int, int);
int GetInputParameter3(FILE*, const char*, const char*, void*, int, int, int);

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

int ns, DoInitUpdateProbs, InterruptRun = 0;
int PlaceDistDistrib[NUM_PLACE_TYPES][MAX_DIST], PlaceSizeDistrib[NUM_PLACE_TYPES][MAX_PLACE_SIZE];

/* int NumPC,NumPCD; */
const int MAXINTFILE = 10;
std::vector<std::string> InterventionFiles;

// default start value for icdf double arrays (was hardcoded as 100)
const double ICDF_START = 100.0;

void GetInverseCdf(FILE* param_file_dat, FILE* preparam_file_dat, const char* icdf_name, InverseCdf* inverseCdf, double start_value = ICDF_START);

int main(int argc, char* argv[])
{
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
		if (sep == std::string::npos) {
			ERR_CRITICAL("Expected argument value to be in the format '<D>,<S>' where <D> is the "
							"timestep interval when a snapshot should be saved and <S> is the "
							"filename in which to save the snapshot");
		}
		parse_double(input.substr(0, sep), P.SnapshotSaveTime);
		parse_read_file(input.substr(sep + 1), snapshot_save_file);
	};

	int cl = clock();

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

	for (int i = 0; i <= 99; i++) {
		std::string param = "CLP" + std::to_string(i);
		std::string description = "Overwrites #" + std::to_string(i) + " wildcard in parameter file";
		args.add_double_option(param, P.clP[i], description);
		if (i < 10) {
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
	fprintf(stderr, "Bitmap Format = *.%s\n", P.BitmapFormat == BitmapFormats::PNG ? "png" : "bmp");
	fprintf(stderr, "sizeof(int)=%i sizeof(long)=%i sizeof(float)=%i sizeof(double)=%i sizeof(unsigned short int)=%i sizeof(int *)=%i\n", (int)sizeof(int), (int)sizeof(long), (int)sizeof(float), (int)sizeof(double), (int)sizeof(unsigned short int), (int)sizeof(int*));

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
	//// **** SET UP OMP / THREADS
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

#ifdef _OPENMP
	P.NumThreads = omp_get_max_threads();
	if ((P.MaxNumThreads > 0) && (P.MaxNumThreads < P.NumThreads)) P.NumThreads = P.MaxNumThreads;
	if (P.NumThreads > MAX_NUM_THREADS)
	{
		fprintf(stderr, "Assigned number of threads (%d) > MAX_NUM_THREADS (%d)\n", P.NumThreads, MAX_NUM_THREADS);
		P.NumThreads = MAX_NUM_THREADS;
	}
	fprintf(stderr, "Using %d threads\n", P.NumThreads);
	omp_set_num_threads(P.NumThreads);

#pragma omp parallel default(shared)
	{
		fprintf(stderr, "Thread %i initialised\n", omp_get_thread_num());
	}
	/* fprintf(stderr,"int=%i\tfloat=%i\tdouble=%i\tint *=%i\n",(int) sizeof(int),(int) sizeof(float),(int) sizeof(double),(int) sizeof(int *));	*/
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

	fprintf(stderr, "Model setup in %lf seconds\n", ((double)(clock() - cl)) / CLOCKS_PER_SEC);



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
			for (int i = 0; (i < P.NumRealisations) && (P.NRactNE < P.NumNonExtinctRealisations); i++)
			{
				if (P.NumRealisations > 1)
				{
					output_file = output_file_base + "." + std::to_string(i);
					fprintf(stderr, "Realisation %i of %i  (time=%lf nr_ne=%i)\n", i + 1, P.NumRealisations, ((double)(clock() - cl)) / CLOCKS_PER_SEC, P.NRactNE);
				}
				///// Set and save seeds
				if (((i == 0) && (P.FitIter == 1)) || (P.ResetSeeds && P.KeepSameSeeds))
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
					if ((P.ModelCalibIteration%14 == 0) && (ModelCalibLoop < 4))
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
					InitModel(i);
					if (!snapshot_load_file.empty()) LoadSnapshot(snapshot_load_file);
					ContCalib = RunModel(i, snapshot_save_file, snapshot_load_file, output_file_base);
				}
				while (ContCalib);
				if (!data_file.empty()) CalcLikelihood(i, data_file, output_file_base);
				if (P.OutputNonSummaryResults)
				{
					if (((!TimeSeries[P.NumSamples - 1].extinct) || (!P.OutputOnlyNonExtinct)) && (P.OutputEveryRealisation))
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
			//sprintf(OutFile, "%s.avE", OutFileBase);
			//SaveSummaryResults();

			Bitmap_Finalise();

			fprintf(stderr, "Extinction in %i out of %i runs\n", P.NRactE, P.NRactNE + P.NRactE);
			fprintf(stderr, "Model ran in %lf seconds\n", ((double)(clock() - cl)) / CLOCKS_PER_SEC);
			fprintf(stderr, "Model finished\n");
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

void parse_intervention_file_option(std::string const& input) {
	std::string output;
	parse_read_file(input, output);
	InterventionFiles.emplace_back(output);
}

void ReadParams(std::string const& ParamFile, std::string const& PreParamFile, std::string const& AdUnitFile)
{
	FILE* ParamFile_dat, * PreParamFile_dat, * AdminFile_dat;
	double s, t, AgeSuscScale;
	int i, j, j1, j2, k, f, nc, na;
	char CountryNameBuf[128 * MAX_COUNTRIES], AdunitListNamesBuf[128 * MAX_ADUNITS];
	char* CountryNames[MAX_COUNTRIES], * AdunitListNames[MAX_ADUNITS];

	AgeSuscScale = 1.0;
	if (!(ParamFile_dat = fopen(ParamFile.c_str(), "rb"))) ERR_CRITICAL("Unable to open parameter file\n");
	PreParamFile_dat = fopen(PreParamFile.c_str(), "rb");
	if (!(AdminFile_dat = fopen(AdUnitFile.c_str(), "rb"))) AdminFile_dat = ParamFile_dat;

	if (P.FitIter == 0)
	{
		for (i = 0; i < MAX_COUNTRIES; i++) {
			CountryNames[i] = CountryNameBuf + 128 * i;
			CountryNames[i][0] = 0;
		}
		for (i = 0; i < MAX_ADUNITS; i++) {
			AdunitListNames[i] = AdunitListNamesBuf + 128 * i;
			AdunitListNames[i][0] = 0;
		}
		for (i = 0; i < 100; i++) P.clP_copies[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, AdminFile_dat, "Longitude cut line", "%lf", (void*)&(P.LongitudeCutLine), 1, 1, 0)) {
			P.LongitudeCutLine = -360.0;
		}
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Update timestep", "%lf", (void*)&(P.TimeStep), 1, 1, 0);
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Sampling timestep", "%lf", (void*)&(P.SampleStep), 1, 1, 0);
		if (P.TimeStep > P.SampleStep) ERR_CRITICAL("Update step must be smaller than sampling step\n");
		t = ceil(P.SampleStep / P.TimeStep - 1e-6);
		P.UpdatesPerSample = (int)t;
		P.TimeStep = P.SampleStep / t;
		P.TimeStepsPerDay = ceil(1.0 / P.TimeStep - 1e-6);
		fprintf(stderr, "Update step = %lf\nSampling step = %lf\nUpdates per sample=%i\nTimeStepsPerDay=%lf\n", P.TimeStep, P.SampleStep, P.UpdatesPerSample, P.TimeStepsPerDay);
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Sampling time", "%lf", (void*)&(P.SampleTime), 1, 1, 0);
		P.NumSamples = 1 + (int)ceil(P.SampleTime / P.SampleStep);
		GetInputParameter(PreParamFile_dat, AdminFile_dat, "Population size", "%i", (void*)&(P.PopSize), 1, 1, 0);
		if (P.NumRealisations == 0)
		{
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Number of realisations", "%i", (void*)&(P.NumRealisations), 1, 1, 0);
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of non-extinct realisations", "%i", (void*)&(P.NumNonExtinctRealisations), 1, 1, 0)) P.NumNonExtinctRealisations = P.NumRealisations;
		}
		else
		{
			P.NumNonExtinctRealisations = P.NumRealisations;
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of cases defining small outbreak", "%i", (void*) & (P.SmallEpidemicCases), 1, 1, 0)) P.SmallEpidemicCases = -1;

		P.NC = -1;
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Number of micro-cells per spatial cell width", "%i", (void*)&(P.NMCL), 1, 1, 0);
		//added parameter to reset seeds after every run
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Reset seeds for every run", "%i", (void*)&(P.ResetSeeds), 1, 1, 0)) P.ResetSeeds = 0;
		if (P.ResetSeeds)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Keep same seeds for every run", "%i", (void*)&(P.KeepSameSeeds), 1, 1, 0)) P.KeepSameSeeds = 0; //added this to control which seeds are used: ggilani 27/11/19
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Reset seeds after intervention", "%i", (void*)&(P.ResetSeedsPostIntervention), 1, 1, 0)) P.ResetSeedsPostIntervention = 0;
		if (P.ResetSeedsPostIntervention)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Time to reset seeds after intervention", "%i", (void*)&(P.TimeToResetSeeds), 1, 1, 0)) P.TimeToResetSeeds = 1000000;
		}
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Include households", "%i", (void*)&(P.DoHouseholds), 1, 1, 0)) P.DoHouseholds = 1;

		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel resolution", "%i", (void*)&P.KernelLookup.size_, 1, 1, 0)) P.KernelLookup.size_ = 4000000;
		if (P.KernelLookup.size_ < 2000000)
		{
			ERR_CRITICAL_FMT("[Kernel resolution] needs to be at least 2000000 - not %d", P.KernelLookup.size_);
		}
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel higher resolution factor", "%i", (void*)&P.KernelLookup.expansion_factor_, 1, 1, 0)) P.KernelLookup.expansion_factor_ = P.KernelLookup.size_ / 1600;
		if (P.KernelLookup.expansion_factor_ < 1 || P.KernelLookup.expansion_factor_ >= P.KernelLookup.size_)
		{
			ERR_CRITICAL_FMT("[Kernel higher resolution factor] needs to be in range [1, P.NKR = %d) - not %d", P.KernelLookup.size_, P.KernelLookup.expansion_factor_);
		}
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputAge", "%i", (void*)&(P.OutputAge), 1, 1, 0)) P.OutputAge = 1;				//// ON  by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputSeverity", "%i", (void*)&(P.OutputSeverity), 1, 1, 0)) P.OutputSeverity = 1;	//// ON  by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputSeverityAdminUnit", "%i", (void*)&(P.OutputSeverityAdminUnit), 1, 1, 0)) P.OutputSeverityAdminUnit = 1;	//// ON  by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputSeverityAge", "%i", (void*)&(P.OutputSeverityAge), 1, 1, 0)) P.OutputSeverityAge = 1;		//// ON  by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputAdUnitAge", "%i", (void*)&(P.OutputAdUnitAge), 1, 1, 0)) P.OutputAdUnitAge = 0;			//// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputR0", "%i", (void*)&(P.OutputR0), 1, 1, 0)) P.OutputR0 = 0;				    //// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputControls", "%i", (void*)&(P.OutputControls), 1, 1, 0)) P.OutputControls = 0;		    //// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputCountry", "%i", (void*)&(P.OutputCountry), 1, 1, 0)) P.OutputCountry = 0;		    //// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputAdUnitVar", "%i", (void*)&(P.OutputAdUnitVar), 1, 1, 0)) P.OutputAdUnitVar = 0;		    //// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputHousehold", "%i", (void*)&(P.OutputHousehold), 1, 1, 0)) P.OutputHousehold = 0;		    //// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputInfType", "%i", (void*)&(P.OutputInfType), 1, 1, 0)) P.OutputInfType = 0;		    //// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputNonSeverity", "%i", (void*)&(P.OutputNonSeverity), 1, 1, 0)) P.OutputNonSeverity = 0;		//// OFF by default.
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "OutputNonSummaryResults", "%i", (void*)&(P.OutputNonSummaryResults), 1, 1, 0)) P.OutputNonSummaryResults = 0;	//// OFF by default.

	if (P.DoHouseholds)
	{
		GetInputParameter(PreParamFile_dat, AdminFile_dat, "Household size distribution", "%lf", (void*)P.HouseholdSizeDistrib[0], MAX_HOUSEHOLD_SIZE, 1, 0);
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Household attack rate", "%lf", (void*)&(P.HouseholdTrans), 1, 1, 0);
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Household transmission denominator power", "%lf", (void*)&(P.HouseholdTransPow), 1, 1, 0);
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Correct age distribution after household allocation to exactly match specified demography", "%i", (void*)&(P.DoCorrectAgeDist), 1, 1, 0)) P.DoCorrectAgeDist = 0;
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
		    P.HouseholdDenomLookup[i] = 1 / pow(((double)(i + 1)), P.HouseholdTransPow);
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Include administrative units within countries", "%i", (void*)&(P.DoAdUnits), 1, 1, 0)) P.DoAdUnits = 1;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Divisor for countries", "%i", (void*)&(P.CountryDivisor), 1, 1, 0)) P.CountryDivisor = 1;
		if (P.DoAdUnits)
		{
			char** AdunitNames, * AdunitNamesBuf;
			if (!(AdunitNames = (char**)malloc(3 * ADUNIT_LOOKUP_SIZE * sizeof(char*)))) ERR_CRITICAL("Unable to allocate temp storage\n");
			if (!(AdunitNamesBuf = (char*)malloc(3 * ADUNIT_LOOKUP_SIZE * 360 * sizeof(char)))) ERR_CRITICAL("Unable to allocate temp storage\n");

			for (i = 0; i < ADUNIT_LOOKUP_SIZE; i++)
			{
				P.AdunitLevel1Lookup[i] = -1;
				AdunitNames[3 * i] = AdunitNamesBuf + 3 * i * 360;
				AdunitNames[3 * i + 1] = AdunitNamesBuf + 3 * i * 360 + 60;
				AdunitNames[3 * i + 2] = AdunitNamesBuf + 3 * i * 360 + 160;
			}
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Divisor for level 1 administrative units", "%i", (void*)&(P.AdunitLevel1Divisor), 1, 1, 0)) P.AdunitLevel1Divisor = 1;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Mask for level 1 administrative units", "%i", (void*)&(P.AdunitLevel1Mask), 1, 1, 0)) P.AdunitLevel1Mask = 1000000000;
			na = (GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Codes and country/province names for admin units", "%s", (void*)AdunitNames, 3 * ADUNIT_LOOKUP_SIZE, 1, 0)) / 3;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Number of countries to include", "%i", (void*)&nc, 1, 1, 0)) nc = 0;
			if ((na > 0) && (nc > 0))
			{
				P.DoAdunitBoundaries = (nc > 0);
				nc = abs(nc);
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "List of names of countries to include", "%s", (nc > 1) ? ((void*)CountryNames) : ((void*)CountryNames[0]), nc, 1, 0);
				P.NumAdunits = 0;
				for (i = 0; i < na; i++)
					for (j = 0; j < nc; j++)
						if ((AdunitNames[3 * i + 1][0]) && (!strcmp(AdunitNames[3 * i + 1], CountryNames[j])) && (atoi(AdunitNames[3 * i]) != 0))
						{
							AdUnits[P.NumAdunits].id = atoi(AdunitNames[3 * i]);
							P.AdunitLevel1Lookup[(AdUnits[P.NumAdunits].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor] = P.NumAdunits;
							if (strlen(AdunitNames[3 * i + 1]) < 100) strcpy(AdUnits[P.NumAdunits].cnt_name, AdunitNames[3 * i + 1]);
							if (strlen(AdunitNames[3 * i + 2]) < 200) strcpy(AdUnits[P.NumAdunits].ad_name, AdunitNames[3 * i + 2]);
							//						fprintf(stderr,"%i %s %s ## ",AdUnits[P.NumAdunits].id,AdUnits[P.NumAdunits].cnt_name,AdUnits[P.NumAdunits].ad_name);
							P.NumAdunits++;
						}
			}
			else
			{
				if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Number of level 1 administrative units to include", "%i", (void*)&(P.NumAdunits), 1, 1, 0)) P.NumAdunits = 0;
				if (P.NumAdunits > 0)
				{
					P.DoAdunitBoundaries = 1;
					if (P.NumAdunits > MAX_ADUNITS) ERR_CRITICAL("MAX_ADUNITS too small.\n");
					GetInputParameter(PreParamFile_dat, AdminFile_dat, "List of level 1 administrative units to include", "%s", (P.NumAdunits > 1) ? ((void*)AdunitListNames) : ((void*)AdunitListNames[0]), P.NumAdunits, 1, 0);
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

			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output incidence by administrative unit", "%i", (void*)&(P.DoAdunitOutput), 1, 1, 0)) P.DoAdunitOutput = 0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Draw administrative unit boundaries on maps", "%i", (void*)&(P.DoAdunitBoundaryOutput), 1, 1, 0)) P.DoAdunitBoundaryOutput = 0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Correct administrative unit populations", "%i", (void*)&(P.DoCorrectAdunitPop), 1, 1, 0)) P.DoCorrectAdunitPop = 0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Fix population size at specified value", "%i", (void*)&(P.DoSpecifyPop), 1, 1, 0)) P.DoSpecifyPop = 0;
			fprintf(stderr, "Using %i administrative units\n", P.NumAdunits);
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Divisor for administrative unit codes for boundary plotting on bitmaps", "%i", (void*)&(P.AdunitBitmapDivisor), 1, 1, 0)) P.AdunitBitmapDivisor = 1;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Only output household to place distance distribution for one administrative unit", "%i", (void*)&(P.DoOutputPlaceDistForOneAdunit), 1, 1, 0)) P.DoOutputPlaceDistForOneAdunit = 0;
			if (P.DoOutputPlaceDistForOneAdunit)
			{
				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Administrative unit for which household to place distance distribution to be output", "%i", (void*)&(P.OutputPlaceDistAdunit), 1, 1, 0)) P.DoOutputPlaceDistForOneAdunit = 0;
			}
		}
		else

		{
			P.DoAdunitBoundaries = P.DoAdunitBoundaryOutput = P.DoAdunitOutput = P.DoCorrectAdunitPop = P.DoSpecifyPop = 0;
			P.AdunitLevel1Divisor = 1; P.AdunitLevel1Mask = 1000000000;
			P.AdunitBitmapDivisor = P.AdunitLevel1Divisor;
		}
	}

	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Include age", "%i", (void*)&(P.DoAge), 1, 1, 0)) P.DoAge = 1;
	if (!P.DoAge)
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

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Initial immunity acts as partial immunity", "%i", (void*)&(P.DoPartialImmunity), 1, 1, 0)) P.DoPartialImmunity = 1;
		if ((P.DoHouseholds) && (!P.DoPartialImmunity))
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Initial immunity applied to all household members", "%i", (void*)&(P.DoWholeHouseholdImmunity), 1, 1, 0)) P.DoWholeHouseholdImmunity = 0;
		}
		else
			P.DoWholeHouseholdImmunity = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Initial immunity profile by age", "%lf", (void*)P.InitialImmunity, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.InitialImmunity[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rates by age", "%lf", (void*)P.RelativeSpatialContact, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.RelativeSpatialContact[i] = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Age-dependent infectiousness", "%lf", (void*)P.AgeInfectiousness, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.AgeInfectiousness[i] = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Age-dependent susceptibility", "%lf", (void*)P.AgeSusceptibility, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.AgeSusceptibility[i] = 1.0;
		GetInputParameter(PreParamFile_dat, AdminFile_dat, "Age distribution of population", "%lf", (void*)P.PropAgeGroup[0], NUM_AGE_GROUPS, 1, 0);
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
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Relative travel rates by age", "%lf", (void*)P.RelativeTravelRate, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.RelativeTravelRate[i] = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "WAIFW matrix", "%lf", (void*)P.WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					P.WAIFW_Matrix[i][j] = 1.0;
		}
		else
		{
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
		P.DoDeath = 0;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)	t += P.AgeInfectiousness[i] * P.PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)	P.AgeInfectiousness[i] /= t;
	}

	if (P.FitIter == 0)
	{
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Include spatial transmission", "%i", (void*)&(P.DoSpatial), 1, 1, 0)) P.DoSpatial = 1;
		GetInputParameter(PreParamFile_dat, AdminFile_dat, "Kernel type", "%i", (void*)&(P.MoveKernel.type_), 1, 1, 0);
		GetInputParameter(PreParamFile_dat, AdminFile_dat, "Kernel scale", "%lf", (void*)&(P.MoveKernel.scale_), 1, 1, 0);
		if (P.KernelOffsetScale != 1)
		{
			P.MoveKernel.scale_ *= P.KernelOffsetScale;
		}
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel 3rd param", "%lf", (void*)&(P.MoveKernel.p3_), 1, 1, 0)) P.MoveKernel.p3_ = 0;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel 4th param", "%lf", (void*)&(P.MoveKernel.p4_), 1, 1, 0)) P.MoveKernel.p4_ = 0;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel Shape", "%lf", (void*)&(P.MoveKernel.shape_), 1, 1, 0)) P.MoveKernel.shape_ = 1.0;
		if (P.KernelPowerScale != 1)
		{
			P.MoveKernel.shape_ *= P.KernelPowerScale;
		}
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Airport Kernel Type", "%i", (void*)&(P.AirportKernel.type_), 1, 1, 0)) P.AirportKernel.type_ = P.MoveKernel.type_;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Airport Kernel Scale", "%lf", (void*)&(P.AirportKernel.scale_), 1, 1, 0)) P.AirportKernel.scale_ = P.MoveKernel.scale_;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Airport Kernel Shape", "%lf", (void*)&(P.AirportKernel.shape_), 1, 1, 0)) P.AirportKernel.shape_ = P.MoveKernel.shape_;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Airport Kernel 3rd param", "%lf", (void*)&(P.AirportKernel.p3_), 1, 1, 0)) P.AirportKernel.p3_ = P.MoveKernel.p3_;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Airport Kernel 4th param", "%lf", (void*)&(P.AirportKernel.p4_), 1, 1, 0)) P.AirportKernel.p4_ = P.MoveKernel.p4_;

		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Include places", "%i", (void*)&(P.DoPlaces), 1, 1, 0)) P.DoPlaces = 1;
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Number of types of places", "%i", (void*)&(P.PlaceTypeNum), 1, 1, 0)) P.PlaceTypeNum = 0;
			if (P.PlaceTypeNum == 0) P.DoPlaces = P.DoAirports = 0;
		}
		else
			P.PlaceTypeNum = P.DoAirports = 0;
	}
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Scaling of household contacts for care home residents", "%lf", (void*)&(P.CareHomeResidentHouseholdScaling), 1, 1, 0)) P.CareHomeResidentHouseholdScaling = 1.0;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Scaling of spatial contacts for care home residents", "%lf", (void*)&(P.CareHomeResidentSpatialScaling), 1, 1, 0)) P.CareHomeResidentSpatialScaling = 1.0;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Scaling of between group (home) contacts for care home residents", "%lf", (void*)&(P.CareHomeResidentPlaceScaling), 1, 1, 0)) P.CareHomeResidentPlaceScaling = 1.0;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Scaling of within group (home) contacts for care home workers", "%lf", (void*)&(P.CareHomeWorkerGroupScaling), 1, 1, 0)) P.CareHomeWorkerGroupScaling = 1.0;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Relative probability that care home residents are hospitalised", "%lf", (void*)&(P.CareHomeRelProbHosp), 1, 1, 0)) P.CareHomeRelProbHosp = 1.0;

		if (P.FitIter == 0)
		{
			if (P.PlaceTypeNum > NUM_PLACE_TYPES) ERR_CRITICAL("Too many place types\n");
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Place type number for care homes", "%i", (void*)&(P.CareHomePlaceType), 1, 1, 0)) P.CareHomePlaceType = -1;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Allow initial infections to be in care homes", "%i", (void*)&(P.CareHomeAllowInitialInfections), 1, 1, 0)) P.CareHomeAllowInitialInfections = 0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Minimum age of care home residents", "%i", (void*)&(P.CareHomeResidentMinimumAge), 1, 1, 0)) P.CareHomeResidentMinimumAge = 1000;

			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Minimum age for age group 1 in place types", "%i", (void*)P.PlaceTypeAgeMin, P.PlaceTypeNum, 1, 0);
			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Maximum age for age group 1 in place types", "%i", (void*)P.PlaceTypeAgeMax, P.PlaceTypeNum, 1, 0);
			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Proportion of age group 1 in place types", "%lf", (void*)&(P.PlaceTypePropAgeGroup), P.PlaceTypeNum, 1, 0);
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Proportion of age group 2 in place types", "%lf", (void*)&(P.PlaceTypePropAgeGroup2), P.PlaceTypeNum, 1, 0))
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
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Minimum age for age group 2 in place types", "%i", (void*)P.PlaceTypeAgeMin2, P.PlaceTypeNum, 1, 0);
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Maximum age for age group 2 in place types", "%i", (void*)P.PlaceTypeAgeMax2, P.PlaceTypeNum, 1, 0);
			}
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Proportion of age group 3 in place types", "%lf", (void*)&(P.PlaceTypePropAgeGroup3), P.PlaceTypeNum, 1, 0))
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
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Minimum age for age group 3 in place types", "%i", (void*)P.PlaceTypeAgeMin3, P.PlaceTypeNum, 1, 0);
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Maximum age for age group 3 in place types", "%i", (void*)P.PlaceTypeAgeMax3, P.PlaceTypeNum, 1, 0);
			}
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel shape params for place types", "%lf", (void*)&(P.PlaceTypeKernelShape), P.PlaceTypeNum, 1, 0))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P.PlaceTypeKernelShape[i] = P.MoveKernel.shape_;
					P.PlaceTypeKernelScale[i] = P.MoveKernel.scale_;
				}
			}
			else
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Kernel scale params for place types", "%lf", (void*)&(P.PlaceTypeKernelScale), P.PlaceTypeNum, 1, 0);
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel 3rd param for place types", "%lf", (void*)&(P.PlaceTypeKernelP3), P.PlaceTypeNum, 1, 0))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P.PlaceTypeKernelP3[i] = P.MoveKernel.p3_;
					P.PlaceTypeKernelP4[i] = P.MoveKernel.p4_;
				}
			}
			else
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Kernel 4th param for place types", "%lf", (void*)&(P.PlaceTypeKernelP4), P.PlaceTypeNum, 1, 0);
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Number of closest places people pick from (0=all) for place types", "%i", (void*)&(P.PlaceTypeNearestNeighb), P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeNearestNeighb[i] = 0;
			if (P.DoAdUnits)
			{
				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Degree to which crossing administrative unit boundaries to go to places is inhibited", "%lf", (void*)&(P.InhibitInterAdunitPlaceAssignment), P.PlaceTypeNum, 1, 0))
					for (i = 0; i < NUM_PLACE_TYPES; i++)
						P.InhibitInterAdunitPlaceAssignment[i] = 0;
			}

			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Include air travel", "%i", (void*)&(P.DoAirports), 1, 1, 0)) P.DoAirports = 0;
			if (!P.DoAirports)
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
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Number of non-airport places", "%i", (void*)&(P.PlaceTypeNoAirNum), 1, 1, 0);
				GetInputParameter(PreParamFile_dat, AdminFile_dat, "Hotel place type", "%i", (void*)&(P.HotelPlaceType), 1, 1, 0);
				if (P.PlaceTypeNoAirNum >= P.PlaceTypeNum) {
					ERR_CRITICAL_FMT("[Number of non-airport places] parameter (%d) is greater than number of places (%d).\n", P.PlaceTypeNoAirNum, P.PlaceTypeNum);
				}
				if (P.HotelPlaceType < P.PlaceTypeNoAirNum || P.HotelPlaceType >= P.PlaceTypeNum) {
					ERR_CRITICAL_FMT("[Hotel place type] parameter (%d) not in the range [%d, %d)\n", P.HotelPlaceType, P.PlaceTypeNoAirNum, P.PlaceTypeNum);
				}

				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Scaling factor for input file to convert to daily traffic", "%lf", (void*)&(P.AirportTrafficScale), 1, 1, 0)) P.AirportTrafficScale = 1.0;
				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of hotel attendees who are local", "%lf", (void*)&(P.HotelPropLocal), 1, 1, 0)) P.HotelPropLocal = 0;
				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Distribution of duration of air journeys", "%lf", (void*)&(P.JourneyDurationDistrib), MAX_TRAVEL_TIME, 1, 0))
				{
					P.JourneyDurationDistrib[0] = 1;
					for (i = 0; i < MAX_TRAVEL_TIME; i++)
						P.JourneyDurationDistrib[i] = 0;
				}
				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Distribution of duration of local journeys", "%lf", (void*)&(P.LocalJourneyDurationDistrib), MAX_TRAVEL_TIME, 1, 0))
				{
					P.LocalJourneyDurationDistrib[0] = 1;
					for (i = 0; i < MAX_TRAVEL_TIME; i++)
						P.LocalJourneyDurationDistrib[i] = 0;
				}
				P.MeanJourneyTime = P.MeanLocalJourneyTime = 0;
				for (i = 0; i < MAX_TRAVEL_TIME; i++)
				{
					P.MeanJourneyTime += ((double)(i)) * P.JourneyDurationDistrib[i];
					P.MeanLocalJourneyTime += ((double)(i)) * P.LocalJourneyDurationDistrib[i];
				}
				fprintf(stderr, "Mean duration of local journeys = %lf days\n", P.MeanLocalJourneyTime);
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
			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Mean size of place types", "%lf", (void*)P.PlaceTypeMeanSize, P.PlaceTypeNum, 1, 0);
			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Param 1 of place group size distribution", "%lf", (void*)P.PlaceTypeGroupSizeParam1, P.PlaceTypeNum, 1, 0);
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Power of place size distribution", "%lf", (void*)P.PlaceTypeSizePower, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeSizePower[i] = 0;
			//added to enable lognormal distribution - ggilani 09/02/17
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Standard deviation of place size distribution", "%lf", (void*)P.PlaceTypeSizeSD, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeSizeSD[i] = 0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Offset of place size distribution", "%lf", (void*)P.PlaceTypeSizeOffset, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeSizeOffset[i] = 0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Maximum of place size distribution", "%lf", (void*)P.PlaceTypeSizeMax, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeSizeMax[i] = 1e20;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Minimum of place size distribution", "%lf", (void*)P.PlaceTypeSizeMin, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeSizeMin[i] = 1.0;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Kernel type for place types", "%i", (void*)P.PlaceTypeKernelType, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.PlaceTypeKernelType[i] = P.MoveKernel.type_;
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Place overlap matrix", "%lf", (void*)P.PlaceExclusivityMatrix, P.PlaceTypeNum * P.PlaceTypeNum, 1, 0))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					for (j = 0; j < NUM_PLACE_TYPES; j++)
						P.PlaceExclusivityMatrix[i * NUM_PLACE_TYPES + j] = (i == j) ? 1 : 0;
			}
		}
		/* Note P.PlaceExclusivityMatrix not used at present - places assumed exclusive (each person belongs to 0 or 1 place) */

		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Proportion of between group place links", "%lf", (void*)P.PlaceTypePropBetweenGroupLinks, P.PlaceTypeNum, 1, 0);
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Relative transmission rates for place types", "%lf", (void*)P.PlaceTypeTrans, P.PlaceTypeNum, 1, 0);
		for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceTypeTrans[i] *= AgeSuscScale;
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Daily seasonality coefficients", "%lf", (void*)P.Seasonality, DAYS_PER_YEAR, 1, 0))
	{
		P.DoSeasonality = 0;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			P.Seasonality[i] = 1;
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
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Number of seed locations", "%i", (void*)&(P.NumSeedLocations), 1, 1, 0)) P.NumSeedLocations = 1;
	if (P.NumSeedLocations > MAX_NUM_SEED_LOCATIONS)
	{
		fprintf(stderr, "Too many seed locations\n");
		P.NumSeedLocations = MAX_NUM_SEED_LOCATIONS;
	}
	GetInputParameter(PreParamFile_dat, AdminFile_dat, "Initial number of infecteds", "%i", (void*)P.NumInitialInfections, P.NumSeedLocations, 1, 0);
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Location of initial infecteds", "%lf", (void*)&(P.LocationInitialInfection[0][0]), P.NumSeedLocations * 2, 1, 0)) P.LocationInitialInfection[0][0] = P.LocationInitialInfection[0][1] = 0.0;
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Minimum population in microcell of initial infection", "%i", (void*)&(P.MinPopDensForInitialInfection), 1, 1, 0)) P.MinPopDensForInitialInfection = 0;
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Maximum population in microcell of initial infection", "%i", (void*)&(P.MaxPopDensForInitialInfection), 1, 1, 0)) P.MaxPopDensForInitialInfection = 10000000;
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Maximum age of initial infections", "%i", (void*)&(P.MaxAgeForInitialInfection), 1, 1, 0)) P.MaxAgeForInitialInfection = 1000;
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Randomise initial infection location", "%i", (void*)&(P.DoRandomInitialInfectionLoc), 1, 1, 0)) P.DoRandomInitialInfectionLoc = 1;
	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "All initial infections located in same microcell", "%i", (void*)&(P.DoAllInitialInfectioninSameLoc), 1, 1, 0)) P.DoAllInitialInfectioninSameLoc = 0;
	if (GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Day of year of start of seeding", "%lf", (void*)&(P.InitialInfectionCalTime), 1, 1, 0))
	{
		if (GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Scaling of infection seeding", "%lf", (void*)&(P.SeedingScaling), 1, 1, 0))
			P.DoNoCalibration = 1;
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
	if (P.FitIter == 0)
	{
		if (P.DoAdUnits)
		{
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Administrative unit to seed initial infection into", "%s", (P.NumSeedLocations > 1) ? ((void*)AdunitListNames) : ((void*)AdunitListNames[0]), P.NumSeedLocations, 1, 0))
				for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
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
			if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Administrative unit seeding weights", "%lf", (void*)&(P.InitialInfectionsAdminUnitWeight[0]), P.NumSeedLocations, 1, 0))
				for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnitWeight[i] = 1.0;
			s = 0;
			for (i = 0; i < P.NumSeedLocations; i++) s += P.InitialInfectionsAdminUnitWeight[i];
			for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnitWeight[i] /= s;
		}
		else
		{
			for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
		}
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Initial rate of importation of infections", "%lf", (void*)&(P.InfectionImportRate1), 1, 1, 0)) P.InfectionImportRate1 = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Changed rate of importation of infections", "%lf", (void*)&(P.InfectionImportRate2), 1, 1, 0)) P.InfectionImportRate2 = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Time when infection rate changes", "%lf", (void*)&(P.InfectionImportChangeTime), 1, 1, 0)) P.InfectionImportChangeTime = 1e10;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Imports via air travel", "%i", (void*)&(P.DoImportsViaAirports), 1, 1, 0)) P.DoImportsViaAirports = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Length of importation time profile provided", "%i", (void*)&(P.DurImportTimeProfile), 1, 1, 0)) P.DurImportTimeProfile = 0;
	if (P.DurImportTimeProfile > 0)
	{
		if (P.DurImportTimeProfile >= MAX_DUR_IMPORT_PROFILE) ERR_CRITICAL("MAX_DUR_IMPORT_PROFILE too small\n");
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Daily importation time profile", "%lf", (void*)P.ImportInfectionTimeProfile, P.DurImportTimeProfile, 1, 0);
	}
	GetInputParameter(ParamFile_dat, PreParamFile_dat, "Reproduction number", "%lf", (void*)&(P.R0), 1, 1, 0);
	if (GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Beta for spatial transmission", "%lf", (void*)&(P.LocalBeta), 1, 1, 0))
		P.FixLocalBeta = 1;
	else
	{
		P.LocalBeta = -1.0;
		P.FixLocalBeta = 0;
	}
	GetInputParameter(ParamFile_dat, PreParamFile_dat, "Infectious period", "%lf", (void*)&(P.InfectiousPeriod), 1, 1, 0);
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "SD of individual variation in susceptibility", "%lf", (void*)&(P.SusceptibilitySD), 1, 1, 0)) P.SusceptibilitySD = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "SD of individual variation in infectiousness", "%lf", (void*)&(P.InfectiousnessSD), 1, 1, 0)) P.InfectiousnessSD = 0;
	if (GetInputParameter2(ParamFile_dat, PreParamFile_dat, "k of individual variation in infectiousness", "%lf", (void*)&s, 1, 1, 0)) P.InfectiousnessSD = 1.0 / sqrt(s);
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "k does not apply in households", "%i", (void*)&P.NoInfectiousnessSDinHH, 1, 1, 0)) P.NoInfectiousnessSDinHH = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Model time varying infectiousness", "%i", (void*)&(P.DoInfectiousnessProfile), 1, 1, 0)) P.DoInfectiousnessProfile = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Power of scaling of spatial R0 with density", "%lf", (void*)&(P.R0DensityScalePower), 1, 1, 0)) P.R0DensityScalePower = 0;
	if (P.DoInfectiousnessProfile)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Infectiousness profile", "%lf", (void*)P.infectious_prof, INFPROF_RES, 1, 0))
		{
			for (i = 0; i < INFPROF_RES; i++)
				P.infectious_prof[i] = 1;
		}
		k = (int)ceil(P.InfectiousPeriod / P.TimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		s = 0;
		P.infectious_prof[INFPROF_RES] = 0;
		for (i = 0; i < MAX_INFECTIOUS_STEPS; i++)	P.infectiousness[i] = 0;
		for (i = 0; i < k; i++)
		{
			t = (((double)i) * P.TimeStep / P.InfectiousPeriod * INFPROF_RES);
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
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Infectious period inverse CDF", "%lf", (void*)P.infectious_icdf.get_values(), CDF_RES + 1, 1, 0))
		{
			P.infectious_icdf.set_neg_log(ICDF_START);
		}
		k = (int)ceil(P.InfectiousPeriod * P.infectious_icdf[CDF_RES] / P.TimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		for (i = 0; i < k; i++) P.infectiousness[i] = 1.0;
		P.infectiousness[k] = 0;
		P.infectious_icdf.assign_exponent();
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Include latent period", "%i", (void*)&(P.DoLatent), 1, 1, 0)) P.DoLatent = 0;
	if (P.DoLatent)
	{
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Latent period", "%lf", (void*) &(P.LatentPeriod), 1, 1, 0);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "Latent period inverse CDF", &P.latent_icdf, 1e10);
	}

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Include symptoms", "%i", (void*)&(P.DoSymptoms), 1, 1, 0)) P.DoSymptoms = 0;
	if (!P.DoSymptoms)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.ProportionSymptomatic[i] = 0;
		P.FalsePositiveRate = 0;
		P.SymptInfectiousness = P.AsymptInfectiousness = 1.0;
		P.LatentToSymptDelay = 0;
	}
	else
	{
		if (P.DoAge)
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Proportion symptomatic by age group", "%lf", (void*)P.ProportionSymptomatic, NUM_AGE_GROUPS, 1, 0);
		else
		{
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Proportion symptomatic", "%lf", (void*)P.ProportionSymptomatic, 1, 1, 0);
			for (i = 1; i < NUM_AGE_GROUPS; i++)
				P.ProportionSymptomatic[i] = P.ProportionSymptomatic[0];
		}
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Delay from end of latent period to start of symptoms", "%lf", (void*)&(P.LatentToSymptDelay), 1, 1, 0);
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Relative rate of random contacts if symptomatic", "%lf", (void*)&(P.SymptSpatialContactRate), 1, 1, 0);
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Symptomatic infectiousness relative to asymptomatic", "%lf", (void*)&(P.SymptInfectiousness), 1, 1, 0)) P.SymptInfectiousness = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Asymptomatic infectiousness relative to symptomatic", "%lf", (void*)&(P.AsymptInfectiousness), 1, 1, 0)) P.AsymptInfectiousness = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Model symptomatic withdrawal to home as true absenteeism", "%i", (void*)&P.DoRealSymptWithdrawal, 1, 1, 0)) P.DoRealSymptWithdrawal = 0;
		if (P.DoPlaces)
		{
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Relative level of place attendance if symptomatic", "%lf", (void*)P.SymptPlaceTypeContactRate, P.PlaceTypeNum, 1, 0);
			if (P.DoRealSymptWithdrawal)
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
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum age of child at home for whom one adult also stays at home", "%i", (void*)&P.CaseAbsentChildAgeCutoff, 1, 1, 0)) P.CaseAbsentChildAgeCutoff = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of children at home for whom one adult also stays at home", "%lf", (void*)&P.CaseAbsentChildPropAdultCarers, 1, 1, 0)) P.CaseAbsentChildPropAdultCarers = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place close round household", "%i", (void*)&P.PlaceCloseRoundHousehold, 1, 1, 0)) P.PlaceCloseRoundHousehold = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Absenteeism place closure", "%i", (void*)&P.AbsenteeismPlaceClosure, 1, 1, 0)) P.AbsenteeismPlaceClosure = 0;
		if (P.AbsenteeismPlaceClosure)
		{
			P.CaseAbsenteeismDelay = 0;  // Set to zero for tracking absenteeism
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Max absent time", "%i", (void*)&P.MaxAbsentTime, 1, 1, 0)) P.MaxAbsentTime = MAX_ABSENT_TIME;
			if (P.MaxAbsentTime > MAX_ABSENT_TIME || P.MaxAbsentTime < 0)
			{
				ERR_CRITICAL_FMT("[Max absent time] out of range (%d), should be in range [0, %d]", P.MaxAbsentTime, MAX_ABSENT_TIME);
			}
		}
		else
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay in starting place absenteeism for cases who withdraw", "%lf", (void*)&P.CaseAbsenteeismDelay, 1, 1, 0)) P.CaseAbsenteeismDelay = 0;
			P.MaxAbsentTime = 0; // Not used when !P.AbsenteeismPlaceClosure
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of place absenteeism for cases who withdraw", "%lf", (void*)&P.CaseAbsenteeismDuration, 1, 1, 0)) P.CaseAbsenteeismDuration = 7;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "False positive rate", "%lf", (void*)&(P.FalsePositiveRate), 1, 1, 0)) P.FalsePositiveRate = 0.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "False positive per capita incidence", "%lf", (void*)&(P.FalsePositivePerCapitaIncidence), 1, 1, 0)) P.FalsePositivePerCapitaIncidence = 0.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "False positive relative incidence by age", "%lf", (void*)P.FalsePositiveAgeRate, NUM_AGE_GROUPS, 1, 0))
			for (j = 0; j < NUM_AGE_GROUPS; j++) P.FalsePositiveAgeRate[j] = 1.0;
	}

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum sensitivity of serology assay", "%lf", (void*)&(P.SeroConvMaxSens), 1, 1, 0)) P.SeroConvMaxSens = 1.0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Seroconversion model parameter 1", "%lf", (void*)&(P.SeroConvP1), 1, 1, 0)) P.SeroConvP1 = 14.0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Seroconversion model parameter 2", "%lf", (void*)&(P.SeroConvP2), 1, 1, 0)) P.SeroConvP2 = 3.0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Specificity of serology assay", "%lf", (void*)&(P.SeroConvSpec), 1, 1, 0)) P.SeroConvSpec = 1.0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Scaling of modelled infection prevalence to match surveys", "%lf", (void*)&(P.InfPrevSurveyScale), 1, 1, 0)) P.InfPrevSurveyScale = 1.0;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Do Severity Analysis", "%i", (void*)&(P.DoSeverity), 1, 1, 0)) P.DoSeverity = 0;
	if (P.DoSeverity)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Factor to scale IFR", "%lf", (void*)&(P.ScaleIFR), 1, 1, 0)) P.ScaleIFR = 1.0;
		//// Means for icdf's.
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "MeanTimeToTest", "%lf", (void*)&(P.Mean_TimeToTest), 1, 1, 0)) P.Mean_TimeToTest = 0.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "MeanTimeToTestOffset", "%lf", (void*)&(P.Mean_TimeToTestOffset), 1, 1, 0)) P.Mean_TimeToTestOffset = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "MeanTimeToTestCriticalOffset", "%lf", (void*)&(P.Mean_TimeToTestCriticalOffset), 1, 1, 0)) P.Mean_TimeToTestCriticalOffset = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "MeanTimeToTestCritRecovOffset", "%lf", (void*)&(P.Mean_TimeToTestCritRecovOffset), 1, 1, 0)) P.Mean_TimeToTestCritRecovOffset = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Age dependent severity delays", "%i", (void*)&i, 1, 1, 0)) i = 0;
		if (!i)
		{
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_MildToRecovery", "%lf", (void*)&(P.Mean_MildToRecovery[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_ILIToRecovery", "%lf", (void*)&(P.Mean_ILIToRecovery[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_SARIToRecovery", "%lf", (void*)&(P.Mean_SARIToRecovery[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_CriticalToCritRecov", "%lf", (void*)&(P.Mean_CriticalToCritRecov[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_CritRecovToRecov", "%lf", (void*)&(P.Mean_CritRecovToRecov[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_ILIToSARI", "%lf", (void*)&(P.Mean_ILIToSARI[0]), 1, 1, 0);
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Mean_ILIToDeath", "%lf", (void*)&(P.Mean_ILIToDeath[0]), 1, 1, 0)) P.Mean_ILIToDeath[0] = 7.0;
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_SARIToCritical", "%lf", (void*)&(P.Mean_SARIToCritical[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_SARIToDeath", "%lf", (void*)&(P.Mean_SARIToDeath[0]), 1, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_CriticalToDeath", "%lf", (void*)&(P.Mean_CriticalToDeath[0]), 1, 1, 0);
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
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_MildToRecovery", "%lf", (void*)(P.Mean_MildToRecovery), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_ILIToRecovery", "%lf", (void*)(P.Mean_ILIToRecovery), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_SARIToRecovery", "%lf", (void*)(P.Mean_SARIToRecovery), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_CriticalToCritRecov", "%lf", (void*)(P.Mean_CriticalToCritRecov), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_CritRecovToRecov", "%lf", (void*)(P.Mean_CritRecovToRecov), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_ILIToSARI", "%lf", (void*)(P.Mean_ILIToSARI), NUM_AGE_GROUPS, 1, 0);
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Mean_ILIToDeath", "%lf", (void*)(P.Mean_ILIToDeath), NUM_AGE_GROUPS, 1, 0))
				for (j = 0; j < NUM_AGE_GROUPS; j++) P.Mean_ILIToDeath[j] = 7.0;
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_SARIToCritical", "%lf", (void*)(P.Mean_SARIToCritical), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_SARIToDeath", "%lf", (void*)(P.Mean_SARIToDeath), NUM_AGE_GROUPS, 1, 0);
			GetInputParameter(ParamFile_dat, PreParamFile_dat, "Mean_CriticalToDeath", "%lf", (void*)(P.Mean_CriticalToDeath), NUM_AGE_GROUPS, 1, 0);
		}

		//// Get InverseCDFs
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "MildToRecovery_icdf", &P.MildToRecovery_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "ILIToRecovery_icdf", &P.ILIToRecovery_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "ILIToDeath_icdf", &P.ILIToDeath_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "SARIToRecovery_icdf", &P.SARIToRecovery_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "CriticalToCritRecov_icdf", &P.CriticalToCritRecov_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "CritRecovToRecov_icdf", &P.CritRecovToRecov_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "ILIToSARI_icdf", &P.ILIToSARI_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "SARIToCritical_icdf", &P.SARIToCritical_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "SARIToDeath_icdf", &P.SARIToDeath_icdf);
		GetInverseCdf(ParamFile_dat, PreParamFile_dat, "CriticalToDeath_icdf", &P.CriticalToDeath_icdf);

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Prop_Mild_ByAge", "%lf", (void*)P.Prop_Mild_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_Mild_ByAge[i] = 0.5;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Prop_ILI_ByAge", "%lf", (void*)P.Prop_ILI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_ILI_ByAge[i] = 0.3;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Prop_SARI_ByAge", "%lf", (void*)P.Prop_SARI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_SARI_ByAge[i] = 0.15;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Prop_Critical_ByAge", "%lf", (void*)P.Prop_Critical_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_Critical_ByAge[i] = 0.05;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "CFR_SARI_ByAge", "%lf", (void*)P.CFR_SARI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.CFR_SARI_ByAge[i] = 0.50;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "CFR_Critical_ByAge", "%lf", (void*)P.CFR_Critical_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.CFR_Critical_ByAge[i] = 0.50;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "CFR_ILI_ByAge", "%lf", (void*)P.CFR_ILI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.CFR_ILI_ByAge[i] = 0.00;

		//Add param to allow severity to be uniformly scaled up or down.
		for (i = 0; i < NUM_AGE_GROUPS; i++)
		{
			P.Prop_SARI_ByAge[i] *= P.ScaleIFR;
			P.Prop_Critical_ByAge[i] *= P.ScaleIFR;
			P.Prop_ILI_ByAge[i] = 1.0 - P.Prop_Mild_ByAge[i] - P.Prop_SARI_ByAge[i] - P.Prop_Critical_ByAge[i];
		}
	}
	if (P.FitIter == 0)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Bounding box for bitmap", "%lf", (void*) &(P.BoundingBox.bottom_left().x), 4, 1, 0))
		{
			P.BoundingBox.bottom_left() = Geometry::Vector2d(0.0, 0.0);
			P.BoundingBox.top_right()   = Geometry::Vector2d(1.0, 1.0);
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Spatial domain for simulation", "%lf", (void*) &(P.SpatialBoundingBox.bottom_left().x), 4, 1, 0))
		{
			P.SpatialBoundingBox.bottom_left() = Geometry::Vector2d(0.0, 0.0);
			P.SpatialBoundingBox.top_right()   = Geometry::Vector2d(1.0, 1.0);
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Grid size", "%lf", (void*)&(P.in_cells_.width), 1, 1, 0)) P.in_cells_.width = 1.0 / 120.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Use long/lat coord system", "%i", (void*)&(P.DoUTM_coords), 1, 1, 0)) P.DoUTM_coords = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Bitmap scale", "%lf", (void*)&(P.BitmapScale), 1, 1, 0)) P.BitmapScale = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Bitmap y:x aspect scaling", "%lf", (void*)&(P.BitmapAspectScale), 1, 1, 0)) P.BitmapAspectScale = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Bitmap movie frame interval", "%i", (void*)&(P.BitmapMovieFrame), 1, 1, 0)) P.BitmapMovieFrame = 250;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output bitmap", "%i", (void*)&(P.OutputBitmap), 1, 1, 0)) P.OutputBitmap = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output bitmap detected", "%i", (void*)&(P.OutputBitmapDetected), 1, 1, 0)) P.OutputBitmapDetected = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output immunity on bitmap", "%i", (void*)&(P.DoImmuneBitmap), 1, 1, 0)) P.DoImmuneBitmap = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output infection tree", "%i", (void*)&(P.DoInfectionTree), 1, 1, 0)) P.DoInfectionTree = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Do one generation", "%i", (void*)&(P.DoOneGen), 1, 1, 0)) P.DoOneGen = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output every realisation", "%i", (void*)&(P.OutputEveryRealisation), 1, 1, 0)) P.OutputEveryRealisation = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number to sample for correlations", "%i", (void*)&(P.MaxCorrSample), 1, 1, 0)) P.MaxCorrSample = 1000000000;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Assume SI model", "%i", (void*)&(P.DoSI), 1, 1, 0)) P.DoSI = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Assume periodic boundary conditions", "%i", (void*)&(P.DoPeriodicBoundaries), 1, 1, 0)) P.DoPeriodicBoundaries = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Only output non-extinct realisations", "%i", (void*)&(P.OutputOnlyNonExtinct), 1, 1, 0)) P.OutputOnlyNonExtinct = 0;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Use cases per thousand threshold for area controls", "%i", (void*)&(P.DoPerCapitaTriggers), 1, 1, 0)) P.DoPerCapitaTriggers = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Use global triggers for interventions", "%i", (void*)&(P.DoGlobalTriggers), 1, 1, 0)) P.DoGlobalTriggers = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Use admin unit triggers for interventions", "%i", (void*)&(P.DoAdminTriggers), 1, 1, 0)) P.DoAdminTriggers = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Use ICU case triggers for interventions", "%i", (void*)&(P.DoICUTriggers), 1, 1, 0)) P.DoICUTriggers = 0;
		if (P.DoGlobalTriggers)  P.DoAdminTriggers = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Divisor for per-capita area threshold (default 1000)", "%i", (void*)&(P.IncThreshPop), 1, 1, 0)) P.IncThreshPop = 1000;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Divisor for per-capita global threshold (default 1000)", "%i", (void*)&(P.GlobalIncThreshPop), 1, 1, 0)) P.GlobalIncThreshPop = 1000;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of sampling intervals over which cumulative incidence measured for global trigger", "%i", (void*)&(P.TriggersSamplingInterval), 1, 1, 0)) P.TriggersSamplingInterval = 10000000;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of cases detected for treatment", "%lf", (void*)&(P.PostAlertControlPropCasesId), 1, 1, 0)) P.PostAlertControlPropCasesId = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of cases detected before outbreak alert", "%lf", (void*)&(P.PreAlertControlPropCasesId), 1, 1, 0)) P.PreAlertControlPropCasesId = 1.0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger alert on deaths", "%i", (void*)&(P.TriggerAlertOnDeaths), 1, 1, 0)) P.TriggerAlertOnDeaths = 0;
	}
	if (P.TriggerAlertOnDeaths)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of deaths accummulated before alert", "%i", (void*)&(P.CaseOrDeathThresholdBeforeAlert), 1, 1, 0)) P.CaseOrDeathThresholdBeforeAlert = 0;
	}
	else
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of detected cases needed before outbreak alert triggered", "%i", (void*)&(P.CaseOrDeathThresholdBeforeAlert), 1, 1, 0)) P.CaseOrDeathThresholdBeforeAlert = 0;
	}

	if (P.CaseOrDeathThresholdBeforeAlert_CommandLine > 0) P.CaseOrDeathThresholdBeforeAlert = P.CaseOrDeathThresholdBeforeAlert_CommandLine;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Alert trigger starts after interventions", "%i", (void*)&(P.DoAlertTriggerAfterInterv), 1, 1, 0)) P.DoAlertTriggerAfterInterv = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Day of year trigger is reached", "%lf", (void*)&(P.DateTriggerReached_CalTime), 1, 1, 0)) P.DateTriggerReached_CalTime = -1;
	if (P.DoAlertTriggerAfterInterv)
	{
		GetInputParameter(ParamFile_dat, PreParamFile_dat, "Day of year interventions start", "%lf", (void*)&(P.Interventions_StartDate_CalTime), 1, 1, 0);
		if (P.DateTriggerReached_CalTime <= P.Interventions_StartDate_CalTime)
			P.DoAlertTriggerAfterInterv = 0;
		else
		{
			P.AlertTriggerAfterIntervThreshold = P.CaseOrDeathThresholdBeforeAlert;
			P.CaseOrDeathThresholdBeforeAlert = 1000;
			fprintf(stderr, "Threshold of %i deaths by day %lg\n", P.AlertTriggerAfterIntervThreshold, P.DateTriggerReached_CalTime);
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

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of days to accummulate cases/deaths before alert", "%i", (void*)&(P.WindowToEvaluateTriggerAlert), 1, 1, 0)) P.WindowToEvaluateTriggerAlert = 1000;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Only treat mixing groups within places", "%i", (void*) & (P.DoPlaceGroupTreat), 1, 1, 0)) P.DoPlaceGroupTreat = 0;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Treatment trigger incidence per cell"				, "%lf", (void*) & (P.TreatCellIncThresh)			, 1, 1, 0)) P.TreatCellIncThresh			= 1000000000;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Case isolation trigger incidence per cell"		, "%lf", (void*) & (P.CaseIsolation_CellIncThresh)	, 1, 1, 0)) P.CaseIsolation_CellIncThresh	= P.TreatCellIncThresh;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Household quarantine trigger incidence per cell"	, "%lf", (void*) & (P.HHQuar_CellIncThresh)			, 1, 1, 0)) P.HHQuar_CellIncThresh			= P.TreatCellIncThresh;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative susceptibility of treated individual", "%lf", (void*) & (P.TreatSuscDrop), 1, 1, 0)) P.TreatSuscDrop = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative infectiousness of treated individual", "%lf", (void*) & (P.TreatInfDrop), 1, 1, 0)) P.TreatInfDrop = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of symptomatic cases resulting in death prevented by treatment", "%lf", (void*) & (P.TreatDeathDrop), 1, 1, 0)) P.TreatDeathDrop = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of symptomatic cases prevented by treatment", "%lf", (void*) & (P.TreatSympDrop), 1, 1, 0)) P.TreatSympDrop = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to treat cell", "%lf", (void*) & (P.TreatDelayMean), 1, 1, 0)) P.TreatDelayMean = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of course of treatment", "%lf", (void*) & (P.TreatCaseCourseLength), 1, 1, 0)) P.TreatCaseCourseLength = 5;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of course of prophylaxis", "%lf", (void*) & (P.TreatProphCourseLength), 1, 1, 0)) P.TreatProphCourseLength = 10;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of detected cases treated", "%lf", (void*) & (P.TreatPropCases), 1, 1, 0)) P.TreatPropCases = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of households of cases treated", "%lf", (void*) & (P.TreatPropCaseHouseholds), 1, 1, 0)) P.TreatPropCaseHouseholds = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of household prophylaxis policy", "%lf", (void*) & (P.TreatHouseholdsDuration), 1, 1, 0)) P.TreatHouseholdsDuration = USHRT_MAX / P.TimeStepsPerDay;
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion treated", "%lf", (void*) & (P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion treated in radial prophylaxis", "%lf", (void*) & (P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Treatment radius", "%lf", (void*) & (P.TreatRadius), 1, 1, 0)) P.TreatRadius = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of place/geographic prophylaxis policy", "%lf", (void*) & (P.TreatPlaceGeogDuration), 1, 1, 0)) P.TreatPlaceGeogDuration = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Treatment start time", "%lf", (void*) & (P.TreatTimeStartBase), 1, 1, 0)) P.TreatTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of places treated after case detected", "%lf", (void*)P.TreatPlaceProbCaseId, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.TreatPlaceProbCaseId[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of people treated in targeted places", "%lf", (void*)P.TreatPlaceTotalProp, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.TreatPlaceTotalProp[i] = 0;
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of doses available", "%lf", (void*) & (P.TreatMaxCoursesBase), 1, 1, 0)) P.TreatMaxCoursesBase = 1e20;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Start time of additional treatment production", "%lf", (void*) & (P.TreatNewCoursesStartTime), 1, 1, 0)) P.TreatNewCoursesStartTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Rate of additional treatment production (courses per day)", "%lf", (void*) & (P.TreatNewCoursesRate), 1, 1, 0)) P.TreatNewCoursesRate = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of people targeted with radial prophylaxis per case", "%i", (void*) & (P.TreatMaxCoursesPerCase), 1, 1, 0)) P.TreatMaxCoursesPerCase = 1000000000;


	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Treat administrative units rather than rings", "%i", (void*) & (P.TreatByAdminUnit), 1, 1, 0)) P.TreatByAdminUnit = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Administrative unit divisor for treatment", "%i", (void*) & (P.TreatAdminUnitDivisor), 1, 1, 0)) P.TreatAdminUnitDivisor = 1;
		if ((P.TreatAdminUnitDivisor == 0) || (P.TreatByAdminUnit == 0)) { P.TreatByAdminUnit = 0; P.TreatAdminUnitDivisor = 1; }
	}
	else
	{
		P.TreatAdminUnitDivisor = 1; P.TreatByAdminUnit = 0;
	}

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Vaccination trigger incidence per cell", "%lf", (void*) & (P.VaccCellIncThresh), 1, 1, 0)) P.VaccCellIncThresh = 1000000000;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative susceptibility of vaccinated individual", "%lf", (void*) & (P.VaccSuscDrop), 1, 1, 0)) P.VaccSuscDrop = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative susceptibility of individual vaccinated after switch time", "%lf", (void*) & (P.VaccSuscDrop2), 1, 1, 0)) P.VaccSuscDrop2 = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Switch time at which vaccine efficacy increases", "%lf", (void*) & (P.VaccTimeEfficacySwitch), 1, 1, 0)) P.VaccTimeEfficacySwitch = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Decay rate of vaccine efficacy (per year)", "%lf", (void*) & (P.VaccEfficacyDecay), 1, 1, 0)) P.VaccEfficacyDecay = 0;
	P.VaccEfficacyDecay /= DAYS_PER_YEAR;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative infectiousness of vaccinated individual", "%lf", (void*) & (P.VaccInfDrop), 1, 1, 0)) P.VaccInfDrop = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of symptomatic cases resulting in death prevented by vaccination", "%lf", (void*) & (P.VaccMortDrop), 1, 1, 0)) P.VaccMortDrop = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of symptomatic cases prevented by vaccination", "%lf", (void*) & (P.VaccSympDrop), 1, 1, 0)) P.VaccSympDrop = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to vaccinate", "%lf", (void*) & (P.VaccDelayMean), 1, 1, 0)) P.VaccDelayMean = 0;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay from vaccination to full protection", "%lf", (void*) & (P.VaccTimeToEfficacy), 1, 1, 0)) P.VaccTimeToEfficacy = 0;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Years between rounds of vaccination", "%lf", (void*) & (P.VaccCampaignInterval), 1, 1, 0)) P.VaccCampaignInterval = 1e10;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Max vaccine doses per day", "%i", (void*) & (P.VaccDosePerDay), 1, 1, 0)) P.VaccDosePerDay = -1;
	P.VaccCampaignInterval *= DAYS_PER_YEAR;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of rounds of vaccination", "%i", (void*) & (P.VaccMaxRounds), 1, 1, 0)) P.VaccMaxRounds = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of households of cases vaccinated", "%lf", (void*) & (P.VaccPropCaseHouseholds), 1, 1, 0)) P.VaccPropCaseHouseholds = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of household vaccination policy", "%lf", (void*) & (P.VaccHouseholdsDuration), 1, 1, 0)) P.VaccHouseholdsDuration = USHRT_MAX / P.TimeStepsPerDay;
	}

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Vaccination start time", "%lf", (void*) & (P.VaccTimeStartBase), 1, 1, 0)) P.VaccTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of population vaccinated", "%lf", (void*) & (P.VaccProp), 1, 1, 0)) P.VaccProp = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Time taken to reach max vaccination coverage (in years)", "%lf", (void*) & (P.VaccCoverageIncreasePeriod), 1, 1, 0)) P.VaccCoverageIncreasePeriod = 0;
	P.VaccCoverageIncreasePeriod *= DAYS_PER_YEAR;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Time to start geographic vaccination", "%lf", (void*) & (P.VaccTimeStartGeo), 1, 1, 0)) P.VaccTimeStartGeo = 1e10;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Vaccination radius", "%lf", (void*) & (P.VaccRadius), 1, 1, 0)) P.VaccRadius = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Minimum radius from case to vaccinate", "%lf", (void*) & (P.VaccMinRadius), 1, 1, 0)) P.VaccMinRadius = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of vaccine courses available", "%lf", (void*) & (P.VaccMaxCoursesBase), 1, 1, 0)) P.VaccMaxCoursesBase = 1e20;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Start time of additional vaccine production", "%lf", (void*) & (P.VaccNewCoursesStartTime), 1, 1, 0)) P.VaccNewCoursesStartTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "End time of additional vaccine production", "%lf", (void*) & (P.VaccNewCoursesEndTime), 1, 1, 0)) P.VaccNewCoursesEndTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Rate of additional vaccine production (courses per day)", "%lf", (void*) & (P.VaccNewCoursesRate), 1, 1, 0)) P.VaccNewCoursesRate = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Apply mass rather than reactive vaccination", "%i", (void*) & (P.DoMassVacc), 1, 1, 0)) P.DoMassVacc = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Priority age range for mass vaccination", "%i", (void*)P.VaccPriorityGroupAge, 2, 1, 0)) { P.VaccPriorityGroupAge[0] = 1; P.VaccPriorityGroupAge[1] = 0; }
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Vaccinate administrative units rather than rings", "%i", (void*) & (P.VaccByAdminUnit), 1, 1, 0)) P.VaccByAdminUnit = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Administrative unit divisor for vaccination", "%i", (void*) & (P.VaccAdminUnitDivisor), 1, 1, 0)) P.VaccAdminUnitDivisor = 1;
		if ((P.VaccAdminUnitDivisor == 0) || (P.VaccByAdminUnit == 0)) P.VaccAdminUnitDivisor = 1;
	}
	else
	{
		P.VaccAdminUnitDivisor = 1; P.VaccByAdminUnit = 0;
	}

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Movement restrictions trigger incidence per cell", "%i", (void*) & (P.MoveRestrCellIncThresh), 1, 1, 0)) P.MoveRestrCellIncThresh = 1000000000;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to start movement restrictions", "%lf", (void*) & (P.MoveDelayMean), 1, 1, 0)) P.MoveDelayMean = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of movement restrictions", "%lf", (void*) & (P.MoveRestrDuration), 1, 1, 0)) P.MoveRestrDuration = 7;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual movements after restrictions", "%lf", (void*) & (P.MoveRestrEffect), 1, 1, 0)) P.MoveRestrEffect = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Minimum radius of movement restrictions", "%lf", (void*) & (P.MoveRestrRadius), 1, 1, 0)) P.MoveRestrRadius = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Movement restrictions start time", "%lf", (void*) & (P.MoveRestrTimeStartBase), 1, 1, 0)) P.MoveRestrTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Impose blanket movement restrictions", "%i", (void*) & (P.DoBlanketMoveRestr), 1, 1, 0)) P.DoBlanketMoveRestr = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Movement restrictions only once", "%i", (void*) & (P.DoMoveRestrOnceOnly), 1, 1, 0)) P.DoMoveRestrOnceOnly = 0;
	if (P.DoMoveRestrOnceOnly) P.DoMoveRestrOnceOnly = 4;
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Movement restrictions in administrative units rather than rings", "%i", (void*) & (P.MoveRestrByAdminUnit), 1, 1, 0)) P.MoveRestrByAdminUnit = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Administrative unit divisor for movement restrictions", "%i", (void*) & (P.MoveRestrAdminUnitDivisor), 1, 1, 0)) P.MoveRestrAdminUnitDivisor = 1;
		if ((P.MoveRestrAdminUnitDivisor == 0) || (P.MoveRestrByAdminUnit == 0)) P.MoveRestrAdminUnitDivisor = 1;
	}
	else
	{
		P.MoveRestrAdminUnitDivisor = 1; P.MoveRestrByAdminUnit = 0;
	}

	//Intervention delays and durations by admin unit: ggilani 16/03/20
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Include intervention delays by admin unit", "%i", (void*) & (P.DoInterventionDelaysByAdUnit), 1, 1, 0)) P.DoInterventionDelaysByAdUnit = 0;
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

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to social distancing by admin unit"			, "%lf", (void*)AdunitDelayToSocialDistance	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToSocialDistance	[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to household quarantine by admin unit"		, "%lf", (void*)AdunitDelayToHQuarantine	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToHQuarantine		[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to case isolation by admin unit"			, "%lf", (void*)AdunitDelayToCaseIsolation	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToCaseIsolation	[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to place closure by admin unit"				, "%lf", (void*)AdunitDelayToPlaceClose		, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToPlaceClose		[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of social distancing by admin unit"		, "%lf", (void*)AdunitDurationSocialDistance, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationSocialDistance	[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of household quarantine by admin unit"	, "%lf", (void*)AdunitDurationHQuarantine	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationHQuarantine		[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of case isolation by admin unit"			, "%lf", (void*)AdunitDurationCaseIsolation	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationCaseIsolation	[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of place closure by admin unit"			, "%lf", (void*)AdunitDurationPlaceClose	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationPlaceClose		[i] = 0;

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
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Include digital contact tracing", "%i", (void*) & (P.DoDigitalContactTracing), 1, 1, 0)) P.DoDigitalContactTracing = 0;
	if (P.DoDigitalContactTracing)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Digital contact tracing trigger incidence per cell", "%lf", (void*) & (P.DigitalContactTracing_CellIncThresh), 1, 1, 0)) P.DigitalContactTracing_CellIncThresh = 1000000000;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of population or households covered by digital contact tracing", "%lf", (void*) & (P.PropPopUsingDigitalContactTracing), 1, 1, 0)) P.PropPopUsingDigitalContactTracing = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of smartphone users by age", "%lf", (void*)P.ProportionSmartphoneUsersByAge, NUM_AGE_GROUPS, 1, 0))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
			{
				P.ProportionSmartphoneUsersByAge[i] = 1;
			}
		}
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Cluster digital app clusters by household", "%i", (void*) & (P.ClusterDigitalContactUsers), 1, 1, 0)) P.ClusterDigitalContactUsers = 0; // by default, don't cluster by location
		}
		else
		{
			P.ClusterDigitalContactUsers = 0;
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of digital contacts who self-isolate", "%lf", (void*) & (P.ProportionDigitalContactsIsolate), 1, 1, 0)) P.ProportionDigitalContactsIsolate = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of contacts to trace per index case", "%i", (void*)&(P.MaxDigitalContactsToTrace), 1, 1, 0)) P.MaxDigitalContactsToTrace = MAX_CONTACTS;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay between isolation of index case and contacts", "%lf", (void*) & (P.DigitalContactTracingDelay), 1, 1, 0)) P.DigitalContactTracingDelay = P.TimeStep;
		//we really need one timestep between to make sure contact is not processed before index
		if (P.DigitalContactTracingDelay == 0) P.DigitalContactTracingDelay = P.TimeStep;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Length of self-isolation for digital contacts", "%lf", (void*) & (P.LengthDigitalContactIsolation), 1, 1, 0)) P.LengthDigitalContactIsolation = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Spatial scaling factor - digital contact tracing", "%lf", (void*) & (P.ScalingFactorSpatialDigitalContacts), 1, 1, 0)) P.ScalingFactorSpatialDigitalContacts = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place scaling factor - digital contact tracing", "%lf", (void*)&(P.ScalingFactorPlaceDigitalContacts), 1, 1, 0)) P.ScalingFactorPlaceDigitalContacts = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Digital contact tracing start time", "%lf", (void*) & (P.DigitalContactTracingTimeStartBase), 1, 1, 0)) P.DigitalContactTracingTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of digital contact tracing policy", "%lf", (void*) & (P.DigitalContactTracingPolicyDuration), 1, 1, 0)) P.DigitalContactTracingPolicyDuration = 7;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output digital contact tracing", "%i", (void*) & (P.OutputDigitalContactTracing), 1, 1, 0)) P.OutputDigitalContactTracing = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output digital contact distribution", "%i", (void*)&(P.OutputDigitalContactDist), 1, 1, 0)) P.OutputDigitalContactDist = 0;

		if (P.DoInterventionDelaysByAdUnit)
		{
			double AdunitDelayToDCT[MAX_ADUNITS];
			double AdunitDurationDCT[MAX_ADUNITS];

			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to digital contact tracing by admin unit", "%lf", (void*)AdunitDelayToDCT, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToDCT[i] = 0;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of digital contact tracing by admin unit", "%lf", (void*)AdunitDurationDCT, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationDCT[i] = 0;
			for (i = 0; i < P.NumAdunits; i++)
			{
				AdUnits[i].DCTDelay = AdunitDelayToDCT[i];
				AdUnits[i].DCTDuration = AdunitDurationDCT[i];
			}
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Isolate index cases in digital contact tracing", "%i", (void*)&(P.DCTIsolateIndexCases), 1, 1, 0)) P.DCTIsolateIndexCases = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual contacts after digital contact tracing isolation", "%lf", (void*)&(P.DCTCaseIsolationEffectiveness), 1, 1, 0)) P.DCTCaseIsolationEffectiveness = P.CaseIsolationEffectiveness;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual household contacts after digital contact tracing isolation", "%lf", (void*)&(P.DCTCaseIsolationHouseEffectiveness), 1, 1, 0)) P.DCTCaseIsolationHouseEffectiveness = P.CaseIsolationHouseEffectiveness;
		//initialise total number of users to 0
		P.NDigitalContactUsers = 0;
		P.NDigitalHouseholdUsers = 0;

		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay between symptom onset and isolation for index case", "%lf", (void*)&(P.DelayFromIndexCaseDetectionToDCTIsolation), 1, 1, 0)) P.DelayFromIndexCaseDetectionToDCTIsolation = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Test index cases and contacts", "%i", (void*)&(P.DoDCTTest), 1, 1, 0)) P.DoDCTTest = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to test index case", "%lf", (void*)&(P.DelayToTestIndexCase), 1, 1, 0)) P.DelayToTestIndexCase = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to test DCT contacts", "%lf", (void*)&(P.DelayToTestDCTContacts), 1, 1, 0)) P.DelayToTestDCTContacts = 7;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Testing specificity - DCT", "%lf", (void*)&(P.SpecificityDCT), 1, 1, 0)) P.SpecificityDCT = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Testing sensitivity - DCT", "%lf", (void*)&(P.SensitivityDCT), 1, 1, 0)) P.SensitivityDCT = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Find contacts of digital contacts", "%i", (void*)&(P.FindContactsOfDCTContacts), 1, 1, 0)) P.FindContactsOfDCTContacts = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Remove contacts of a negative index case", "%i", (void*)&(P.RemoveContactsOfNegativeIndexCase), 1, 1, 0)) P.RemoveContactsOfNegativeIndexCase = 0;
	}
	else
	{
		//Set these to 1 so it doesn't interfere with code if we aren't using digital contact tracing.

		P.ScalingFactorSpatialDigitalContacts = 1;
		P.ScalingFactorPlaceDigitalContacts = 1;
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** PLACE CLOSURE
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****


	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for place closure", "%i", (void*) & (P.PlaceCloseCellIncThresh1), 1, 1, 0)) P.PlaceCloseCellIncThresh1 = 1000000000;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for second place closure", "%i", (void*)&(P.PlaceCloseCellIncThresh2), 1, 1, 0)) P.PlaceCloseCellIncThresh2 = 1000000000;
	if (P.PlaceCloseCellIncThresh1 < 0) P.PlaceCloseCellIncThresh1 = 1000000000;
	if (P.PlaceCloseCellIncThresh2 < 0) P.PlaceCloseCellIncThresh2 = 1000000000;
	if(!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for end of place closure", "%i", (void*) & (P.PlaceCloseCellIncStopThresh), 1, 1, 0)) P.PlaceCloseCellIncStopThresh = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to start place closure", "%lf", (void*) & (P.PlaceCloseDelayMean), 1, 1, 0)) P.PlaceCloseDelayMean = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of place closure", "%lf", (void*) & (P.PlaceCloseDurationBase), 1, 1, 0)) P.PlaceCloseDurationBase = 7;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of second place closure", "%lf", (void*) & (P.PlaceCloseDuration2), 1, 1, 0)) P.PlaceCloseDuration2 = 7;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of places remaining open after closure by place type", "%lf", (void*)P.PlaceCloseEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.PlaceCloseEffect[i] = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportional attendance after closure by place type", "%lf", (void*)P.PlaceClosePropAttending, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.PlaceClosePropAttending[i] = 0;
	}
	if (P.DoHouseholds)
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rate after closure", "%lf", (void*)& P.PlaceCloseHouseholdRelContact, 1, 1, 0)) P.PlaceCloseHouseholdRelContact = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rate after closure", "%lf", (void*)& P.PlaceCloseSpatialRelContact, 1, 1, 0)) P.PlaceCloseSpatialRelContact = 1;

	if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Include holidays", "%i", (void*) & (P.DoHolidays), 1, 1, 0)) P.DoHolidays = 0;
	if (P.DoHolidays)
	{
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Proportion of places remaining open during holidays by place type", "%lf", (void*)P.HolidayEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.HolidayEffect[i] = 1;
		if (!GetInputParameter2(PreParamFile_dat, AdminFile_dat, "Number of holidays", "%i", (void*) & (P.NumHolidays), 1, 1, 0)) P.NumHolidays = 0;
		if (P.NumHolidays > DAYS_PER_YEAR) P.NumHolidays = DAYS_PER_YEAR;
		if (P.NumHolidays > 0)
		{
			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Holiday start times", "%lf", (void*)P.HolidayStartTime, P.NumHolidays, 1, 0);
			GetInputParameter(PreParamFile_dat, AdminFile_dat, "Holiday durations", "%lf", (void*)P.HolidayDuration, P.NumHolidays, 1, 0);
		}
	}
	else
		P.NumHolidays = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Minimum radius for place closure", "%lf", (void*) & (P.PlaceCloseRadius), 1, 1, 0)) P.PlaceCloseRadius = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure start time", "%lf", (void*) & (P.PlaceCloseTimeStartBase), 1, 1, 0)) P.PlaceCloseTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure second start time", "%lf", (void*) & (P.PlaceCloseTimeStartBase2), 1, 1, 0)) P.PlaceCloseTimeStartBase2 = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Places close only once", "%i", (void*) & (P.DoPlaceCloseOnceOnly), 1, 1, 0)) P.DoPlaceCloseOnceOnly = 0;
	if (P.DoPlaceCloseOnceOnly) P.DoPlaceCloseOnceOnly = 4;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure incidence threshold", "%i", (void*) & (P.PlaceCloseIncTrig1), 1, 1, 0)) P.PlaceCloseIncTrig1 = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure second incidence threshold", "%i", (void*)&(P.PlaceCloseIncTrig2), 1, 1, 0)) P.PlaceCloseIncTrig2 = P.PlaceCloseIncTrig1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure fractional incidence threshold", "%lf", (void*) & (P.PlaceCloseFracIncTrig), 1, 1, 0)) P.PlaceCloseFracIncTrig = 0;
	if ((P.DoAdUnits) && (P.DoPlaces))
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure in administrative units rather than rings", "%i", (void*) & (P.PlaceCloseByAdminUnit), 1, 1, 0)) P.PlaceCloseByAdminUnit = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Administrative unit divisor for place closure", "%i", (void*) & (P.PlaceCloseAdminUnitDivisor), 1, 1, 0)) P.PlaceCloseAdminUnitDivisor = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place types to close for admin unit closure (0/1 array)", "%i", (void*) & (P.PlaceCloseAdunitPlaceTypes), P.PlaceTypeNum, 1, 0))
			for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceCloseAdunitPlaceTypes[i] = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Cumulative proportion of place members needing to become sick for admin unit closure", "%lf", (void*) & (P.PlaceCloseCasePropThresh), 1, 1, 0)) P.PlaceCloseCasePropThresh = 2;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of places in admin unit needing to pass threshold for place closure", "%lf", (void*) & (P.PlaceCloseAdunitPropThresh), 1, 1, 0)) P.PlaceCloseAdunitPropThresh = 2;
		if ((P.PlaceCloseAdminUnitDivisor < 1) || (P.PlaceCloseByAdminUnit == 0)) P.PlaceCloseAdminUnitDivisor = 1;
	}
	else
	{
		P.PlaceCloseAdminUnitDivisor = 1; P.PlaceCloseByAdminUnit = 0;
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** SOCIAL DISTANCING
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for social distancing", "%i", (void*) & (P.SocDistCellIncThresh), 1, 1, 0)) P.SocDistCellIncThresh = 1000000000;
	if(!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for end of social distancing", "%i", (void*) & (P.SocDistCellIncStopThresh), 1, 1, 0)) P.SocDistCellIncStopThresh = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of social distancing", "%lf", (void*) & (P.SocDistDuration), 1, 1, 0)) P.SocDistDuration = 7;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of social distancing after change", "%lf", (void*) & (P.SocDistDuration2), 1, 1, 0)) P.SocDistDuration2 = 7;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative place contact rate given social distancing by place type", "%lf", (void*)P.SocDistPlaceEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.SocDistPlaceEffect[i] = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative place contact rate given enhanced social distancing by place type", "%lf", (void*)P.EnhancedSocDistPlaceEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.EnhancedSocDistPlaceEffect[i] = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative place contact rate given social distancing by place type after change", "%lf", (void*)P.SocDistPlaceEffect2, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.SocDistPlaceEffect2[i] = P.SocDistPlaceEffect[i];
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative place contact rate given enhanced social distancing by place type after change", "%lf", (void*)P.EnhancedSocDistPlaceEffect2, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.EnhancedSocDistPlaceEffect2[i] = P.EnhancedSocDistPlaceEffect[i];
	}
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rate given social distancing", "%lf", (void*)&P.SocDistHouseholdEffect, 1, 1, 0)) P.SocDistHouseholdEffect = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rate given enhanced social distancing", "%lf", (void*)&P.EnhancedSocDistHouseholdEffect, 1, 1, 0)) P.EnhancedSocDistHouseholdEffect = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rate given social distancing  after change", "%lf", (void*)&P.SocDistHouseholdEffect2, 1, 1, 0)) P.SocDistHouseholdEffect2 = P.SocDistHouseholdEffect;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rate given enhanced social distancing after change", "%lf", (void*)&P.EnhancedSocDistHouseholdEffect2, 1, 1, 0)) P.EnhancedSocDistHouseholdEffect2 = P.EnhancedSocDistHouseholdEffect;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Cluster compliance with enhanced social distancing by household", "%i", (void*)&P.EnhancedSocDistClusterByHousehold, 1, 1, 0)) P.EnhancedSocDistClusterByHousehold = 0;
	}
	else
		P.EnhancedSocDistClusterByHousehold = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rate given social distancing", "%lf", (void*)& P.SocDistSpatialEffect, 1, 1, 0)) P.SocDistSpatialEffect = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rate given social distancing after change", "%lf", (void*)&P.SocDistSpatialEffect2, 1, 1, 0)) P.SocDistSpatialEffect2 = P.SocDistSpatialEffect;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Minimum radius for social distancing", "%lf", (void*) & (P.SocDistRadius), 1, 1, 0)) P.SocDistRadius = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Social distancing start time", "%lf", (void*) & (P.SocDistTimeStartBase), 1, 1, 0)) P.SocDistTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay for change in effectiveness of social distancing", "%lf", (void*)&(P.SocDistChangeDelay), 1, 1, 0)) P.SocDistChangeDelay = USHRT_MAX / P.TimeStepsPerDay;
	if(!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion compliant with enhanced social distancing by age group", "%lf", (void*)P.EnhancedSocDistProportionCompliant, NUM_AGE_GROUPS, 1, 0))
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion compliant with enhanced social distancing", "%lf", (void*)&t, 1, 1, 0)) t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.EnhancedSocDistProportionCompliant[i] = t;
	}
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rate given enhanced social distancing", "%lf", (void*)& P.EnhancedSocDistSpatialEffect, 1, 1, 0)) P.EnhancedSocDistSpatialEffect = 1;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rate given enhanced social distancing after change", "%lf", (void*)&P.EnhancedSocDistSpatialEffect2, 1, 1, 0)) P.EnhancedSocDistSpatialEffect2 = P.EnhancedSocDistSpatialEffect;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Social distancing only once", "%i", (void*) & (P.DoSocDistOnceOnly), 1, 1, 0)) P.DoSocDistOnceOnly = 0;
	if (P.DoSocDistOnceOnly) P.DoSocDistOnceOnly = 4;

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Airport closure effectiveness", "%lf", (void*) & (P.AirportCloseEffectiveness), 1, 1, 0)) P.AirportCloseEffectiveness = 0;
	P.AirportCloseEffectiveness = 1.0 - P.AirportCloseEffectiveness;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Airport closure start time", "%lf", (void*) & (P.AirportCloseTimeStartBase), 1, 1, 0)) P.AirportCloseTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Airport closure duration", "%lf", (void*) & (P.AirportCloseDuration), 1, 1, 0)) P.AirportCloseDuration = USHRT_MAX / P.TimeStepsPerDay;

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** HOUSEHOLD QUARANTINE
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Retrigger household quarantine with each new case in quarantine window", "%i", (void*) & (P.DoHQretrigger), 1, 1, 0)) P.DoHQretrigger =0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Household quarantine start time", "%lf", (void*) & (P.HQuarantineTimeStartBase), 1, 1, 0)) P.HQuarantineTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to start household quarantine", "%lf", (void*) & (P.HQuarantineDelay), 1, 1, 0)) P.HQuarantineDelay = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Length of time households are quarantined", "%lf", (void*) & (P.HQuarantineHouseDuration), 1, 1, 0)) P.HQuarantineHouseDuration = 0;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of household quarantine policy", "%lf", (void*) & (P.HQuarantinePolicyDuration), 1, 1, 0)) P.HQuarantinePolicyDuration = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rate after quarantine", "%lf", (void*) & (P.HQuarantineHouseEffect), 1, 1, 0)) P.HQuarantineHouseEffect = 1;
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual place contacts after household quarantine by place type", "%lf", (void*)P.HQuarantinePlaceEffect, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++) P.HQuarantinePlaceEffect[i] = 1;
		}
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual spatial contacts after household quarantine", "%lf", (void*) & (P.HQuarantineSpatialEffect), 1, 1, 0)) P.HQuarantineSpatialEffect = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Household level compliance with quarantine", "%lf", (void*) & (P.HQuarantinePropHouseCompliant), 1, 1, 0)) P.HQuarantinePropHouseCompliant = 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Individual level compliance with quarantine", "%lf", (void*) & (P.HQuarantinePropIndivCompliant), 1, 1, 0)) P.HQuarantinePropIndivCompliant = 1;
	}
	else
		P.HQuarantineTimeStartBase = 1e10;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Case isolation start time", "%lf", (void*) & (P.CaseIsolationTimeStartBase), 1, 1, 0)) P.CaseIsolationTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of detected cases isolated", "%lf", (void*) & (P.CaseIsolationProp), 1, 1, 0)) P.CaseIsolationProp = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Delay to start case isolation", "%lf", (void*) & (P.CaseIsolationDelay), 1, 1, 0)) P.CaseIsolationDelay = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of case isolation", "%lf", (void*) & (P.CaseIsolationDuration), 1, 1, 0)) P.CaseIsolationDuration = 0;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of case isolation policy", "%lf", (void*) & (P.CaseIsolationPolicyDuration), 1, 1, 0)) P.CaseIsolationPolicyDuration = 1e10;
	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual contacts after case isolation", "%lf", (void*) & (P.CaseIsolationEffectiveness), 1, 1, 0)) P.CaseIsolationEffectiveness = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual household contacts after case isolation", "%lf", (void*) & (P.CaseIsolationHouseEffectiveness), 1, 1, 0))
			P.CaseIsolationHouseEffectiveness = P.CaseIsolationEffectiveness;
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** VARIABLE EFFICACIES OVER TIME
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Vary efficacies over time", "%i", (void*) & (P.VaryEfficaciesOverTime), 1, 1, 0)) P.VaryEfficaciesOverTime = 0;
	//// **** number of change times
	if (!P.VaryEfficaciesOverTime)
	{
		P.Num_SD_ChangeTimes = 1;
		P.Num_CI_ChangeTimes = 1;
		P.Num_HQ_ChangeTimes = 1;
		P.Num_PC_ChangeTimes = 1;
		P.Num_DCT_ChangeTimes = 1;
	}
	else
	{
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of change times for levels of social distancing"		, "%i", (void*) & (P.Num_SD_ChangeTimes)	, 1, 1, 0)) P.Num_SD_ChangeTimes	= 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of change times for levels of case isolation"			, "%i", (void*) & (P.Num_CI_ChangeTimes)	, 1, 1, 0)) P.Num_CI_ChangeTimes	= 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of change times for levels of household quarantine"	, "%i", (void*) & (P.Num_HQ_ChangeTimes)	, 1, 1, 0)) P.Num_HQ_ChangeTimes	= 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of change times for levels of place closure"			, "%i", (void*) & (P.Num_PC_ChangeTimes)	, 1, 1, 0)) P.Num_PC_ChangeTimes	= 1;
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of change times for levels of digital contact tracing"	, "%i", (void*) & (P.Num_DCT_ChangeTimes)	, 1, 1, 0)) P.Num_DCT_ChangeTimes	= 1;
	}

	//// **** change times:
	//// By default, initialize first change time to zero and all subsequent change times to occur after simulation time, i.e. single value of efficacy for social distancing.
	P.SD_ChangeTimes	[0] = 0;
	P.CI_ChangeTimes	[0] = 0;
	P.HQ_ChangeTimes	[0] = 0;
	P.PC_ChangeTimes	[0] = 0;
	P.DCT_ChangeTimes	[0] = 0;
	for (int ChangeTime = 1; ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES; ChangeTime++)
	{
		P.SD_ChangeTimes	[ChangeTime] = 1e10;
		P.CI_ChangeTimes	[ChangeTime] = 1e10;
		P.HQ_ChangeTimes	[ChangeTime] = 1e10;
		P.PC_ChangeTimes	[ChangeTime] = 1e10;
		P.DCT_ChangeTimes	[ChangeTime] = 1e10;
	}
	//// Get real values from (pre)param file
	GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Change times for levels of social distancing"		, "%lf", (void*)P.SD_ChangeTimes	, P.Num_SD_ChangeTimes	, 1, 0);
	GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Change times for levels of case isolation"			, "%lf", (void*)P.CI_ChangeTimes	, P.Num_CI_ChangeTimes	, 1, 0);
	GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Change times for levels of household quarantine"	, "%lf", (void*)P.HQ_ChangeTimes	, P.Num_HQ_ChangeTimes	, 1, 0);
	GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Change times for levels of place closure"			, "%lf", (void*)P.PC_ChangeTimes	, P.Num_PC_ChangeTimes	, 1, 0);
	GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Change times for levels of digital contact tracing", "%lf", (void*)P.DCT_ChangeTimes, P.Num_DCT_ChangeTimes	, 1, 0);

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
	//// soc dist
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rates over time given social distancing"			, "%lf", (void*)P.SD_SpatialEffects_OverTime, P.Num_SD_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) P.SD_SpatialEffects_OverTime[ChangeTime] = P.SocDistSpatialEffect; //// by default, initialize to Relative spatial contact rate given social distancing
	//// enhanced soc dist
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rates over time given enhanced social distancing"	, "%lf", (void*)P.Enhanced_SD_SpatialEffects_OverTime, P.Num_SD_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) P.Enhanced_SD_SpatialEffects_OverTime[ChangeTime] = P.EnhancedSocDistSpatialEffect; //// by default, initialize to Relative spatial contact rate given enhanced social distancing
	//// case isolation
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual contacts after case isolation over time"							, "%lf", (void*)P.CI_SpatialAndPlaceEffects_OverTime, P.Num_CI_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_CI_ChangeTimes; ChangeTime++) P.CI_SpatialAndPlaceEffects_OverTime[ChangeTime] = P.CaseIsolationEffectiveness;
	//// household quarantine
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual spatial contacts over time after household quarantine"				, "%lf", (void*)P.HQ_SpatialEffects_OverTime, P.Num_HQ_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) P.HQ_SpatialEffects_OverTime[ChangeTime] = P.HQuarantineSpatialEffect;
	//// place closure
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative spatial contact rates over time after place closure"				, "%lf", (void*)P.PC_SpatialEffects_OverTime, P.Num_PC_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) P.PC_SpatialEffects_OverTime[ChangeTime] = P.PlaceCloseSpatialRelContact;
	//// digital contact tracing
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual contacts after digital contact tracing isolation over time"			, "%lf", (void*)P.DCT_SpatialAndPlaceEffects_OverTime, P.Num_DCT_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_DCT_ChangeTimes; ChangeTime++) P.DCT_SpatialAndPlaceEffects_OverTime[ChangeTime] = P.DCTCaseIsolationEffectiveness;

	///// Household contact rates over time
	if (P.DoHouseholds)
	{
		//// soc dist
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rates over time given social distancing"			, "%lf", (void*)P.SD_HouseholdEffects_OverTime, P.Num_SD_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) P.SD_HouseholdEffects_OverTime[ChangeTime] = P.SocDistHouseholdEffect;
		//// enhanced soc dist
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rates over time given enhanced social distancing"	, "%lf", (void*)P.Enhanced_SD_HouseholdEffects_OverTime, P.Num_SD_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) P.Enhanced_SD_HouseholdEffects_OverTime[ChangeTime] = P.EnhancedSocDistHouseholdEffect;
		//// case isolation
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual household contacts after case isolation over time"					, "%lf", (void*)P.CI_HouseholdEffects_OverTime, P.Num_CI_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_CI_ChangeTimes; ChangeTime++) P.CI_HouseholdEffects_OverTime[ChangeTime] = P.CaseIsolationHouseEffectiveness;
		//// household quarantine
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rates over time after quarantine"					, "%lf", (void*)P.HQ_HouseholdEffects_OverTime, P.Num_HQ_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) P.HQ_HouseholdEffects_OverTime[ChangeTime] = P.HQuarantineHouseEffect;
		//// place closure
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative household contact rates over time after place closure"				, "%lf", (void*)P.PC_HouseholdEffects_OverTime, P.Num_PC_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) P.PC_HouseholdEffects_OverTime[ChangeTime] = P.PlaceCloseHouseholdRelContact;
		//// digital contact tracing
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual household contacts after digital contact tracing isolation over time", "%lf", (void*)P.DCT_HouseholdEffects_OverTime, P.Num_DCT_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_DCT_ChangeTimes; ChangeTime++) P.DCT_HouseholdEffects_OverTime[ChangeTime] = P.DCTCaseIsolationHouseEffectiveness;
	}

	///// place contact rates over time
	if (P.DoPlaces)
	{
		//// soc dist
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative place contact rates over time given social distancing by place type", "%lf", (void*) &P.SD_PlaceEffects_OverTime[0][0], P.Num_SD_ChangeTimes * P.PlaceTypeNum, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) //// by default populate to values of P.SocDistPlaceEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.SocDistPlaceEffect[PlaceType];

		//// enhanced soc dist
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Relative place contact rates over time given enhanced social distancing by place type", "%lf", (void*) &P.Enhanced_SD_PlaceEffects_OverTime[0][0], P.Num_SD_ChangeTimes * P.PlaceTypeNum, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) //// by default populate to values of P.EnhancedSocDistPlaceEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.Enhanced_SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.EnhancedSocDistPlaceEffect[PlaceType];

		//// household quarantine
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Residual place contacts over time after household quarantine by place type", "%lf", (void*) &P.HQ_PlaceEffects_OverTime[0][0], P.Num_HQ_ChangeTimes * P.PlaceTypeNum, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) //// by default populate to values of P.HQuarantinePlaceEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.HQ_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.HQuarantinePlaceEffect[PlaceType];

		//// place closure
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of places remaining open after closure by place type over time", "%lf", (void*) &P.PC_PlaceEffects_OverTime[0][0], P.Num_PC_ChangeTimes * P.PlaceTypeNum, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) //// by default populate to values of P.PlaceCloseEffect
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.PC_PlaceEffects_OverTime[ChangeTime][PlaceType] = P.PlaceCloseEffect[PlaceType];

		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportional attendance after closure by place type over time", "%lf", (void*) &P.PC_PropAttending_OverTime[0][0], P.Num_PC_ChangeTimes * P.PlaceTypeNum, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) //// by default populate to values of P.PlaceClosePropAttending
				for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
					P.PC_PropAttending_OverTime[ChangeTime][PlaceType] = P.PlaceClosePropAttending[PlaceType];
	}


	//// ****  compliance
	//// case isolation
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of detected cases isolated over time", "%lf", (void*)P.CI_Prop_OverTime, P.Num_CI_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_CI_ChangeTimes; ChangeTime++) P.CI_Prop_OverTime[ChangeTime] = P.CaseIsolationProp;
	//// household quarantine (individual level)
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Individual level compliance with quarantine over time"	, "%lf", (void*)P.HQ_Individual_PropComply_OverTime, P.Num_HQ_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) P.HQ_Individual_PropComply_OverTime[ChangeTime] = P.HQuarantinePropIndivCompliant;
	//// household quarantine (Household level)
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Household level compliance with quarantine over time"	, "%lf", (void*)P.HQ_Household_PropComply_OverTime, P.Num_HQ_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) P.HQ_Household_PropComply_OverTime[ChangeTime] = P.HQuarantinePropHouseCompliant;
	//// digital contact tracing
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of digital contacts who self-isolate over time", "%lf", (void*)P.DCT_Prop_OverTime, P.Num_DCT_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_DCT_ChangeTimes; ChangeTime++) P.DCT_Prop_OverTime[ChangeTime] = P.ProportionDigitalContactsIsolate;
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Maximum number of contacts to trace per index case over time", "%i", (void*)P.DCT_MaxToTrace_OverTime, P.Num_DCT_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_DCT_ChangeTimes; ChangeTime++) P.DCT_MaxToTrace_OverTime[ChangeTime] = P.MaxDigitalContactsToTrace;
	if (P.DoPlaces)
	{
		//// ****  thresholds
		//// place closure (global threshold)
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure incidence threshold over time", "%lf", (void*)P.PC_IncThresh_OverTime, P.Num_PC_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) P.PC_IncThresh_OverTime[ChangeTime] = P.PlaceCloseIncTrig1;
		//// place closure (fractional global threshold)
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Place closure fractional incidence threshold over time", "%lf", (void*)P.PC_FracIncThresh_OverTime, P.Num_PC_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) P.PC_FracIncThresh_OverTime[ChangeTime] = P.PlaceCloseFracIncTrig;
		//// place closure (cell incidence threshold)
		if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for place closure over time", "%i", (void*)P.PC_CellIncThresh_OverTime, P.Num_PC_ChangeTimes, 1, 0))
			for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) P.PC_CellIncThresh_OverTime[ChangeTime] = P.PlaceCloseCellIncThresh1;
		for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) if(P.PC_CellIncThresh_OverTime[ChangeTime]<0) P.PC_CellIncThresh_OverTime[ChangeTime] = 1000000000; // allows -1 to be used as a proxy for no cell-based triggering
	}
	//// household quarantine
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Household quarantine trigger incidence per cell over time", "%lf", (void*)P.HQ_CellIncThresh_OverTime, P.Num_HQ_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++) P.HQ_CellIncThresh_OverTime[ChangeTime] = P.HHQuar_CellIncThresh;
	//// case isolation
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Case isolation trigger incidence per cell over time", "%lf", (void*)P.CI_CellIncThresh_OverTime, P.Num_CI_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_CI_ChangeTimes; ChangeTime++) P.CI_CellIncThresh_OverTime[ChangeTime] = P.CaseIsolation_CellIncThresh;
	//// soc dists
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for social distancing over time", "%i", (void*)P.SD_CellIncThresh_OverTime, P.Num_SD_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++) P.SD_CellIncThresh_OverTime[ChangeTime] = P.SocDistCellIncThresh;

	//// **** Durations (later add Case isolation and Household quarantine)
	// place closure
	if (!P.VaryEfficaciesOverTime || !GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of place closure over time", "%lf", (void*)P.PC_Durs_OverTime, P.Num_PC_ChangeTimes, 1, 0))
		for (int ChangeTime = 0; ChangeTime < P.Num_PC_ChangeTimes; ChangeTime++) P.PC_Durs_OverTime[ChangeTime] = P.PlaceCloseDurationBase;

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
	if(P.FitIter==0)
	{
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of key workers randomly distributed in the population", "%i", (void*) & (P.KeyWorkerPopNum), 1, 1, 0)) P.KeyWorkerPopNum = 0;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Number of key workers in different places by place type", "%i", (void*)P.KeyWorkerPlaceNum, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++) P.KeyWorkerPlaceNum[i] = 0;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of staff who are key workers per chosen place by place type", "%lf", (void*)P.KeyWorkerPropInKeyPlaces, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++) P.KeyWorkerPropInKeyPlaces[i] = 1.0;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Trigger incidence per cell for key worker prophylaxis", "%i", (void*) & (P.KeyWorkerProphCellIncThresh), 1, 1, 0)) P.KeyWorkerProphCellIncThresh = 1000000000;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Key worker prophylaxis start time", "%lf", (void*) & (P.KeyWorkerProphTimeStartBase), 1, 1, 0)) P.KeyWorkerProphTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Duration of key worker prophylaxis", "%lf", (void*) & (P.KeyWorkerProphDuration), 1, 1, 0)) P.KeyWorkerProphDuration = 0;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Time interval from start of key worker prophylaxis before policy restarted", "%lf", (void*) & (P.KeyWorkerProphRenewalDuration), 1, 1, 0)) P.KeyWorkerProphRenewalDuration = P.KeyWorkerProphDuration;
			if (P.DoHouseholds)
			{
				if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Proportion of key workers whose households are also treated as key workers", "%lf", (void*) & (P.KeyWorkerHouseProp), 1, 1, 0)) P.KeyWorkerHouseProp = 0;
			}
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Minimum radius for key worker prophylaxis", "%lf", (void*) & (P.KeyWorkerProphRadius), 1, 1, 0)) P.KeyWorkerProphRadius = 0;
		}
		else
		{
			P.KeyWorkerPopNum = 0;
			P.KeyWorkerProphTimeStartBase = 1e10;
		}

		//Added this to parameter list so that recording infection events (and the number to record) can easily be turned off and on: ggilani - 10/10/2014
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Record infection events", "%i", (void*) & (P.DoRecordInfEvents), 1, 1, 0)) P.DoRecordInfEvents = 0;
		if (P.DoRecordInfEvents)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Max number of infection events to record", "%i", (void*) & (P.MaxInfEvents), 1, 1, 0)) P.MaxInfEvents = 1000;
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Record infection events per run", "%i", (void*) & (P.RecordInfEventsPerRun), 1, 1, 0)) P.RecordInfEventsPerRun = 0;
		}
		else
		{
			P.MaxInfEvents = 0;
		}
		//Include a limit to the number of infections to simulate, if this happens before time runs out
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Limit number of infections", "%i", (void*) & (P.LimitNumInfections), 1, 1, 0)) P.LimitNumInfections = 0;
		if (P.LimitNumInfections)
		{
			if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Max number of infections", "%i", (void*) & (P.MaxNumInfections), 1, 1, 0)) P.MaxNumInfections = 60000;
		}
		//Add origin-destination matrix parameter
		if (!GetInputParameter2(ParamFile_dat, PreParamFile_dat, "Output origin destination matrix", "%i", (void*) & (P.DoOriginDestinationMatrix), 1, 1, 0)) P.DoOriginDestinationMatrix = 0;

		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Mean child age gap", "%i", (void*) & (P.MeanChildAgeGap), 1, 1, 0)) P.MeanChildAgeGap=2;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Min adult age", "%i", (void*)&(P.MinAdultAge), 1, 1, 0)) P.MinAdultAge = 19;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Max MF partner age gap", "%i", (void*) & (P.MaxMFPartnerAgeGap), 1, 1, 0)) P.MaxMFPartnerAgeGap = 5;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Max FM partner age gap", "%i", (void*) & (P.MaxFMPartnerAgeGap), 1, 1, 0)) P.MaxFMPartnerAgeGap = 5;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Min parent age gap", "%i", (void*) & (P.MinParentAgeGap), 1, 1, 0)) P.MinParentAgeGap = 19;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Max parent age gap", "%i", (void*) & (P.MaxParentAgeGap), 1, 1, 0)) P.MaxParentAgeGap = 44;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Max child age", "%i", (void*) & (P.MaxChildAge), 1, 1, 0)) P.MaxChildAge = 20;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "One Child Two Pers Prob", "%lf", (void*) & (P.OneChildTwoPersProb), 1, 1, 0)) P.OneChildTwoPersProb = 0.08;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Two Child Three Pers Prob", "%lf", (void*) & (P.TwoChildThreePersProb), 1, 1, 0)) P.TwoChildThreePersProb = 0.11;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "One Pers House Prob Old", "%lf", (void*) & (P.OnePersHouseProbOld), 1, 1, 0)) P.OnePersHouseProbOld = 0.5;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Two Pers House Prob Old", "%lf", (void*) & (P.TwoPersHouseProbOld), 1, 1, 0)) P.TwoPersHouseProbOld = 0.5;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "One Pers House Prob Young", "%lf", (void*) & (P.OnePersHouseProbYoung), 1, 1, 0)) P.OnePersHouseProbYoung = 0.23;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Two Pers House Prob Young", "%lf", (void*) & (P.TwoPersHouseProbYoung), 1, 1, 0)) P.TwoPersHouseProbYoung = 0.23;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "One Child Prob Youngest Child Under Five", "%lf", (void*) & (P.OneChildProbYoungestChildUnderFive), 1, 1, 0)) P.OneChildProbYoungestChildUnderFive = 0.5;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Two Children Prob Youngest Under Five", "%lf", (void*) & (P.TwoChildrenProbYoungestUnderFive), 1, 1, 0)) P.TwoChildrenProbYoungestUnderFive = 0.0;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Prob Youngest Child Under Five", "%lf", (void*) & (P.ProbYoungestChildUnderFive), 1, 1, 0)) P.ProbYoungestChildUnderFive = 0;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Zero Child Three Pers Prob", "%lf", (void*) & (P.ZeroChildThreePersProb), 1, 1, 0)) P.ZeroChildThreePersProb = 0.25;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "One Child Four Pers Prob", "%lf", (void*) & (P.OneChildFourPersProb), 1, 1, 0)) P.OneChildFourPersProb = 0.2;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Young And Single Slope", "%lf", (void*) & (P.YoungAndSingleSlope), 1, 1, 0)) P.YoungAndSingleSlope = 0.7;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Young And Single", "%i", (void*) & (P.YoungAndSingle), 1, 1, 0)) P.YoungAndSingle = 36;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "No Child Pers Age", "%i", (void*) & (P.NoChildPersAge), 1, 1, 0)) P.NoChildPersAge = 44;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Old Pers Age", "%i", (void*) & (P.OldPersAge), 1, 1, 0)) P.OldPersAge = 60;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Three Child Five Pers Prob", "%lf", (void*) & (P.ThreeChildFivePersProb), 1, 1, 0)) P.ThreeChildFivePersProb = 0.5;
		if (!GetInputParameter2all(ParamFile_dat, PreParamFile_dat,AdminFile_dat, "Older Gen Gap", "%i", (void*) & (P.OlderGenGap), 1, 1, 0)) P.OlderGenGap = 19;
	}
	// Close input files.
	fclose(ParamFile_dat);
	fclose(PreParamFile_dat);
	fclose(AdminFile_dat);

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
		for (int j = 0; j < P.PlaceTypeNum; j++) P.PlaceTypeTrans[j] *= P.R0scale;
		fprintf(stderr, "Rescaled transmission coefficients by factor of %lg\n", P.R0scale);
	}
	fprintf(stderr, "Parameters read\n");
}
void ReadInterventions(std::string const& IntFile)
{
	FILE* dat;
	double r, s, startt, stopt;
	int j, k, au, ni, f, nsr;
	char buf[65536], txt[65536];
	Intervention CurInterv;

	fprintf(stderr, "Reading intervention file.\n");
	if (!(dat = fopen(IntFile.c_str(), "rb"))) ERR_CRITICAL("Unable to open intervention file\n");
	if(fscanf(dat, "%*[^<]") != 0) { // needs to be separate line because start of file
        ERR_CRITICAL("fscanf failed in ReadInterventions\n");
    }
	if(fscanf(dat, "<%[^>]", txt) != 1) {
        ERR_CRITICAL("fscanf failed in ReadInterventions\n");
    }
	if (strcmp(txt, "\?xml version=\"1.0\" encoding=\"ISO-8859-1\"\?") != 0) ERR_CRITICAL("Intervention file not XML.\n");
	if(fscanf(dat, "%*[^<]<%[^>]", txt) != 1) {
        ERR_CRITICAL("fscanf failed in ReadInterventions\n");
    }
	if (strcmp(txt, "InterventionSettings") != 0) ERR_CRITICAL("Intervention has no top level.\n");
	ni = 0;
	while (!feof(dat))
	{
		if(fscanf(dat, "%*[^<]<%[^>]", txt) != 1) {
            ERR_CRITICAL("fscanf failed in ReadInterventions\n");
        }
		if (strcmp(txt, "intervention") == 0)
		{
			ni++;
			if(fscanf(dat, "%*[^<]<%[^>]", txt) != 1) {
                ERR_CRITICAL("fscanf failed in ReadInterventions\n");
            }
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
				sscanf(txt, "%i", &CurInterv.InterventionType);
			if (!GetXMLNode(dat, "AUThresh", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%i", &CurInterv.DoAUThresh);
			if (!GetXMLNode(dat, "StartTime", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.StartTime);
			startt = CurInterv.StartTime;
			if (!GetXMLNode(dat, "StopTime", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.StopTime);
			stopt = CurInterv.StopTime;
			if (!GetXMLNode(dat, "MinDuration", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.MinDuration);
			CurInterv.MinDuration *= DAYS_PER_YEAR;
			if (!GetXMLNode(dat, "RepeatInterval", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.RepeatInterval);
			CurInterv.RepeatInterval *= DAYS_PER_YEAR;
			if (!GetXMLNode(dat, "MaxPrevAtStart", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.StartThresholdHigh);
			if (!GetXMLNode(dat, "MinPrevAtStart", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.StartThresholdLow);
			if (!GetXMLNode(dat, "MaxPrevAtStop", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.StopThreshold);
			if (GetXMLNode(dat, "NoStartAfterMinDur", "parameters", txt, 1))
				sscanf(txt, "%i", &CurInterv.NoStartAfterMin);
			else
				CurInterv.NoStartAfterMin = 0;
			if (!GetXMLNode(dat, "Level", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%lf", &CurInterv.Level);
			if (GetXMLNode(dat, "LevelCellVar", "parameters", txt, 1))
				sscanf(txt, "%lf", &CurInterv.LevelCellVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelAUVar", "parameters", txt, 1))
				sscanf(txt, "%lf", &CurInterv.LevelAUVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelCountryVar", "parameters", txt, 1))
				sscanf(txt, "%lf", &CurInterv.LevelCountryVar);
			else
				CurInterv.LevelCellVar = 0;
			if (GetXMLNode(dat, "LevelClustering", "parameters", txt, 1))
				sscanf(txt, "%lf", &CurInterv.LevelClustering);
			else
				CurInterv.LevelClustering = 0;
			if (GetXMLNode(dat, "ControlParam", "parameters", txt, 1))
				sscanf(txt, "%lf", &CurInterv.ControlParam);
			else
				CurInterv.ControlParam = 0;
			if (GetXMLNode(dat, "TimeOffset", "parameters", txt, 1))
				sscanf(txt, "%lf", &CurInterv.TimeOffset);
			else
				CurInterv.TimeOffset = 0;

			if (!GetXMLNode(dat, "MaxRounds", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%u", &CurInterv.MaxRounds);
			if (!GetXMLNode(dat, "MaxResource", "parameters", txt, 1)) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			sscanf(txt, "%u", &CurInterv.MaxResource);
			if (GetXMLNode(dat, "NumSequentialReplicas", "parameters", txt, 1))
				sscanf(txt, "%i", &nsr);
			else
				nsr = 0;
			do {
                if(fscanf(dat, "%*[^<]<%[^>]", txt) != 1) {
                    ERR_CRITICAL("fscanf failed in ReadInterventions\n");
                }
            } while ((strcmp(txt, "/intervention") != 0) && (strcmp(txt, "/parameters") != 0) && (!feof(dat)));
			if (strcmp(txt, "/parameters") != 0) ERR_CRITICAL("Incomplete intervention parameter specification in intervention file\n");
			if(fscanf(dat, "%*[^<]<%[^>]", txt) != 1) {
                ERR_CRITICAL("fscanf failed in ReadInterventions\n");
            }
			if ((strcmp(txt, "adunits") != 0) && (strcmp(txt, "countries") != 0)) ERR_CRITICAL("Incomplete adunits/countries specification in intervention file\n");
			if (strcmp(txt, "adunits") == 0)
			{
				while (GetXMLNode(dat, "A", "adunits", buf, 0))
				{
					sscanf(buf, "%s", txt);
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
					sscanf(buf, "%s", txt);
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
			if(fscanf(dat, "%*[^<]<%[^>]", txt) != 1) {
                ERR_CRITICAL("fscanf failed in ReadInterventions\n");
            }
			if (strcmp(txt, "/intervention") != 0) ERR_CRITICAL("Incorrect intervention specification in intervention file\n");
		}
	}
	if (strcmp(txt, "/InterventionSettings") != 0) ERR_CRITICAL("Intervention has no top level closure.\n");
	fprintf(stderr, "%i interventions read\n", ni);
	fclose(dat);
}
int GetXMLNode(FILE* dat, const char* NodeName, const char* ParentName, char* Value, int ResetFilePos)
{
	// ResetFilePos=1 leaves dat cursor in same position as when function was called. 0 leaves it at end of NodeName closure
	// GetXMLNode returns 1 if NodeName found, 0 otherwise. If NodeName not found, ParentName closure must be

	char buf[65536], CloseNode[2048], CloseParent[2048];
	int CurPos, ret;

	sprintf(CloseParent, "/%s", ParentName);
	CurPos = ftell(dat);
	do
	{
		if(fscanf(dat, "%*[^<]<%[^>]", buf) != 1) {
            ERR_CRITICAL("fscanf failed in GetXMLNode");
        }
	} while ((strcmp(buf, CloseParent) != 0) && (strcmp(buf, NodeName) != 0) && (!feof(dat)));
	if (strcmp(buf, CloseParent) == 0)
		ret = 0;
	else
	{
		if (strcmp(buf, NodeName) != 0) ERR_CRITICAL("Incomplete node specification in XML file\n");
		if(fscanf(dat, ">%[^<]", buf) != 1) {
            ERR_CRITICAL("fscanf failed in GetXMLNode");
        }
		if (strlen(buf) < 2048) strcpy(Value, buf);
		//		fprintf(stderr,"# %s=%s\n",NodeName,Value);
		if(fscanf(dat, "<%[^>]", buf) != 1) {
            ERR_CRITICAL("fscanf failed in GetXMLNode");
        }
		sprintf(CloseNode, "/%s", NodeName);
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

	fprintf(stderr, "Reading airport data...\nAirports with no connections = ");
	if (!(dat = fopen(air_travel_file.c_str(), "rb"))) ERR_CRITICAL("Unable to open airport file\n");
	if(fscanf(dat, "%i %i", &P.Nairports, &P.Air_popscale) != 2) {
        ERR_CRITICAL("fscanf failed in void ReadAirTravel\n");
    }
	sc = (float)((double)P.PopSize / (double)P.Air_popscale);
	if (P.Nairports > MAX_AIRPORTS) ERR_CRITICAL("Too many airports\n");
	if (P.Nairports < 2) ERR_CRITICAL("Too few airports\n");
	buf = (float*)Memory::xcalloc(P.Nairports + 1, sizeof(float));
	Airports = (Airport*)Memory::xcalloc(P.Nairports, sizeof(Airport));
	for (i = 0; i < P.Nairports; i++)
	{
		if(fscanf(dat, "%f %f %lf", &(Airports[i].loc.x), &(Airports[i].loc.y), &traf) != 3) {
            ERR_CRITICAL("fscanf failed in void ReadAirTravel\n");
        }
		traf *= (P.AirportTrafficScale * sc);
		if (!P.SpatialBoundingBox.inside(Geometry::Vector2d(Airports[i].loc)))
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
			if(fscanf(dat, "%f", buf + j) != 1) {
                ERR_CRITICAL("fscanf failed in void ReadAirTravel\n");
            }
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
				fprintf(stderr, "#%i# ", i);
			else
				fprintf(stderr, "%i ", i);
		}
	}
	fclose(dat);
	Memory::xfree(buf);
	fprintf(stderr, "\nAirport data read OK.\n");
	for (i = 0; i < P.Nairports; i++)
	{
		/*		fprintf(stderr,"(%f %i|",Airports[i].total_traffic,Airports[i].num_connected);
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
	/*		fprintf(stderr,"%f %i ",t,k);
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
						fprintf(stderr,"<%f> ",Airports[i].prop_traffic[Airports[i].num_connected-1]);
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
	/*		fprintf(stderr,"%f) ",Airports[i].total_traffic);
	*/
	}
	fprintf(stderr, "Airport data clipped OK.\n");
	for (i = 0; i < MAX_DIST; i++) AirTravelDist[i] = 0;
	for (i = 0; i < P.Nairports; i++)
		if (Airports[i].total_traffic > 0)
		{
			for (j = 0; j < Airports[i].num_connected; j++)
			{
				k = (int)Airports[i].conn_airports[j];
				traf = floor(sqrt(dist2_raw(Airports[i].loc.x, Airports[i].loc.y, Airports[k].loc.x, Airports[k].loc.y)) / OUTPUT_DIST_SCALE);
				l = (int)traf;
				//fprintf(stderr,"%(%i) ",l);
				if (l < MAX_DIST)
					AirTravelDist[l] += Airports[i].total_traffic * Airports[i].prop_traffic[j];
			}
		}
	outname = output_file_base + ".airdist.xls";
	if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open air travel output file\n");
	fprintf(dat, "dist\tfreq\n");
	for (i = 0; i < MAX_DIST; i++)
		fprintf(dat, "%i\t%.10f\n", i, AirTravelDist[i]);
	fclose(dat);
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
		//	sprintf(buf, "%s.ge" DIRECTORY_SEPARATOR "%s.avi", OutFile, OutFile);
		//	avi = CreateAvi(buf, P.BitmapMovieFrame, NULL);
		//}
#endif
		for (unsigned p = 0; p < bmh->imagesize; p++)
		{
			bmInfected[p] = bmRecovered[p] = bmTreated[p] = 0;
		}
	}
	ns = 0;
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
			Hosts[k].inf = InfStat_Susceptible;
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
				Hosts[k].RecoveringFromCritical_time = USHRT_MAX - 1;
				Hosts[k].Severity_Current = Severity::Asymptomatic;
				Hosts[k].Severity_Final = Severity::Asymptomatic;
				Hosts[k].inf = InfStat_Susceptible;
			}
		}

#pragma omp parallel for reduction(+:nim) schedule(static,1) default(none) \
		shared(P, Cells, Hosts, Households)
	for (int tn = 0; tn < P.NumThreads; tn++)
	{
		for (int i = tn; i < P.NC; i+=P.NumThreads)
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
	for (int l = 0; l < P.NMCP; l++)
	{
		int i = (int)(McellLookup[l] - Mcells);
		Mcells[i].vacc_start_time = Mcells[i].treat_start_time = USHRT_MAX - 1;
		Mcells[i].treat_end_time = 0;
		Mcells[i].treat_trig = Mcells[i].vacc_trig = Mcells[i].vacc = Mcells[i].treat = 0;
		Mcells[i].place_trig = Mcells[i].move_trig = Mcells[i].socdist_trig = Mcells[i].keyworkerproph_trig =
			Mcells[i].placeclose = Mcells[i].moverest = Mcells[i].socdist = Mcells[i].keyworkerproph = 0;
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
	fprintf(stderr, "Finished InitModel.\n");
}

void SeedInfection(double t, int* NumSeedingInfections_byLocation, int rf, int run) //adding run number to pass it to event log
{
	/* *NumSeedingInfections_byLocation is an array of the number of seeding infections by location. During runtime, usually just a single int (given by a poisson distribution)*/
	/*rf set to 0 when initializing model, otherwise set to 1 during runtime. */

	int i /*seed location index*/;
	int j /*microcell number*/;
	int k, l /*k,l are grid coords at first, then l changed to be person within Microcell j, then k changed to be index of new infection*/;
	int m = 0/*guard against too many infections and infinite loop*/;
	int f /*range = {0, 1000}*/;
	int n /*number of seed locations?*/;

	n = ((rf == 0) ? P.NumSeedLocations : 1);
	for (i = 0; i < n; i++)
	{
		if ((!P.DoRandomInitialInfectionLoc) || ((P.DoAllInitialInfectioninSameLoc) && (rf))) //// either non-random locations, doing all initial infections in same location, and not initializing.
		{
			k = (int)(P.LocationInitialInfection[i][0] / P.in_microcells_.width);
			l = (int)(P.LocationInitialInfection[i][1] / P.in_microcells_.height);
			j = k * P.total_microcells_high_ + l;
			m = 0;
			for (k = 0; (k < NumSeedingInfections_byLocation[i]) && (m < 10000); k++)
			{
				l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))]; //// randomly choose member of microcell j. Name this member l
				if (Hosts[l].inf == InfStat_Susceptible) //// If Host l is uninfected.
				{
					if ((CalcPersonSusc(l, 0, 0, 0) > 0) && (Hosts[l].age <= P.MaxAgeForInitialInfection) && (P.CareHomeAllowInitialInfections || P.CareHomePlaceType<0 || Hosts[l].PlaceLinks[P.CareHomePlaceType]<0))
					{
						//only reset the initial location if rf==0, i.e. when initial seeds are being set, not when imported cases are being set
						if (rf == 0)
						{
							P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc.x;
							P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc.y;
						}
						Hosts[l].infector = -2;
						Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
						DoInfect(l, t, 0, run); ///// guessing this updates a number of things about person l at time t in thread 0 for this run.
						m = 0;
					}
				}
				else { k--; m++; } //// k-- means if person l chosen is already infected, go again. The m < 10000 is a guard against a) too many infections; b) an infinite loop if no more uninfected people left.
			}
		}
		else if (P.DoAllInitialInfectioninSameLoc)
		{
			f = 0;
			do
			{
				m = 0;
				do
				{
					l = (int)(ranf() * ((double)P.PopSize));
					j = Hosts[l].mcell;
					//fprintf(stderr,"%i ",AdUnits[Mcells[j].adunit].id);
				} while ((Mcells[j].n < NumSeedingInfections_byLocation[i]) || (Mcells[j].n > P.MaxPopDensForInitialInfection)
					|| (Mcells[j].n < P.MinPopDensForInitialInfection)
					|| ((P.InitialInfectionsAdminUnit[i] > 0) && ((AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[i] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));
				for (k = 0; (k < NumSeedingInfections_byLocation[i]) && (m < 10000); k++)
				{
					l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
					if (Hosts[l].inf == InfStat_Susceptible)
					{
						if ((CalcPersonSusc(l, 0, 0, 0) > 0) && (Hosts[l].age <= P.MaxAgeForInitialInfection) && (P.CareHomeAllowInitialInfections || P.CareHomePlaceType < 0 || Hosts[l].PlaceLinks[P.CareHomePlaceType] < 0))
						{
							P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc.x;
							P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc.y;
							Hosts[l].infector = -2; Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
							DoInfect(l, t, 0, run);
							m = 0;
						}
					}
					else
					{
						k--; m++;
					}
				}
				if (m)
					f++;
				else
					f = 0;
			} while ((f > 0) && (f < 1000));
		}
		else
		{
			m = 0;
			for (k = 0; (k < NumSeedingInfections_byLocation[i]) && (m < 10000); k++)
			{
				do
				{
					l = (int)(ranf() * ((double)P.PopSize));
					j = Hosts[l].mcell;
					//fprintf(stderr,"@@ %i %i ",AdUnits[Mcells[j].adunit].id, (int)(AdUnits[Mcells[j].adunit].id / P.CountryDivisor));
				} while ((Mcells[j].n == 0) || (Mcells[j].n > P.MaxPopDensForInitialInfection)
					|| (Mcells[j].n < P.MinPopDensForInitialInfection)
					|| ((P.InitialInfectionsAdminUnit[i] > 0) && ((AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[i] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));
				l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
				if (Hosts[l].inf == InfStat_Susceptible)
				{
					if ((CalcPersonSusc(l, 0, 0, 0) > 0) && (Hosts[l].age<=P.MaxAgeForInitialInfection) && (P.CareHomeAllowInitialInfections || P.CareHomePlaceType < 0 || Hosts[l].PlaceLinks[P.CareHomePlaceType] < 0))
					{
						P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc.x;
						P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc.y;
						Hosts[l].infector = -2; Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
						DoInfect(l, t, 0, run);
						m = 0;
					}
					else
					{
						k--; m++;
					}
				}
				else
				{
					k--; m++;
				}
			}
		}
	}
	if (m > 0) fprintf(stderr, "### Seeding error ###\n");
}


int RunModel(int run, std::string const& snapshot_save_file, std::string const& snapshot_load_file, std::string const& output_file_base)
{
	int j, k, l, fs, fs2, NumSeedingInfections, NumSeedingInfections_byLocation[MAX_NUM_SEED_LOCATIONS] /*Denotes either Num imported Infections given rate ir, or number false positive "infections"*/;
	double ir; // infection import rate?;
	double t, cI, lcI, t2;
	unsigned short int ts;
	int continueEvents = 1;


/*	fprintf(stderr, "Checking consistency of initial state...\n");
	int i, i2, k2;
	for (i = j = k = ni = fs2 = i2 = 0; i < P.N; i++)
	{
		if (i % 1000 == 0) fprintf(stderr, "\r*** %i              ", i);
		if (Hosts[i].inf == 0) j++;
		if ((Hosts[i].pcell < P.NC) && (Hosts[i].pcell >= 0))
		{
			if (Cells[Hosts[i].pcell].susceptible[Hosts[i].listpos] != i)
			{
				k++;
				for (l = fs = 0; (l < Cells[Hosts[i].pcell].n) && (!fs); l++)
					fs = (Cells[Hosts[i].pcell].susceptible[l] == i);
				if (!fs) NumSeedingInfections++;
			}
			else
			{
				if ((Hosts[i].listpos > Cells[Hosts[i].pcell].S - 1) && (Hosts[i].inf == InfStat_Susceptible)) i2++;
				if ((Hosts[i].listpos < Cells[Hosts[i].pcell].S + Cells[Hosts[i].pcell].L + Cells[Hosts[i].pcell].I - 1) && (abs(Hosts[i].inf) == InfStat_Recovered)) i2++;
			}
			if ((Cells[Hosts[i].pcell].S + Cells[Hosts[i].pcell].L + Cells[Hosts[i].pcell].I + Cells[Hosts[i].pcell].R + Cells[Hosts[i].pcell].D) != Cells[Hosts[i].pcell].n)
			{
				k2++;
			}
		}
		else
			fs2++;
	}
	fprintf(stderr, "\n*** susceptibles=%i\nincorrect listpos=%i\nhosts not found in cell list=%i\nincorrect cell refs=%i\nincorrect positioning in cell susc list=%i\nwrong cell totals=%i\n", j, k, NumSeedingInfections, fs2, i2, k2);
*/
	InterruptRun = 0;
	lcI = 1;
	if (!snapshot_load_file.empty())
	{
		P.ts_age = (int)(P.SnapshotLoadTime * P.TimeStepsPerDay);
		t = ((double)P.ts_age) * P.TimeStep;
	}
	else
	{
		t = 0;
		P.ts_age = 0;
	}
	fs = 1;
	fs2 = 0;

	for (ns = 1; ((ns < P.NumSamples) && (!InterruptRun)); ns++) //&&(continueEvents) <-removed this
	{

		RecordSample(t, ns - 1, output_file_base);
		CalibrationThresholdCheck(t, ns - 1);

		fprintf(stderr, "\r    t=%lg   %i    %i|%i    %i     %i [=%i]  %i (%lg %lg %lg)   %lg    ", t,
			State.S, State.L, State.I, State.R, State.D, State.S + State.L + State.I + State.R + State.D, State.cumD, State.cumT, State.cumV, State.cumVG, sqrt(State.maxRad2) / 1000); //added State.cumVG
		if (!InterruptRun)
		{
			//Only run to a certain number of infections: ggilani 28/10/14
			if (P.LimitNumInfections) continueEvents = (State.cumI < P.MaxNumInfections);

			for (j = 0; ((j < P.UpdatesPerSample) && (!InterruptRun) && (continueEvents)); j++)
			{
				ts = (unsigned short int) (P.TimeStepsPerDay * t);

				//if we are to reset random numbers after an intervention event, specific time
				if (P.ResetSeedsPostIntervention)
					if ((P.ResetSeedsFlag == 0) && (ts >= (P.TimeToResetSeeds * P.TimeStepsPerDay)))
					{
						setall(&P.nextRunSeed1, &P.nextRunSeed2);
						P.ResetSeedsFlag = 1;
					}

				if (fs)
				{
					if (P.DoAirports) TravelDepartSweep(t);
					k = (int)t;
					if (P.DurImportTimeProfile > 0)
					{
						if (k < P.DurImportTimeProfile)
							ir = P.ImportInfectionTimeProfile[k] * ((t > P.InfectionImportChangeTime) ? (P.InfectionImportRate2 / P.InfectionImportRate1) : 1.0);
						else
							ir = 0;
					}
					else	ir = (t > P.InfectionImportChangeTime) ? P.InfectionImportRate2 : P.InfectionImportRate1;
					if (ir > 0) //// if infection import rate > 0, seed some infections
					{
						for (k = NumSeedingInfections = 0; k < P.NumSeedLocations; k++) NumSeedingInfections += (NumSeedingInfections_byLocation[k] = (int)ignpoi(P.TimeStep * ir * P.InitialInfectionsAdminUnitWeight[k] * P.SeedingScaling)); //// sample number imported infections from Poisson distribution.
						if (NumSeedingInfections > 0)		SeedInfection(t, NumSeedingInfections_byLocation, 1, run);
					}
					if (P.FalsePositivePerCapitaIncidence > 0)
					{
						NumSeedingInfections = (int)ignpoi(P.TimeStep * P.FalsePositivePerCapitaIncidence * ((double)P.PopSize));
						if (NumSeedingInfections > 0)
						{
							for (k = 0; k < NumSeedingInfections; k++)
							{
								do
								{
									l = (int)(((double)P.PopSize) * ranf()); //// choose person l randomly from entire population. (but change l if while condition not satisfied?)
								} while ((abs(Hosts[l].inf) == InfStat_Dead) || (ranf() > P.FalsePositiveAgeRate[HOST_AGE_GROUP(l)]));
								DoFalseCase(l, t, ts, 0);
							}
						}
					}
					InfectSweep(t, run);  // loops over all infectious people and decides which susceptible people to infect (at household, place and spatial level), and adds them to queue. Then changes each person's various characteristics using DoInfect function.  adding run number as a parameter to infect sweep so we can track run number: ggilani - 15/10/14
					//// IncubRecoverySweep loops over all infecteds (either latent or infectious). If t is the right time, latent people moved to being infected, and infectious people moved to being clinical cases. Possibly also add them to recoveries or deaths. Add them to hospitalisation & hospitalisation discharge queues.
					if (!P.DoSI) IncubRecoverySweep(t, run);
					// If doing new contact tracing, update numbers of people under contact tracing after each time step

					if (P.DoDigitalContactTracing) DigitalContactTracingSweep(t);

					fs2 = ((P.DoDeath) || (State.L + State.I > 0) || (ir > 0) || (P.FalsePositivePerCapitaIncidence > 0));
					///// TreatSweep loops over microcells to decide which cells are treated (either with treatment, vaccine, social distancing, movement restrictions etc.). Calls DoVacc, DoPlaceClose, DoProphNoDelay etc. to change (threaded) State variables
					if (!TreatSweep(t))
					{
						if ((!fs2) && (State.L + State.I == 0) && (P.FalsePositivePerCapitaIncidence == 0)) { if ((ir == 0) && (((int)t) > P.DurImportTimeProfile)) fs = 0; }
					}
					if (P.DoAirports) TravelReturnSweep(t);
					UpdateHostClosure();
				}
				t += P.TimeStep;
				if (P.DoDeath) P.ts_age++;
				if (!snapshot_save_file.empty() && (t <= P.SnapshotSaveTime) && (t + P.TimeStep > P.SnapshotSaveTime)) SaveSnapshot(snapshot_save_file);
				if (t > P.TreatNewCoursesStartTime) P.TreatMaxCourses += P.TimeStep * P.TreatNewCoursesRate;
				if ((t > P.VaccNewCoursesStartTime) && (t < P.VaccNewCoursesEndTime)) P.VaccMaxCourses += P.TimeStep * P.VaccNewCoursesRate;
				cI = ((double)(State.S)) / ((double)P.PopSize);
				if ((lcI - cI) > 0.2)
				{
					lcI = cI;
					UpdateProbs(0);
					DoInitUpdateProbs = 1;
				}
			}
		}
	}
	if (!InterruptRun) RecordSample(t, P.NumSamples - 1, output_file_base);
	fprintf(stderr, "\nEnd of run\n");
	t2 = t + P.SampleTime;
//	if(!InterruptRun)
	while (fs)
	{
		fs = TreatSweep(t2);
		t2 += P.SampleStep;
	}
	//	fprintf(stderr,"End RunModel\n");
	if (P.DoAirports)
	{
		t2 = t;
		for (t2 = t; t2 <= t + MAX_TRAVEL_TIME; t2 += P.TimeStep)
			TravelReturnSweep(t2);
	}
/*	if (!InterruptRun)
	{
		fprintf(stderr, "Checking consistency of final state...\n");
		int i, i2, k2;
		for (i = j = k = NumSeedingInfections = fs2 = i2 = 0; i < P.PopSize; i++)
		{
			if (i % 1000 == 0) fprintf(stderr, "\r*** %i              ", i);
			if (Hosts[i].inf == 0) j++;
			if ((Hosts[i].pcell < P.NC) && (Hosts[i].pcell >= 0))
			{
				if (Cells[Hosts[i].pcell].susceptible[Hosts[i].listpos] != i)
				{
					k++;
					for (l = fs = 0; (l < Cells[Hosts[i].pcell].n) && (!fs); l++)
						fs = (Cells[Hosts[i].pcell].susceptible[l] == i);
					if (!fs) NumSeedingInfections++;
				}
				else
				{
					if ((Hosts[i].listpos > Cells[Hosts[i].pcell].S - 1) && (Hosts[i].inf == InfStat_Susceptible)) i2++;
					if ((Hosts[i].listpos < Cells[Hosts[i].pcell].S + Cells[Hosts[i].pcell].L + Cells[Hosts[i].pcell].I - 1) && (abs(Hosts[i].inf) == InfStat_Recovered)) i2++;
				}
				if ((Cells[Hosts[i].pcell].S + Cells[Hosts[i].pcell].L + Cells[Hosts[i].pcell].I + Cells[Hosts[i].pcell].R + Cells[Hosts[i].pcell].D) != Cells[Hosts[i].pcell].n)
				{
					k2++;
				}
			}
			else
				fs2++;
		}
		fprintf(stderr, "\n*** susceptibles=%i\nincorrect listpos=%i\nhosts not found in cell list=%i\nincorrect cell refs=%i\nincorrect positioning in cell susc list=%i\nwrong cell totals=%i\n", j, k, NumSeedingInfections, fs2, i2, k2);
	}
*/
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
						fprintf(stderr, "*%i %i: %i %i", i, j, Hosts[i].PlaceLinks[j], P.Nplace[j]);
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
						fprintf(stderr, "*%i %i: %i ", i, j, Hosts[i].PlaceLinks[j]);
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
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "dist");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != P.HotelPlaceType)
				fprintf(dat, "\tfreq_p%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < MAX_DIST; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != P.HotelPlaceType)
					fprintf(dat, "\t%i", PlaceDistDistrib[j][i]);
			fprintf(dat, "\n");
		}
		fclose(dat);
		outname = output_file_base + ".placesize.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "size");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != P.HotelPlaceType)
				fprintf(dat, "\tfreq_p%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < MAX_PLACE_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != P.HotelPlaceType)
					fprintf(dat, "\t%i", PlaceSizeDistrib[j][i]);
			fprintf(dat, "\n");
		}
		fclose(dat);
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
	if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "0,");
	for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "%i,", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
	fprintf(dat, "\n");
	for (i = 0; i < P.NumAdunits; i++)
	{
		fprintf(dat, "%i,", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
		for (j = 0; j < P.NumAdunits; j++)
		{
			fprintf(dat, "%.10f,", AdUnits[i].origin_dest[j]);
		}
		fprintf(dat, "\n");
	}
	fclose(dat);
}

void SaveResults(std::string const& output_file_base)
{
	int i, j;
	FILE* dat;
	std::string outname;

	if (P.OutputNonSeverity)
	{
		outname = output_file_base + ".xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tL\tI\tR\tD\tincI\tincR\tincFC\tincC\tincDC\tincTC\tincCT\tincCC\tcumT\tcumTP\tcumV\tcumVG\tExtinct\trmsRad\tmaxRad\n");//\t\t%.10f\t%.10f\t%.10f\n",P.R0household,P.R0places,P.R0spatial);
		for(i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10ft%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].L, TimeSeries[i].I,
				TimeSeries[i].R, TimeSeries[i].D, TimeSeries[i].incI,
				TimeSeries[i].incR, TimeSeries[i].incFC, TimeSeries[i].incC, TimeSeries[i].incDC, TimeSeries[i].incTC, TimeSeries[i].incCT, TimeSeries[i].incCC,
				TimeSeries[i].cumT, TimeSeries[i].cumTP, TimeSeries[i].cumV, TimeSeries[i].cumVG, TimeSeries[i].extinct, TimeSeries[i].rmsRad, TimeSeries[i].maxRad);
		}
		fclose(dat);
	}

	if ((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		outname = output_file_base + ".adunit.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tI_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tC_%s", AdUnits[i].ad_name);
		for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tDC_%s", AdUnits[i].ad_name);

		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f", TimeSeries[i].t);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", TimeSeries[i].incI_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", TimeSeries[i].incC_adunit[j]);
			for (j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", TimeSeries[i].incDC_adunit[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}

	if ((P.DoDigitalContactTracing) && (P.DoAdUnits) && (P.OutputDigitalContactTracing))
	{
		outname = output_file_base + ".digitalcontacttracing.xls"; //modifying to csv file
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
    		fprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++)
		{
			fprintf(dat, "\tincDCT_%s", AdUnits[i].ad_name);
		}
		for (i = 0; i < P.NumAdunits; i++)
		{
			fprintf(dat, "\tDCT_%s", AdUnits[i].ad_name);
		}
		fprintf(dat, "\n");
		//print actual output
		for(i=0; i<P.NumSamples; i++)
		{
			fprintf(dat, "%.10lf", TimeSeries[i].t);
			for (j = 0; j < P.NumAdunits; j++)
			{
				fprintf(dat, "\t%.10lf", TimeSeries[i].incDCT_adunit[j]);
			}
			for (j = 0; j < P.NumAdunits; j++)
			{
				fprintf(dat, "\t%.10lf", TimeSeries[i].DCT_adunit[j]);
			}
		fprintf(dat, "\n");
		}
		fclose(dat);

	}

	if (P.DoDigitalContactTracing && P.OutputDigitalContactDist)
	{
		outname = output_file_base + ".digitalcontactdist.xls"; //modifying to csv file
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		//print headers
		fprintf(dat, "nContacts\tFrequency\n");
		for (i = 0; i < (MAX_CONTACTS + 1); i++)
		{
			fprintf(dat, "%i\t%i\n", i, State.contact_dist[i]);
		}
		fclose(dat);
	}

	if(P.KeyWorkerProphTimeStartBase < P.SampleTime)
	{
		outname = output_file_base + ".keyworker.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for(i = 0; i < 2; i++) fprintf(dat, "\tI%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tC%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tT%i", i);
		fprintf(dat, "\t%i\t%i\n", P.KeyWorkerNum, P.KeyWorkerIncHouseNum);
		for(i = 0; i < P.NumSamples; i++)
			{
			fprintf(dat, "%.10f", TimeSeries[i].t);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", TimeSeries[i].incI_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", TimeSeries[i].incC_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", TimeSeries[i].cumT_keyworker[j]);
			fprintf(dat, "\n");
			}
		fclose(dat);
		}

	if(P.DoInfectionTree)
	{
		outname = output_file_base + "%s.tree.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		for(i = 0; i < P.PopSize; i++)
			if(Hosts[i].infect_type % INFECT_TYPE_MASK > 0)
				fprintf(dat, "%i\t%i\t%i\t%i\n", i, Hosts[i].infector, Hosts[i].infect_type % INFECT_TYPE_MASK, (int)HOST_AGE_YEAR(i));
		fclose(dat);
		}
#if defined(_WIN32) || defined(IMAGE_MAGICK)
	static int dm[13] ={0,31,28,31,30,31,30,31,31,30,31,30,31};
	int d, m, y, dml, f;
#ifdef _WIN32
	//if(P.OutputBitmap == 1) CloseAvi(avi);
	//if((TimeSeries[P.NumSamples - 1].extinct) && (P.OutputOnlyNonExtinct))
	//	{
	//	outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base + ".avi";
	//	DeleteFile(outname);
	//	}
#endif
	if(P.OutputBitmap >= 1 && P.BitmapFormat == BitmapFormats::PNG)
		{
		// Generate Google Earth .kml file
		outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base + ".ge.kml"; // outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base ".kml";
		if(!(dat = fopen(outname.c_str(), "wb")))
			{
			ERR_CRITICAL("Unable to open output kml file\n");
			}
		fprintf(dat, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.2\">\n<Document>\n");
		fprintf(dat, "<name>%s</name>\n", output_file_base.c_str());
		y = 2009;
		m = 1;
		d = 1;
		for(i = 0; i < P.NumSamples; i++)
			{
			fprintf(dat, "<GroundOverlay>\n<name>Snapshot %i</name>\n", i + 1);
			fprintf(dat, "<TimeSpan>\n<begin>%i-%02i-%02iT00:00:00Z</begin>\n", y, m, d);
			d += (int)P.SampleStep; // SampleStep has to be an integer here.
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
				fprintf(dat, "<end>%i-%02i-%02iT00:00:00Z</end>\n</TimeSpan>\n", y, m, d);
				outname = output_file_base + ".ge" DIRECTORY_SEPARATOR + output_file_base + "." + std::to_string(i + 1) + ".png";
				fprintf(dat, "<Icon>\n<href>%s</href>\n</Icon>\n", outname.c_str());
				fprintf(dat, "<LatLonBox>\n<north>%.10f</north>\n<south>%.10f</south>\n<east>%.10f</east>\n<west>%.10f</west>\n</LatLonBox>\n",
					P.SpatialBoundingBox.top_right().y, P.SpatialBoundingBox.bottom_left().y, P.SpatialBoundingBox.top_right().x, P.SpatialBoundingBox.bottom_left().x);
				fprintf(dat, "</GroundOverlay>\n");
			}
		fprintf(dat, "</Document>\n</kml>\n");
		fclose(dat);
		}
#endif


	if((P.DoSeverity)&&(P.OutputSeverity))
	{
		outname = output_file_base + ".severity.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open severity output file\n");
		fprintf(dat, "t\tRt\tTG\tSI\tS\tI\tR\tincI\tMild\tILI\tSARI\tCritical\tCritRecov\tincMild\tincILI\tincSARI\tincCritical\tincCritRecov\tincDeath\tincDeath_ILI\tincDeath_SARI\tincDeath_Critical\tcumMild\tcumILI\tcumSARI\tcumCritical\tcumCritRecov\tcumDeath\tcumDeath_ILI\tcumDeath_SARI\tcumDeath_Critical\n");//\t\t%.10f\t%.10f\t%.10f\n",P.R0household,P.R0places,P.R0spatial);
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				TimeSeries[i].t, TimeSeries[i].Rdenom, TimeSeries[i].meanTG, TimeSeries[i].meanSI, TimeSeries[i].S, TimeSeries[i].I, TimeSeries[i].R, TimeSeries[i].incI,
				TimeSeries[i].Mild		, TimeSeries[i].ILI		, TimeSeries[i].SARI	, TimeSeries[i].Critical	, TimeSeries[i].CritRecov	,
				TimeSeries[i].incMild	, TimeSeries[i].incILI	, TimeSeries[i].incSARI	, TimeSeries[i].incCritical	, TimeSeries[i].incCritRecov,
				TimeSeries[i].incD,	TimeSeries[i].incDeath_ILI, TimeSeries[i].incDeath_SARI, TimeSeries[i].incDeath_Critical,
				TimeSeries[i].cumMild	, TimeSeries[i].cumILI	, TimeSeries[i].cumSARI	, TimeSeries[i].cumCritical	, TimeSeries[i].cumCritRecov, TimeSeries[i].D	,
				TimeSeries[i].cumDeath_ILI, TimeSeries[i].cumDeath_SARI, TimeSeries[i].cumDeath_Critical);
		}
		fclose(dat);

		if((P.DoAdUnits) && (P.OutputSeverityAdminUnit))
		{
			//// output severity results by admin unit
			outname = output_file_base + ".severity.adunit.xls";
			if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
			fprintf(dat, "t");

			/////// ****** /////// ****** /////// ****** COLNAMES
			std::string colnames[] = {
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
					fprintf(dat, "\t%s%s", colname.c_str(), AdUnits[i].ad_name);
				}
			}

			fprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table.
			for(i = 0; i < P.NumSamples; i++)
			{
				fprintf(dat, "%.10f", TimeSeries[i].t);

				//// prevalence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].Mild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].Critical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].CritRecov_adunit[j]);

				//// incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].incDeath_Critical_adunit[j]);

				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumDeath_Critical_adunit[j]);

				if(i != P.NumSamples - 1) fprintf(dat, "\n");
			}
			fclose(dat);
		}
	}

	if (P.DoAdUnits && P.OutputAdUnitAge)
	{
		//// output infections by age and admin unit
		outname = output_file_base + ".age.adunit.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");

		// colnames
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				fprintf(dat, "\tincInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// incidence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				fprintf(dat, "\tprevInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// prevalence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				fprintf(dat, "\tcumInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// cumulative incidence
		fprintf(dat, "\n");

		// Populate
		for (int Time = 0; Time < P.NumSamples; Time++)
		{
			fprintf(dat, "%.10f", TSMean[Time].t);
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					fprintf(dat, "\t%.10f", TimeSeries[Time].incInf_age_adunit[AgeGroup][AdUnit]);	// incidence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					fprintf(dat, "\t%.10f", TimeSeries[Time].prevInf_age_adunit[AgeGroup][AdUnit]);	// prevalence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					fprintf(dat, "\t%.10f", TimeSeries[Time].cumInf_age_adunit[AgeGroup][AdUnit]);	// cumulative incidence
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
}

void SaveSummaryResults(std::string const& output_file_base) //// calculates and saves summary results (called for average of extinct and non-extinct realisation time series - look in main)
{
	int i, j;
	double c, t;
	FILE* dat;
	std::string outname;

	c = 1 / ((double)(P.NRactE + P.NRactNE));

	if (P.OutputNonSeverity)
	{
		outname = output_file_base + ".xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		//// set colnames
		fprintf(dat, "t\tS\tL\tI\tR\tD\tincI\tincR\tincD\tincC\tincDC\tincTC\tcumT\tcumTmax\tcumTP\tcumV\tcumVmax\tExtinct\trmsRad\tmaxRad\tvS\tvI\tvR\tvD\tvincI\tvincR\tvincFC\tvincC\tvincDC\tvincTC\tvrmsRad\tvmaxRad\t\t%i\t%i\t%.10f\t%.10f\t%.10f\t\t%.10f\t%.10f\t%.10f\t%.10f\n",
			P.NRactNE, P.NRactE, P.R0household, P.R0places, P.R0spatial, c * PeakHeightSum, c * PeakHeightSS - c * c * PeakHeightSum * PeakHeightSum, c * PeakTimeSum, c * PeakTimeSS - c * c * PeakTimeSum * PeakTimeSum);
		c = 1 / ((double)P.NRactual);

		//// populate table
		for(i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].L, c * TSMean[i].I, c * TSMean[i].R,
				c * TSMean[i].D, c * TSMean[i].incI, c * TSMean[i].incR, c * TSMean[i].incFC, c * TSMean[i].incC, c * TSMean[i].incDC, c * TSMean[i].incTC,
				c * TSMean[i].cumT, TSMean[i].cumTmax, c * TSMean[i].cumTP, c * TSMean[i].cumV, TSMean[i].cumVmax, c * TSMean[i].extinct, c * TSMean[i].rmsRad, c * TSMean[i].maxRad);
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
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
		fclose(dat);
	}

	if (P.OutputControls)
	{
		outname = output_file_base + ".controls.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tincC\tincTC\tincFC\tcumT\tcumUT\tcumTP\tcumV\tincHQ\tincAC\tincAH\tincAA\tincACS\tincAPC\tincAPA\tincAPCS\tpropSocDist");
		for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tprClosed_%i", j);
		fprintf(dat, "t\tvS\tvincC\tvincTC\tvincFC\tvcumT\tvcumUT\tvcumTP\tvcumV");
		for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tvprClosed_%i", j);
		fprintf(dat, "\n");
		for(i = 0; i < P.NumSamples; i++)
			{
			fprintf(dat, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].incC, c * TSMean[i].incTC, c * TSMean[i].incFC,
				c * TSMean[i].cumT, c * TSMean[i].cumUT, c * TSMean[i].cumTP, c * TSMean[i].cumV, c * TSMean[i].incHQ,
				c * TSMean[i].incAC, c * TSMean[i].incAH, c * TSMean[i].incAA, c * TSMean[i].incACS,
				c * TSMean[i].incAPC, c * TSMean[i].incAPA, c * TSMean[i].incAPCS,c*TSMean[i].PropSocDist);
			for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\t%lf", c * TSMean[i].PropPlacesClosed[j]);
			fprintf(dat, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				c * TSVar[i].S - c * c * TSMean[i].S * TSMean[i].S,
				c * TSVar[i].incC - c * c * TSMean[i].incC * TSMean[i].incC,
				c * TSVar[i].incTC - c * c * TSMean[i].incTC * TSMean[i].incTC,
				c * TSVar[i].incFC - c * c * TSMean[i].incFC * TSMean[i].incFC,
				c * TSVar[i].cumT - c * c * TSMean[i].cumT * TSMean[i].cumT,
				c * TSVar[i].cumUT - c * c * TSMean[i].cumUT * TSMean[i].cumUT,
				c * TSVar[i].cumTP - c * c * TSMean[i].cumTP * TSMean[i].cumTP,
				c * TSVar[i].cumV - c * c * TSMean[i].cumV * TSMean[i].cumV);
			for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\t%lf", TSVar[i].PropPlacesClosed[j]);
			fprintf(dat, "\n");
			}
		fclose(dat);

	}

	if (P.OutputAge)
	{
		outname = output_file_base + ".age.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for(i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "\tI%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for(i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "\tC%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		for(i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "\tD%i-%i", AGE_GROUP_WIDTH * i, AGE_GROUP_WIDTH * (i + 1));
		fprintf(dat, "\n");
		for(i = 0; i < P.NumSamples; i++)
			{
			fprintf(dat, "%.10f", c * TSMean[i].t);
			for(j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incIa[j]);
			for(j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incCa[j]);
			for(j = 0; j < NUM_AGE_GROUPS; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incDa[j]);
			fprintf(dat, "\n");
			}
		fprintf(dat, "dist");
		for(j = 0; j < NUM_AGE_GROUPS; j++)
			fprintf(dat, "\t%.10f", AgeDist[j]);
		fprintf(dat, "\n");
		fclose(dat);
	}

	if((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		outname = output_file_base + ".adunit.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tI_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tC_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tDC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tT_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\t%i", AdUnits[i].n);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\t%.10f", P.PopByAdunit[i][1]);
		fprintf(dat, "\n");
		for(i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f", c * TSMean[i].t);
			for(j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incI_adunit[j]);
			for(j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incC_adunit[j]);
			for(j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
			for(j = 0; j < P.NumAdunits; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].cumT_adunit[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);

		if (P.OutputAdUnitVar)
		{
			outname = output_file_base + ".adunitVar.xls";
			if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
			fprintf(dat, "t");
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tI_%s", AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tC_%s", AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tDC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tT_%s", AdUnits[i].ad_name);
			fprintf(dat, "\n");
			for (i = 0; i < P.NumSamples; i++)
			{
				fprintf(dat, "%.10f", c * TSMean[i].t);
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "\t%.10f", c * TSVar[i].incI_adunit[j] - c * c * TSMean[i].incI_adunit[j] * TSMean[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "\t%.10f", c * TSVar[i].incC_adunit[j] - c * c * TSMean[i].incC_adunit[j] * TSMean[i].incC_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "\t%.10f", c * TSVar[i].incDC_adunit[j] - c * c * TSMean[i].incDC_adunit[j] * TSMean[i].incDC_adunit[j]); //added detected cases: ggilani 03/02/15
				for (j = 0; j < P.NumAdunits; j++)
					fprintf(dat, "\t%.10f", c * TSVar[i].cumT_adunit[j] - c * c * TSMean[i].cumT_adunit[j] * TSMean[i].cumT_adunit[j]);
				fprintf(dat, "\n");
			}
			fclose(dat);
		}
	}

	if ((P.DoDigitalContactTracing) && (P.DoAdUnits) && (P.OutputDigitalContactTracing))
	{
		outname = output_file_base + ".digitalcontacttracing.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for (i = 0; i < P.NumAdunits; i++)
		{
			fprintf(dat, "\tincDCT_%s", AdUnits[i].ad_name); // //printing headers for inc per admin unit
		}
		for (i = 0; i < P.NumAdunits; i++)
		{
			fprintf(dat, "\tDCT_%s", AdUnits[i].ad_name); // //printing headers for prevalence of digital contact tracing per admin unit
		}
		fprintf(dat, "\n");
		//print actual output
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10lf", c* TSMean[i].t);
			for (j = 0; j < P.NumAdunits; j++)
			{
				fprintf(dat, "\t%.10lf", c * TSMean[i].incDCT_adunit[j]);
			}
			for (j = 0; j < P.NumAdunits; j++)
			{
				fprintf(dat, "\t%.10lf", c * TSMean[i].DCT_adunit[j]);
			}
			fprintf(dat, "\n");
		}

		fclose(dat);

	}

	if(P.KeyWorkerProphTimeStartBase < P.SampleTime)
	{
		outname = output_file_base + ".keyworker.xls";
		if(!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for(i = 0; i < 2; i++) fprintf(dat, "\tI%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tC%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tT%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tvI%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tvC%i", i);
		for(i = 0; i < 2; i++) fprintf(dat, "\tvT%i", i);
		fprintf(dat, "\t%i\t%i\n", P.KeyWorkerNum, P.KeyWorkerIncHouseNum);
		for(i = 0; i < P.NumSamples; i++)
			{
			fprintf(dat, "%.10f", c * TSMean[i].t);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incI_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].incC_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", c * TSMean[i].cumT_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", c * TSVar[i].incI_keyworker[j] - c * c * TSMean[i].incI_keyworker[j] * TSMean[i].incI_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", c * TSVar[i].incC_keyworker[j] - c * c * TSMean[i].incC_keyworker[j] * TSMean[i].incC_keyworker[j]);
			for(j = 0; j < 2; j++)
				fprintf(dat, "\t%.10f", c * TSVar[i].cumT_keyworker[j] - c * c * TSMean[i].cumT_keyworker[j] * TSMean[i].cumT_keyworker[j]);
			fprintf(dat, "\n");
			}
		fclose(dat);
	}

	if (P.OutputInfType)
	{
		outname = output_file_base + ".inftype.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tR\tTG\tSI");
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\tRtype_%i", j);
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\tincItype_%i", j);
		for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\tRage_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lf\t%lf\t%lf\t%lf", c * TSMean[i].t, c * TSMean[i].Rdenom, c* TSMean[i].meanTG, c* TSMean[i].meanSI);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\t%lf", c * TSMean[i].Rtype[j]);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\t%lf", c * TSMean[i].incItype[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%lf", c * TSMean[i].Rage[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}

	if (P.OutputR0)
	{
		outname = output_file_base + ".R0.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_SEC_REC; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < MAX_GEN_REC; j++)
				fprintf(dat, "\t%.10f", c * indivR0_av[i][j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
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
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%i", i);
		fprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				fprintf(dat, "\t%.10f", inf_household_av[j][i] * c);
			fprintf(dat, "\n");
		}
		fprintf(dat, "\n");
		for (i = 1; i <= MAX_HOUSEHOLD_SIZE; i++)
			fprintf(dat, "\t%i", i);
		fprintf(dat, "\n");
		for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 1; j <= MAX_HOUSEHOLD_SIZE; j++)
				fprintf(dat, "\t%.10f", case_household_av[j][i] * c);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}

	if (P.OutputCountry)
	{
		outname = output_file_base + ".country.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_COUNTRIES; i++)
			fprintf(dat, "%i\t%.10f\t%.10f\n", i, infcountry_av[i] * c, infcountry_num[i] * c);
		fclose(dat);
	}

	if ((P.DoSeverity)&&(P.OutputSeverity))
	{
		//// output separate severity file (can integrate with main if need be)
		outname = output_file_base + ".severity.xls";

		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open severity output file\n");
		fprintf(dat, "t\tPropSocDist\tRt\tTG\tSI\tS\tI\tR\tincI\tincC\tMild\tILI\tSARI\tCritical\tCritRecov\tSARIP\tCriticalP\tCritRecovP\tprevQuarNotInfected\tprevQuarNotSymptomatic\tincMild\tincILI\tincSARI\tincCritical\tincCritRecov\tincSARIP\tincCriticalP\tincCritRecovP\tincDeath\tincDeath_ILI\tincDeath_SARI\tincDeath_Critical\tcumMild\tcumILI\tcumSARI\tcumCritical\tcumCritRecov\tcumDeath\tcumDeath_ILI\tcumDeath_SARI\tcumDeath_Critical\t");
		fprintf(dat, "PropSocDist_v\tRt_v\tTG_v\tSI_v\tS_v\tI_v\tR_v\tincI_v\tincC_v\tMild_v\tILI_v\tSARI_v\tCritical_v\tCritRecov_v\tincMild_v\tincILI_v\tincSARI_v\tincCritical_v\tincCritRecov_v\tincDeath_v\tincDeath_ILI_v\tincDeath_SARI_v\tincDeath_Critical_v\tcumMild_v\tcumILI_v\tcumSARI_v\tcumCritical_v\tcumCritRecov_v\tcumDeath_v\tcumDeath_ILI_v\tcumDeath_SARI_v\tcumDeath_Critical_v\n");
		double SARI, Critical, CritRecov, incSARI, incCritical, incCritRecov, sc1, sc2,sc3,sc4; //this stuff corrects bed prevalence for exponentially distributed time to test results in hospital
		sc1 = (P.Mean_TimeToTest > 0) ? exp(-1.0 / P.Mean_TimeToTest) : 0.0;
		sc2 = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestOffset / P.Mean_TimeToTest) : 0.0;
		sc3 = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCriticalOffset / P.Mean_TimeToTest) : 0.0;
		sc4 = (P.Mean_TimeToTest > 0) ? exp(-P.Mean_TimeToTestCritRecovOffset / P.Mean_TimeToTest) : 0.0;
		incSARI = incCritical = incCritRecov = 0;
		for (i = 0; i < P.NumSamples; i++)
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

			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.17f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t",
				c* TSMean[i].t, c* TSMean[i].PropSocDist, c* TSMean[i].Rdenom, c* TSMean[i].meanTG, c* TSMean[i].meanSI, c* TSMean[i].S, c* TSMean[i].I, c* TSMean[i].R, c* TSMean[i].incI, c* TSMean[i].incC,
				c* TSMean[i].Mild, c* TSMean[i].ILI, c* TSMean[i].SARI,c* TSMean[i].Critical, c* TSMean[i].CritRecov,c* (TSMean[i].SARI - SARI), c* (TSMean[i].Critical - Critical), c* (TSMean[i].CritRecov - CritRecov),
				c * TSMean[i].prevQuarNotInfected, c * TSMean[i].prevQuarNotSymptomatic,
				c * TSMean[i].incMild, c * TSMean[i].incILI, c * TSMean[i].incSARI, c * TSMean[i].incCritical, c * TSMean[i].incCritRecov, c * incSARI, c * incCritical, c * incCritRecov, c * TSMean[i].incD,
				c * TSMean[i].incDeath_ILI, c * TSMean[i].incDeath_SARI, c * TSMean[i].incDeath_Critical,
				c * TSMean[i].cumMild, c * TSMean[i].cumILI, c * TSMean[i].cumSARI, c * TSMean[i].cumCritical, c * TSMean[i].cumCritRecov, c*TSMean[i].D,
				c * TSMean[i].cumDeath_ILI, c * TSMean[i].cumDeath_SARI, c * TSMean[i].cumDeath_Critical);
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
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
		fclose(dat);
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
			if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
			fprintf(dat, "t");

			/////// ****** /////// ****** /////// ****** COLNAMES
			//// prevalance
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tMild_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tILI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tSARI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tCritical_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tCritRecov_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tSARIP_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tCriticalP_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tCritRecovP_%i", i);

			//// incidence
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincMild_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincILI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincSARI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincCritical_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincCritRecov_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincSARIP_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincCriticalP_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincCritRecovP_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincDeath_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincDeath_ILI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincDeath_SARI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tincDeath__Critical_%i", i);

			//// cumulative incidence
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumMild_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumILI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumSARI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumCritical_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumCritRecov_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumDeaths_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumDeaths_ILI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumDeaths_SARI_%i", i);
			for (i = 0; i < NUM_AGE_GROUPS; i++) fprintf(dat, "\tcumDeaths_Critical_%i", i);

			fprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table.
			for (i = 0; i < P.NumSamples; i++)
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
				fprintf(dat, "%.10f", c * TSMean[i].t);
				//// prevalance
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].Mild_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].ILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].SARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].Critical_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].CritRecov_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * (TSMean[i].SARI_age[j] - SARI_a[j]));
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * (TSMean[i].Critical_age[j] - Critical_a[j]));
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * (TSMean[i].CritRecov_age[j] - CritRecov_a[j]));

				//// incidence
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incIa[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incMild_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incSARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incCritical_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incCritRecov_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * incSARI_a[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * incCritical_a[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * incCritRecov_a[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incDa[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incDeath_ILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incDeath_SARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].incDeath_Critical_age[j]);

				//// cumulative incidence
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumMild_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumSARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumCritical_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumCritRecov_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * (TSMean[i].cumDeath_ILI_age[j] + TSMean[i].cumDeath_SARI_age[j] + TSMean[i].cumDeath_Critical_age[j]));
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_ILI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_SARI_age[j]);
				for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_Critical_age[j]);
				fprintf(dat, "\n");
			}
			fclose(dat);
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
			if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
			fprintf(dat, "t");

			/////// ****** /////// ****** /////// ****** COLNAMES
			//// prevalance
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tMild_%s"					, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tILI_%s"					, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tSARI_%s"					, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritical_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritRecov_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tSARIP_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCriticalP_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritRecovP_%s"			, AdUnits[i].ad_name);

			//// incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincI_%s"					, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincMild_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincILI_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincSARI_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritical_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritRecov_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincSARIP_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCriticalP_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritRecovP_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincDeath_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincDeath_ILI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincDeath_SARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincDeath__Critical_%s"	, AdUnits[i].ad_name);

			//// cumulative incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumMild_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumILI_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumSARI_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumCritical_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumCritRecov_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumDeaths_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumDeaths_ILI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumDeaths_SARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumDeaths_Critical_%s"	, AdUnits[i].ad_name);

			fprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table.
			for (i = 0; i < P.NumSamples; i++)
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
				fprintf(dat, "%.10f", c*TSMean[i].t);
				//// prevalance
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].Mild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].Critical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].CritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * (TSMean[i].SARI_adunit[j] - SARI_a[j]));
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * (TSMean[i].Critical_adunit[j] - Critical_a[j]));
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * (TSMean[i].CritRecov_adunit[j] - CritRecov_a[j]));

				//// incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * incSARI_a[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * incCritical_a[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * incCritRecov_a[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incDeath_Critical_adunit[j]);

				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumD_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumDeath_Critical_adunit[j]);

				fprintf(dat, "\n");
			}
			fclose(dat);
			Memory::xfree(SARI_a); Memory::xfree(Critical_a); Memory::xfree(CritRecov_a);
			Memory::xfree(incSARI_a); Memory::xfree(incCritical_a); Memory::xfree(incCritRecov_a);
		}
	}

	if (P.DoAdUnits && P.OutputAdUnitAge)
	{
		//// output infections by age and admin unit
		outname = output_file_base + ".age.adunit.xls";
		if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");

		// colnames
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				fprintf(dat, "\tincInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// incidence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				fprintf(dat, "\tprevInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// prevalence
		for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				fprintf(dat, "\tcumInf_AG_%i_%s", AgeGroup, AdUnits[AdUnit].ad_name);	// cumulative incidence
		fprintf(dat, "\n");

		// Populate
		for (int Time = 0; Time < P.NumSamples; Time++)
		{
			fprintf(dat, "%.10f", c * TSMean[Time].t);
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					fprintf(dat, "\t%.10f", c * TSMean[Time].incInf_age_adunit[AgeGroup][AdUnit]);	// incidence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					fprintf(dat, "\t%.10f", c * TSMean[Time].prevInf_age_adunit[AgeGroup][AdUnit]);	// prevalence
			for (int AdUnit = 0; AdUnit < P.NumAdunits; AdUnit++)
				for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
					fprintf(dat, "\t%.10f", c * TSMean[Time].cumInf_age_adunit[AgeGroup][AdUnit]);	// cumulative incidence
			fprintf(dat, "\n");
		}
		fclose(dat);
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
	FILE* dat;

	std::string outname = output_file_base + ".seeds.xls";
	if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "%i\t%i\n", P.nextRunSeed1, P.nextRunSeed2);
	fclose(dat);
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
	FILE* dat;

	std::string outname = output_file_base + ".infevents.xls";
	if (!(dat = fopen(outname.c_str(), "wb"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "type,t,thread,ind_infectee,cell_infectee,listpos_infectee,adunit_infectee,x_infectee,y_infectee,t_infector,ind_infector,cell_infector\n");
	for (i = 0; i < nEvents; i++)
	{
		fprintf(dat, "%i\t%.10f\t%i\t%i\t%i\t%i\t%i\t%.10f\t%.10f\t%.10f\t%i\t%i\n",
			InfEventLog[i].type, InfEventLog[i].t, InfEventLog[i].thread, InfEventLog[i].infectee_ind, InfEventLog[i].infectee_cell, InfEventLog[i].listpos, InfEventLog[i].infectee_adunit, InfEventLog[i].infectee_x, InfEventLog[i].infectee_y, InfEventLog[i].t_infector, InfEventLog[i].infector_ind, InfEventLog[i].infector_cell);
	}
	fclose(dat);
}

void LoadSnapshot(std::string const& snapshot_load_file)
{
	FILE* dat;
	int i, j, * CellMemberArray, * CellSuscMemberArray;
	int32_t l;
	long long CM_offset, CSM_offset;
	double t;
	int** Array_InvCDF;
	float* Array_tot_prob, ** Array_cum_trans, ** Array_max_trans;

	if (!(dat = fopen(snapshot_load_file.c_str(), "rb"))) ERR_CRITICAL("Unable to open snapshot file\n");
	fprintf(stderr, "Loading snapshot.");
	Array_InvCDF = (int**)Memory::xcalloc(P.NCP, sizeof(int*));
	Array_max_trans = (float**)Memory::xcalloc(P.NCP, sizeof(float*));
	Array_cum_trans = (float**)Memory::xcalloc(P.NCP, sizeof(float*));
	Array_tot_prob = (float*)Memory::xcalloc(P.NCP, sizeof(float));
	for (i = 0; i < P.NCP; i++)
	{
		Array_InvCDF[i] = Cells[i].InvCDF;
		Array_max_trans[i] = Cells[i].max_trans;
		Array_cum_trans[i] = Cells[i].cum_trans;
		Array_tot_prob[i] = Cells[i].tot_prob;
	}

	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.PopSize) ERR_CRITICAL_FMT("Incorrect N (%i %i) in snapshot file.\n", P.PopSize, i);
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.NH) ERR_CRITICAL("Incorrect NH in snapshot file.\n");
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.NC) ERR_CRITICAL_FMT("## %i neq %i\nIncorrect NC in snapshot file.", i, P.NC);
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.NCP) ERR_CRITICAL("Incorrect NCP in snapshot file.\n");
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.ncw) ERR_CRITICAL("Incorrect ncw in snapshot file.\n");
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.nch) ERR_CRITICAL("Incorrect nch in snapshot file.\n");
	fread_big((void*)& l, sizeof(int32_t), 1, dat); if (l != P.setupSeed1) ERR_CRITICAL("Incorrect setupSeed1 in snapshot file.\n");
	fread_big((void*)& l, sizeof(int32_t), 1, dat); if (l != P.setupSeed2) ERR_CRITICAL("Incorrect setupSeed2 in snapshot file.\n");
	fread_big((void*)& t, sizeof(double), 1, dat); if (t != P.TimeStep) ERR_CRITICAL("Incorrect TimeStep in snapshot file.\n");
	fread_big((void*) & (P.SnapshotLoadTime), sizeof(double), 1, dat);
	P.NumSamples = 1 + (int)ceil((P.SampleTime - P.SnapshotLoadTime) / P.SampleStep);
	fprintf(stderr, ".");
	fread_big((void*)& CellMemberArray, sizeof(int*), 1, dat);
	fprintf(stderr, ".");
	fread_big((void*)& CellSuscMemberArray, sizeof(int*), 1, dat);
	fprintf(stderr, ".");
	CM_offset = State.CellMemberArray - CellMemberArray;
	CSM_offset = State.CellSuscMemberArray - CellSuscMemberArray;

	fread_big((void*)Hosts, sizeof(Person), (size_t)P.PopSize, dat);
	fprintf(stderr, ".");
	fread_big((void*)Households, sizeof(Household), (size_t)P.NH, dat);
	fprintf(stderr, ".");
	fread_big((void*)Cells, sizeof(Cell), (size_t)P.NC, dat);
	fprintf(stderr, ".");
	fread_big((void*)Mcells, sizeof(Microcell), (size_t)P.NMC, dat);
	fprintf(stderr, ".");
	fread_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	fprintf(stderr, ".");
	fread_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	fprintf(stderr, ".");
	for (i = 0; i < P.NC; i++)
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
	for (i = 0; i < P.NMC; i++)
		if (Mcells[i].n > 0)
			Mcells[i].members += CM_offset;

	for (i = 0; i < P.NCP; i++)
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
	fprintf(stderr, "\n");
	fclose(dat);
}

void SaveSnapshot(std::string const& snapshot_save_file)
{
	FILE* dat;
	int i = 1;

	if (!(dat = fopen(snapshot_save_file.c_str(), "wb"))) ERR_CRITICAL("Unable to open snapshot file\n");

	fwrite_big((void*) & (P.PopSize), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.NH), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.NC), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.NCP), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.ncw), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.nch), sizeof(int), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.setupSeed1), sizeof(int32_t), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.setupSeed2), sizeof(int32_t), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.TimeStep), sizeof(double), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.SnapshotSaveTime), sizeof(double), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (State.CellMemberArray), sizeof(int*), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (State.CellSuscMemberArray), sizeof(int*), 1, dat);
	fprintf(stderr, "## %i\n", i++);

	fwrite_big((void*)Hosts, sizeof(Person), (size_t)P.PopSize, dat);

	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)Households, sizeof(Household), (size_t)P.NH, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)Cells, sizeof(Cell), (size_t)P.NC, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)Mcells, sizeof(Microcell), (size_t)P.NMC, dat);
	fprintf(stderr, "## %i\n", i++);

	fwrite_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.PopSize, dat);
	fprintf(stderr, "## %i\n", i++);

	fclose(dat);
}

void UpdateProbs(int DoPlace)
{
	if (!DoPlace)
	{
#pragma omp parallel for schedule(static,500) default(none) \
			shared(P, CellLookup)
		for (int j = 0; j < P.NCP; j++)
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
		for (int j = 0; j < P.NCP; j++)
		{
			CellLookup[j]->S0 = CellLookup[j]->S;
			CellLookup[j]->tot_prob = 0;
		}
	}
#pragma omp parallel for schedule(static,500) default(none) \
		shared(P, CellLookup)
	for (int j = 0; j < P.NCP; j++)
	{
		int m, k;
		float t;
		CellLookup[j]->cum_trans[0] = ((float)(CellLookup[0]->S0)) * CellLookup[j]->max_trans[0];
		t = ((float)CellLookup[0]->n) * CellLookup[j]->max_trans[0];
		for (m = 1; m < P.NCP; m++)
		{
				CellLookup[j]->cum_trans[m] = CellLookup[j]->cum_trans[m - 1] + ((float)(CellLookup[m]->S0)) * CellLookup[j]->max_trans[m];
				t += ((float)CellLookup[m]->n) * CellLookup[j]->max_trans[m];
		}
		CellLookup[j]->tot_prob = CellLookup[j]->cum_trans[P.NCP - 1];
		for (m = 0; m < P.NCP; m++)
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
	else VariableAndValue = INT_MAX; //// i.e. if not doing triggering (at either admin or global level) then set value to be arbitrarily large so that it will surpass any trigger threshold. Probably other ways around this if anybody wants to correct?

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


void UpdateEfficaciesAndComplianceProportions(double t)
{
	//// **** social distancing
	for (int ChangeTime = 0; ChangeTime < P.Num_SD_ChangeTimes; ChangeTime++)
		if (t == P.SD_ChangeTimes[ChangeTime])
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
		if (t == P.CI_ChangeTimes[ChangeTime])
		{
			P.CaseIsolationEffectiveness		= P.CI_SpatialAndPlaceEffects_OverTime	[ChangeTime]; //// spatial / place
			P.CaseIsolationHouseEffectiveness	= P.CI_HouseholdEffects_OverTime		[ChangeTime]; //// household

			P.CaseIsolationProp					= P.CI_Prop_OverTime					[ChangeTime]; //// compliance
			P.CaseIsolation_CellIncThresh		= P.CI_CellIncThresh_OverTime			[ChangeTime]; //// cell incidence threshold
		}

	////// **** household quarantine
	if (P.DoHouseholds)
		for (int ChangeTime = 0; ChangeTime < P.Num_HQ_ChangeTimes; ChangeTime++)
			if (t == P.HQ_ChangeTimes[ChangeTime])
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
			if (t == P.PC_ChangeTimes[ChangeTime])
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
				if(P.PlaceCloseTimeStart<1e10) P.PlaceCloseTimeStart = t;

				// ensure that new duration doesn't go over next change time. Judgement call here - talk to Neil if this is what he wants.
				if ((ChangeTime < P.Num_PC_ChangeTimes - 1) && (P.PlaceCloseTimeStart + P.PlaceCloseDuration >= P.PC_ChangeTimes[ChangeTime + 1]))
					P.PlaceCloseDuration = P.PC_ChangeTimes[ChangeTime + 1] - P.PC_ChangeTimes[ChangeTime]; // -1;
				//fprintf(stderr, "\nt=%lf, n=%i (%i)  PlaceCloseDuration = %lf  (%lf) \n", t, ChangeTime, P.Num_PC_ChangeTimes, P.PlaceCloseDuration, P.PC_ChangeTimes[ChangeTime+1]);
			}
	}

	//// **** digital contact tracing
	for (int ChangeTime = 0; ChangeTime < P.Num_DCT_ChangeTimes; ChangeTime++)
		if (t == P.DCT_ChangeTimes[ChangeTime])
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
				if (Hosts[Person].inf == InfStat_Susceptible || abs(Hosts[Person].inf) == InfStat_Recovered) QuarNotInfected++;
				if (Hosts[Person].inf > 0) QuarNotSymptomatic++;
			}

	TimeSeries[n].prevQuarNotInfected		= (double) QuarNotInfected;
	TimeSeries[n].prevQuarNotSymptomatic	= (double) QuarNotSymptomatic;
}

void RecordSample(double t, int n, std::string const& output_file_base)
{
	int j, k, S, L, I, R, D, N, cumC, cumTC, cumI, cumR, cumD, cumDC, cumFC, cumTG, cumSI, nTG;
	int cumCT; //added cumulative number of contact traced: ggilani 15/06/17
	int cumCC; //added cumulative number of cases who are contacts: ggilani 28/05/2019
	int cumDCT; //added cumulative number of cases who are digitally contact traced: ggilani 11/03/20
	int cumHQ, cumAC, cumAH, cumAA, cumACS, cumAPC, cumAPA, cumAPCS, numPC, trigDetectedCases;
	int cumC_country[MAX_COUNTRIES]; //add cumulative cases per country
	unsigned short int ts;

	//// Severity quantities
	int Mild, ILI, SARI, Critical, CritRecov, cumMild, cumILI, cumSARI, cumCritical, cumCritRecov, cumDeath_ILI, cumDeath_SARI, cumDeath_Critical;

	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	//// initialize to zero
	S = L = I = R = D = cumI = cumC = cumDC = cumTC = cumFC = cumHQ = cumAC = cumAA = cumAH = cumACS
		= cumAPC = cumAPA = cumAPCS = cumD = cumCT = cumCC = cumDCT = cumTG = cumSI = nTG = 0;
	for (int i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] = 0;
	if (P.DoSeverity)
		Mild = ILI = SARI = Critical = CritRecov = cumMild = cumILI = cumSARI = cumCritical = cumCritRecov = cumDeath_ILI = cumDeath_SARI = cumDeath_Critical = 0;

#pragma omp parallel for schedule(static,10000) reduction(+:S,L,I,R,D,cumTC) default(none) \
		shared(P, CellLookup)
	for (int i = 0; i < P.NCP; i++)
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
	if (N != P.PopSize) fprintf(stderr, "## %i #\n", P.PopSize - N);
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
	TimeSeries[n].incI = (double)(cumI - State.cumI);
	TimeSeries[n].incC = (double)(cumC - State.cumC);
	TimeSeries[n].incFC = (double)(cumFC - State.cumFC);
	TimeSeries[n].incCT = (double)(cumCT - State.cumCT); // added contact tracing
	TimeSeries[n].incCC = (double)(cumCC - State.cumCC); // added cases who are contacts
	TimeSeries[n].incDCT = (double)(cumDCT - State.cumDCT); //added cases who are digitally contact traced
	TimeSeries[n].incDC = (double)(cumDC - State.cumDC); //added incidence of detected cases
	TimeSeries[n].incTC = (double)(cumTC - State.cumTC);
	TimeSeries[n].incR = (double)(cumR - State.cumR);
	TimeSeries[n].incD = (double)(cumD - State.cumD);
	TimeSeries[n].incHQ = (double)(cumHQ - State.cumHQ);
	TimeSeries[n].incAC = (double)(cumAC - State.cumAC);
	TimeSeries[n].incAH = (double)(cumAH - State.cumAH);
	TimeSeries[n].incAA = (double)(cumAA - State.cumAA);
	TimeSeries[n].incACS = (double)(cumACS - State.cumACS);
	TimeSeries[n].incAPC = (double)(cumAPC - State.cumAPC);
	TimeSeries[n].incAPA = (double)(cumAPA - State.cumAPA);
	TimeSeries[n].incAPCS = (double)(cumAPCS - State.cumAPCS);
	TimeSeries[n].cumT = State.cumT;
	TimeSeries[n].cumUT = State.cumUT;
	TimeSeries[n].cumTP = State.cumTP;
	TimeSeries[n].cumV = State.cumV;
	TimeSeries[n].cumVG = State.cumVG; //added VG;
	TimeSeries[n].cumDC = cumDC;
	TimeSeries[n].meanTG = TimeSeries[n].meanSI = 0;
	if (nTG > 0) // record mean generation time and serial interval in timestep
	{
		TimeSeries[n].meanTG = P.TimeStep * ((double)cumTG) / ((double)nTG);
		TimeSeries[n].meanSI = P.TimeStep * ((double)cumSI) / ((double)nTG);
	}
	//fprintf(stderr, "\ncumD=%i last_cumD=%i incD=%lg\n ", cumD, State.cumD, TimeSeries[n].incD);
	//incidence per country
	for (int i = 0; i < MAX_COUNTRIES; i++) TimeSeries[n].incC_country[i] = (double)(cumC_country[i] - State.cumC_country[i]);
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
		TimeSeries[n].incMild = (double)(cumMild - State.cumMild);
		TimeSeries[n].incILI = (double)(cumILI - State.cumILI);
		TimeSeries[n].incSARI = (double)(cumSARI - State.cumSARI);
		TimeSeries[n].incCritical = (double)(cumCritical - State.cumCritical);
		TimeSeries[n].incCritRecov = (double)(cumCritRecov - State.cumCritRecov);
		TimeSeries[n].incDeath_ILI = (double)(cumDeath_ILI - State.cumDeath_ILI);
		TimeSeries[n].incDeath_SARI = (double)(cumDeath_SARI - State.cumDeath_SARI);
		TimeSeries[n].incDeath_Critical = (double)(cumDeath_Critical - State.cumDeath_Critical);

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
	for (int i = k = 0; i < P.NMC; i++) if (Mcells[i].socdist == 2) k++;
	TimeSeries[n].PropSocDist = ((double)k) / ((double)P.NMC);

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
	unsigned short int ts;

	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	trigAlertCases = State.cumDC;
	if (n >= P.WindowToEvaluateTriggerAlert) {
		trigAlertCases -= (int)TimeSeries[n - P.WindowToEvaluateTriggerAlert].cumDC;
	}

	if (P.TriggerAlertOnDeaths)  //// if using deaths as trigger (as opposed to detected cases)
	{
		trigAlert = (int)TimeSeries[n].D;
		if (n >= P.WindowToEvaluateTriggerAlert) trigAlert -= (int) TimeSeries[n - P.WindowToEvaluateTriggerAlert].D;
	}
	else
	{
		trigAlert = trigAlertCases;
	}

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
//					fprintf(stderr, "@@## trigAlertCases=%i P.HolidaysStartDay_SimTime=%lg \n",trigAlertCases, P.HolidaysStartDay_SimTime);
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
					fprintf(stderr, "\n** %i %lf %lf | %lg / %lg \t", P.ModelCalibIteration, t, P.DateTriggerReached_SimTime + P.DateTriggerReached_CalTime - P.Interventions_StartDate_CalTime, P.HolidaysStartDay_SimTime, RatioPredictedObserved);
					fprintf(stderr, "| %i %i %i %i -> ", trigAlert, trigAlertCases, P.AlertTriggerAfterIntervThreshold, P.CaseOrDeathThresholdBeforeAlert);

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
							P.SeedingScaling /= pow(RatioPredictedObserved, 0.3 + 0.2 * ranf()); // include random number to prevent loops
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
					fprintf(stderr, "%i : %lg\n", P.CaseOrDeathThresholdBeforeAlert, P.SeedingScaling);

					if(P.StopCalibration)
						fprintf(stderr, "Calibration ended.\n");
					else
						InterruptRun = 1;
				}
				else
					P.StopCalibration = 1;
			}
		}
		P.ControlPropCasesId = P.PostAlertControlPropCasesId;

		if (P.VaryEfficaciesOverTime)
			UpdateEfficaciesAndComplianceProportions(t - P.Epidemic_StartDate_CalTime);

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
			fprintf(stderr, "\nSecond place closure period (t=%lg)\n", t);
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

	static int DataAlreadyRead = 0, ncols,nrows, *ColTypes;
	static double **Data,NegBinK,sumL;

	if (!DataAlreadyRead)
	{
		char FieldName[1024];
		if (!(dat = fopen(DataFile.c_str(), "r"))) ERR_CRITICAL("Unable to open data file\n");
		fscanf(dat, "%i %i %lg", &nrows, &ncols, &NegBinK);
		if (!(ColTypes = (int*)calloc(ncols, sizeof(int)))) ERR_CRITICAL("Unable to allocate data file storage\n");
		if (!(Data = (double**)calloc(nrows, sizeof(double *)))) ERR_CRITICAL("Unable to allocate data file storage\n");
		for(int i=0;i<nrows;i++)
			if (!(Data[i] = (double*)calloc(ncols, sizeof(double)))) ERR_CRITICAL("Unable to allocate data file storage\n");
		for (int i = 0; i < ncols; i++)
		{
			ColTypes[i] = -100;
			fscanf(dat, "%s", FieldName);
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
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				fscanf(dat, "%lg", &(Data[i][j]));
		fclose(dat);
		DataAlreadyRead = 1;
	}

	// calculate likelihood function
	double c, LL=0.0;
	double kp = (P.clP[99] > 0) ? P.clP[99] : NegBinK; //clP[99] reserved for fitting overdispersion
	c = 1.0; // 1 / ((double)(P.NRactE + P.NRactNE));
	int offset= (P.Interventions_StartDate_CalTime > 0) ? ((int)(P.Interventions_StartDate_CalTime - P.DateTriggerReached_SimTime)) : 0;
	for (int i = 1; i < ncols;i++)
	{
		if ((ColTypes[i] >= 0)&&(ColTypes[i] <= 2))
		{
			double ModelValueSum = 0.0;
			for (int j = 0; j < nrows; j++)
			{
				int day = (int)Data[j][0]; // day is day of year - directly indexes TimeSeries[]
				if ((Data[j][i]>=-1) && (day < P.NumSamples)) // data is not NA (-ve) and within time range of model run
				{
					double ModelValue;
					if (ColTypes[i]==0)
						ModelValue=c*TimeSeries[day-offset].incD; // all deaths by date of death
					else if (ColTypes[i] == 1)
						ModelValue=c*(TimeSeries[day-offset].incDeath_Critical+ TimeSeries[day-offset].incDeath_SARI); // hospital deaths (SARI and Critical) by date of death
					else if (ColTypes[i] == 2)
						ModelValue = c * TimeSeries[day-offset].incDeath_ILI; // care home deaths (ILI) by date of death
					ModelValueSum += ModelValue;
					if (Data[j][i] >= 0)
					{
						if ((j > 0) && (Data[j - 1][i] == -1)) // cumulative column: -1 means sum column up to first >=0 value
						{
							ModelValue = ModelValueSum;
							ModelValueSum = 0.0;  // reset cumulative sum
						}
						if (NegBinK >= 10000)
							//prob model and data from same underlying poisson
							LL += lgamma(2 * (Data[j][i] + ModelValue) + 1) - lgamma(Data[j][i] + ModelValue + 1) - lgamma(Data[j][i] + 1) - lgamma(ModelValue + 1) - (3 * (Data[j][i] + ModelValue) + 1) * log(2);
						else
						{
							//neg bin LL (NegBinK=1 implies no over-dispersion. >1 implies more)
							double knb = 1.0 + ModelValue / kp;
							double pnb = kp / (1.0 + kp);
							LL += lgamma(Data[j][i] + knb) - lgamma(Data[j][i] + 1) - lgamma(knb) + knb * log(1.0 - pnb) + Data[j][i] * log(pnb);
						}
					}
				}
			}
		}
		else if (ColTypes[i] == 3) //seroprevalence by date of sample
		{
			for (int j = 0; j < nrows; j++)
			{
				int day = (int)Data[j][0]; // day is day of year - directly indexes TimeSeries[]
				if ((Data[j][i] >= 0) && (day < P.NumSamples)) // data is not NA (-ve) and within time range of model run
				{
					double m = Data[j][i]; // numerator
					double N = Data[j][i + 1]; // denominator
					double ModelValue;
					for (int k = offset; k < day; k++) // loop over all days of infection up to day of sample
					{
						double prob_seroconvert = P.SeroConvMaxSens*(1.0-0.5*((exp(-((double)(day - k))*P.SeroConvP1) + 1.0)*exp(-((double)(day - k))*P.SeroConvP2))); // add P1 to P2 to prevent degeneracy
						ModelValue += c * TimeSeries[k - offset].incI * prob_seroconvert;
					}
					ModelValue += c * TimeSeries[day-offset].S * (1.0 - P.SeroConvSpec);
					ModelValue /= ((double)P.PopSize);
					LL += m * log((ModelValue + 1e-20)/(m/N+1e-20)) + (N - m) * log((1.0 - ModelValue + 1e-20)/(1.0-m/N+1e-20));  // subtract saturated likelihood
				}
			}
		}
		else if (ColTypes[i] == 5) // infection prevalence by date of sample
		{
			for (int j = 0; j < nrows; j++)
			{
				int day = (int)Data[j][0]; // day is day of year - directly indexes TimeSeries[]
				if ((Data[j][i] >= 0) && (day < P.NumSamples)) // data is not NA (-ve) and within time range of model run
				{
					double m = Data[j][i]; // numerator
					double N = Data[j][i + 1]; // denominator
					double ModelValue = P.InfPrevSurveyScale * c * TimeSeries[day-offset].I / ((double)P.PopSize);
					LL += m * log(ModelValue + 1e-20) + (N - m) * log(1.0 - ModelValue);
				}
			}
		}
	}
	fprintf(stderr, "Log-likelihood = %lg\n", LL);
	if (run == 0)
		sumL = LL;
	else
	{
		double maxLL = LL;
		if (sumL > maxLL) maxLL = sumL;
		sumL = maxLL + log(exp(sumL - maxLL) + exp(LL - maxLL));
	}

	if (run + 1 == P.NumRealisations)
	{
		LL = sumL - log((double)P.NumRealisations);
		std::string TmpFile = OutFileBase + ".ll.tmp";
		std::string OutFile = OutFileBase + ".ll.txt";
		if (!(dat = fopen(TmpFile.c_str(), "w"))) ERR_CRITICAL("Unable to open likelihood file\n");
		fprintf(dat, "%i\t%.8lg\n", P.FitIter, LL);
		fclose(dat);
		rename(TmpFile.c_str(), OutFile.c_str()); // rename only when file is complete and closed
	}
}

void RecordInfTypes(void)
{
	int i, j, k, l, lc, lc2, b, c, n, nf, i2;
	double* res, * res_av, * res_var, t, s;

	for (n = 0; n < P.NumSamples; n++)
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
	for (b = 0; b < P.NC; b++)
		if ((Cells[b].S != Cells[b].n) || (Cells[b].R > 0))
			for (c = 0; c < Cells[b].n; c++)
				Hosts[Cells[b].members[c]].listpos = 0;
	//	for(b=0;b<P.NC;b++)
	//		if((Cells[b].S!=Cells[b].n)||(Cells[b].R>0))
	{
		j = k = l = lc = lc2 = 0; t = 1e10;
		//			for(c=0;c<Cells[b].n;c++)
		for (i = 0; i < P.PopSize; i++)
		{
			//				i=Cells[b].members[c];
			if (j == 0) j = k = Households[Hosts[i].hh].nh;
			if ((Hosts[i].inf != InfStat_Susceptible) && (Hosts[i].inf != InfStat_ImmuneAtStart))
			{
				if (Hosts[i].latent_time * P.TimeStep <= P.SampleTime)
					TimeSeries[(int)(Hosts[i].latent_time * P.TimeStep / P.SampleStep)].Rdenom++;
				infcountry[mcell_country[Hosts[i].mcell]]++;
				if (abs(Hosts[i].inf) < InfStat_Recovered)
					l = -1;
				else if (l >= 0)
					l++;
				if ((l >= 0) && ((Hosts[i].inf == InfStat_RecoveredFromSymp) || (Hosts[i].inf == InfStat_Dead_WasSymp)))
				{
					lc2++;
					if (Hosts[i].latent_time * P.TimeStep <= t) // This convoluted logic is to pick up households where the index is symptomatic
					{
						lc = 1; t = Hosts[i].latent_time * P.TimeStep;
					}
				}
				else if ((l > 0) && (Hosts[i].latent_time * P.TimeStep < t))
				{
					lc = 0; t = Hosts[i].latent_time * P.TimeStep;
				}
				i2 = Hosts[i].infector;
				if (i2 >= 0)
				{
					Hosts[i2].listpos++;
					if (Hosts[i2].latent_time * P.TimeStep <= P.SampleTime)
					{
						TimeSeries[(int)(Hosts[i2].latent_time * P.TimeStep / P.SampleStep)].Rtype[Hosts[i].infect_type % INFECT_TYPE_MASK]++;
						TimeSeries[(int)(Hosts[i2].latent_time * P.TimeStep / P.SampleStep)].Rage[HOST_AGE_GROUP(i)]++;
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
	for (b = 0; b < P.NC; b++)
		if ((Cells[b].S != Cells[b].n) || (Cells[b].R > 0))
			for (c = 0; c < Cells[b].n; c++)
			{
				i = Cells[b].members[c];
				if ((abs(Hosts[i].inf) == InfStat_Recovered) || (abs(Hosts[i].inf) == InfStat_Dead))
				{
					l = Hosts[i].infect_type / INFECT_TYPE_MASK;
					if ((l < MAX_GEN_REC) && (Hosts[i].listpos < MAX_SEC_REC)) indivR0[Hosts[i].listpos][l]++;
				}
			}
	/* 	if(!TimeSeries[P.NumSamples-1].extinct) */
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
	for (n = 0; n < P.NumSamples; n++)
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

	if (TimeSeries[P.NumSamples - 1].extinct)
	{
		TSMean = TSMeanE; TSVar = TSVarE; P.NRactE++;
	}
	else
	{
		TSMean = TSMeanNE; TSVar = TSVarNE; P.NRactNE++;
	}
	lc = -k;

	// This calculates sum and sum of squares of entire TimeSeries array
	for (n = 0; n < P.NumSamples; n++)
	{
		if ((n + lc >= 0) && (n + lc < P.NumSamples))
		{
			if (s < TimeSeries[n + lc].incC) { s = TimeSeries[n + lc].incC; t = P.SampleStep * ((double)(n + lc)); }
			res = (double*)&TimeSeries[n + lc];
			res_av = (double*)&TSMean[n];
			res_var = (double*)&TSVar[n];
			for (std::size_t i = ResultsDoubleOffsetStart /* skip over initial fields */; i < nf; i++)
			{
				res_av[i] += res[i];
				res_var[i] += res[i] * res[i];
			}
			if (P.DoAdUnits && P.OutputAdUnitAge)
				for (std::size_t age = 0; age < NUM_AGE_GROUPS; ++age)
					for (std::size_t adunit = 0; adunit < P.NumAdunits; ++adunit)
					{
						TSMean[n].prevInf_age_adunit[age][adunit] = TimeSeries[n + lc].prevInf_age_adunit[age][adunit];
						TSMean[n].cumInf_age_adunit [age][adunit] = TimeSeries[n + lc].cumInf_age_adunit [age][adunit];
						TSMean[n].incInf_age_adunit [age][adunit] = TimeSeries[n + lc].incInf_age_adunit [age][adunit];
					}

			if (TSMean[n].cumTmax < TimeSeries[n + lc].cumT) TSMean[n].cumTmax = TimeSeries[n + lc].cumT;
			if (TSMean[n].cumVmax < TimeSeries[n + lc].cumV) TSMean[n].cumVmax = TimeSeries[n + lc].cumV;
		}
		TSMean[n].t += ((double) n )* P.SampleStep;
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
		for (int i = tn; i < P.NCP; i += P.NumThreads)
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
					ptrdiff_t mcl_from = cl_from_mcl + l + k * P.total_microcells_high_;
					if (Mcells[mcl_from].n > 0)
					{
						//get proportion of each population of cell that exists in each admin unit
						pop_dens_from[Mcells[mcl_from].adunit] += (((double)Mcells[mcl_from].n) / ((double)Cells[cl_from].n));
					}
				}
			}

			for (int j = i; j < P.NCP; j++)
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
					total_flow = Cells[cl_from].cum_trans[j] * Cells[cl_from].n;
				}
				else
				{
					total_flow = (Cells[cl_from].cum_trans[j] - Cells[cl_from].cum_trans[j - 1]) * Cells[cl_from].n;
				}

				//loop over microcells within destination cell
				for (int m = 0; m < P.NMCL; m++)
				{
					for (int p = 0; p < P.NMCL; p++)
					{
						//get index of microcell
						ptrdiff_t mcl_to = cl_to_mcl + p + m * P.total_microcells_high_;
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

//// Get parameters code (called by ReadParams function)
int GetInputParameter(FILE* dat, FILE* dat2, const char* SItemName, const char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	int FindFlag;

	FindFlag = GetInputParameter2(dat, dat2, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	if (!FindFlag)
	{
		ERR_CRITICAL_FMT("\nUnable to find parameter `%s' in input file. Aborting program...\n", SItemName);
	}
	return FindFlag;
}
int GetInputParameter2(FILE* dat, FILE* dat2, const char* SItemName, const char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	int FindFlag = 0;

	if (dat2) FindFlag = GetInputParameter3(dat2, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	if (!FindFlag)
		FindFlag = GetInputParameter3(dat, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	return FindFlag;
}

int GetInputParameter2all(FILE* dat, FILE* dat2, FILE* dat3, const char* SItemName, const char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	int FindFlag = 0;

	if (dat3) FindFlag = GetInputParameter3(dat3, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	if (!FindFlag)
	{
		if (dat2) FindFlag = GetInputParameter3(dat2, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
		if ((!FindFlag)&&(dat))
			FindFlag = GetInputParameter3(dat, SItemName, ItemType, ItemPtr, NumItem, NumItem2, Offset);
	}
	return FindFlag;
}

/*
    Reads a string (as per fscanf %s).
    Returns true if it succeeds, false on EOF, and does not return on error.
*/
bool readString(const char* SItemName, FILE* dat, char *buf) {
    int r = fscanf(dat, "%s", buf);
    if(r == 1) {
        return true;
    } else if (r == EOF) {
        if(ferror(dat)) {
            ERR_CRITICAL_FMT("fscanf failed for %s: %s.\n", SItemName, strerror(errno));
        } else {
            // EOF
            return false;
        }
    } else {
        ERR_CRITICAL_FMT("Unexpected fscanf result %d for %s.\n", r, SItemName);
    }
}

// #define to map *nix strncasecmp to windows equivalent
#ifdef WIN32
#define strncasecmp _strnicmp
#endif

int GetInputParameter3(FILE* dat, const char* SItemName, const char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
{
	char match[10000] = "", ReadItemName[10000] = "", ItemName[10000];
	int FindFlag = 0, EndString, CurPos, i, j, n;

	n = 0;
	fseek(dat, 0, 0);
	sprintf(ItemName, "[%s]", SItemName);
	while (!FindFlag)
	{
		if(!readString(SItemName, dat, match)) return 0;
		FindFlag = (!strncasecmp(match, ItemName, strlen(match)));
		if (FindFlag)
		{
			CurPos = ftell(dat);
			strcpy(ReadItemName, match);
			EndString = (match[strlen(match) - 1] == ']');
			while ((!EndString) && (FindFlag))
			{
				if(!readString(SItemName, dat, match)) return 0;
				strcat(ReadItemName, " ");
				strcat(ReadItemName, match);
				FindFlag = (!strncmp(ReadItemName, ItemName, strlen(ReadItemName)));
				EndString = (ReadItemName[strlen(ReadItemName) - 1] == ']');
			}
			if (!EndString)
			{
				fseek(dat, CurPos, 0);
				FindFlag = 0;
			}
		}
	}
	if (FindFlag)
	{
		FindFlag = 0;
		if (!strcmp(ItemType, "%lf"))	n = 1;
		else if (!strcmp(ItemType, "%i"))	n = 2;
		else if (!strcmp(ItemType, "%s"))	n = 3;
		if (NumItem2 < 2)
		{
			if (NumItem == 1)
			{
				if(fscanf(dat, "%s", match) != 1) { ERR_CRITICAL_FMT("fscanf failed for %s\n", SItemName); }
				if ((match[0] == '#') && (match[1] >= '0') && (match[1] <= '9'))
				{
					int cln;
					if ((match[2] >= '0') && (match[2] <= '9'))
						cln = ((int)(match[1] - '0')) * 10 + ((int)(match[2] - '0'));
					else
						cln = (int)(match[1] - '0');
					if ((n<3)&&(P.clP_copies[cln] < MAX_CLP_COPIES))
					{
						P.clP_type[cln][P.clP_copies[cln]] = n;
						P.clP_ptr[cln][P.clP_copies[cln]] = ItemPtr;
						P.clP_copies[cln]++;
					}
					FindFlag++;
					if (n == 1)
						*((double*)ItemPtr) = P.clP[cln];
					else if (n == 2)
						*((int*)ItemPtr) = (int)P.clP[cln];
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] != '[') && (!feof(dat)))
				{
					FindFlag++;
					if (n == 1)
						sscanf(match, "%lf", (double*)ItemPtr);
					else if (n == 2)
						sscanf(match, "%i", (int*)ItemPtr);
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
			}
			else
			{
				for (CurPos = 0; CurPos < NumItem; CurPos++)
				{
					if(fscanf(dat, "%s", match) != 1) { ERR_CRITICAL_FMT("fscanf failed for %s\n", SItemName); }
					if ((match[0] == '#') && (match[1] >= '0') && (match[1] <= '9'))
					{
						int cln;
						if ((match[2] >= '0') && (match[2] <= '9'))
							cln = ((int)(match[1] - '0')) * 10 + ((int)(match[2] - '0'));
						else
							cln = (int)(match[1] - '0');
						FindFlag++;
						if (n == 1)
						{
							double* ipt = (((double*)ItemPtr) + CurPos + Offset);
							*ipt = P.clP[cln];
							if (P.clP_copies[cln] < MAX_CLP_COPIES)
							{
								P.clP_type[cln][P.clP_copies[cln]] = 1;
								P.clP_ptr[cln][P.clP_copies[cln]] = (void *)ipt;
								P.clP_copies[cln]++;
							}
						}
						else if (n == 2)
						{
							int* ipt = (((int*)ItemPtr) + CurPos + Offset);
							*ipt = (int)P.clP[cln];
							if (P.clP_copies[cln] < MAX_CLP_COPIES)
							{
								P.clP_type[cln][P.clP_copies[cln]] = 2;
								P.clP_ptr[cln][P.clP_copies[cln]] = (void*)ipt;
								P.clP_copies[cln]++;
							}
						}
						else if (n == 3)
							sscanf(match, "%s", *(((char**)ItemPtr) + CurPos + Offset));
					}
					else if ((match[0] != '[') && (!feof(dat)))
					{
						FindFlag++;
						if (n == 1)
							sscanf(match, "%lf", ((double*)ItemPtr) + CurPos + Offset);
						else if (n == 2)
							sscanf(match, "%i", ((int*)ItemPtr) + CurPos + Offset);
						else if (n == 3)
							sscanf(match, "%s", *(((char**)ItemPtr) + CurPos + Offset));
					}
					else
						CurPos = NumItem;
				}
			}
		}
		else
		{
			for (j = 0; j < NumItem; j++)
			{ //added these braces
				for (i = 0; i < NumItem2; i++)
				{
					if(fscanf(dat, "%s", match) != 1) { ERR_CRITICAL_FMT("fscanf failed for %s\n", SItemName); }
					if ((match[0] == '#') && (match[1] >= '0') && (match[1] <= '9'))
					{
						FindFlag++;
						int cln;
						if ((match[2] >= '0') && (match[2] <= '9'))
							cln = ((int)(match[1] - '0')) * 10 + ((int)(match[2] - '0'));
						else
							cln = (int)(match[1] - '0');
						if (n == 1)
						{
							double* ipt = (((double**)ItemPtr)[j + Offset] + i + Offset);
							*ipt = P.clP[cln];
							if (P.clP_copies[cln] < MAX_CLP_COPIES)
							{
								P.clP_type[cln][P.clP_copies[cln]] = 1;
								P.clP_ptr[cln][P.clP_copies[cln]] = (void*)ipt;
								P.clP_copies[cln]++;
							}
						}
						else if (n == 2)
						{
							int* ipt = (((int**)ItemPtr)[j + Offset] + i + Offset);
							*ipt = (int)P.clP[cln];
							if (P.clP_copies[cln] < MAX_CLP_COPIES)
							{
								P.clP_type[cln][P.clP_copies[cln]] = 2;
								P.clP_ptr[cln][P.clP_copies[cln]] = (void*)ipt;
								P.clP_copies[cln]++;
							}
						}
						else if (n == 3)
							sscanf(match, "%s", *(((char**)ItemPtr) + CurPos + Offset));
					}
					else if((match[0] != '[') && (!feof(dat)))
					{
						FindFlag++;
						if (n == 1)
							sscanf(match, "%lf", ((double**)ItemPtr)[j + Offset] + i + Offset);
						else
							sscanf(match, "%i", ((int**)ItemPtr)[j + Offset] + i + Offset);
					}
					else
					{
						i = NumItem2;
						j = NumItem;
					}
				}

			}
		}
	}
	//	fprintf(stderr,"%s\n",SItemName);
	return FindFlag;
}


void GetInverseCdf(FILE* param_file_dat, FILE* preparam_file_dat, const char* icdf_name, InverseCdf* inverseCdf,
	double start_value)
{
	if (!GetInputParameter2(param_file_dat, preparam_file_dat, icdf_name, "%lf", (void*)inverseCdf->get_values(), CDF_RES + 1, 1, 0))
	{
		inverseCdf->set_neg_log(start_value);
	}
	inverseCdf->assign_exponent();
}
