/*
(c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission.
*/

#include <errno.h>
#include <stddef.h>

#include "SpatialSim.h"
#include "binio.h"
#include "Rand.h"
#include "Error.h"
#include "Dist.h"
#include "Kernels.h"
#include "Bitmap.h"
#include "Model.h"
#include "Param.h"
#include "SetupModel.h"
#include "SharedFuncs.h"
#include "ModelMacros.h"
#include "InfStat.h"
#include "CalcInfSusc.h"
#include "Update.h"
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

void ReadParams(char*, char*);
void ReadInterventions(char*);
int GetXMLNode(FILE*, const char*, const char*, char*, int);
void ReadAirTravel(char*);
void InitModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void SeedInfection(double, int*, int, int); //adding run number as a parameter for event log: ggilani - 15/10/2014
int RunModel(int); //adding run number as a parameter for event log: ggilani - 15/10/2014
void TravelReturnSweep(double);
void TravelDepartSweep(double);
void InfectSweep(double, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14
void IncubRecoverySweep(double, int); //added int as argument to record run number: ggilani - 15/10/14
int TreatSweep(double);
//void HospitalSweep(double); //added hospital sweep function: ggilani - 10/11/14
void DigitalContactTracingSweep(double); // added function to update contact tracing number
void SaveDistribs(void);
void SaveOriginDestMatrix(void); //added function to save origin destination matrix so it can be done separately to the main results: ggilani - 13/02/15
void SaveResults(void);
void SaveSummaryResults(void);
void SaveRandomSeeds(void); //added this function to save random seeds for each run: ggilani - 09/03/17
void SaveEvents(void); //added this function to save infection events from all realisations: ggilani - 15/10/14
void LoadSnapshot(void);
void SaveSnapshot(void);
void RecordInfTypes(void);
void RecordSample(double, int);

double ContactTracingEffect(int, unsigned short int); //added function to account for effect of contact tracing: ggilani 13/11/14
void CalcOriginDestMatrix_adunit(void); //added function to calculate origin destination matrix: ggilani 28/01/15
int BedsAvailablePerAdUnit(double, int); //added function to check current number of beds available in an admin unit

int GetInputParameter(FILE*, FILE*, const char*, const char*, void*, int, int, int);
int GetInputParameter2(FILE*, FILE*, const char*, const char*, void*, int, int, int);
int GetInputParameter3(FILE*, const char*, const char*, void*, int, int, int);


///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// 
///// ***** ///// ***** ///// ***** ///// ***** ///// ***** GLOBAL VARIABLES (some structures in SpatialSim.h file and some containers) - memory allocated later. 
///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// ***** ///// 

param P;
person* Hosts;
household* Households;
popvar State, StateT[MAX_NUM_THREADS];
cell* Cells, ** CellLookup;
microcell* Mcells, ** McellLookup;
place** Places;
adminunit AdUnits[MAX_ADUNITS];
//// Time Series defs: 
//// TimeSeries is an array of type results, used to store (unsurprisingly) a time series of every quantity in results. Mostly used in RecordSample.
//// TSMeanNE and TSVarNE are the mean and variance of non-extinct time series. TSMeanE and TSVarE are the mean and variance of extinct time series. TSMean and TSVar are pointers that point to either extinct or non-extinct. 
results* TimeSeries, * TSMean, * TSVar, * TSMeanNE, * TSVarNE, * TSMeanE, * TSVarE; //// TimeSeries used in RecordSample, RecordInfTypes, SaveResults. TSMean and TSVar
airport* Airports;
bitmap_header* bmh;
//added declaration of pointer to events log: ggilani - 10/10/2014
events* InfEventLog;
int* nEvents;
int* RevCellLookup;

//added counter for number of types ranf and ranf_mt are called
int count_ranf, count_ranfmt;
double ** PopDensity, * mcell_dens;
int* mcell_adunits, * mcell_num, * mcell_country;
double inftype[INFECT_TYPE_MASK], inftype_av[INFECT_TYPE_MASK], infcountry[MAX_COUNTRIES], infcountry_av[MAX_COUNTRIES], infcountry_num[MAX_COUNTRIES];
double indivR0[MAX_SEC_REC][MAX_GEN_REC], indivR0_av[MAX_SEC_REC][MAX_GEN_REC];
double inf_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], denom_household[MAX_HOUSEHOLD_SIZE + 1];
double inf_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], AgeDist[NUM_AGE_GROUPS], AgeDist2[NUM_AGE_GROUPS];
double case_household[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1], case_household_av[MAX_HOUSEHOLD_SIZE + 1][MAX_HOUSEHOLD_SIZE + 1];
double PropPlaces[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES];
double PropPlacesC[NUM_AGE_GROUPS * AGE_GROUP_WIDTH][NUM_PLACE_TYPES], AirTravelDist[MAX_DIST];
double PeakHeightSum, PeakHeightSS, PeakTimeSum, PeakTimeSS;
float* bmi, * bmi2, * bmi3, * bmi4;

char OutFile[1024], OutFileBase[1024], OutDensFile[1024], SnapshotLoadFile[1024], SnapshotSaveFile[1024], AdunitFile[1024];

int ns, DoInitUpdateProbs, InterruptRun = 0;
unsigned int cntr = 0;
int PlaceDistDistrib[NUM_PLACE_TYPES][MAX_DIST], PlaceSizeDistrib[NUM_PLACE_TYPES][MAX_PLACE_SIZE];


FILE* KMLFile, * KMLFile2;

/* int NumPC,NumPCD; */
#define MAXINTFILE 10

int main(int argc, char* argv[])
{
	char ParamFile[1024]{}, DensityFile[1024]{}, NetworkFile[1024]{}, AirTravelFile[1024]{}, SchoolFile[1024]{}, RegDemogFile[1024]{}, InterventionFile[MAXINTFILE][1024]{}, PreParamFile[1024]{}, buf[2048]{}, * sep;
	int i, GotP, GotPP, GotO, GotD, GotL, GotS, GotA, GotAP, GotScF, GotIF, Perr, cl;

	///// Flags to ensure various parameters have been read; set to false as default. 
	GotP = GotO = GotD = GotL = GotS = GotA = GotAP = GotScF = GotIF = GotPP = 0;

	Perr = 0;
	fprintf(stderr, "sizeof(int)=%i sizeof(long)=%i sizeof(float)=%i sizeof(double)=%i sizeof(unsigned short int)=%i sizeof(int *)=%i\n", (int)sizeof(int), (int)sizeof(long), (int)sizeof(float), (int)sizeof(double), (int)sizeof(unsigned short int), (int)sizeof(int*));
	cl = clock();

	///// Read in command line arguments - lots of things, e.g. random number seeds; (pre)parameter files; binary files; population data; output directory? etc. 

	if (argc < 7)	Perr = 1;
	else
	{
		///// Get seeds. 
		i = argc - 4;
		sscanf(argv[i], "%li", &P.seed1);
		sscanf(argv[i + 1], "%li", &P.seed2);
		sscanf(argv[i + 2], "%li", &P.seed3);
		sscanf(argv[i + 3], "%li", &P.seed4);

		///// Set paramter defaults - read them in after
		P.PlaceCloseIndepThresh = P.LoadSaveNetwork = P.DoHeteroDensity = P.DoPeriodicBoundaries = P.DoSchoolFile = P.DoAdunitDemog = P.OutputDensFile = P.MaxNumThreads = P.DoInterventionFile = 0;
		P.PreControlClusterIdCaseThreshold = 0;
		P.R0scale = 1.0;
		P.KernelOffsetScale = P.KernelPowerScale = 1.0; //added this so that kernel parameters are only changed if input from the command line: ggilani - 15/10/2014
		P.DoSaveSnapshot = P.DoLoadSnapshot  = 0; 

		//// scroll through command line arguments, anticipating what they can be using various if statements. 
		for (i = 1; i < argc - 4; i++)
		{
			if ((argv[i][0] != '/') && ((argv[i][2] != ':') && (argv[i][3] != ':'))) Perr = 1;
			if (argv[i][1] == 'P' && argv[i][2] == ':')
			{
				GotP = 1;
				sscanf(&argv[i][3], "%s", ParamFile);
			}
			else if (argv[i][1] == 'O' && argv[i][2] == ':')
			{
				GotO = 1;
				sscanf(&argv[i][3], "%s", OutFileBase);
			}
			else if (argv[i][1] == 'D' && argv[i][2] == ':')
			{
				GotD = 1;
				sscanf(&argv[i][3], "%s", DensityFile);
				P.DoHeteroDensity = 1;
				P.DoPeriodicBoundaries = 0;
			}
			else if (argv[i][1] == 'A' && argv[i][2] == ':')
			{
				GotA = 1;
				sscanf(&argv[i][3], "%s", AdunitFile);
			}
			else if (argv[i][1] == 'L' && argv[i][2] == ':')
			{
				GotL = 1;
				P.LoadSaveNetwork = 1;
				sscanf(&argv[i][3], "%s", NetworkFile);
			}
			else if (argv[i][1] == 'S' && argv[i][2] == ':')
			{
				P.LoadSaveNetwork = 2;
				GotS = 1;
				sscanf(&argv[i][3], "%s", NetworkFile);
			}
			else if (argv[i][1] == 'R' && argv[i][2] == ':')
			{
				sscanf(&argv[i][3], "%lf", &P.R0scale);
			}
			else if (argv[i][1] == 'K' && argv[i][2] == 'P' && argv[i][3] == ':') //added Kernel Power and Offset scaling so that it can easily be altered from the command line in order to vary the kernel quickly: ggilani - 15/10/14
			{
				sscanf(&argv[i][4], "%lf", &P.KernelPowerScale);
			}
			else if (argv[i][1] == 'K' && argv[i][2] == 'O' && argv[i][3] == ':')
			{
				sscanf(&argv[i][4], "%lf", &P.KernelOffsetScale);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '1' && argv[i][5] == ':') // generic command line specified param - matched to #1 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP1);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '2' && argv[i][5] == ':') // generic command line specified param - matched to #2 in param file
			{
				sscanf(&argv[i][6], "%lf", &P.clP2);
			}
			else if(argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '3' && argv[i][5] == ':') // generic command line specified param - matched to #3 in param file
				{
				sscanf(&argv[i][6], "%lf", &P.clP3);
				}
			else if(argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '4' && argv[i][5] == ':') // generic command line specified param - matched to #4 in param file
				{
				sscanf(&argv[i][6], "%lf", &P.clP4);
				}
			else if(argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '5' && argv[i][5] == ':') // generic command line specified param - matched to #5 in param file
				{
				sscanf(&argv[i][6], "%lf", &P.clP5);
				}
			else if(argv[i][1] == 'C' && argv[i][2] == 'L' && argv[i][3] == 'P' && argv[i][4] == '6' && argv[i][5] == ':') // generic command line specified param - matched to #6 in param file
				{
				sscanf(&argv[i][6], "%lf", &P.clP6);
				}
			else if (argv[i][1] == 'A' && argv[i][2] == 'P' && argv[i][3] == ':')
			{
				GotAP = 1;
				sscanf(&argv[i][3], "%s", AirTravelFile);
			}
			else if (argv[i][1] == 's' && argv[i][2] == ':')
			{
				GotScF = 1;
				sscanf(&argv[i][3], "%s", SchoolFile);
			}
			else if (argv[i][1] == 'T' && argv[i][2] == ':')
			{
				sscanf(&argv[i][3], "%i", &P.PreControlClusterIdCaseThreshold);
			}
			else if (argv[i][1] == 'C' && argv[i][2] == ':')
			{
				sscanf(&argv[i][3], "%i", &P.PlaceCloseIndepThresh);
			}
			else if (argv[i][1] == 'd' && argv[i][2] == ':')
			{
				P.DoAdunitDemog = 1;
				sscanf(&argv[i][3], "%s", RegDemogFile);
			}
			else if (argv[i][1] == 'c' && argv[i][2] == ':')
			{
				sscanf(&argv[i][3], "%i", &P.MaxNumThreads);
			}
			else if (argv[i][1] == 'M' && argv[i][2] == ':')
			{
				P.OutputDensFile = 1;
				sscanf(&argv[i][3], "%s", OutDensFile);
			}
			else if (argv[i][1] == 'I' && argv[i][2] == ':')
			{
				GotIF = 1;
				sscanf(&argv[i][3], "%s", InterventionFile[P.DoInterventionFile]);
				P.DoInterventionFile++;
			}
			else if (argv[i][1] == 'L' && argv[i][2] == 'S' && argv[i][3] == ':')
			{
				sscanf(&argv[i][4], "%s", SnapshotLoadFile);
				P.DoLoadSnapshot = 1;
			}
			else if (argv[i][1] == 'P' && argv[i][2] == 'P' && argv[i][3] == ':')
			{
				sscanf(&argv[i][4], "%s", PreParamFile);
				GotPP = 1;
			}
			else if (argv[i][1] == 'S' && argv[i][2] == 'S' && argv[i][3] == ':')
			{
				sscanf(&argv[i][4], "%s", buf);
				fprintf(stderr, "### %s\n", buf);
				sep = strchr(buf, ',');
				if (!sep)
					Perr = 1;
				else
				{
					P.DoSaveSnapshot = 1;
					*sep = ' ';
					sscanf(buf, "%lf %s", &(P.SnapshotSaveTime), SnapshotSaveFile);
				}
			}
		}
		if (((GotS) && (GotL)) || (!GotP) || (!GotO)) Perr = 1;
	}

	///// END Read in command line arguments

	sprintf(OutFile, "%s", OutFileBase);

	fprintf(stderr, "Param=%s\nOut=%s\nDens=%s\n", ParamFile, OutFile, DensityFile);
	if (Perr) ERR_CRITICAL("Syntax:\nSpatialSim /P:ParamFile /O:OutputFile [/A:AirTravelFile] [/s:SchoolFile] [/D:DensityFile] [/L:NetworkFileToLoad | /S:NetworkFileToSave] [/R:R0scaling] Seed1 Seed2 Seed3 Seed4\n");

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
	if (!GotPP)
	{
		sprintf(PreParamFile, ".." DIRECTORY_SEPARATOR "Pre_%s", ParamFile);
	}

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** READ IN PARAMETERS, DATA ETC.
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 


	ReadParams(ParamFile, PreParamFile);
	if (GotScF) P.DoSchoolFile = 1;
	if (P.DoAirports)
	{
		if (!GotAP) ERR_CRITICAL("Syntax:\nSpatialSim /P:ParamFile /O:OutputFile /AP:AirTravelFile [/s:SchoolFile] [/D:DensityFile] [/L:NetworkFileToLoad | /S:NetworkFileToSave] [/R:R0scaling] Seed1 Seed2 Seed3 Seed4\n");
		ReadAirTravel(AirTravelFile);
	}

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** INITIALIZE
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 

	//initialise random number counts to zero
	count_ranf = 0;
	count_ranfmt = 0;

	///// initialize model (for all realisations). 
	SetupModel(DensityFile, NetworkFile, SchoolFile, RegDemogFile);

	for (i = 0; i < MAX_ADUNITS; i++) AdUnits[i].NI = 0;
	if (P.DoInterventionFile > 0)
		for (i = 0; i < P.DoInterventionFile; i++)
			ReadInterventions(InterventionFile[i]);

	fprintf(stderr, "Model setup in %lf seconds\n", ((double)(clock() - cl)) / CLOCKS_PER_SEC);

	//print out number of calls to random number generator

	///// Set seeds
	if (!P.ResetSeeds)	setall(P.seed3, P.seed4);

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** RUN MODEL
	//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** 


	P.NRactE = P.NRactNE = 0;
	for (i = 0; (i < P.NR) && (P.NRactNE < P.NRN) && (!InterruptRun); i++)
	{
		if (P.NR > 1)
		{
			sprintf(OutFile, "%s.%i", OutFileBase, i);
			fprintf(stderr, "Realisation %i   (time=%lf nr_ne=%i)\n", i + 1, ((double)(clock() - cl)) / CLOCKS_PER_SEC, P.NRactNE);
		}

		///// Set and save seeds
		if (P.ResetSeeds)
		{
			if (P.KeepSameSeeds)
			{
				P.newseed1 = P.seed3;
				P.newseed2 = P.seed4;
			}
			else
			{
				if (i == 0)
				{
					P.newseed1 = P.seed3;
					P.newseed2 = P.seed4;
				}
				else
				{
					//sample 4 new seeds to use as random number generators
					P.newseed1 = (int)(ranf() * 1e8);
					P.newseed2 = (int)(ranf() * 1e8);
				}
			}
			//save these seeds to file
			SaveRandomSeeds();
			//reset seeds
			setall(P.newseed1, P.newseed2);
			//fprintf(stderr, "%i, %i\n", P.newseed1,P.newseed2);
			//fprintf(stderr, "%f\n", ranf());
		}

		///// initialize model (for this realisation). 
		InitModel(i); //passing run number into RunModel so we can save run number in the infection event log: ggilani - 15/10/2014
		//fprintf(stderr, "%f\n", ranf());
		//fprintf(stderr, "%f\n", ranf_mt(5));
		if (P.DoLoadSnapshot) LoadSnapshot();
		if (RunModel(i))
		{  // has been interrupted to reset holiday time
			InitModel(i);
			RunModel(i);
		}
		if (P.OuptutNonSummaryResults)
		{
			if (((!TimeSeries[P.NumSamples - 1].extinct) || (!P.OutputOnlyNonExtinct)) && (P.OutputEveryRealisation))
			{
				SaveResults();
			}
		}
		if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun == 1))
		{
			SaveEvents();
		}
	}
	sprintf(OutFile, "%s", OutFileBase);

	//Calculate origin destination matrix if needed
	if ((P.DoAdUnits) && (P.DoOriginDestinationMatrix))
	{
		CalcOriginDestMatrix_adunit();
		SaveOriginDestMatrix();
	}

	P.NRactual = P.NRactNE;
	TSMean = TSMeanNE; TSVar = TSVarNE;
	if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun == 0))
	{
		SaveEvents();
	}
	sprintf(OutFile, "%s.avNE", OutFileBase);
	SaveSummaryResults();
	P.NRactual = P.NRactE;
	TSMean = TSMeanE; TSVar = TSVarE;
	sprintf(OutFile, "%s.avE", OutFileBase);
	//SaveSummaryResults();


#ifdef WIN32_BM
	Gdiplus::GdiplusShutdown(m_gdiplusToken);
#endif
	fprintf(stderr, "Extinction in %i out of %i runs\n", P.NRactE, P.NRactNE + P.NRactE);
	fprintf(stderr, "Model ran in %lf seconds\n", ((double)(clock() - cl)) / CLOCKS_PER_SEC);
	fprintf(stderr, "Model finished\n");
}


void ReadParams(char* ParamFile, char* PreParamFile)
{
	FILE* dat, * dat2, *dat3;
	double s, t, AgeSuscScale;
	int i, j, k, f, nc, na;
	char CountryNameBuf[128 * MAX_COUNTRIES], AdunitListNamesBuf[128 * MAX_ADUNITS];
	char** AdunitNames, * AdunitNamesBuf;

	char* CountryNames[MAX_COUNTRIES];
	for (i = 0; i < MAX_COUNTRIES; i++) { CountryNames[i] = CountryNameBuf + 128 * i; CountryNames[i][0] = 0; }
	char* AdunitListNames[MAX_ADUNITS];
	for (i = 0; i < MAX_ADUNITS; i++) { AdunitListNames[i] = AdunitListNamesBuf + 128 * i; AdunitListNames[i][0] = 0; }
	if (!(dat = fopen(ParamFile, "r"))) ERR_CRITICAL("Unable to open parameter file\n");
	dat2 = fopen(PreParamFile, "r");
	if (!(dat3 = fopen(AdunitFile, "r"))) dat3 = dat2;
	AgeSuscScale = 1.0;
		GetInputParameter(dat, dat2, "Update timestep", "%lf", (void*) & (P.TimeStep), 1, 1, 0);
	GetInputParameter(dat, dat2, "Sampling timestep", "%lf", (void*) & (P.SampleStep), 1, 1, 0);
	if (P.TimeStep > P.SampleStep) ERR_CRITICAL("Update step must be smaller than sampling step\n");
	t = ceil(P.SampleStep / P.TimeStep - 1e-6);
	P.UpdatesPerSample = (int)t;
	P.TimeStep = P.SampleStep / t;
	P.TimeStepsPerDay = ceil(1.0 / P.TimeStep - 1e-6);
	P.TimeStepsPerDayInt = (int)P.TimeStepsPerDay;
	P.TimeStepsPerYear = P.TimeStepsPerDay * DAYS_PER_YEAR;
	fprintf(stderr, "Update step = %lf\nSampling step = %lf\nUpdates per sample=%i\nTimeStepsPerDay=%lf\n", P.TimeStep, P.SampleStep, P.UpdatesPerSample, P.TimeStepsPerDay);
	GetInputParameter(dat, dat2, "Sampling time", "%lf", (void*) & (P.SampleTime), 1, 1, 0);
	P.NumSamples = 1 + (int)ceil(P.SampleTime / P.SampleStep);
	GetInputParameter(dat, dat3, "Population size", "%i", (void*) & (P.N), 1, 1, 0);
	GetInputParameter(dat, dat2, "Number of realisations", "%i", (void*) & (P.NR), 1, 1, 0);
	if (!GetInputParameter2(dat, dat2, "Number of non-extinct realisations", "%i", (void*) & (P.NRN), 1, 1, 0)) P.NRN = P.NR;
	if (!GetInputParameter2(dat, dat2, "Maximum number of cases defining small outbreak", "%i", (void*) & (P.SmallEpidemicCases), 1, 1, 0)) P.SmallEpidemicCases = -1;
	P.NC = -1;
	GetInputParameter(dat, dat2, "Number of micro-cells per spatial cell width", "%i", (void*) & (P.NMCL), 1, 1, 0);
	//added parameter to reset seeds after every run
	if (!GetInputParameter2(dat, dat2, "Reset seeds for every run", "%i", (void*) & (P.ResetSeeds), 1, 1, 0)) P.ResetSeeds = 0;
	if (P.ResetSeeds)
	{
		if (!GetInputParameter2(dat, dat2, "Keep same seeds for every run", "%i", (void*) & (P.KeepSameSeeds), 1, 1, 0)) P.KeepSameSeeds = 0; //added this to control which seeds are used: ggilani 27/11/19
	}
	if (!GetInputParameter2(dat, dat2, "Reset seeds after intervention", "%i", (void*) & (P.ResetSeedsPostIntervention), 1, 1, 0)) P.ResetSeedsPostIntervention = 1;
	if (P.ResetSeedsPostIntervention)
	{
		if (!GetInputParameter2(dat, dat2, "Time to reset seeds after intervention", "%i", (void*) & (P.TimeToResetSeeds), 1, 1, 0)) P.TimeToResetSeeds = 1000000;
	}
	if (!GetInputParameter2(dat, dat3, "Include households", "%i", (void*) & (P.DoHouseholds), 1, 1, 0)) P.DoHouseholds = 1;

	if (!GetInputParameter2(dat, dat3, "OutputAge"					, "%i", (void*) & (P.OutputAge)					, 1, 1, 0)) P.OutputAge = 1;				//// ON  by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputSeverityAdminUnit"	, "%i", (void*) & (P.OutputSeverityAdminUnit)	, 1, 1, 0)) P.OutputSeverityAdminUnit = 1;	//// ON  by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputR0"					, "%i", (void*) & (P.OutputR0)					, 1, 1, 0)) P.OutputR0 = 0;				    //// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputControls"				, "%i", (void*) & (P.OutputControls)			, 1, 1, 0)) P.OutputControls = 0;		    //// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputCountry"				, "%i", (void*) & (P.OutputCountry)				, 1, 1, 0)) P.OutputCountry = 0;		    //// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputAdUnitVar"			, "%i", (void*) & (P.OutputAdUnitVar)			, 1, 1, 0)) P.OutputAdUnitVar = 0;		    //// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputHousehold"			, "%i", (void*) & (P.OutputHousehold)			, 1, 1, 0)) P.OutputHousehold = 0;		    //// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputInfType"				, "%i", (void*) & (P.OutputInfType)				, 1, 1, 0)) P.OutputInfType = 0;		    //// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OutputNonSeverity"			, "%i", (void*) & (P.OutputNonSeverity)			, 1, 1, 0)) P.OutputNonSeverity = 0;		//// OFF by defualt. 
	if (!GetInputParameter2(dat, dat3, "OuptutNonSummaryResults"	, "%i", (void*) & (P.OuptutNonSummaryResults)	, 1, 1, 0)) P.OuptutNonSummaryResults = 0;	//// OFF by defualt. 

	if (P.DoHouseholds)
	{
		GetInputParameter(dat, dat3, "Household size distribution", "%lf", (void*)P.HouseholdSizeDistrib[0], MAX_HOUSEHOLD_SIZE, 1, 0);
		GetInputParameter(dat, dat2, "Household attack rate", "%lf", (void*) & (P.HouseholdTrans), 1, 1, 0);
		GetInputParameter(dat, dat2, "Household transmission denominator power", "%lf", (void*) & (P.HouseholdTransPow), 1, 1, 0);
		if (!GetInputParameter2(dat, dat3, "Correct age distribution after household allocation to exactly match specified demography", "%i", (void*)&(P.DoCorrectAgeDist), 1, 1, 0)) P.DoCorrectAgeDist = 0;	
	}
	else
	{
		P.HouseholdTrans = 0.0;
		P.HouseholdTransPow = 1.0;
		P.HouseholdSizeDistrib[0][0] = 1.0;
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P.HouseholdSizeDistrib[0][i] = 0;
	}
	for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
		P.HouseholdSizeDistrib[0][i] = P.HouseholdSizeDistrib[0][i] + P.HouseholdSizeDistrib[0][i - 1];
	for (i = 0; i < MAX_HOUSEHOLD_SIZE; i++)
		P.HouseholdDenomLookup[i] = 1 / pow(((double)(i + 1)), P.HouseholdTransPow);
	if (!GetInputParameter2(dat, dat3, "Include administrative units within countries", "%i", (void*) & (P.DoAdUnits), 1, 1, 0)) P.DoAdUnits = 1;
	if (!GetInputParameter2(dat, dat3, "Divisor for countries", "%i", (void*) & (P.CountryDivisor), 1, 1, 0)) P.CountryDivisor = 1;
	if (P.DoAdUnits)
	{
		if (!(AdunitNames = (char**)malloc(3 * ADUNIT_LOOKUP_SIZE * sizeof(char*)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		if (!(AdunitNamesBuf = (char*)malloc(3 * ADUNIT_LOOKUP_SIZE * 360 * sizeof(char)))) ERR_CRITICAL("Unable to allocate temp storage\n");

		for (i = 0; i < ADUNIT_LOOKUP_SIZE; i++)
		{
			P.AdunitLevel1Lookup[i] = -1;
			AdunitNames[3 * i] = AdunitNamesBuf + 3 * i * 360;
			AdunitNames[3 * i + 1] = AdunitNamesBuf + 3 * i * 360 + 60;
			AdunitNames[3 * i + 2] = AdunitNamesBuf + 3 * i * 360 + 160;
		}
		if (!GetInputParameter2(dat, dat3, "Divisor for level 1 administrative units", "%i", (void*)&(P.AdunitLevel1Divisor), 1, 1, 0)) P.AdunitLevel1Divisor = 1;
		if (!GetInputParameter2(dat, dat3, "Mask for level 1 administrative units", "%i", (void*)&(P.AdunitLevel1Mask), 1, 1, 0)) P.AdunitLevel1Mask = 1000000000;
		na = (GetInputParameter2(dat, dat3, "Codes and country/province names for admin units", "%s", (void*)AdunitNames, 3 * ADUNIT_LOOKUP_SIZE, 1, 0)) / 3;
		if (!GetInputParameter2(dat, dat3, "Number of countries to include", "%i", (void*)&nc, 1, 1, 0)) nc = 0;
		if ((na > 0) && (nc>0))
		{
			P.DoAdunitBoundaries = (nc > 0);
			nc = abs(nc);
			GetInputParameter(dat, dat3, "List of names of countries to include", "%s", (nc > 1) ? ((void*)CountryNames) : ((void*)CountryNames[0]), nc, 1, 0);
			P.NumAdunits = 0;
			for (i = 0; i < na; i++)
				for (j = 0; j < nc; j++)
				{
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
		}
		else
		{
			if (!GetInputParameter2(dat, dat3, "Number of level 1 administrative units to include", "%i", (void*) & (P.NumAdunits), 1, 1, 0)) P.NumAdunits = 0;
			if (P.NumAdunits > 0)
			{
				P.DoAdunitBoundaries = 1;
				if (P.NumAdunits > MAX_ADUNITS) ERR_CRITICAL("MAX_ADUNITS too small.\n");
				GetInputParameter(dat, dat3, "List of level 1 administrative units to include", "%s", (P.NumAdunits > 1) ? ((void*)AdunitListNames) : ((void*)AdunitListNames[0]), P.NumAdunits, 1, 0);
				na = P.NumAdunits;
				for (i = 0; i < P.NumAdunits; i++)
				{
					f = 0;
					if (na > 0)
					{
						for (j = 0; (j < na) && (!f); j++) f = (!strcmp(AdunitNames[3 * j + 2], AdunitListNames[i]));
						if(f) k = atoi(AdunitNames[3 * (j-1)]);
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
			free(AdunitNames);
			free(AdunitNamesBuf);
		}
		
		if (!GetInputParameter2(dat, dat3, "Output incidence by administrative unit", "%i", (void*) & (P.DoAdunitOutput), 1, 1, 0)) P.DoAdunitOutput = 0;
		if (!GetInputParameter2(dat, dat3, "Draw administrative unit boundaries on maps", "%i", (void*) & (P.DoAdunitBoundaryOutput), 1, 1, 0)) P.DoAdunitBoundaryOutput = 0;
		if (!GetInputParameter2(dat, dat3, "Correct administrative unit populations", "%i", (void*) & (P.DoCorrectAdunitPop), 1, 1, 0)) P.DoCorrectAdunitPop = 0;
		if (!GetInputParameter2(dat, dat3, "Fix population size at specified value", "%i", (void*) & (P.DoSpecifyPop), 1, 1, 0)) P.DoSpecifyPop = 0;
		fprintf(stderr, "Using %i administrative units\n", P.NumAdunits);
		if (!GetInputParameter2(dat, dat3, "Divisor for administrative unit codes for boundary plotting on bitmaps", "%i", (void*) & (P.AdunitBitmapDivisor), 1, 1, 0)) P.AdunitBitmapDivisor = 1;
		if (!GetInputParameter2(dat, dat2, "Only output household to place distance distribution for one administrative unit", "%i", (void*) & (P.DoOutputPlaceDistForOneAdunit), 1, 1, 0)) P.DoOutputPlaceDistForOneAdunit = 0;
		if (P.DoOutputPlaceDistForOneAdunit)
		{
			if (!GetInputParameter2(dat, dat2, "Administrative unit for which household to place distance distribution to be output", "%i", (void*) & (P.OutputPlaceDistAdunit), 1, 1, 0)) P.DoOutputPlaceDistForOneAdunit = 0;
		}
	}
	else
	{
		P.DoAdunitBoundaries = P.DoAdunitBoundaryOutput = P.DoAdunitOutput = P.DoCorrectAdunitPop = P.DoSpecifyPop = 0; 
		P.AdunitLevel1Divisor = 1; P.AdunitLevel1Mask = 1000000000; 
		P.AdunitBitmapDivisor = P.AdunitLevel1Divisor;
	}

	if (!GetInputParameter2(dat, dat3, "Include age", "%i", (void*) & (P.DoAge), 1, 1, 0)) P.DoAge = 1;
	if (!P.DoAge)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.PropAgeGroup[0][i] = 1.0 / NUM_AGE_GROUPS;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
		{
			P.InitialImmunity[i] = 0;
			P.AgeInfectiousness[i] = P.AgeSusceptibility[i] = 1;
			P.RelativeSpatialContact[i] = P.RelativeTravelRate[i] = 1.0;
		}
	}
	else
	{

		if (P.DoHouseholds)
		{
			if (!GetInputParameter2(dat, dat2, "Initial immunity applied to all household members", "%i", (void*) & (P.DoWholeHouseholdImmunity), 1, 1, 0)) P.DoWholeHouseholdImmunity = 0;
		}
		else
			P.DoWholeHouseholdImmunity = 0;
		if (!GetInputParameter2(dat, dat2, "Initial immunity profile by age", "%lf", (void*)P.InitialImmunity, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.InitialImmunity[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Relative spatial contact rates by age", "%lf", (void*)P.RelativeSpatialContact, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.RelativeSpatialContact[i] = 1;
		if (!GetInputParameter2(dat, dat2, "Age-dependent infectiousness", "%lf", (void*)P.AgeInfectiousness, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.AgeInfectiousness[i] = 1.0;
		if (!GetInputParameter2(dat, dat2, "Age-dependent susceptibility", "%lf", (void*)P.AgeSusceptibility, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.AgeSusceptibility[i] = 1.0;
		GetInputParameter(dat, dat3, "Age distribution of population", "%lf", (void*)P.PropAgeGroup[0], NUM_AGE_GROUPS, 1, 0);
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
		if (!GetInputParameter2(dat, dat3, "Relative travel rates by age", "%lf", (void*)P.RelativeTravelRate, NUM_AGE_GROUPS, 1, 0))
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.RelativeTravelRate[i] = 1;
		if (!GetInputParameter2(dat, dat2, "WAIFW matrix", "%lf", (void*)P.WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0))
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
	if (!GetInputParameter2(dat, dat3, "Include spatial transmission", "%i", (void*) & (P.DoSpatial), 1, 1, 0)) P.DoSpatial = 1;
	GetInputParameter(dat, dat3, "Kernel type", "%i", (void*) & (P.MoveKernelType), 1, 1, 0);
	GetInputParameter(dat, dat3, "Kernel scale", "%lf", (void*) & (P.MoveKernelScale), 1, 1, 0);
	if (P.KernelOffsetScale != 1)
	{
		P.MoveKernelScale *= P.KernelOffsetScale;
	}
	if (!GetInputParameter2(dat, dat3, "Kernel 3rd param", "%lf", (void*) & (P.MoveKernelP3), 1, 1, 0)) P.MoveKernelP3 = 0;
	if (!GetInputParameter2(dat, dat3, "Kernel 4th param", "%lf", (void*) & (P.MoveKernelP4), 1, 1, 0)) P.MoveKernelP4 = 0;
	if (!GetInputParameter2(dat, dat3, "Kernel Shape", "%lf", (void*) & (P.MoveKernelShape), 1, 1, 0)) P.MoveKernelShape = 1.0;
	if (P.KernelPowerScale != 1)
	{
		P.MoveKernelShape *= P.KernelPowerScale;
	}
	if (!GetInputParameter2(dat, dat3, "Airport Kernel Type", "%i", (void*) & (P.AirportKernelType), 1, 1, 0)) P.AirportKernelType = P.MoveKernelType;
	if (!GetInputParameter2(dat, dat3, "Airport Kernel Scale", "%lf", (void*) & (P.AirportKernelScale), 1, 1, 0)) P.AirportKernelScale = P.MoveKernelScale;
	if (!GetInputParameter2(dat, dat3, "Airport Kernel Shape", "%lf", (void*) & (P.AirportKernelShape), 1, 1, 0)) P.AirportKernelShape = P.MoveKernelShape;
	if (!GetInputParameter2(dat, dat3, "Airport Kernel 3rd param", "%lf", (void*) & (P.AirportKernelP3), 1, 1, 0)) P.AirportKernelP3 = P.MoveKernelP3;
	if (!GetInputParameter2(dat, dat3, "Airport Kernel 4th param", "%lf", (void*) & (P.AirportKernelP4), 1, 1, 0)) P.AirportKernelP4 = P.MoveKernelP4;

	if (!GetInputParameter2(dat, dat3, "Include places", "%i", (void*)&(P.DoPlaces), 1, 1, 0)) P.DoPlaces = 1; 
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat3, "Number of types of places", "%i", (void*)&(P.PlaceTypeNum), 1, 1, 0)) P.PlaceTypeNum = 0;
		if (P.PlaceTypeNum == 0) P.DoPlaces = P.DoAirports = 0;
	}
	else
		P.PlaceTypeNum = P.DoAirports = 0;
	if (P.DoPlaces)
	{
		if (P.PlaceTypeNum > NUM_PLACE_TYPES) ERR_CRITICAL("Too many place types\n");
		GetInputParameter(dat, dat3, "Minimum age for age group 1 in place types", "%i", (void*)P.PlaceTypeAgeMin, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat3, "Maximum age for age group 1 in place types", "%i", (void*)P.PlaceTypeAgeMax, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat3, "Proportion of age group 1 in place types", "%lf", (void*) & (P.PlaceTypePropAgeGroup), P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat3, "Proportion of age group 2 in place types", "%lf", (void*) & (P.PlaceTypePropAgeGroup2), P.PlaceTypeNum, 1, 0))
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
			GetInputParameter(dat, dat3, "Minimum age for age group 2 in place types", "%i", (void*)P.PlaceTypeAgeMin2, P.PlaceTypeNum, 1, 0);
			GetInputParameter(dat, dat3, "Maximum age for age group 2 in place types", "%i", (void*)P.PlaceTypeAgeMax2, P.PlaceTypeNum, 1, 0);
		}
		if (!GetInputParameter2(dat, dat3, "Proportion of age group 3 in place types", "%lf", (void*) & (P.PlaceTypePropAgeGroup3), P.PlaceTypeNum, 1, 0))
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
			GetInputParameter(dat, dat3, "Minimum age for age group 3 in place types", "%i", (void*)P.PlaceTypeAgeMin3, P.PlaceTypeNum, 1, 0);
			GetInputParameter(dat, dat3, "Maximum age for age group 3 in place types", "%i", (void*)P.PlaceTypeAgeMax3, P.PlaceTypeNum, 1, 0);
		}
		if (!GetInputParameter2(dat, dat3, "Kernel shape params for place types", "%lf", (void*) & (P.PlaceTypeKernelShape), P.PlaceTypeNum, 1, 0))
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++)
			{
				P.PlaceTypeKernelShape[i] = P.MoveKernelShape;
				P.PlaceTypeKernelScale[i] = P.MoveKernelScale;
			}
		}
		else
			GetInputParameter(dat, dat3, "Kernel scale params for place types", "%lf", (void*) & (P.PlaceTypeKernelScale), P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat3, "Kernel 3rd param for place types", "%lf", (void*) & (P.PlaceTypeKernelP3), P.PlaceTypeNum, 1, 0))
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++)
			{
				P.PlaceTypeKernelP3[i] = P.MoveKernelP3;
				P.PlaceTypeKernelP4[i] = P.MoveKernelP4;
			}
		}
		else
			GetInputParameter(dat, dat3, "Kernel 4th param for place types", "%lf", (void*) & (P.PlaceTypeKernelP4), P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat3, "Number of closest places people pick from (0=all) for place types", "%i", (void*) & (P.PlaceTypeNearestNeighb), P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeNearestNeighb[i] = 0;
		if (P.DoAdUnits)
		{
			if (!GetInputParameter2(dat, dat2, "Degree to which crossing administrative unit boundaries to go to places is inhibited", "%lf", (void*) & (P.InhibitInterAdunitPlaceAssignment), P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++)
					P.InhibitInterAdunitPlaceAssignment[i] = 0;
		}

#ifdef FAST_US
		P.DoAirports = 0;
#else
		if (!GetInputParameter2(dat, dat2, "Include air travel", "%i", (void*)&(P.DoAirports), 1, 1, 0)) P.DoAirports = 0;
#endif
		if (P.DoAirports)
		{
			if (!GetInputParameter2(dat, dat2, "Scaling factor for input file to convert to daily traffic", "%lf", (void*) & (P.AirportTrafficScale), 1, 1, 0)) P.AirportTrafficScale = 1.0;
			if (!GetInputParameter2(dat, dat2, "Proportion of hotel attendees who are local", "%lf", (void*) & (P.HotelPropLocal), 1, 1, 0)) P.HotelPropLocal = 0;
			if (!GetInputParameter2(dat, dat2, "Distribution of duration of air journeys", "%lf", (void*) & (P.JourneyDurationDistrib), MAX_TRAVEL_TIME, 1, 0))
			{
				P.JourneyDurationDistrib[0] = 1;
				for (i = 0; i < MAX_TRAVEL_TIME; i++)
					P.JourneyDurationDistrib[i] = 0;
			}
			if (!GetInputParameter2(dat, dat2, "Distribution of duration of local journeys", "%lf", (void*) & (P.LocalJourneyDurationDistrib), MAX_TRAVEL_TIME, 1, 0))
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
			for (i = j = 0; i <= 1024; i++)
			{
				s = ((double)i) / 1024;
				while (P.JourneyDurationDistrib[j] < s)j++;
				P.InvJourneyDurationDistrib[i] = j;
			}
			for (i = j = 0; i <= 1024; i++)
			{
				s = ((double)i) / 1024;
				while (P.LocalJourneyDurationDistrib[j] < s)j++;
				P.InvLocalJourneyDurationDistrib[i] = j;
			}
		}
		GetInputParameter(dat, dat3, "Mean size of place types", "%lf", (void*)P.PlaceTypeMeanSize, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat3, "Param 1 of place group size distribution", "%lf", (void*)P.PlaceTypeGroupSizeParam1, P.PlaceTypeNum, 1, 0);
		if (!GetInputParameter2(dat, dat3, "Power of place size distribution", "%lf", (void*)P.PlaceTypeSizePower, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizePower[i] = 0;
		//added to enable lognormal distribution - ggilani 09/02/17
		if (!GetInputParameter2(dat, dat3, "Standard deviation of place size distribution", "%lf", (void*)P.PlaceTypeSizeSD, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizeSD[i] = 0;
		if (!GetInputParameter2(dat, dat3, "Offset of place size distribution", "%lf", (void*)P.PlaceTypeSizeOffset, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizeOffset[i] = 0;
		if (!GetInputParameter2(dat, dat3, "Maximum of place size distribution", "%lf", (void*)P.PlaceTypeSizeMax, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeSizeMax[i] = 1e20;
		if (!GetInputParameter2(dat, dat3, "Kernel type for place types", "%i", (void*)P.PlaceTypeKernelType, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++)
				P.PlaceTypeKernelType[i] = P.MoveKernelType;
		GetInputParameter(dat, dat3, "Place overlap matrix", "%lf", (void*)P.PlaceExclusivityMatrix, P.PlaceTypeNum * P.PlaceTypeNum, 1, 0); //changed from P.PlaceTypeNum,P.PlaceTypeNum,0);
/* Note P.PlaceExclusivityMatrix not used at present - places assumed exclusive (each person belongs to 0 or 1 place) */

		GetInputParameter(dat, dat2, "Proportion of between group place links", "%lf", (void*)P.PlaceTypePropBetweenGroupLinks, P.PlaceTypeNum, 1, 0);
		GetInputParameter(dat, dat2, "Relative transmission rates for place types", "%lf", (void*)P.PlaceTypeTrans, P.PlaceTypeNum, 1, 0);
		for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceTypeTrans[i] *= AgeSuscScale;
	}
	if (!GetInputParameter2(dat, dat2, "Daily seasonality coefficients", "%lf", (void*)P.Seasonality, DAYS_PER_YEAR, 1, 0))
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
	if (!GetInputParameter2(dat, dat3, "Number of seed locations", "%i", (void*) & (P.NumSeedLocations), 1, 1, 0)) P.NumSeedLocations = 1;
	if (P.NumSeedLocations > MAX_NUM_SEED_LOCATIONS)
	{
		fprintf(stderr, "Too many seed locations\n");
		P.NumSeedLocations = MAX_NUM_SEED_LOCATIONS;
	}
	GetInputParameter(dat, dat3, "Initial number of infecteds", "%i", (void*)P.NumInitialInfections, P.NumSeedLocations, 1, 0);
	if (!GetInputParameter2(dat, dat3, "Location of initial infecteds", "%lf", (void*)&(P.LocationInitialInfection[0][0]), P.NumSeedLocations * 2, 1, 0)) P.LocationInitialInfection[0][0] = P.LocationInitialInfection[0][1] = 0.0;
	if (!GetInputParameter2(dat, dat3, "Minimum population in microcell of initial infection", "%i", (void*) & (P.MinPopDensForInitialInfection), 1, 1, 0)) P.MinPopDensForInitialInfection = 0;
	if (!GetInputParameter2(dat, dat3, "Maximum population in microcell of initial infection", "%i", (void*)&(P.MaxPopDensForInitialInfection), 1, 1, 0)) P.MaxPopDensForInitialInfection = 10000000;
	if (!GetInputParameter2(dat, dat3, "Randomise initial infection location", "%i", (void*) & (P.DoRandomInitialInfectionLoc), 1, 1, 0)) P.DoRandomInitialInfectionLoc=1;
	if (!GetInputParameter2(dat, dat3, "All initial infections located in same microcell", "%i", (void*) & (P.DoAllInitialInfectioninSameLoc), 1, 1, 0)) P.DoAllInitialInfectioninSameLoc=0;
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat3, "Administrative unit to seed initial infection into", "%s", (P.NumSeedLocations > 1) ? ((void*)AdunitListNames) : ((void*)AdunitListNames[0]), P.NumSeedLocations, 1, 0))
			for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
		else
			for (i = 0; i < P.NumSeedLocations; i++)
			{
				f = 0;
				if (P.NumAdunits > 0)
				{
					for (j = 0; (j < P.NumAdunits) && (!f); j++) f = (!strcmp(AdUnits[j].ad_name, AdunitListNames[i]));
					if (f) k = AdUnits[j-1].id;
				}
				if (!f) k = atoi(AdunitListNames[i]);
				P.InitialInfectionsAdminUnit[i] = k;
				P.InitialInfectionsAdminUnitId[i]=P.AdunitLevel1Lookup[(k % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor];
			}
		if(!GetInputParameter2(dat, dat3, "Administrative unit seeding weights", "%lf", (void*) & (P.InitialInfectionsAdminUnitWeight[0]), P.NumSeedLocations, 1, 0))
			for(i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnitWeight[i] = 1.0;
		s=0;
		for(i = 0; i < P.NumSeedLocations; i++) s+=P.InitialInfectionsAdminUnitWeight[i];
		for(i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnitWeight[i]/=s;
	//	for (i = 0; i < P.NumSeedLocations; i++) fprintf(stderr, "## %i %s %i %i %lf\n",i, AdUnits[P.InitialInfectionsAdminUnitId[i]].ad_name, P.InitialInfectionsAdminUnitId[i], P.InitialInfectionsAdminUnit[i], P.InitialInfectionsAdminUnitWeight[i]);
	}
	else
	{
		for (i = 0; i < P.NumSeedLocations; i++) P.InitialInfectionsAdminUnit[i] = 0;
	}
	GetInputParameter(dat, dat2, "Reproduction number", "%lf", (void*) & (P.R0), 1, 1, 0);
	GetInputParameter(dat, dat2, "Infectious period", "%lf", (void*) & (P.InfectiousPeriod), 1, 1, 0);
	if (!GetInputParameter2(dat, dat2, "Assume SIS model", "%i", (void*) & (P.DoSIS), 1, 1, 0)) P.DoSIS = 0;
	if (!GetInputParameter2(dat, dat2, "SD of individual variation in infectiousness", "%lf", (void*) & (P.InfectiousnessSD), 1, 1, 0)) P.InfectiousnessSD = 0;
	if (GetInputParameter2(dat, dat2, "k of individual variation in infectiousness", "%lf", (void*)& s, 1, 1, 0)) P.InfectiousnessSD = 1.0 / sqrt(s);
	if (P.InfectiousnessSD > 0)
	{
		P.InfectiousnessGamA = P.InfectiousnessGamR = 1 / (P.InfectiousnessSD * P.InfectiousnessSD);
	}
	if (P.DoSIS)
	{
		if (!GetInputParameter2(dat, dat2, "Factor by which susceptibility is multiplied per infection", "%lf", (void*) & (P.SuscReductionFactorPerInfection), 1, 1, 0)) P.SuscReductionFactorPerInfection = 1;
	}
	else
		P.SuscReductionFactorPerInfection = 0;
	if (!GetInputParameter2(dat, dat2, "Model time varying infectiousness", "%i", (void*) & (P.DoInfectiousnessProfile), 1, 1, 0)) P.DoInfectiousnessProfile = 0;
	if (!GetInputParameter2(dat, dat2, "Power of scaling of spatial R0 with density", "%lf", (void*) & (P.R0DensityScalePower), 1, 1, 0)) P.R0DensityScalePower = 0;
	if (P.DoInfectiousnessProfile)
	{
		if (!GetInputParameter2(dat, dat2, "Infectiousness profile", "%lf", (void*)P.infectious_prof, INFPROF_RES, 1, 0))
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
		for (i = 0; i <= CDF_RES; i++) P.infectious_icdf[i] = exp(-1.0);
	}
	else
	{
		if (!GetInputParameter2(dat, dat2, "Infectious period inverse CDF", "%lf", (void*)P.infectious_icdf, CDF_RES + 1, 1, 0))
		{
			P.infectious_icdf[CDF_RES] = 100;
			for (i = 0; i < CDF_RES; i++)
				P.infectious_icdf[i] = -log(1 - ((double)i) / CDF_RES);
		}
		k = (int)ceil(P.InfectiousPeriod * P.infectious_icdf[CDF_RES] / P.TimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		for (i = 0; i < k; i++) P.infectiousness[i] = 1.0;
		P.infectiousness[k] = 0;
		for (i = 0; i <= CDF_RES; i++) P.infectious_icdf[i] = exp(-P.infectious_icdf[i]);
	}
	if (!GetInputParameter2(dat, dat2, "Include latent period", "%i", (void*) & (P.DoLatent), 1, 1, 0)) P.DoLatent = 0;
	if (P.DoLatent)
	{
		GetInputParameter(dat, dat2, "Latent period", "%lf", (void*) & (P.LatentPeriod), 1, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Latent period inverse CDF", "%lf", (void*)P.latent_icdf, CDF_RES + 1, 1, 0))
		{
			P.latent_icdf[CDF_RES] = 1e10;
			for (i = 0; i < CDF_RES; i++)
				P.latent_icdf[i] = -log(1 - ((double)i) / CDF_RES);
		}
		for (i = 0; i <= CDF_RES; i++)
			P.latent_icdf[i] = exp(-P.latent_icdf[i]);
	}

	if (!GetInputParameter2(dat, dat2, "Include symptoms", "%i", (void*) & (P.DoSymptoms), 1, 1, 0)) P.DoSymptoms = 0;
	if (!P.DoSymptoms)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.ProportionSymptomatic[i] = 0;
		P.FalsePositiveRate = 0;
		P.SymptInfectiousness = 1.0;
		P.LatentToSymptDelay = 0;
	}
	else
	{
		if (P.DoAge)
			GetInputParameter(dat, dat2, "Proportion symptomatic by age group", "%lf", (void*)P.ProportionSymptomatic, NUM_AGE_GROUPS, 1, 0);
		else
		{
			GetInputParameter(dat, dat2, "Proportion symptomatic", "%lf", (void*)P.ProportionSymptomatic, 1, 1, 0);
			for (i = 1; i < NUM_AGE_GROUPS; i++)
				P.ProportionSymptomatic[i] = P.ProportionSymptomatic[0];
		}
		GetInputParameter(dat, dat2, "Delay from end of latent period to start of symptoms", "%lf", (void*) & (P.LatentToSymptDelay), 1, 1, 0);
		GetInputParameter(dat, dat2, "Relative rate of random contacts if symptomatic", "%lf", (void*) & (P.SymptSpatialContactRate), 1, 1, 0);
		GetInputParameter(dat, dat2, "Symptomatic infectiousness relative to asymptomatic", "%lf", (void*) & (P.SymptInfectiousness), 1, 1, 0);
		if (!GetInputParameter2(dat, dat2, "Model symptomatic withdrawal to home as true absenteeism", "%i", (void*)& P.DoRealSymptWithdrawal, 1, 1, 0)) P.DoRealSymptWithdrawal = 0;
		if (P.DoPlaces)
		{
			GetInputParameter(dat, dat2, "Relative level of place attendance if symptomatic", "%lf", (void*)P.SymptPlaceTypeContactRate, P.PlaceTypeNum, 1, 0);
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
		if (!GetInputParameter2(dat, dat2, "Maximum age of child at home for whom one adult also stays at home", "%i", (void*)& P.CaseAbsentChildAgeCutoff, 1, 1, 0)) P.CaseAbsentChildAgeCutoff = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of children at home for whom one adult also stays at home", "%lf", (void*)& P.CaseAbsentChildPropAdultCarers, 1, 1, 0)) P.CaseAbsentChildPropAdultCarers = 0;
#ifdef ABSENTEEISM_PLACE_CLOSURE
		P.CaseAbsenteeismDelay = 0;  // Set to zero for tracking absenteeism
#else
		if (!GetInputParameter2(dat, dat2, "Delay in starting place absenteeism for cases who withdraw", "%lf", (void*)& P.CaseAbsenteeismDelay, 1, 1, 0)) P.CaseAbsenteeismDelay = 0;
#endif
		if (!GetInputParameter2(dat, dat2, "Duration of place absenteeism for cases who withdraw", "%lf", (void*)& P.CaseAbsenteeismDuration, 1, 1, 0)) P.CaseAbsenteeismDuration = 7;

		if (!GetInputParameter2(dat, dat2, "False positive rate", "%lf", (void*) & (P.FalsePositiveRate), 1, 1, 0)) P.FalsePositiveRate = 0.0;
		if (!GetInputParameter2(dat, dat2, "False positive per capita incidence", "%lf", (void*) & (P.FalsePositivePerCapitaIncidence), 1, 1, 0)) P.FalsePositivePerCapitaIncidence = 0.0;
		if (!GetInputParameter2(dat, dat2, "False positive relative incidence by age", "%lf", (void*)P.FalsePositiveAgeRate, NUM_AGE_GROUPS, 1, 0))
			for (j = 0; j < NUM_AGE_GROUPS; j++) P.FalsePositiveAgeRate[j] = 1.0;
	}

	if(!GetInputParameter2(dat, dat2, "Do Severity Analysis", "%i", (void*) & (P.DoSeverity), 1, 1, 0)) P.DoSeverity = 0;
	if(P.DoSeverity == 1)
	{

		//// Means for icdf's.
		GetInputParameter(dat, dat2, "Mean_MildToRecovery", "%lf", (void*) & (P.Mean_MildToRecovery), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_ILIToRecovery", "%lf", (void*) & (P.Mean_ILIToRecovery), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_SARIToRecovery", "%lf", (void*) & (P.Mean_SARIToRecovery), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_CriticalToCritRecov", "%lf", (void*) & (P.Mean_CriticalToCritRecov), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_CritRecovToRecov", "%lf", (void*) & (P.Mean_CritRecovToRecov), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_ILIToSARI", "%lf", (void*) & (P.Mean_ILIToSARI), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_SARIToCritical", "%lf", (void*) & (P.Mean_SARIToCritical), 1, 1, 0);
		GetInputParameter(dat, dat2, "Mean_CriticalToDeath", "%lf", (void*) & (P.Mean_CriticalToDeath), 1, 1, 0);

		//// Get ICDFs
		if(!GetInputParameter2(dat, dat2, "MildToRecovery_icdf", "%lf", (void*)P.MildToRecovery_icdf, CDF_RES + 1, 1, 0))
			{
			P.MildToRecovery_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.MildToRecovery_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.MildToRecovery_icdf[i] = exp(-P.MildToRecovery_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "ILIToRecovery_icdf", "%lf", (void*)P.ILIToRecovery_icdf, CDF_RES + 1, 1, 0))
			{
			P.ILIToRecovery_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.ILIToRecovery_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.ILIToRecovery_icdf[i] = exp(-P.ILIToRecovery_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "SARIToRecovery_icdf", "%lf", (void*)P.SARIToRecovery_icdf, CDF_RES + 1, 1, 0))
			{
			P.SARIToRecovery_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.SARIToRecovery_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.SARIToRecovery_icdf[i] = exp(-P.SARIToRecovery_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "CriticalToCritRecov_icdf", "%lf", (void*)P.CriticalToCritRecov_icdf, CDF_RES + 1, 1, 0))
			{
			P.CriticalToCritRecov_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.CriticalToCritRecov_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.CriticalToCritRecov_icdf[i] = exp(-P.CriticalToCritRecov_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "CritRecovToRecov_icdf", "%lf", (void*)P.CritRecovToRecov_icdf, CDF_RES + 1, 1, 0))
			{
			P.CritRecovToRecov_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.CritRecovToRecov_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.CritRecovToRecov_icdf[i] = exp(-P.CritRecovToRecov_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "ILIToSARI_icdf", "%lf", (void*)P.ILIToSARI_icdf, CDF_RES + 1, 1, 0))
			{
			P.ILIToSARI_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.ILIToSARI_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.ILIToSARI_icdf[i] = exp(-P.ILIToSARI_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "SARIToCritical_icdf", "%lf", (void*)P.SARIToCritical_icdf, CDF_RES + 1, 1, 0))
			{
			P.SARIToCritical_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.SARIToCritical_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.SARIToCritical_icdf[i] = exp(-P.SARIToCritical_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "CriticalToDeath_icdf", "%lf", (void*)P.CriticalToDeath_icdf, CDF_RES + 1, 1, 0))
			{
			P.CriticalToDeath_icdf[CDF_RES] = 100;
			for(i = 0; i < CDF_RES; i++)
				P.CriticalToDeath_icdf[i] = -log(1 - ((double)i) / CDF_RES);
			}
		for(i = 0; i <= CDF_RES; i++) P.CriticalToDeath_icdf[i] = exp(-P.CriticalToDeath_icdf[i]);

		if(!GetInputParameter2(dat, dat2, "Prop_Mild_ByAge", "%lf", (void*)P.Prop_Mild_ByAge, NUM_AGE_GROUPS, 1, 0))
			for(i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_Mild_ByAge[i] = 0.5;

		if(!GetInputParameter2(dat, dat2, "Prop_ILI_ByAge", "%lf", (void*)P.Prop_ILI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for(i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_ILI_ByAge[i] = 0.3;

		if(!GetInputParameter2(dat, dat2, "Prop_SARI_ByAge", "%lf", (void*)P.Prop_SARI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for(i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_SARI_ByAge[i] = 0.15;

		if(!GetInputParameter2(dat, dat2, "Prop_Critical_ByAge", "%lf", (void*)P.Prop_Critical_ByAge, NUM_AGE_GROUPS, 1, 0))
			for(i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_Critical_ByAge[i] = 0.05;

		if(!GetInputParameter2(dat, dat2, "CFR_SARI_ByAge", "%lf", (void*)P.CFR_SARI_ByAge, NUM_AGE_GROUPS, 1, 0))
			for(i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_Critical_ByAge[i] = 0.50;

		if(!GetInputParameter2(dat, dat2, "CFR_Critical_ByAge", "%lf", (void*)P.CFR_Critical_ByAge, NUM_AGE_GROUPS, 1, 0))
			for(i = 0; i < NUM_AGE_GROUPS; i++)
				P.Prop_Critical_ByAge[i] = 0.50;

	}

	if (!GetInputParameter2(dat, dat2, "Bounding box for bitmap", "%lf", (void*) & (P.BoundingBox[0]), 4, 1, 0))
	{
		P.BoundingBox[0] = P.BoundingBox[1] = 0.0;
		P.BoundingBox[2] = P.BoundingBox[3] = 1.0;
	}
	if (!GetInputParameter2(dat, dat2, "Spatial domain for simulation", "%lf", (void*) & (P.SpatialBoundingBox[0]), 4, 1, 0))
	{
		P.SpatialBoundingBox[0] = P.SpatialBoundingBox[1] = 0.0;
		P.SpatialBoundingBox[2] = P.SpatialBoundingBox[3] = 1.0;
	}
	if (!GetInputParameter2(dat, dat2, "Grid size", "%lf", (void*) & (P.cwidth), 1, 1, 0)) P.cwidth = 1.0 / 120.0;
	if (!GetInputParameter2(dat, dat2, "Use long/lat coord system", "%i", (void*) & (P.DoUTM_coords), 1, 1, 0)) P.DoUTM_coords = 1;
	if (!GetInputParameter2(dat, dat2, "Bitmap scale", "%lf", (void*) & (P.BitmapScale), 1, 1, 0)) P.BitmapScale = 1.0;
	if (!GetInputParameter2(dat, dat2, "Bitmap y:x aspect scaling", "%lf", (void*) & (P.BitmapAspectScale), 1, 1, 0)) P.BitmapAspectScale = 1.0;
	if (!GetInputParameter2(dat, dat2, "Bitmap movie frame interval", "%i", (void*) & (P.BitmapMovieFrame), 1, 1, 0)) P.BitmapMovieFrame = 250;
	if (!GetInputParameter2(dat, dat2, "Output bitmap", "%i", (void*) & (P.OutputBitmap), 1, 1, 0)) P.OutputBitmap = 0;
	if (!GetInputParameter2(dat, dat2, "Output bitmap detected", "%i", (void*) & (P.OutputBitmapDetected), 1, 1, 0)) P.OutputBitmapDetected = 0;
	if (!GetInputParameter2(dat, dat2, "Output immunity on bitmap", "%i", (void*) & (P.DoImmuneBitmap), 1, 1, 0)) P.DoImmuneBitmap = 0;
	if (!GetInputParameter2(dat, dat2, "Output infection tree", "%i", (void*) & (P.DoInfectionTree), 1, 1, 0)) P.DoInfectionTree = 0;
	if (!GetInputParameter2(dat, dat2, "Do one generation", "%i", (void*) & (P.DoOneGen), 1, 1, 0)) P.DoOneGen = 0;
	if (!GetInputParameter2(dat, dat2, "Output every realisation", "%i", (void*) & (P.OutputEveryRealisation), 1, 1, 0)) P.OutputEveryRealisation = 0;
	if (!GetInputParameter2(dat, dat2, "Maximum number to sample for correlations", "%i", (void*) & (P.MaxCorrSample), 1, 1, 0)) P.MaxCorrSample = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Assume SI model", "%i", (void*) & (P.DoSI), 1, 1, 0)) P.DoSI = 0;
	if (!GetInputParameter2(dat, dat2, "Assume periodic boundary conditions", "%i", (void*) & (P.DoPeriodicBoundaries), 1, 1, 0)) P.DoPeriodicBoundaries = 0;
	if (!GetInputParameter2(dat, dat2, "Only output non-extinct realisations", "%i", (void*) & (P.OutputOnlyNonExtinct), 1, 1, 0)) P.OutputOnlyNonExtinct = 0;

	if (!GetInputParameter2(dat, dat2, "Use cases per thousand threshold for area controls", "%i", (void*) & (P.DoPerCapitaTriggers), 1, 1, 0)) P.DoPerCapitaTriggers = 0;
	if (!GetInputParameter2(dat, dat2, "Use global triggers for interventions", "%i", (void*) & (P.DoGlobalTriggers), 1, 1, 0)) P.DoGlobalTriggers = 0;
	if (!GetInputParameter2(dat, dat2, "Use admin unit triggers for interventions", "%i", (void*) & (P.DoAdminTriggers), 1, 1, 0)) P.DoAdminTriggers = 0;
	if(!GetInputParameter2(dat, dat2, "Use ICU case triggers for interventions", "%i", (void*) & (P.DoICUTriggers), 1, 1, 0)) P.DoICUTriggers = 0;
	if (P.DoGlobalTriggers)  P.DoAdminTriggers = 0;
	if (!GetInputParameter2(dat, dat2, "Divisor for per-capita area threshold (default 1000)", "%i", (void*) & (P.IncThreshPop), 1, 1, 0)) P.IncThreshPop = 1000;
	if (!GetInputParameter2(dat, dat2, "Divisor for per-capita global threshold (default 1000)", "%i", (void*) & (P.GlobalIncThreshPop), 1, 1, 0)) P.GlobalIncThreshPop = 1000;


	if (!GetInputParameter2(dat, dat2, "Number of sampling intervals over which cumulative incidence measured for global trigger", "%i", (void*) & (P.TriggersSamplingInterval), 1, 1, 0)) P.TriggersSamplingInterval = 10000000;
	if (!GetInputParameter2(dat, dat2, "Proportion of cases detected for treatment", "%lf", (void*) & (P.PostAlertControlPropCasesId), 1, 1, 0)) P.PostAlertControlPropCasesId = 1;
	if (!GetInputParameter2(dat, dat2, "Proportion of cases detected before outbreak alert", "%lf", (void*) & (P.PreAlertControlPropCasesId), 1, 1, 0)) P.PreAlertControlPropCasesId = 1.0;
	if (!GetInputParameter2(dat, dat2, "Trigger alert on deaths", "%i", (void*)&(P.PreControlClusterIdUseDeaths), 1, 1, 0)) P.PreControlClusterIdUseDeaths = 0;
	if (P.PreControlClusterIdUseDeaths)
	{
		if (P.PreControlClusterIdCaseThreshold == 0)
			if (!GetInputParameter2(dat, dat2, "Number of deaths accummulated before alert", "%i", (void*)&(P.PreControlClusterIdCaseThreshold), 1, 1, 0)) P.PreControlClusterIdCaseThreshold = 0;
	}
	else if (P.PreControlClusterIdCaseThreshold == 0)
	{
		if (!GetInputParameter2(dat, dat2, "Number of detected cases needed before outbreak alert triggered", "%i", (void*) & (P.PreControlClusterIdCaseThreshold), 1, 1, 0)) P.PreControlClusterIdCaseThreshold = 0;
	}
	if (!GetInputParameter2(dat, dat2, "Number of days to accummulate cases/deaths before alert", "%i", (void*)&(P.PreControlClusterIdDuration), 1, 1, 0)) P.PreControlClusterIdDuration = 1000;
	if (!GetInputParameter2(dat, dat2, "Day of year trigger is reached", "%lf", (void*)&(P.PreControlClusterIdCalTime), 1, 1, 0)) P.PreControlClusterIdCalTime = -1;
	P.PreControlClusterIdHolOffset = P.PreControlClusterIdTime = 0;
	if (!GetInputParameter2(dat, dat2, "Only use confirmed cases to trigger alert", "%i", (void*) & (P.DoEarlyCaseDiagnosis), 1, 1, 0)) P.DoEarlyCaseDiagnosis = 0;
	if (!GetInputParameter2(dat, dat2, "Only treat mixing groups within places", "%i", (void*) & (P.DoPlaceGroupTreat), 1, 1, 0)) P.DoPlaceGroupTreat = 0;

	if (!GetInputParameter2(dat, dat2, "Treatment trigger incidence per cell"				, "%lf", (void*) & (P.TreatCellIncThresh)			, 1, 1, 0)) P.TreatCellIncThresh			= 1000000000;
	if (!GetInputParameter2(dat, dat2, "Case isolation trigger incidence per cell"			, "%lf", (void*) & (P.CaseIsolation_CellIncThresh)	, 1, 1, 0)) P.CaseIsolation_CellIncThresh	= 1000000000;
	if (!GetInputParameter2(dat, dat2, "Household quarantine trigger incidence per cell"			, "%lf", (void*) & (P.HHQuar_CellIncThresh)			, 1, 1, 0)) P.HHQuar_CellIncThresh			= 1000000000;

	if (!GetInputParameter2(dat, dat2, "Relative susceptibility of treated individual", "%lf", (void*) & (P.TreatSuscDrop), 1, 1, 0)) P.TreatSuscDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Relative infectiousness of treated individual", "%lf", (void*) & (P.TreatInfDrop), 1, 1, 0)) P.TreatInfDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases resulting in death prevented by treatment", "%lf", (void*) & (P.TreatDeathDrop), 1, 1, 0)) P.TreatDeathDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases prevented by treatment", "%lf", (void*) & (P.TreatSympDrop), 1, 1, 0)) P.TreatSympDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to treat cell", "%lf", (void*) & (P.TreatDelayMean), 1, 1, 0)) P.TreatDelayMean = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of course of treatment", "%lf", (void*) & (P.TreatCaseCourseLength), 1, 1, 0)) P.TreatCaseCourseLength = 5;
	if (!GetInputParameter2(dat, dat2, "Duration of course of prophylaxis", "%lf", (void*) & (P.TreatProphCourseLength), 1, 1, 0)) P.TreatProphCourseLength = 10;
	if (!GetInputParameter2(dat, dat2, "Proportion of detected cases treated", "%lf", (void*) & (P.TreatPropCases), 1, 1, 0)) P.TreatPropCases = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of households of cases treated", "%lf", (void*) & (P.TreatPropCaseHouseholds), 1, 1, 0)) P.TreatPropCaseHouseholds = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household prophylaxis policy", "%lf", (void*) & (P.TreatHouseholdsDuration), 1, 1, 0)) P.TreatHouseholdsDuration = USHRT_MAX / P.TimeStepsPerDay;
	}
	if (!GetInputParameter2(dat, dat2, "Proportion treated", "%lf", (void*) & (P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	if (!GetInputParameter2(dat, dat2, "Proportion treated in radial prophylaxis", "%lf", (void*) & (P.TreatPropRadial), 1, 1, 0)) P.TreatPropRadial = 1.0;
	if (!GetInputParameter2(dat, dat2, "Treatment radius", "%lf", (void*) & (P.TreatRadius), 1, 1, 0)) P.TreatRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of place/geographic prophylaxis policy", "%lf", (void*) & (P.TreatPlaceGeogDuration), 1, 1, 0)) P.TreatPlaceGeogDuration = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Treatment start time", "%lf", (void*) & (P.TreatTimeStartBase), 1, 1, 0)) P.TreatTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of places treated after case detected", "%lf", (void*)P.TreatPlaceProbCaseId, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.TreatPlaceProbCaseId[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of people treated in targeted places", "%lf", (void*)P.TreatPlaceTotalProp, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.TreatPlaceTotalProp[i] = 0;
	}
	if (!GetInputParameter2(dat, dat2, "Maximum number of doses available", "%lf", (void*) & (P.TreatMaxCoursesBase), 1, 1, 0)) P.TreatMaxCoursesBase = 1e20;
	if (!GetInputParameter2(dat, dat2, "Start time of additional treatment production", "%lf", (void*) & (P.TreatNewCoursesStartTime), 1, 1, 0)) P.TreatNewCoursesStartTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Rate of additional treatment production (courses per day)", "%lf", (void*) & (P.TreatNewCoursesRate), 1, 1, 0)) P.TreatNewCoursesRate = 0;
	if (!GetInputParameter2(dat, dat2, "Maximum number of people targetted with radial prophylaxis per case", "%i", (void*) & (P.TreatMaxCoursesPerCase), 1, 1, 0)) P.TreatMaxCoursesPerCase = 1000000000;


	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Treat administrative units rather than rings", "%i", (void*) & (P.TreatByAdminUnit), 1, 1, 0)) P.TreatByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for treatment", "%i", (void*) & (P.TreatAdminUnitDivisor), 1, 1, 0)) P.TreatAdminUnitDivisor = 1;
		if ((P.TreatAdminUnitDivisor == 0) || (P.TreatByAdminUnit == 0)) { P.TreatByAdminUnit = 0; P.TreatAdminUnitDivisor = 1; }
	}
	else
	{
		P.TreatAdminUnitDivisor = 1; P.TreatByAdminUnit = 0;
	}

	if (!GetInputParameter2(dat, dat2, "Vaccination trigger incidence per cell", "%lf", (void*) & (P.VaccCellIncThresh), 1, 1, 0)) P.VaccCellIncThresh = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Relative susceptibility of vaccinated individual", "%lf", (void*) & (P.VaccSuscDrop), 1, 1, 0)) P.VaccSuscDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Relative susceptibility of individual vaccinated after switch time", "%lf", (void*) & (P.VaccSuscDrop2), 1, 1, 0)) P.VaccSuscDrop2 = 1;
	if (!GetInputParameter2(dat, dat2, "Switch time at which vaccine efficacy increases", "%lf", (void*) & (P.VaccTimeEfficacySwitch), 1, 1, 0)) P.VaccTimeEfficacySwitch = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Decay rate of vaccine efficacy (per year)", "%lf", (void*) & (P.VaccEfficacyDecay), 1, 1, 0)) P.VaccEfficacyDecay = 0;
	P.VaccEfficacyDecay /= DAYS_PER_YEAR;
	if (!GetInputParameter2(dat, dat2, "Relative infectiousness of vaccinated individual", "%lf", (void*) & (P.VaccInfDrop), 1, 1, 0)) P.VaccInfDrop = 1;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases resulting in death prevented by vaccination", "%lf", (void*) & (P.VaccMortDrop), 1, 1, 0)) P.VaccMortDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Proportion of symptomatic cases prevented by vaccination", "%lf", (void*) & (P.VaccSympDrop), 1, 1, 0)) P.VaccSympDrop = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to vaccinate", "%lf", (void*) & (P.VaccDelayMean), 1, 1, 0)) P.VaccDelayMean = 0;

	if (!GetInputParameter2(dat, dat2, "Delay from vaccination to full protection", "%lf", (void*) & (P.VaccTimeToEfficacy), 1, 1, 0)) P.VaccTimeToEfficacy = 0;

	if (!GetInputParameter2(dat, dat2, "Years between rounds of vaccination", "%lf", (void*) & (P.VaccCampaignInterval), 1, 1, 0)) P.VaccCampaignInterval = 1e10;
	if (!GetInputParameter2(dat, dat2, "Max vaccine doses per day", "%i", (void*) & (P.VaccDosePerDay), 1, 1, 0)) P.VaccDosePerDay = -1;
	P.VaccCampaignInterval *= DAYS_PER_YEAR;
	if (!GetInputParameter2(dat, dat2, "Maximum number of rounds of vaccination", "%i", (void*) & (P.VaccMaxRounds), 1, 1, 0)) P.VaccMaxRounds = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of households of cases vaccinated", "%lf", (void*) & (P.VaccPropCaseHouseholds), 1, 1, 0)) P.VaccPropCaseHouseholds = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household vaccination policy", "%lf", (void*) & (P.VaccHouseholdsDuration), 1, 1, 0)) P.VaccHouseholdsDuration = USHRT_MAX / P.TimeStepsPerDay;
	}

	if (!GetInputParameter2(dat, dat2, "Vaccination start time", "%lf", (void*) & (P.VaccTimeStartBase), 1, 1, 0)) P.VaccTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Proportion of population vaccinated", "%lf", (void*) & (P.VaccProp), 1, 1, 0)) P.VaccProp = 0;
	if (!GetInputParameter2(dat, dat2, "Time taken to reach max vaccination coverage (in years)", "%lf", (void*) & (P.VaccCoverageIncreasePeriod), 1, 1, 0)) P.VaccCoverageIncreasePeriod = 0;
	P.VaccCoverageIncreasePeriod *= DAYS_PER_YEAR;
	if (!GetInputParameter2(dat, dat2, "Time to start geographic vaccination", "%lf", (void*) & (P.VaccTimeStartGeo), 1, 1, 0)) P.VaccTimeStartGeo = 1e10;
	if (!GetInputParameter2(dat, dat2, "Vaccination radius", "%lf", (void*) & (P.VaccRadius), 1, 1, 0)) P.VaccRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Minimum radius from case to vaccinate", "%lf", (void*) & (P.VaccMinRadius), 1, 1, 0)) P.VaccMinRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Maximum number of vaccine courses available", "%lf", (void*) & (P.VaccMaxCoursesBase), 1, 1, 0)) P.VaccMaxCoursesBase = 1e20;
	if (!GetInputParameter2(dat, dat2, "Start time of additional vaccine production", "%lf", (void*) & (P.VaccNewCoursesStartTime), 1, 1, 0)) P.VaccNewCoursesStartTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "End time of additional vaccine production", "%lf", (void*) & (P.VaccNewCoursesEndTime), 1, 1, 0)) P.VaccNewCoursesEndTime = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Rate of additional vaccine production (courses per day)", "%lf", (void*) & (P.VaccNewCoursesRate), 1, 1, 0)) P.VaccNewCoursesRate = 0;
	if (!GetInputParameter2(dat, dat2, "Apply mass rather than reactive vaccination", "%i", (void*) & (P.DoMassVacc), 1, 1, 0)) P.DoMassVacc = 0;
	if (!GetInputParameter2(dat, dat2, "Priority age range for mass vaccination", "%i", (void*)P.VaccPriorityGroupAge, 2, 1, 0)) { P.VaccPriorityGroupAge[0] = 1; P.VaccPriorityGroupAge[1] = 0; }
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Vaccinate administrative units rather than rings", "%i", (void*) & (P.VaccByAdminUnit), 1, 1, 0)) P.VaccByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for vaccination", "%i", (void*) & (P.VaccAdminUnitDivisor), 1, 1, 0)) P.VaccAdminUnitDivisor = 1;
		if ((P.VaccAdminUnitDivisor == 0) || (P.VaccByAdminUnit == 0)) P.VaccAdminUnitDivisor = 1;
	}
	else
	{
		P.VaccAdminUnitDivisor = 1; P.VaccByAdminUnit = 0;
	}

	if (!GetInputParameter2(dat, dat2, "Movement restrictions trigger incidence per cell", "%i", (void*) & (P.MoveRestrCellIncThresh), 1, 1, 0)) P.MoveRestrCellIncThresh = 1000000000;
	if (!GetInputParameter2(dat, dat2, "Delay to start movement restrictions", "%lf", (void*) & (P.MoveDelayMean), 1, 1, 0)) P.MoveDelayMean = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of movement restrictions", "%lf", (void*) & (P.MoveRestrDuration), 1, 1, 0)) P.MoveRestrDuration = 7;
	if (!GetInputParameter2(dat, dat2, "Residual movements after restrictions", "%lf", (void*) & (P.MoveRestrEffect), 1, 1, 0)) P.MoveRestrEffect = 0;
	if (!GetInputParameter2(dat, dat2, "Minimum radius of movement restrictions", "%lf", (void*) & (P.MoveRestrRadius), 1, 1, 0)) P.MoveRestrRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Movement restrictions start time", "%lf", (void*) & (P.MoveRestrTimeStartBase), 1, 1, 0)) P.MoveRestrTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Impose blanket movement restrictions", "%i", (void*) & (P.DoBlanketMoveRestr), 1, 1, 0)) P.DoBlanketMoveRestr = 0;
	if (!GetInputParameter2(dat, dat2, "Movement restrictions only once", "%i", (void*) & (P.DoMoveRestrOnceOnly), 1, 1, 0)) P.DoMoveRestrOnceOnly = 0;
	if (P.DoMoveRestrOnceOnly) P.DoMoveRestrOnceOnly = 4;
	if (P.DoAdUnits)
	{
		if (!GetInputParameter2(dat, dat2, "Movement restrictions in administrative units rather than rings", "%i", (void*) & (P.MoveRestrByAdminUnit), 1, 1, 0)) P.MoveRestrByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for movement restrictions", "%i", (void*) & (P.MoveRestrAdminUnitDivisor), 1, 1, 0)) P.MoveRestrAdminUnitDivisor = 1;
		if ((P.MoveRestrAdminUnitDivisor == 0) || (P.MoveRestrByAdminUnit == 0)) P.MoveRestrAdminUnitDivisor = 1;
	}
	else
	{
		P.MoveRestrAdminUnitDivisor = 1; P.MoveRestrByAdminUnit = 0;
	}

	//Intervention delays and durations by admin unit: ggilani 16/03/20
	if (!GetInputParameter2(dat, dat2, "Include intervention delays by admin unit", "%i", (void*) & (P.DoInterventionDelaysByAdUnit), 1, 1, 0)) P.DoInterventionDelaysByAdUnit = 0;
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

		if (!GetInputParameter2(dat, dat2, "Delay to social distancing by admin unit"		, "%lg", (void*)AdunitDelayToSocialDistance	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToSocialDistance	[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Delay to household quarantine by admin unit"		, "%lg", (void*)AdunitDelayToHQuarantine	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToHQuarantine		[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Delay to case isolation by admin unit"			, "%lg", (void*)AdunitDelayToCaseIsolation	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToCaseIsolation	[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Delay to place closure by admin unit"			, "%lg", (void*)AdunitDelayToPlaceClose		, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDelayToPlaceClose		[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of social distancing by admin unit"	, "%lg", (void*)AdunitDurationSocialDistance, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationSocialDistance	[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household quarantine by admin unit"	, "%lg", (void*)AdunitDurationHQuarantine	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationHQuarantine		[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of case isolation by admin unit"		, "%lg", (void*)AdunitDurationCaseIsolation	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationCaseIsolation	[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of place closure by admin unit"		, "%lg", (void*)AdunitDurationPlaceClose	, P.NumAdunits, 1, 0)) for (i = 0; i < P.NumAdunits; i++) AdunitDurationPlaceClose		[i] = 0;
		
		for (i = 0; i < P.NumAdunits; i++)
		{
			AdUnits[i].SocialDistanceDelay		= AdunitDelayToSocialDistance	[i];
			AdUnits[i].SocialDistanceDuration	= AdunitDurationSocialDistance	[i];
			AdUnits[i].HQuarantineDelay			= AdunitDelayToHQuarantine		[i];
			AdUnits[i].HQuarantineDuration		= AdunitDurationHQuarantine		[i];
			AdUnits[i].CaseIsolationDelay		= AdunitDelayToCaseIsolation	[i];
			AdUnits[i].CaseIsolationDuration	= AdunitDurationCaseIsolation	[i];
			AdUnits[i].PlaceCloseDelay			= AdunitDelayToPlaceClose		[i];
			AdUnits[i].PlaceCloseDuration		= AdunitDurationPlaceClose		[i];
		}
	}

	//New code for digital contact tracing - ggilani: 09/03/20
	if (!GetInputParameter2(dat, dat2, "Include digital contact tracing", "%i", (void*) & (P.DoDigitalContactTracing), 1, 1, 0)) P.DoDigitalContactTracing = 0;
	if (P.DoDigitalContactTracing)
	{
		if (!GetInputParameter2(dat, dat2, "Digital contact tracing trigger incidence per cell", "%lf", (void*) & (P.DigitalContactTracing_CellIncThresh), 1, 1, 0)) P.DigitalContactTracing_CellIncThresh = 1000000000;

		if (!GetInputParameter2(dat, dat2, "Proportion of population or households covered by digital contact tracing", "%lf", (void*) & (P.PropPopUsingDigitalContactTracing), 1, 1, 0)) P.PropPopUsingDigitalContactTracing = 1;
		if (!GetInputParameter2(dat, dat2, "Proportion of smartphone users by age", "%lf", (void*)P.ProportionSmartphoneUsersByAge, NUM_AGE_GROUPS, 1, 0))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
			{
				P.ProportionSmartphoneUsersByAge[i] = 1;
			}
		}
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(dat, dat2, "Cluster digital app clusters by household", "%i", (void*) & (P.ClusterDigitalContactUsers), 1, 1, 0)) P.ClusterDigitalContactUsers = 0; // by default, don't cluster by location
		}
		else
		{
			P.ClusterDigitalContactUsers = 0;
		}
		if (!GetInputParameter2(dat, dat2, "Proportion of digital contacts who self-isolate", "%lf", (void*) & (P.ProportionDigitalContactsIsolate), 1, 1, 0)) P.ProportionDigitalContactsIsolate = 0;
		if (!GetInputParameter2(dat, dat2, "Delay before digital contacts self-isolate", "%lf", (void*) & (P.DigitalContactTracingDelay), 1, 1, 0)) P.DigitalContactTracingDelay = 0;
		if (!GetInputParameter2(dat, dat2, "Length of self-isolation for digital contacts", "%lf", (void*) & (P.LengthDigitalContactIsolation), 1, 1, 0)) P.LengthDigitalContactIsolation = 0;
		if (!GetInputParameter2(dat, dat2, "Place scaling factor - digital contact tracing", "%lf", (void*) & (P.ScalingFactorPlaceDigitalContacts), 1, 1, 0)) P.ScalingFactorPlaceDigitalContacts = 1;
		if (!GetInputParameter2(dat, dat2, "Spatial scaling factor - digital contact tracing", "%lf", (void*) & (P.ScalingFactorSpatialDigitalContacts), 1, 1, 0)) P.ScalingFactorSpatialDigitalContacts = 1;
		if (!GetInputParameter2(dat, dat2, "Digital contact tracing start time", "%lf", (void*) & (P.DigitalContactTracingTimeStartBase), 1, 1, 0)) P.DigitalContactTracingTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Duration of digital contact tracing policy", "%lf", (void*) & (P.DigitalContactTracingPolicyDuration), 1, 1, 0)) P.DigitalContactTracingPolicyDuration = 7;
		if (!GetInputParameter2(dat, dat2, "Output digital contact tracing", "%i", (void*) & (P.OutputDigitalContactTracing), 1, 1, 0)) P.OutputDigitalContactTracing = 0;

		//initialise total number of users to 0
		P.NDigitalContactUsers = 0;
		P.NDigitalHouseholdUsers = 0;
	}
	else
	{
		//Set these to 1 so it doesn't interfere with code if we aren't using digital contact tracing.
		P.ScalingFactorPlaceDigitalContacts = 1;
		P.ScalingFactorSpatialDigitalContacts = 1;
	}


	if (!GetInputParameter2(dat, dat2, "Trigger incidence per cell for place closure", "%i", (void*) & (P.PlaceCloseCellIncThresh), 1, 1, 0)) P.PlaceCloseCellIncThresh = 1000000000;
	if(!GetInputParameter2(dat, dat2, "Trigger incidence per cell for end of place closure", "%i", (void*) & (P.PlaceCloseCellIncStopThresh), 1, 1, 0)) P.PlaceCloseCellIncStopThresh = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to start place closure", "%lf", (void*) & (P.PlaceCloseDelayMean), 1, 1, 0)) P.PlaceCloseDelayMean = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of place closure", "%lf", (void*) & (P.PlaceCloseDurationBase), 1, 1, 0)) P.PlaceCloseDurationBase = 7;
	if (!GetInputParameter2(dat, dat2, "Duration of second place closure", "%lf", (void*) & (P.PlaceCloseDuration2), 1, 1, 0)) P.PlaceCloseDuration2 = 7;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Proportion of places remaining open after closure by place type", "%lf", (void*)P.PlaceCloseEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.PlaceCloseEffect[i] = 1;
	}
	if (P.DoHouseholds)
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate after closure", "%lf", (void*)& P.PlaceCloseHouseholdRelContact, 1, 1, 0)) P.PlaceCloseHouseholdRelContact = 1;
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate after closure", "%lf", (void*)& P.PlaceCloseSpatialRelContact, 1, 1, 0))
	{
		P.PlaceCloseSpatialRelContact = 1;
	}

	if (!GetInputParameter2(dat, dat3, "Include holidays", "%i", (void*) & (P.DoHolidays), 1, 1, 0)) P.DoHolidays = 0;
	if (P.DoHolidays)
	{
		if (!GetInputParameter2(dat, dat3, "Proportion of places remaining open during holidays by place type", "%lf", (void*)P.HolidayEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.HolidayEffect[i] = 1;
		if (!GetInputParameter2(dat, dat3, "Number of holidays", "%i", (void*) & (P.NumHolidays), 1, 1, 0)) P.NumHolidays = 0;
		if (P.NumHolidays > DAYS_PER_YEAR) P.NumHolidays = DAYS_PER_YEAR;
		if (P.NumHolidays > 0)
		{
			GetInputParameter(dat, dat3, "Holiday start times", "%lf", (void*)P.HolidayStartTime, P.NumHolidays, 1, 0);
			GetInputParameter(dat, dat3, "Holiday durations", "%lf", (void*)P.HolidayDuration, P.NumHolidays, 1, 0);
		}
	}
	else
		P.NumHolidays = 0;
	if (!GetInputParameter2(dat, dat2, "Minimum radius for place closure", "%lf", (void*) & (P.PlaceCloseRadius), 1, 1, 0)) P.PlaceCloseRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Place closure start time", "%lf", (void*) & (P.PlaceCloseTimeStartBase), 1, 1, 0)) P.PlaceCloseTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Place closure second start time", "%lf", (void*) & (P.PlaceCloseTimeStartBase2), 1, 1, 0)) P.PlaceCloseTimeStartBase2 = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Places close only once", "%i", (void*) & (P.DoPlaceCloseOnceOnly), 1, 1, 0)) P.DoPlaceCloseOnceOnly = 0;
	if (P.DoPlaceCloseOnceOnly) P.DoPlaceCloseOnceOnly = 4;
	if (!GetInputParameter2(dat, dat2, "Place closure incidence threshold", "%i", (void*) & (P.PlaceCloseIncTrig), 1, 1, 0)) P.PlaceCloseIncTrig = 1;
	if (!GetInputParameter2(dat, dat2, "Place closure fractional incidence threshold", "%lf", (void*) & (P.PlaceCloseFracIncTrig), 1, 1, 0)) P.PlaceCloseFracIncTrig = 0;
	if ((P.DoAdUnits) && (P.DoPlaces))
	{
		if (!GetInputParameter2(dat, dat2, "Place closure in administrative units rather than rings", "%i", (void*) & (P.PlaceCloseByAdminUnit), 1, 1, 0)) P.PlaceCloseByAdminUnit = 0;
		if (!GetInputParameter2(dat, dat2, "Administrative unit divisor for place closure", "%i", (void*) & (P.PlaceCloseAdminUnitDivisor), 1, 1, 0)) P.PlaceCloseAdminUnitDivisor = 1;
		if (!GetInputParameter2(dat, dat2, "Place types to close for admin unit closure (0/1 array)", "%i", (void*) & (P.PlaceCloseAdunitPlaceTypes), P.PlaceTypeNum, 1, 0))
			for (i = 0; i < P.PlaceTypeNum; i++) P.PlaceCloseAdunitPlaceTypes[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Cumulative proportion of place members needing to become sick for admin unit closure", "%lf", (void*) & (P.PlaceCloseCasePropThresh), 1, 1, 0)) P.PlaceCloseCasePropThresh = 2;
		if (!GetInputParameter2(dat, dat2, "Proportion of places in admin unit needing to pass threshold for place closure", "%lf", (void*) & (P.PlaceCloseAdunitPropThresh), 1, 1, 0)) P.PlaceCloseAdunitPropThresh = 2;
		if ((P.PlaceCloseAdminUnitDivisor < 1) || (P.PlaceCloseByAdminUnit == 0)) P.PlaceCloseAdminUnitDivisor = 1;
	}
	else
	{
		P.PlaceCloseAdminUnitDivisor = 1; P.PlaceCloseByAdminUnit = 0;
	}

	if (!GetInputParameter2(dat, dat2, "Trigger incidence per cell for social distancing", "%i", (void*) & (P.SocDistCellIncThresh), 1, 1, 0)) P.SocDistCellIncThresh = 1000000000;
	if(!GetInputParameter2(dat, dat2, "Trigger incidence per cell for end of social distancing", "%i", (void*) & (P.SocDistCellIncStopThresh), 1, 1, 0)) P.SocDistCellIncStopThresh = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of social distancing", "%lf", (void*) & (P.SocDistDuration), 1, 1, 0)) P.SocDistDuration = 7;
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Relative place contact rate given social distancing by place type", "%lf", (void*)P.SocDistPlaceEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.SocDistPlaceEffect[i] = 1;
		if (!GetInputParameter2(dat, dat2, "Relative place contact rate given enhanced social distancing by place type", "%lf", (void*)P.ESocDistPlaceEffect, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.ESocDistPlaceEffect[i] = 1;
		if (!GetInputParameter2(dat, dat2, "Relative place contact rate given social distancing by place type after change", "%lf", (void*)P.SocDistPlaceEffect2, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.SocDistPlaceEffect2[i] = P.SocDistPlaceEffect[i];
		if (!GetInputParameter2(dat, dat2, "Relative place contact rate given enhanced social distancing by place type after change", "%lf", (void*)P.ESocDistPlaceEffect2, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.ESocDistPlaceEffect2[i] = P.ESocDistPlaceEffect[i];
	}
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate given social distancing", "%lf", (void*)& P.SocDistHouseholdEffect, 1, 1, 0)) P.SocDistHouseholdEffect = 1;
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate given enhanced social distancing", "%lf", (void*)& P.ESocDistHouseholdEffect, 1, 1, 0)) P.ESocDistHouseholdEffect = 1;
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate given social distancing  after change", "%lf", (void*)&P.SocDistHouseholdEffect2, 1, 1, 0)) P.SocDistHouseholdEffect2 = P.SocDistHouseholdEffect;
		if (!GetInputParameter2(dat, dat2, "Relative houshold contact rate given enhanced social distancing after change", "%lf", (void*)&P.ESocDistHouseholdEffect2, 1, 1, 0)) P.ESocDistHouseholdEffect2 = P.ESocDistHouseholdEffect;
	}
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate given social distancing", "%lf", (void*)& P.SocDistSpatialEffect, 1, 1, 0)) P.SocDistSpatialEffect = 1;
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate given social distancing after change", "%lf", (void*)&P.SocDistSpatialEffect2, 1, 1, 0)) P.SocDistSpatialEffect2 = P.SocDistSpatialEffect;
	if (!GetInputParameter2(dat, dat2, "Minimum radius for social distancing", "%lf", (void*) & (P.SocDistRadius), 1, 1, 0)) P.SocDistRadius = 0;
	if (!GetInputParameter2(dat, dat2, "Social distancing start time", "%lf", (void*) & (P.SocDistTimeStartBase), 1, 1, 0)) P.SocDistTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Delay for change in effectiveness of social distancing", "%lf", (void*)&(P.SocDistChangeDelay), 1, 1, 0)) P.SocDistChangeDelay = USHRT_MAX / P.TimeStepsPerDay;
	if(!GetInputParameter2(dat, dat2, "Proportion compliant with enhanced social distancing by age group", "%lf", (void*)P.ESocProportionCompliant, NUM_AGE_GROUPS, 1, 0))
	{
		if (!GetInputParameter2(dat, dat2, "Proportion compliant with enhanced social distancing", "%lf", (void*)&t, 1, 1, 0)) t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P.ESocProportionCompliant[i] = t;
	}
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate given enhanced social distancing", "%lf", (void*)& P.ESocDistSpatialEffect, 1, 1, 0)) P.ESocDistSpatialEffect = 1;
	if (!GetInputParameter2(dat, dat2, "Relative spatial contact rate given enhanced social distancing after change", "%lf", (void*)&P.ESocDistSpatialEffect2, 1, 1, 0)) P.ESocDistSpatialEffect2 = P.ESocDistSpatialEffect;

	if (!GetInputParameter2(dat, dat2, "Social distancing only once", "%i", (void*) & (P.DoSocDistOnceOnly), 1, 1, 0)) P.DoSocDistOnceOnly = 0;
	if (P.DoSocDistOnceOnly) P.DoSocDistOnceOnly = 4;

	if (!GetInputParameter2(dat, dat2, "Airport closure effectiveness", "%lf", (void*) & (P.AirportCloseEffectiveness), 1, 1, 0)) P.AirportCloseEffectiveness = 0;
	P.AirportCloseEffectiveness = 1.0 - P.AirportCloseEffectiveness;
	if (!GetInputParameter2(dat, dat2, "Airport closure start time", "%lf", (void*) & (P.AirportCloseTimeStartBase), 1, 1, 0)) P.AirportCloseTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Airport closure duration", "%lf", (void*) & (P.AirportCloseDuration), 1, 1, 0)) P.AirportCloseDuration = USHRT_MAX / P.TimeStepsPerDay;

	if (P.DoHouseholds)
	{
		if(!GetInputParameter2(dat, dat2, "Retrigger household quarantine with each new case in quarantine window", "%i", (void*) & (P.DoHQretrigger), 1, 1, 0)) P.DoHQretrigger =0;
		if (!GetInputParameter2(dat, dat2, "Household quarantine start time", "%lf", (void*) & (P.HQuarantineTimeStartBase), 1, 1, 0)) P.HQuarantineTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Delay to start household quarantine", "%lf", (void*) & (P.HQuarantineHouseDelay), 1, 1, 0)) P.HQuarantineHouseDelay = 0;
		if (!GetInputParameter2(dat, dat2, "Length of time households are quarantined", "%lf", (void*) & (P.HQuarantineHouseDuration), 1, 1, 0)) P.HQuarantineHouseDuration = 0;
		if (!GetInputParameter2(dat, dat2, "Duration of household quarantine policy", "%lf", (void*) & (P.HQuarantinePolicyDuration), 1, 1, 0)) P.HQuarantinePolicyDuration = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Relative household contact rate after quarantine", "%lf", (void*) & (P.HQuarantineHouseEffect), 1, 1, 0)) P.HQuarantineHouseEffect = 1;
		if (P.DoPlaces)
		{
			if (!GetInputParameter2(dat, dat2, "Residual place contacts after household quarantine by place type", "%lf", (void*)P.HQuarantinePlaceEffect, P.PlaceTypeNum, 1, 0))
				for (i = 0; i < NUM_PLACE_TYPES; i++) P.HQuarantinePlaceEffect[i] = 1;
		}
		if (!GetInputParameter2(dat, dat2, "Residual spatial contacts after household quarantine", "%lf", (void*) & (P.HQuarantineSpatialEffect), 1, 1, 0)) P.HQuarantineSpatialEffect = 1;
		if (!GetInputParameter2(dat, dat2, "Household level compliance with quarantine", "%lf", (void*) & (P.HQuarantinePropHouseCompliant), 1, 1, 0)) P.HQuarantinePropHouseCompliant = 1;
		if (!GetInputParameter2(dat, dat2, "Individual level compliance with quarantine", "%lf", (void*) & (P.HQuarantinePropIndivCompliant), 1, 1, 0)) P.HQuarantinePropIndivCompliant = 1;
	}
	else
		P.HQuarantineTimeStartBase = 1e10;
	if (!GetInputParameter2(dat, dat2, "Case isolation start time", "%lf", (void*) & (P.CaseIsolationTimeStartBase), 1, 1, 0)) P.CaseIsolationTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
	if (!GetInputParameter2(dat, dat2, "Proportion of detected cases isolated", "%lf", (void*) & (P.CaseIsolationProp), 1, 1, 0)) P.CaseIsolationProp = 0;
	if (!GetInputParameter2(dat, dat2, "Delay to start case isolation", "%lf", (void*) & (P.CaseIsolationDelay), 1, 1, 0)) P.CaseIsolationDelay = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of case isolation", "%lf", (void*) & (P.CaseIsolationDuration), 1, 1, 0)) P.CaseIsolationDuration = 0;
	if (!GetInputParameter2(dat, dat2, "Duration of case isolation policy", "%lf", (void*) & (P.CaseIsolationPolicyDuration), 1, 1, 0)) P.CaseIsolationPolicyDuration = 1e10;
	if (!GetInputParameter2(dat, dat2, "Residual contacts after case isolation", "%lf", (void*) & (P.CaseIsolationEffectiveness), 1, 1, 0)) P.CaseIsolationEffectiveness = 1;
	if (P.DoHouseholds)
	{
		if (!GetInputParameter2(dat, dat2, "Residual household contacts after case isolation", "%lf", (void*) & (P.CaseIsolationHouseEffectiveness), 1, 1, 0))
			P.CaseIsolationHouseEffectiveness = P.CaseIsolationEffectiveness;
	}
	if (P.DoPlaces)
	{
		if (!GetInputParameter2(dat, dat2, "Number of key workers randomly distributed in the population", "%i", (void*) & (P.KeyWorkerPopNum), 1, 1, 0)) P.KeyWorkerPopNum = 0;
		if (!GetInputParameter2(dat, dat2, "Number of key workers in different places by place type", "%i", (void*)P.KeyWorkerPlaceNum, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.KeyWorkerPlaceNum[i] = 0;
		if (!GetInputParameter2(dat, dat2, "Proportion of staff who are key workers per chosen place by place type", "%lf", (void*)P.KeyWorkerPropInKeyPlaces, P.PlaceTypeNum, 1, 0))
			for (i = 0; i < NUM_PLACE_TYPES; i++) P.KeyWorkerPropInKeyPlaces[i] = 1.0;
		if (!GetInputParameter2(dat, dat2, "Trigger incidence per cell for key worker prophylaxis", "%i", (void*) & (P.KeyWorkerProphCellIncThresh), 1, 1, 0)) P.KeyWorkerProphCellIncThresh = 1000000000;
		if (!GetInputParameter2(dat, dat2, "Key worker prophylaxis start time", "%lf", (void*) & (P.KeyWorkerProphTimeStartBase), 1, 1, 0)) P.KeyWorkerProphTimeStartBase = USHRT_MAX / P.TimeStepsPerDay;
		if (!GetInputParameter2(dat, dat2, "Duration of key worker prophylaxis", "%lf", (void*) & (P.KeyWorkerProphDuration), 1, 1, 0)) P.KeyWorkerProphDuration = 0;
		if (!GetInputParameter2(dat, dat2, "Time interval from start of key worker prophylaxis before policy restarted", "%lf", (void*) & (P.KeyWorkerProphRenewalDuration), 1, 1, 0)) P.KeyWorkerProphRenewalDuration = P.KeyWorkerProphDuration;
		if (P.DoHouseholds)
		{
			if (!GetInputParameter2(dat, dat2, "Proportion of key workers whose households are also treated as key workers", "%lf", (void*) & (P.KeyWorkerHouseProp), 1, 1, 0)) P.KeyWorkerHouseProp = 0;
		}
		if (!GetInputParameter2(dat, dat2, "Minimum radius for key worker prophylaxis", "%lf", (void*) & (P.KeyWorkerProphRadius), 1, 1, 0)) P.KeyWorkerProphRadius = 0;
	}
	else
	{
		P.KeyWorkerPopNum = 0;
		P.KeyWorkerProphTimeStartBase = 1e10;
	}
	if (!GetInputParameter2(dat, dat2, "Initial rate of importation of infections", "%lf", (void*) & (P.InfectionImportRate1), 1, 1, 0)) P.InfectionImportRate1 = 0;
	if (!GetInputParameter2(dat, dat2, "Changed rate of importation of infections", "%lf", (void*) & (P.InfectionImportRate2), 1, 1, 0)) P.InfectionImportRate2 = 0;
	if (!GetInputParameter2(dat, dat2, "Time when infection rate changes", "%lf", (void*) & (P.InfectionImportChangeTime), 1, 1, 0)) P.InfectionImportChangeTime = 1e10;
	if (!GetInputParameter2(dat, dat2, "Imports via air travel", "%i", (void*) & (P.DoImportsViaAirports), 1, 1, 0)) P.DoImportsViaAirports = 0;
	if (!GetInputParameter2(dat, dat2, "Length of importation time profile provided", "%i", (void*) & (P.DurImportTimeProfile), 1, 1, 0)) P.DurImportTimeProfile = 0;
	if (P.DurImportTimeProfile > 0)
	{
		if (P.DurImportTimeProfile >= MAX_DUR_IMPORT_PROFILE) ERR_CRITICAL("MAX_DUR_IMPORT_PROFILE too small\n");
		GetInputParameter(dat, dat2, "Daily importation time profile", "%lf", (void*)P.ImportInfectionTimeProfile, P.DurImportTimeProfile, 1, 0);
	}
	//Added this to parameter list so that recording infection events (and the number to record) can easily be turned off and on: ggilani - 10/10/2014
	if (!GetInputParameter2(dat, dat2, "Record infection events", "%i", (void*) & (P.DoRecordInfEvents), 1, 1, 0)) P.DoRecordInfEvents = 0;
	if (P.DoRecordInfEvents)
	{
		if (!GetInputParameter2(dat, dat2, "Max number of infection events to record", "%i", (void*) & (P.MaxInfEvents), 1, 1, 0)) P.MaxInfEvents = 1000;
		if (!GetInputParameter2(dat, dat2, "Record infection events per run", "%i", (void*) & (P.RecordInfEventsPerRun), 1, 1, 0)) P.RecordInfEventsPerRun = 0;
	}
	else
	{
		P.MaxInfEvents = 0;
	}
	//Include a limit to the number of infections to simulate, if this happens before time runs out
	if (!GetInputParameter2(dat, dat2, "Limit number of infections", "%i", (void*) & (P.LimitNumInfections), 1, 1, 0)) P.LimitNumInfections = 0;
	if (P.LimitNumInfections)
	{
		if (!GetInputParameter2(dat, dat2, "Max number of infections", "%i", (void*) & (P.MaxNumInfections), 1, 1, 0)) P.MaxNumInfections = 60000;
	}
	//Add origin-destination matrix parameter
	if (!GetInputParameter2(dat, dat2, "Output origin destination matrix", "%i", (void*) & (P.DoOriginDestinationMatrix), 1, 1, 0)) P.DoOriginDestinationMatrix = 0;

	fclose(dat);

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
	if (P.TreatCellIncThresh < 1) P.TreatCellIncThresh = 1;
	if (P.CaseIsolation_CellIncThresh < 1) P.CaseIsolation_CellIncThresh = 1;
	if (P.DigitalContactTracing_CellIncThresh < 1) P.DigitalContactTracing_CellIncThresh = 1;
	if (P.HHQuar_CellIncThresh < 1) P.HHQuar_CellIncThresh = 1;
	if (P.MoveRestrCellIncThresh < 1) P.MoveRestrCellIncThresh = 1;
	if (P.PlaceCloseCellIncThresh < 1) P.PlaceCloseCellIncThresh = 1;
	if (P.KeyWorkerProphCellIncThresh < 1) P.KeyWorkerProphCellIncThresh = 1;


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
			asin2sqx[i] = asin(sqrt(((double)(i)) / 1000));
			asin2sqx[i] = asin2sqx[i] * asin2sqx[i];
		}
		for (t = 0; t <= 360; t++)
		{
			sinx[(int)t] = sin(PI * t / 180);
			cosx[(int)t] = cos(PI * t / 180);
		}
	}
	fprintf(stderr, "Parameters read\n");
}
void ReadInterventions(char* IntFile)
{
	FILE* dat;
	double r, s, startt, stopt;
	int j, k, au, ni, f, nsr;
	char buf[65536], txt[65536];
	intervention CurInterv;

	fprintf(stderr, "Reading intervention file.\n");
	if (!(dat = fopen(IntFile, "r"))) ERR_CRITICAL("Unable to open intervention file\n");
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
void ReadAirTravel(char* AirTravelFile)
{
	int i, j, k, l;
	float sc, t, t2;
	float* buf;
	double traf;
	char outname[1024];
	FILE* dat;

	fprintf(stderr, "Reading airport data...\nAirports with no connections = ");
	if (!(dat = fopen(AirTravelFile, "r"))) ERR_CRITICAL("Unable to open airport file\n");
	if(fscanf(dat, "%i %i", &P.Nairports, &P.Air_popscale) != 2) {
        ERR_CRITICAL("fscanf failed in void ReadAirTravel\n");
    }
	sc = (float)((double)P.N / (double)P.Air_popscale);
	if (P.Nairports > MAX_AIRPORTS) ERR_CRITICAL("Too many airports\n");
	if (P.Nairports < 2) ERR_CRITICAL("Too few airports\n");
	if (!(buf = (float*)calloc(P.Nairports + 1, sizeof(float)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	if (!(Airports = (airport*)calloc(P.Nairports, sizeof(airport)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	for (i = 0; i < P.Nairports; i++)
	{
		if(fscanf(dat, "%f %f %lf", &(Airports[i].loc_x), &(Airports[i].loc_y), &traf) != 3) {
            ERR_CRITICAL("fscanf failed in void ReadAirTravel\n");
        }
		traf *= (P.AirportTrafficScale * sc);
		if ((Airports[i].loc_x < P.SpatialBoundingBox[0]) || (Airports[i].loc_x > P.SpatialBoundingBox[2])
			|| (Airports[i].loc_y < P.SpatialBoundingBox[1]) || (Airports[i].loc_y > P.SpatialBoundingBox[3]))
		{
			Airports[i].loc_x = Airports[i].loc_y = -1;
			Airports[i].total_traffic = 0;
		}
		else
		{
			//fprintf(stderr,"(%f\t%f) ",Airports[i].loc_x,Airports[i].loc_y);
			Airports[i].loc_x -= (float)P.SpatialBoundingBox[0];
			Airports[i].loc_y -= (float)P.SpatialBoundingBox[1];
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
			if (!(Airports[i].prop_traffic = (float*)calloc(Airports[i].num_connected, sizeof(float)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			if (!(Airports[i].conn_airports = (unsigned short int*) calloc(Airports[i].num_connected, sizeof(unsigned short int)))) ERR_CRITICAL("Unable to allocate airport storage\n");
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
	free(buf);
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
				traf = floor(sqrt(dist2_raw(Airports[i].loc_x, Airports[i].loc_y, Airports[k].loc_x, Airports[k].loc_y)) / OUTPUT_DIST_SCALE);
				l = (int)traf;
				//fprintf(stderr,"%(%i) ",l);
				if (l < MAX_DIST)
					AirTravelDist[l] += Airports[i].total_traffic * Airports[i].prop_traffic[j];
			}
		}
	sprintf(outname, "%s.airdist.xls", OutFile);
	if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open air travel output file\n");
	fprintf(dat, "dist\tfreq\n");
	for (i = 0; i < MAX_DIST; i++)
		fprintf(dat, "%i\t%.10f\n", i, AirTravelDist[i]);
	fclose(dat);
}

void InitModel(int run) // passing run number so we can save run number in the infection event log: ggilani - 15/10/2014
{
	int i, j, k, l, m, j2, tn, nim, stt, stp;
	int nsi[MAX_NUM_SEED_LOCATIONS];

	if (P.OutputBitmap)
	{
#ifdef WIN32_BM
		//if (P.OutputBitmap == 1)
		//{
		//	char buf[200];
		//	sprintf(buf, "%s.ge" DIRECTORY_SEPARATOR "%s.avi", OutFile, OutFile);
		//	avi = CreateAvi(buf, P.BitmapMovieFrame, NULL);
		//}
#endif
		for (unsigned p = 0; p < bmh->imagesize; p++)
		{
			bmi2[p] = bmi3[p] = bmi4[p] = 0;
		}
	}
	ns = 0;
	State.S = P.N;
	State.L = State.I = State.R = State.D = 0;
	State.cumI = State.cumR = State.cumC = State.cumFC = State.cumH = State.cumCT = State.cumCC = State.cumTC = State.cumD = State.cumDC = State.trigDC = State.DCT = State.cumDCT
		= State.cumInf_h = State.cumInf_n = State.cumInf_s = State.cumHQ
		= State.cumAC = State.cumAH = State.cumAA = State.cumACS
		= State.cumAPC = State.cumAPA = State.cumAPCS = 0;
	State.cumT = State.cumUT = State.cumTP = State.cumV = State.sumRad2 = State.maxRad2 = State.cumV_daily = State.cumVG = 0; //added State.cumVG
	State.mvacc_cum = 0;

	if (P.DoSeverity)
	{
		State.Mild		= State.ILI			= State.SARI	= State.Critical	= State.CritRecov		= 0;
		State.cumMild	= State.cumILI		= State.cumSARI = State.cumCritical = State.cumCritRecov	= 0;
		for (int AdminUnit = 0; AdminUnit <= P.NumAdunits; AdminUnit++)
			State.Mild_adunit[AdminUnit] = State.ILI_adunit[AdminUnit] =
			State.SARI_adunit[AdminUnit] = State.Critical_adunit[AdminUnit] = State.CritRecov_adunit[AdminUnit] = 0;
		for (int AdminUnit = 0; AdminUnit <= P.NumAdunits; AdminUnit++)
			State.cumMild_adunit[AdminUnit] = State.cumILI_adunit[AdminUnit] =
			State.cumSARI_adunit[AdminUnit] = State.cumCritical_adunit[AdminUnit] = State.cumCritRecov_adunit[AdminUnit] = State.cumD_adunit[AdminUnit] = 0;
	}
	
	for (i = 0; i < NUM_AGE_GROUPS; i++) State.cumCa[i] = State.cumIa[i] = State.cumDa[i] = 0;
	for (i = 0; i < 2; i++) State.cumC_keyworker[i] = State.cumI_keyworker[i] = State.cumT_keyworker[i] = 0;
	for (i = 0; i < NUM_PLACE_TYPES; i++) State.NumPlacesClosed[i] = 0;
	for (i = 0; i < INFECT_TYPE_MASK; i++) State.cumItype[i] = 0;
	//initialise cumulative case counts per country to zero: ggilani 12/11/14
	for (i = 0; i < MAX_COUNTRIES; i++) State.cumC_country[i] = 0;
	if (P.DoAdUnits)
		for (i = 0; i <= P.NumAdunits; i++)
		{
			State.cumI_adunit[i] = State.cumC_adunit[i] = State.cumD_adunit[i] = State.cumT_adunit[i] = State.cumH_adunit[i] = 
				State.cumDC_adunit[i] = State.cumCT_adunit[i] = State.cumCC_adunit[i] = State.trigDC_adunit[i] = State.DCT_adunit[i] = State.cumDCT_adunit[i] = 0; //added hospitalisation, added detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
			AdUnits[i].place_close_trig = 0;
			AdUnits[i].CaseIsolationTimeStart = AdUnits[i].HQuarantineTimeStart = AdUnits[i].DigitalContactTracingTimeStart = AdUnits[i].SocialDistanceTimeStart = AdUnits[i].PlaceCloseTimeStart = 1e10;
			if (P.DoUpdateCaseDetection)
				AdUnits[i].caseDetectRate = AdUnits[i].caseDetectInit;
			AdUnits[i].ndct = 0; //noone being digitally contact traced at beginning of run
			AdUnits[i].ndct_queue = 0;

		}
	for (j = 0; j < MAX_NUM_THREADS; j++)
	{
		StateT[j].L = StateT[j].I = StateT[j].R = StateT[j].D = 0;
		StateT[j].cumI = StateT[j].cumR = StateT[j].cumC = StateT[j].cumFC = StateT[j].cumH = StateT[j].cumCT = StateT[j].cumCC = StateT[j].DCT = StateT[j].cumDCT = StateT[j].cumTC = StateT[j].cumD = StateT[j].cumDC
			= StateT[j].cumInf_h = StateT[j].cumInf_n = StateT[j].cumInf_s = StateT[j].cumHQ = StateT[j].cumAC = StateT[j].cumACS
			= StateT[j].cumAH = StateT[j].cumAA = StateT[j].cumAPC = StateT[j].cumAPA = StateT[j].cumAPCS = 0;
		StateT[j].cumT = StateT[j].cumUT = StateT[j].cumTP = StateT[j].cumV = StateT[j].sumRad2 = StateT[j].maxRad2 = StateT[j].cumV_daily = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++) StateT[j].cumCa[i] = StateT[j].cumIa[i] = StateT[j].cumDa[i] = 0;
		for (i = 0; i < 2; i++) StateT[j].cumC_keyworker[i] = StateT[j].cumI_keyworker[i] = StateT[j].cumT_keyworker[i] = 0;
		for (i = 0; i < NUM_PLACE_TYPES; i++) StateT[j].NumPlacesClosed[i] = 0;
		for (i = 0; i < INFECT_TYPE_MASK; i++) StateT[j].cumItype[i] = 0;
		//initialise cumulative case counts per country per thread to zero: ggilani 12/11/14
		for (i = 0; i < MAX_COUNTRIES; i++) StateT[j].cumC_country[i] = 0;
		if (P.DoAdUnits)
			for (i = 0; i <= P.NumAdunits; i++)
				StateT[j].cumI_adunit[i] = StateT[j].cumC_adunit[i] = StateT[j].cumD_adunit[i] = StateT[j].cumT_adunit[i] = StateT[j].cumH_adunit[i] = StateT[j].cumDC_adunit[i] =
				StateT[j].cumCT_adunit[i] = StateT[j].cumCC_adunit[i] = StateT[j].nct_queue[i] = StateT[j].ndct_queue[i] = StateT[j].cumDCT_adunit[i] = StateT[j].DCT_adunit[i] = 0; //added hospitalisation, detected cases, contact tracing per adunit, cases who are contacts: ggilani 03/02/15, 15/06/17
	
		if (P.DoSeverity)
		{
			StateT[j].Mild		= StateT[j].ILI		= StateT[j].SARI	= StateT[j].Critical	= StateT[j].CritRecov		= 0;
			StateT[j].cumMild	= StateT[j].cumILI	= StateT[j].cumSARI = StateT[j].cumCritical = StateT[j].cumCritRecov	= 0;
			for (int AdminUnit = 0; AdminUnit <= P.NumAdunits; AdminUnit++)
				StateT[j].Mild_adunit[AdminUnit] = StateT[j].ILI_adunit[AdminUnit] =
				StateT[j].SARI_adunit[AdminUnit] = StateT[j].Critical_adunit[AdminUnit] = StateT[j].CritRecov_adunit[AdminUnit] = 0;
			for (int AdminUnit = 0; AdminUnit <= P.NumAdunits; AdminUnit++)
				StateT[j].cumMild_adunit[AdminUnit] = StateT[j].cumILI_adunit[AdminUnit] =
				StateT[j].cumSARI_adunit[AdminUnit] = StateT[j].cumCritical_adunit[AdminUnit] = StateT[j].cumCritRecov_adunit[AdminUnit] = StateT[j].cumD_adunit[AdminUnit] = 0;
		}

	}
	nim = 0;


#pragma omp parallel for private(i,j,k,l,m,tn,stt,stp) reduction(+:nim) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)
	{
		stp = P.NC / P.NumThreads + 1;
		stt = (tn + 1) * stp;
		if (stt > P.NC) stt = P.NC;
		for (i = tn * stp; i < stt; i++)
		{
			if ((Cells[i].tot_treat != 0) || (Cells[i].tot_vacc != 0) || (Cells[i].S != Cells[i].n) || (Cells[i].D > 0) || (Cells[i].R > 0))
			{
				Cells[i].S = Cells[i].n;
				Cells[i].L = Cells[i].I = Cells[i].R = Cells[i].cumTC = Cells[i].D = 0;
				Cells[i].infected = Cells[i].latent = Cells[i].susceptible + Cells[i].S;
				Cells[i].tot_treat = Cells[i].tot_vacc = 0;

				for (l = 0; l < MAX_INTERVENTION_TYPES; l++) Cells[i].CurInterv[l] = -1;
				for (j = 0; j < Cells[i].n; j++)
				{
					k = Cells[i].members[j];
					Cells[i].susceptible[j] = k; //added this in here instead
					if (P.DoAirports) Hosts[k].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
					Hosts[k].quar_start_time = Hosts[k].isolation_start_time = Hosts[k].absent_start_time = Hosts[k].dct_start_time = USHRT_MAX - 1;
					Hosts[k].absent_stop_time = 0;
					Hosts[k].quar_comply = 2;
					Hosts[k].susc = 1.0;
					Hosts[k].to_die = 0;
					Hosts[k].Travelling = 0;
					Hosts[k].detected = 0; //set detected to zero initially: ggilani - 19/02/15
					Hosts[k].inf = InfStat_Susceptible;
					Hosts[k].listpos = j;
					Hosts[k].treat_stop_time = Hosts[k].num_treats = 0;
					Hosts[k].vacc_start_time = Hosts[k].treat_start_time = USHRT_MAX - 1;
					Hosts[k].latent_time = Hosts[k].recovery_time = 0; //also set hospitalisation time to zero: ggilani 28/10/2014
					Hosts[k].infector = -1;
					Hosts[k].infect_type = 0;
					Hosts[k].ncontacts = 0; //reset to 0 for digital contact tracing

					if (P.DoSeverity)
					{
						Hosts[k].SARI_time						= USHRT_MAX - 1; //// think better to set to initialize to maximum possible value, but keep this way for now. 
						Hosts[k].Critical_time					= USHRT_MAX - 1; 
						Hosts[k].RecoveringFromCritical_time	= USHRT_MAX - 1;  
						Hosts[k].Severity_Current				= Severity_Asymptomatic;
						Hosts[k].Severity_Final					= Severity_Asymptomatic;
					}
				}
				// Next loop needs to count down for DoImmune host list reordering to work
				for (j = Cells[i].n - 1; j >= 0; j--)
				{
					k = Cells[i].members[j];
					if (P.DoWholeHouseholdImmunity)
					{
						if (P.InitialImmunity[0] != 0)
						{
							if (Households[Hosts[k].hh].FirstPerson == k)
							{
								for (m = 0; m < Households[Hosts[k].hh].nh; m++)
									Hosts[k + m].inf = InfStat_Susceptible;
								if ((P.InitialImmunity[0] == 1) || (ranf_mt(tn) < P.InitialImmunity[0]))
									for (m = Households[Hosts[k].hh].nh - 1; m >= 0; m--)
										DoImmune(k + m);
							}
						}
					}
					Hosts[k].inf = InfStat_Susceptible;
					if (P.DoInitEquilib)
					{
						if (P.DoSIS)
						{
							if (P.SuscReductionFactorPerInfection > 0)
								Hosts[k].susc = (float)exp(-P.MeanAnnualDeathRate * HOST_AGE_YEAR(k));
							else
								Hosts[k].susc = (float)( (ranf_mt(tn) > exp(-P.MeanAnnualDeathRate * ((double)HOST_AGE_YEAR(k)))) ? 0.0 : 1.0);
						}
						else if (ranf_mt(tn) > exp(-P.MeanAnnualDeathRate * ((double)HOST_AGE_YEAR(k))))
							DoImmune(k);
					}
					else
					{
						m = HOST_AGE_GROUP(k);
						if (P.DoSIS)
						{
							if (P.SuscReductionFactorPerInfection > 0)
								Hosts[k].susc = (float)(1.0 - P.InitialImmunity[m]);
							else
								nim += (int)(Hosts[k].susc = (ranf_mt(tn) < P.InitialImmunity[m]) ? 0.0f : 1.0f);
						}
						else
						{
							if ((P.InitialImmunity[m] == 1) || ((P.InitialImmunity[m] > 0) && (ranf_mt(tn) < P.InitialImmunity[m]))) { DoImmune(k); nim += 1; }
						}
					}
				}
			}
		}
	}
	fprintf(stderr, "Finished cell init - %i people assigned as immune.\n", nim);


#pragma omp parallel for private(i,j,k,j2,l) schedule(static,500)
	for (l = 0; l < P.NMCP; l++)
	{
		i = (int)(McellLookup[l] - Mcells);
		Mcells[i].vacc_start_time = Mcells[i].treat_start_time = USHRT_MAX - 1;
		Mcells[i].treat_end_time = 0;
		Mcells[i].treat_trig = Mcells[i].vacc_trig = Mcells[i].vacc = Mcells[i].treat = 0;
		Mcells[i].place_trig = Mcells[i].move_trig = Mcells[i].socdist_trig = Mcells[i].keyworkerproph_trig =
			Mcells[i].placeclose = Mcells[i].moverest = Mcells[i].socdist = Mcells[i].keyworkerproph = 0;
		Mcells[i].move_start_time = USHRT_MAX - 1;
		Mcells[i].place_end_time = Mcells[i].move_end_time =
			Mcells[i].socdist_end_time = Mcells[i].keyworkerproph_end_time = 0;
		if (P.DoPlaces)
			for (j = 0; j < P.PlaceTypeNum; j++)
				for (k = 0; k < Mcells[i].np[j]; k++)
				{
					j2 = Mcells[i].places[j][k];
					Places[j][j2].treat = Places[j][j2].control_trig = 0;
					Places[j][j2].treat_end_time = Places[j][j2].close_end_time = 0;
					Places[j][j2].close_start_time = USHRT_MAX - 1;
#ifdef ABSENTEEISM_PLACE_CLOSURE
					Places[j][j2].AbsentLastUpdateTime = 0;
					for (int i2 = 0; i2 < MAX_ABSENT_TIME; i2++) Places[j][j2].Absent[i2] = 0;
#endif
				}
	}

	P.SocDistDurationC = P.SocDistDuration;
	P.SocDistHouseholdEffectC = P.SocDistHouseholdEffect;
	P.SocDistSpatialEffectC = P.SocDistSpatialEffect;
	P.ESocDistHouseholdEffectC = P.ESocDistHouseholdEffect;
	P.ESocDistSpatialEffectC = P.ESocDistSpatialEffect;
	for (i = 0; i < P.PlaceTypeNum; i++)
	{
		P.SocDistPlaceEffectC[i] = P.SocDistPlaceEffect[i];
		P.ESocDistPlaceEffectC[i] = P.ESocDistPlaceEffect[i];
	}

	for (i = 0; i < MAX_NUM_THREADS; i++)
	{
		for (j = 0; j < MAX_NUM_THREADS; j++)
			StateT[i].n_queue[j] = 0;
		for (j = 0; j < P.PlaceTypeNum; j++)
			StateT[i].np_queue[j] = 0;
	}
	if (DoInitUpdateProbs)
	{
		UpdateProbs(0);
		DoInitUpdateProbs = 0;
	}
	//initialise event log to zero at the beginning of every run: ggilani - 10/10/2014. UPDATE: 15/10/14 - we are now going to store all events from all realisations in one file
	if ((P.DoRecordInfEvents) && (P.RecordInfEventsPerRun))
	{
		*nEvents = 0;
		for (i = 0; i < P.MaxInfEvents; i++)
		{
			InfEventLog[i].t = InfEventLog[i].infectee_x = InfEventLog[i].infectee_y = InfEventLog[i].t_infector = 0.0;
			InfEventLog[i].infectee_ind = InfEventLog[i].infector_ind = 0;
			InfEventLog[i].infectee_adunit = InfEventLog[i].listpos = InfEventLog[i].infectee_cell = InfEventLog[i].infector_cell = InfEventLog[i].thread = 0;
		}
	}

	for (i = 0; i < P.NumSeedLocations; i++) nsi[i] = (int) (((double) P.NumInitialInfections[i]) * P.InitialInfectionsAdminUnitWeight[i]+0.5);
	SeedInfection(0, nsi, 0, run);
	P.ControlPropCasesId = P.PreAlertControlPropCasesId;
	P.TreatTimeStart = 1e10;

	//  Mass Vacc now starts after outbreak alert trigger
	//	if(P.DoMassVacc)
	//		P.VaccTimeStart=P.VaccTimeStartBase;
	//	else
	P.VaccTimeStart = 1e10;
	P.MoveRestrTimeStart = 1e10;
	P.PlaceCloseTimeStart = 1e10;
	P.PlaceCloseTimeStart2 = 1e10;
	P.SocDistTimeStart = 1e10;
	P.AirportCloseTimeStart = 1e10;
	P.CaseIsolationTimeStart = 1e10;
	P.HQuarantineTimeStart = 1e10;
	P.KeyWorkerProphTimeStart = 1e10;
	P.TreatMaxCourses = P.TreatMaxCoursesBase;
	P.VaccMaxCourses = P.VaccMaxCoursesBase;
	P.PlaceCloseDuration = P.PlaceCloseDurationBase;
	P.ResetSeedsFlag = 0; //added this to allow resetting seeds part way through run: ggilani 27/11/2019

	fprintf(stderr, "Finished InitModel.\n");
}

void SeedInfection(double t, int* nsi, int rf, int run) //adding run number to pass it to event log
{
	/* *nsi is an array of the number of seeding infections by location (I think). During runtime, usually just a single int (given by a poisson distribution)*/
	/*rf set to 0 when initializing model, otherwise set to 1 during runtime. */

	int i /*seed location index*/;
	int j /*microcell number*/;
	int k, l /*k,l are grid coords at first, then l changed to be person within Microcell j, then k changed to be index of new infection*/;
	int m /*guard against too many infections and infinite loop*/;
	int f /*range = {0, 1000}*/, f2 /*not used, think extraneous*/;
	int n /*number of seed locations?*/;

	f2 = ((t >= 0) && (P.DoImportsViaAirports));
	n = ((rf == 0) ? P.NumSeedLocations : 1);
	for (i = 0; i < n; i++)
	{
		if ((!P.DoRandomInitialInfectionLoc) || ((P.DoAllInitialInfectioninSameLoc) && (rf))) //// either non-random locations, doing all initial infections in same location, and not initializing.
		{
			k = (int)(P.LocationInitialInfection[i][0] / P.mcwidth);
			l = (int)(P.LocationInitialInfection[i][1] / P.mcheight);
			j = k * P.nmch + l;
			m = 0;
			for (k = 0; (k < nsi[i]) && (m < 10000); k++)
			{
				l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))]; //// randomly choose member of microcell j. Name this member l
				if (Hosts[l].inf == InfStat_Susceptible) //// If Host l is uninfected. 
				{
					if (CalcPersonSusc(l, 0, 0, 0) > 0)
					{
						//only reset the initial location if rf==0, i.e. when initial seeds are being set, not when imported cases are being set
						if (rf == 0)
						{
							P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc_x;
							P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc_y;
						}
						Hosts[l].infector = -2;
						Hosts[l].infect_type = INFECT_TYPE_MASK - 1;
						DoInfect(l, t, 0, run); ///// guessing this updates a number of things about person l at time t in thread 0 for this run. 
						m = 0;
					}
				}
				else { k--; m++; } //// think k-- means if person l chosen is already infected, go again. The m < 10000 is a guard against a) too many infections; b) an infinite loop if no more uninfected people left. 
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
					l = (int)(ranf() * ((double)P.N));
					j = Hosts[l].mcell;
					//fprintf(stderr,"%i ",AdUnits[Mcells[j].adunit].id);
				} while ((Mcells[j].n < nsi[i]) || (Mcells[j].n > P.MaxPopDensForInitialInfection)
					|| (Mcells[j].n < P.MinPopDensForInitialInfection) 
					|| ((P.InitialInfectionsAdminUnit[i] > 0) && ((AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[i] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));
				for (k = 0; (k < nsi[i]) && (m < 10000); k++)
				{
					l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
					if (Hosts[l].inf == InfStat_Susceptible)
					{
						if (CalcPersonSusc(l, 0, 0, 0) > 0)
						{
							P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc_x;
							P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc_y;
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
			for (k = 0; (k < nsi[i]) && (m < 10000); k++)
			{
				do
				{
					l = (int)(ranf() * ((double)P.N));
					j = Hosts[l].mcell;
					//fprintf(stderr,"@@ %i %i ",AdUnits[Mcells[j].adunit].id, (int)(AdUnits[Mcells[j].adunit].id / P.CountryDivisor));
				} while ((Mcells[j].n == 0) || (Mcells[j].n > P.MaxPopDensForInitialInfection)
					|| (Mcells[j].n < P.MinPopDensForInitialInfection) 
					|| ((P.InitialInfectionsAdminUnit[i] > 0) && ((AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor != (P.InitialInfectionsAdminUnit[i] % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor)));
				l = Mcells[j].members[(int)(ranf() * ((double)Mcells[j].n))];
				if (Hosts[l].inf == InfStat_Susceptible)
				{
					if (CalcPersonSusc(l, 0, 0, 0) > 0)
					{
						P.LocationInitialInfection[i][0] = Households[Hosts[l].hh].loc_x;
						P.LocationInitialInfection[i][1] = Households[Hosts[l].hh].loc_y;
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


int RunModel(int run) //added run number as parameter
{
	int j, k, l, fs, fs2, nu, ni,nsi[MAX_NUM_SEED_LOCATIONS] /*Denotes either Num imported Infections given rate ir, or number false positive "infections"*/;
	double ir; // infection import rate?; 
	double t, cI, lcI, t2;
	unsigned short int ts;
	int continueEvents = 1;

	InterruptRun = 0;
	lcI = 1;
	if (P.DoLoadSnapshot)
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
	nu = 0;

	for (ns = 1; ((ns < P.NumSamples) && (!InterruptRun)); ns++) //&&(continueEvents) <-removed this
	{
		RecordSample(t, ns - 1);

		//Only run to a certain number of infections: ggilani 28/10/14
		if (P.LimitNumInfections) continueEvents = (State.cumI < P.MaxNumInfections);

		fprintf(stderr, "\r    t=%lg   %i    %i|%i    %i     %i [=%i]  %i (%lg %lg %lg)   %lg    ", t, 
			State.S, State.L, State.I, State.R, State.D, State.S+State.L+State.I+State.R+State.D, State.cumD, State.cumT, State.cumV, State.cumVG, sqrt(State.maxRad2) / 1000); //added State.cumVG

		for (j = 0; ((j < P.UpdatesPerSample) && (!InterruptRun) && (continueEvents)); j++)
		{
			ts = (unsigned short int) (P.TimeStepsPerDay * t);

			//if we are to reset random numbers after an intervention event, specific time
			if (P.ResetSeedsPostIntervention)
				if ((P.ResetSeedsFlag == 0) && (ts >= (P.TimeToResetSeeds * P.TimeStepsPerDay)))
				{
					P.newseed1 = (int)(ranf() * 1e8);
					P.newseed2 = (int)(ranf() * 1e8);
					setall(P.newseed1, P.newseed2);
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
#ifdef IMPORT_POP_SIZE
				ir *= ((double)P.N) / IMPORT_POP_SIZE; //// change infection import rate to account for differences between imported population size IMPORT_POP_SIZE and runtime population size P.N. 
#endif
				if (ir > 0) //// if infection import rate > 0, seed some infections
				{
					for(k=ni=0;k<P.NumSeedLocations;k++) ni+=(nsi[k] = (int)ignpoi(P.TimeStep * ir*P.InitialInfectionsAdminUnitWeight[k])); //// sample number imported infections from Poisson distribution. 
					if (ni > 0)		SeedInfection(t, nsi, 1, run);
				}
				if (P.FalsePositivePerCapitaIncidence > 0)
				{
					ni = (int)ignpoi(P.TimeStep * P.FalsePositivePerCapitaIncidence * ((double)P.N));
					if (ni > 0)
					{
						for (k = 0; k < ni; k++)
						{
							do
							{
								l = (int)(((double)P.N) * ranf()); //// choose person l randomly from entire population. (but change l if while condition not satisfied?)
							} while ((abs(Hosts[l].inf) == InfStat_Dead) || (ranf() > P.FalsePositiveAgeRate[HOST_AGE_GROUP(l)]));
							DoFalseCase(l, t, ts, 0);
						}
					}
				}
				InfectSweep(t, run);  // loops over all infectious people and decides which susceptible people to infect (at household, place and spatial level), and adds them to queue. Then changes each person's various characteristics using DoInfect function.  adding run number as a parameter to infect sweep so we can track run number: ggilani - 15/10/14
				//// IncubRecoverySweep loops over all infecteds (either latent or infectious). If t is the right time, latent people moved to being infected, and infectious people moved to being clinical cases. Possibly also add them to recoveries or deaths. Add them to hospitalisation & hospitalisation discharge queues. 
				if (!P.DoSI) IncubRecoverySweep(t, run);
				// If doing new contact tracing, update numbers of people under contact tracing after each time step

				if (P.DoDigitalContactTracing)
				{
					DigitalContactTracingSweep(t);
				}

				nu++;
				fs2 = ((P.DoDeath) || (P.DoSIS) || (State.L + State.I > 0) || (ir > 0) || (P.FalsePositivePerCapitaIncidence > 0));
				///// TreatSweep loops over microcells to decide which cells are treated (either with treatment, vaccine, social distancing, movement restrictions etc.). Calls DoVacc, DoPlaceClose, DoProphNoDelay etc. to change (threaded) State variables
				if (!TreatSweep(t))
				{
					if ((!fs2) && (State.L + State.I == 0) && (P.FalsePositivePerCapitaIncidence == 0)) { if ((ir == 0) && (((int)t) > P.DurImportTimeProfile)) fs = 0; }
				}
				if (P.DoAirports) TravelReturnSweep(t);
			}
			t += P.TimeStep;
			if (P.DoDeath) P.ts_age++;
			if ((P.DoSaveSnapshot) && (t <= P.SnapshotSaveTime) && (t + P.TimeStep > P.SnapshotSaveTime)) SaveSnapshot();
			if (t > P.TreatNewCoursesStartTime) P.TreatMaxCourses += P.TimeStep * P.TreatNewCoursesRate;
			if ((t > P.VaccNewCoursesStartTime) && (t < P.VaccNewCoursesEndTime)) P.VaccMaxCourses += P.TimeStep * P.VaccNewCoursesRate;
			cI = ((double)(State.S)) / ((double)P.N);
			if (((lcI - cI) > 0.2) && (!P.DoSIS))
			{
				lcI = cI;
				UpdateProbs(0);
				DoInitUpdateProbs = 1;
			}
		}
	}
	RecordSample(t, P.NumSamples - 1);
	fprintf(stderr, "\nEnd of run\n");
	t2 = t + P.SampleTime;
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
/*		fprintf(stderr,"Checking consistency of final state...\n");
        int i, i2, k2;
		for(i=j=k=ni=fs2=i2=0;i<P.N;i++)
			{
			if(i%1000==0) fprintf(stderr,"\r*** %i              ",i);
			if(Hosts[i].inf==0) j++;
			if((Hosts[i].pcell<P.NC)&&(Hosts[i].pcell>=0))
				{
				if(Cells[Hosts[i].pcell].susceptible[Hosts[i].listpos]!=i)
					{
					k++;
					for(l=fs=0;(l<Cells[Hosts[i].pcell].n)&&(!fs);l++)
						fs=(Cells[Hosts[i].pcell].susceptible[l]==i);
					if(!fs) ni++;
					}
				else
					{
					if((Hosts[i].listpos>Cells[Hosts[i].pcell].S-1)&&(Hosts[i].inf== InfStat_Susceptible)) i2++;
					if ((Hosts[i].listpos < Cells[Hosts[i].pcell].S+ Cells[Hosts[i].pcell].L+ Cells[Hosts[i].pcell].I-1) && (abs(Hosts[i].inf) == InfStat_Recovered)) i2++;
					}
				if ((Cells[Hosts[i].pcell].S + Cells[Hosts[i].pcell].L + Cells[Hosts[i].pcell].I + Cells[Hosts[i].pcell].R + Cells[Hosts[i].pcell].D) != Cells[Hosts[i].pcell].n)
				{
					k2++;
				}
				}
			else
				fs2++;
			}
		fprintf(stderr,"\n*** susceptibles=%i\nincorrect listpos=%i\nhosts not found in cell list=%i\nincorrect cell refs=%i\nincorrect positioning in cell susc list=%i\nwrong cell totals=%i\n",j,k,ni,fs2,i2,k2);
*/		
	if(!InterruptRun) RecordInfTypes();
	return (InterruptRun);
}

void SaveDistribs(void)
{
	int i, j, k;
	FILE* dat;
	char outname[1024];
	double s;

	if (P.DoPlaces)
	{
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
			{
				for (i = 0; i < P.Nplace[j]; i++)
					Places[j][i].n = 0;
				for (i = 0; i < P.N; i++)
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
		for (i = 0; i < P.N; i++)
			for (j = 0; j < P.PlaceTypeNum; j++)
				if ((j != HOTEL_PLACE_TYPE) && (Hosts[i].PlaceLinks[j] >= 0))
				{
					if (Hosts[i].PlaceLinks[j] >= P.Nplace[j])
						fprintf(stderr, "*%i %i: %i ", i, j, Hosts[i].PlaceLinks[j]);
					else if ((!P.DoOutputPlaceDistForOneAdunit) ||
						((AdUnits[Mcells[Hosts[i].mcell].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor == (P.OutputPlaceDistAdunit % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor))
					{
						k = Hosts[i].PlaceLinks[j];
						s = sqrt(dist2_raw(Households[Hosts[i].hh].loc_x, Households[Hosts[i].hh].loc_y, Places[j][k].loc_x, Places[j][k].loc_y)) / OUTPUT_DIST_SCALE;
						k = (int)s;
						if (k < MAX_DIST) PlaceDistDistrib[j][k]++;
					}
				}
		for (j = 0; j < P.PlaceTypeNum; j++)
			for (i = 0; i < MAX_PLACE_SIZE; i++)
				PlaceSizeDistrib[j][i] = 0;
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
				for (i = 0; i < P.Nplace[j]; i++)
					if (Places[j][i].n < MAX_PLACE_SIZE)
						PlaceSizeDistrib[j][Places[j][i].n]++;
		sprintf(outname, "%s.placedist.xls", OutFile);
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "dist");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
				fprintf(dat, "\tfreq_p%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < MAX_DIST; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != HOTEL_PLACE_TYPE)
					fprintf(dat, "\t%i", PlaceDistDistrib[j][i]);
			fprintf(dat, "\n");
		}
		fclose(dat);
		sprintf(outname, "%s.placesize.xls", OutFile);
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "size");
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
				fprintf(dat, "\tfreq_p%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < MAX_PLACE_SIZE; i++)
		{
			fprintf(dat, "%i", i);
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (j != HOTEL_PLACE_TYPE)
					fprintf(dat, "\t%i", PlaceSizeDistrib[j][i]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}
}
void SaveOriginDestMatrix(void)
{
	/** function: SaveOriginDestMatrix
	 *
	 * purpose: to save the calculated origin destination matrix to file
	 * parameters: none
	 * returns: none
	 *
	 * author: ggilani, 13/02/15
	 */
	int i, j;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.origdestmat.xls", OutFile);
	if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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

void SaveResults(void)
{
	int i, j;
	FILE* dat;
	char outname[1024];

	if (P.OutputNonSeverity)
	{
		sprintf(outname, "%s.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tL\tI\tR\tD\tincI\tincR\tincFC\tincC\tincDC\tincTC\tincH\tincCT\tincCC\tcumT\tcumTP\tcumV\tcumVG\tExtinct\trmsRad\tmaxRad\n");//\t\t%.10f\t%.10f\t%.10f\n",P.R0household,P.R0places,P.R0spatial);
		for(i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].L, TimeSeries[i].I,
				TimeSeries[i].R, TimeSeries[i].D, TimeSeries[i].incI,
				TimeSeries[i].incR, TimeSeries[i].incFC, TimeSeries[i].incC, TimeSeries[i].incDC, TimeSeries[i].incTC, TimeSeries[i].incH, TimeSeries[i].incCT, TimeSeries[i].incCC, 
				TimeSeries[i].cumT, TimeSeries[i].cumTP, TimeSeries[i].cumV, TimeSeries[i].cumVG, TimeSeries[i].extinct, TimeSeries[i].rmsRad, TimeSeries[i].maxRad);
		}
		fclose(dat);
	}

	if ((P.DoAdUnits) && (P.DoAdunitOutput))
	{
		sprintf(outname, "%s.adunit.xls", OutFile); 
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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
		sprintf(outname, "%s.digitalcontacttracing.xls", OutFile); //modifying to csv file
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\t");
		for (i = 0; i < P.NumAdunits; i++)
		{
			fprintf(dat, "incDCT%i\t", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor); //"\tT%i" //printing headers for inc per admin unit
		}
		for (i = 0; i < P.NumAdunits; i++)
		{
			fprintf(dat, "DCT%i\t", (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor); //"\tT%i" //printing headers for prevalence of digital contact tracing per admin unit
		}
		//print actual output
		for (j = 0; j < P.NumAdunits; j++)
		{
			fprintf(dat, "%.10f,", TimeSeries[i].incDCT_adunit[j]); //"\t%lf" //added contact tracing
		}
		for (j = 0; j < P.NumAdunits; j++)
		{
			fprintf(dat, "%.10f,", TimeSeries[i].DCT_adunit[j]); //"\t%lf" //added contact tracing
		}
		fprintf(dat, "\n");
		fclose(dat);

	}

	if(P.KeyWorkerProphTimeStartBase < P.SampleTime)
		{
		sprintf(outname, "%s.keyworker.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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
		sprintf(outname, "%s.tree.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		for(i = 0; i < P.N; i++)
			if(Hosts[i].infect_type % INFECT_TYPE_MASK > 0)
				fprintf(dat, "%i\t%i\t%i\t%i\n", i, Hosts[i].infector, Hosts[i].infect_type % INFECT_TYPE_MASK, (int)HOST_AGE_YEAR(i));
		fclose(dat);
		}
#if defined(WIN32_BM) || defined(IMAGE_MAGICK)
	static int dm[13] ={0,31,28,31,30,31,30,31,31,30,31,30,31};
	int d, m, y, dml, f;
#ifdef WIN32_BM
	//if(P.OutputBitmap == 1) CloseAvi(avi);
	//if((TimeSeries[P.NumSamples - 1].extinct) && (P.OutputOnlyNonExtinct))
	//	{
	//	sprintf(outname, "%s.ge" DIRECTORY_SEPARATOR "%s.avi", OutFile, OutFile);
	//	DeleteFile(outname);
	//	}
#endif
	if(P.OutputBitmap >= 1)
		{
		// Generate Google Earth .kml file
		sprintf(outname, "%s.ge" DIRECTORY_SEPARATOR "%s.ge.kml", OutFile, OutFile); //sprintf(outname,"%s.ge" DIRECTORY_SEPARATOR "%s.kml",OutFileBase,OutFile);
		if(!(dat = fopen(outname, "wb")))
			{
			ERR_CRITICAL("Unable to open output kml file\n");
			}
		fprintf(dat, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.2\">\n<Document>\n");
		fprintf(dat, "<name>%s</name>\n", OutFile);
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
				sprintf(outname, "%s.ge" DIRECTORY_SEPARATOR "%s.%i.png", OutFile, OutFile, i + 1);
				fprintf(dat, "<Icon>\n<href>%s</href>\n</Icon>\n", outname);
				fprintf(dat, "<LatLonBox>\n<north>%.10f</north>\n<south>%.10f</south>\n<east>%.10f</east>\n<west>%.10f</west>\n</LatLonBox>\n",
					P.SpatialBoundingBox[3], P.SpatialBoundingBox[1], P.SpatialBoundingBox[2], P.SpatialBoundingBox[0]);
				fprintf(dat, "</GroundOverlay>\n");
			}
		fprintf(dat, "</Document>\n</kml>\n");
		fclose(dat);
		}
#endif


	if(P.DoSeverity)
	{
		printf("\nSaveResults: writing Severity");
		//// output separate severity file (can integrate with main if need be)
		sprintf(outname, "%s.severity.xls", OutFile);
		printf("\nFilename: %s", outname);

		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open severity output file\n");
		fprintf(dat, "t\tS\tI\tR\tincI\tMild\tILI\tSARI\tCritical\tCritRecov\tincMild\tincILI\tincSARI\tincCritical\tincCritRecov\tincDeath\tcumMild\tcumILI\tcumSARI\tcumCritical\tcumCritRecov\tcumDeath\n");//\t\t%.10f\t%.10f\t%.10f\n",P.R0household,P.R0places,P.R0spatial);
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				TimeSeries[i].t, TimeSeries[i].S, TimeSeries[i].I, TimeSeries[i].R, TimeSeries[i].incI,
				TimeSeries[i].Mild		, TimeSeries[i].ILI		, TimeSeries[i].SARI	, TimeSeries[i].Critical	, TimeSeries[i].CritRecov	,
				TimeSeries[i].incMild	, TimeSeries[i].incILI	, TimeSeries[i].incSARI	, TimeSeries[i].incCritical	, TimeSeries[i].incCritRecov, TimeSeries[i].incD,
				TimeSeries[i].cumMild	, TimeSeries[i].cumILI	, TimeSeries[i].cumSARI	, TimeSeries[i].cumCritical	, TimeSeries[i].cumCritRecov, TimeSeries[i].D);
		}
		fclose(dat);

		if((P.DoAdUnits) && (P.OutputSeverityAdminUnit))
		{
			//// output severity results by admin unit
			sprintf(outname, "%s.severity.adunit.xls", OutFile);
			if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
			fprintf(dat, "t");

			/////// ****** /////// ****** /////// ****** COLNAMES
			//// prevalance
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tMild_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tILI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tSARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritical_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritRecov_%s", AdUnits[i].ad_name);

			//// incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincMild_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincILI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincSARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritical_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritRecov_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincDeath_adu%s"	, AdUnits[i].ad_name);

			//// cumulative incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumMild_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumILI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumSARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumCritical_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumCritRecov_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumDeaths_%s"	, AdUnits[i].ad_name);

			fprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table. 
			for(i = 0; i < P.NumSamples; i++)
			{
				fprintf(dat, "%.10f", TimeSeries[i].t);

				//// prevalance
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
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%lg", TimeSeries[i].incD_adunit[j]);

				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", TimeSeries[i].cumD_adunit[j]);

				if(i != P.NumSamples - 1) fprintf(dat, "\n");
			}
			fclose(dat);
		}
	}
}

void SaveSummaryResults(void) //// calculates and saves summary results (called for average of extinct and non-extinct realisation time series - look in main) 
{
	int i, j;
	double c, t,s;
	FILE* dat;
	char outname[1024];

	c = 1 / ((double)(P.NRactE + P.NRactNE));

	if (P.OutputNonSeverity)
	{
		sprintf(outname, "%s.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		//// set colnames
		fprintf(dat, "t\tS\tL\tI\tR\tD\tincI\tincR\tincD\tincC\tincDC\tincTC\tincH\tcumT\tcumTmax\tcumTP\tcumV\tcumVmax\tExtinct\trmsRad\tmaxRad\tvS\tvI\tvR\tvD\tvincI\tvincR\tvincFC\tvincC\tvincDC\tvincTC\tvincH\tvrmsRad\tvmaxRad\t\t%i\t%i\t%.10f\t%.10f\t%.10f\t\t%.10f\t%.10f\t%.10f\t%.10f\n",
			P.NRactNE, P.NRactE, P.R0household, P.R0places, P.R0spatial, c * PeakHeightSum, c * PeakHeightSS - c * c * PeakHeightSum * PeakHeightSum, c * PeakTimeSum, c * PeakTimeSS - c * c * PeakTimeSum * PeakTimeSum);
		c = 1 / ((double)P.NRactual);
		//added this for sake of test border control
		//c=1;
		//// populate table
		for(i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%.10f\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t%10lf\t",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].L, c * TSMean[i].I, c * TSMean[i].R,
				c * TSMean[i].D, c * TSMean[i].incI, c * TSMean[i].incR, c * TSMean[i].incFC, c * TSMean[i].incC, c * TSMean[i].incDC, c * TSMean[i].incTC, c * TSMean[i].incH, 
				c * TSMean[i].cumT, TSMean[i].cumTmax, c * TSMean[i].cumTP, c * TSMean[i].cumV, TSMean[i].cumVmax, c * TSMean[i].extinct, c * TSMean[i].rmsRad, c * TSMean[i].maxRad);
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
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
				c * TSVar[i].incH	- c * c * TSMean[i].incH	* TSMean[i].incH, //added hospitalisation
				c * TSVar[i].rmsRad - c * c * TSMean[i].rmsRad	* TSMean[i].rmsRad,
				c * TSVar[i].maxRad - c * c * TSMean[i].maxRad	* TSMean[i].maxRad);
		}
		fclose(dat);
	}

	if (P.OutputControls)
	{
		sprintf(outname, "%s.controls.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tS\tincC\tincTC\tincFC\tincH\tcumT\tcumUT\tcumTP\tcumV\tincHQ\tincAC\tincAH\tincAA\tincACS\tincAPC\tincAPA\tincAPCS\tpropSocDist");
		for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tprClosed_%i", j);
		fprintf(dat, "t\tvS\tvincC\tvincTC\tvincFC\tvincH\tvcumT\tvcumUT\tvcumTP\tvcumV");
		for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\tvprClosed_%i", j);
		fprintf(dat, "\n");
		for(i = 0; i < P.NumSamples; i++)
			{
			fprintf(dat, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				c * TSMean[i].t, c * TSMean[i].S, c * TSMean[i].incC, c * TSMean[i].incTC, c * TSMean[i].incFC, c * TSMean[i].incH,
				c * TSMean[i].cumT, c * TSMean[i].cumUT, c * TSMean[i].cumTP, c * TSMean[i].cumV, c * TSMean[i].incHQ,
				c * TSMean[i].incAC, c * TSMean[i].incAH, c * TSMean[i].incAA, c * TSMean[i].incACS,
				c * TSMean[i].incAPC, c * TSMean[i].incAPA, c * TSMean[i].incAPCS,c*TSMean[i].PropSocDist);
			for(j = 0; j < NUM_PLACE_TYPES; j++) fprintf(dat, "\t%lf", c * TSMean[i].PropPlacesClosed[j]);
			fprintf(dat, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				c * TSVar[i].S - c * c * TSMean[i].S * TSMean[i].S,
				c * TSVar[i].incC - c * c * TSMean[i].incC * TSMean[i].incC,
				c * TSVar[i].incTC - c * c * TSMean[i].incTC * TSMean[i].incTC,
				c * TSVar[i].incFC - c * c * TSMean[i].incFC * TSMean[i].incFC,
				c * TSVar[i].incH - c * c * TSMean[i].incH * TSMean[i].incH,
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
		sprintf(outname, "%s.age.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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
		sprintf(outname, "%s.adunit.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t");
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tI_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tC_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tDC_%s", AdUnits[i].ad_name); //added detected cases: ggilani 03/02/15
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tT_%s", AdUnits[i].ad_name);
		for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\t%.10f", P.PopByAdunit[i][0]);
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
			sprintf(outname, "%s.adunitVar.xls", OutFile);
			if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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

	if(P.KeyWorkerProphTimeStartBase < P.SampleTime)
	{
		sprintf(outname, "%s.keyworker.xls", OutFile);
		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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
		sprintf(outname, "%s.inftype.xls", OutFile);
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		fprintf(dat, "t\tR");
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\tRtype_%i", j);
		for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\tincItype_%i", j);
		for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\tRage_%i", j);
		fprintf(dat, "\n");
		for (i = 0; i < P.NumSamples; i++)
		{
			fprintf(dat, "%lf\t%lf", c * TSMean[i].t, c * TSMean[i].Rdenom);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\t%lf", c * TSMean[i].Rtype[j]);
			for (j = 0; j < INFECT_TYPE_MASK; j++) fprintf(dat, "\t%lf", c * TSMean[i].incItype[j]);
			for (j = 0; j < NUM_AGE_GROUPS; j++) fprintf(dat, "\t%lf", c * TSMean[i].Rage[j]);
			fprintf(dat, "\n");
		}
		fclose(dat);
	}

	if (P.OutputR0)
	{
		sprintf(outname, "%s.R0.xls", OutFile);
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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
		sprintf(outname, "%s.household.xls", OutFile);
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
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
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
		sprintf(outname, "%s.country.xls", OutFile);
		if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
		for (i = 0; i < MAX_COUNTRIES; i++)
			fprintf(dat, "%i\t%.10f\t%.10f\n", i, infcountry_av[i] * c, infcountry_num[i] * c);
		fclose(dat);
	}

	if(P.DoSeverity)
	{
		printf("\nSaveSummaryResults: writing Severity");
		//// output separate severity file (can integrate with main if need be)
		sprintf(outname, "%s.severity.xls", OutFile);
		printf("\nFilename: %s", outname);

		if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open severity output file\n");
		//// set colnames
		fprintf(dat, "t\tPropSchClosed\tPropSocDist\tS\tL\tI\tR\tD\tincI\tincC\tincDC\tMild\tILI\tSARI\tCritical\tCritRecov\tincMild\tincILI\tincSARI\tincCritical\tincCritRecov\tincDeath\tcumMild\tcumILI\tcumSARI\tcumCritical\tcumCritRecov\tvarS\tvarI\tvarR\tvarD\tvarMild\tvarILI\tvarSARI\tvarCritical\tvarCritRecov\tvarincI\tvarincC\tvarincDC\tvarincMild\tvarincILI\tvarincSARI\tvarincCritical\tvarincCritRecov\tvarincDeath\tvarcumMild\tvarcumILI\tvarcumSARI\tvarcumCritical\tvarcumCritRecov\n");//\t\t%lf\t%lf\t%lf\n",P.R0household,P.R0places,P.R0spatial);
		c = 1 / ((double)P.NRactual);
		//// populate table
		s = 0.0;
		for(i = 0; i < P.NumSamples; i++)
		{
			t=0.0;
			//// means
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t",
				c * TSMean[i].t			, c * TSMean[i].PropPlacesClosed[0], c * TSMean[i].PropSocDist, 
				c * TSMean[i].S			, c * TSMean[i].L		, c * TSMean[i].I		, c * TSMean[i].R			, c * TSMean[i].D, c * TSMean[i].incI, c * TSMean[i].incC, c * TSMean[i].incDC,
				c * TSMean[i].Mild		, c * TSMean[i].ILI		, c * TSMean[i].SARI	, c * TSMean[i].Critical	, c * TSMean[i].CritRecov	, 
				c * TSMean[i].incMild	, c * TSMean[i].incILI	, c * TSMean[i].incSARI	, c * TSMean[i].incCritical	, c * TSMean[i].incCritRecov, c * TSMean[i].incD, 
				c * TSMean[i].cumMild	, c * TSMean[i].cumILI	, c * TSMean[i].cumSARI	, c * TSMean[i].cumCritical	, c * TSMean[i].cumCritRecov);
			//// variances
			fprintf(dat, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				c * TSVar[i].S				- c * c * TSMean[i].S				* TSMean[i].S,
				c * TSVar[i].I				- c * c * TSMean[i].I				* TSMean[i].I,
				c * TSVar[i].R				- c * c * TSMean[i].R				* TSMean[i].R,
				c * TSVar[i].D				- c * c * TSMean[i].D				* TSMean[i].D,
				c * TSVar[i].Mild			- c * c * TSMean[i].Mild			* TSMean[i].Mild,
				c * TSVar[i].ILI			- c * c * TSMean[i].ILI				* TSMean[i].ILI,
				c * TSVar[i].SARI			- c * c * TSMean[i].SARI			* TSMean[i].SARI,
				c * TSVar[i].Critical		- c * c * TSMean[i].Critical		* TSMean[i].Critical,
				c * TSVar[i].CritRecov		- c * c * TSMean[i].CritRecov		* TSMean[i].CritRecov,
				c * TSVar[i].incI			- c * c * TSMean[i].incI			* TSMean[i].incI,
				c * TSVar[i].incC			- c * c * TSMean[i].incC			* TSMean[i].incC,
				c * TSVar[i].incDC			- c * c * TSMean[i].incDC			* TSMean[i].incDC,
				c * TSVar[i].incMild		- c * c * TSMean[i].incMild			* TSMean[i].incMild,
				c * TSVar[i].incILI			- c * c * TSMean[i].incILI			* TSMean[i].incILI,
				c * TSVar[i].incSARI		- c * c * TSMean[i].incSARI			* TSMean[i].incSARI,
				c * TSVar[i].incCritical	- c * c * TSMean[i].incCritical		* TSMean[i].incCritical,
				c * TSVar[i].incCritRecov	- c * c * TSMean[i].incCritRecov	* TSMean[i].incCritRecov,
				c * TSVar[i].incD			- c * c * TSMean[i].incD			* TSMean[i].incD,
				c * TSVar[i].cumMild		- c * c * TSMean[i].cumMild			* TSMean[i].cumMild,
				c * TSVar[i].cumILI			- c * c * TSMean[i].cumILI			* TSMean[i].cumILI,
				c * TSVar[i].cumSARI		- c * c * TSMean[i].cumSARI			* TSMean[i].cumSARI,
				c * TSVar[i].cumCritical	- c * c * TSMean[i].cumCritical		* TSMean[i].cumCritical,
				c * TSVar[i].cumCritRecov	- c * c * TSMean[i].cumCritRecov	* TSMean[i].cumCritRecov);
		}
		fclose(dat);

		if((P.DoAdUnits) && (P.OutputSeverityAdminUnit))
		{
			//// output severity results by admin unit
			sprintf(outname, "%s.severity.adunit.xls", OutFile);
			if(!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
			fprintf(dat, "t");

			/////// ****** /////// ****** /////// ****** COLNAMES
			//// **** means
			//// prevalance
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tMild_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tILI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tSARI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritical_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tCritRecov_%s"	, AdUnits[i].ad_name);

			//// incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincI_%s"				, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincMild_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincILI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincSARI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritical_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincCritRecov_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tincDeath_%s"			, AdUnits[i].ad_name);

			//// cumulative incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumMild_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumILI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumSARI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumCritical_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumCritRecov_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tcumDeath_%s"			, AdUnits[i].ad_name);

			//// **** variances
			//// prevalance
			for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarMild_%s"			, AdUnits[i].ad_name);
			for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarILI_%s"			, AdUnits[i].ad_name);
			for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarSARI_%s"			, AdUnits[i].ad_name);
			for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarCritical_%s"		, AdUnits[i].ad_name);
			for(i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarCritRecov_%s"		, AdUnits[i].ad_name);

			//// incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincI_%s"			, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincMild_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincILI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincSARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincCritical_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincCritRecov_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarincDeath_%s"		, AdUnits[i].ad_name);

			//// cumulative incidence
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarcumMild_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarcumILI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarcumSARI_%s"		, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarcumCritical_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarcumCritRecov_%s"	, AdUnits[i].ad_name);
			for (i = 0; i < P.NumAdunits; i++) fprintf(dat, "\tvarcumDeath_%s"		, AdUnits[i].ad_name);

			fprintf(dat, "\n");

			/////// ****** /////// ****** /////// ****** Populate table. 
			for(i = 0; i < P.NumSamples; i++)
			{
				fprintf(dat, "%.10f", c * TSMean[i].t);
				//// **** means
				//// prevalance
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].Mild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].Critical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].CritRecov_adunit[j]);
				    
				//// incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].incD_adunit[j]);

				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSMean[i].cumD_adunit[j]);

				//// **** variances
				//// prevalance
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].Mild_adunit[j]		 - c * c * TSMean[i].Mild_adunit[j]			* TSMean[i].Mild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].ILI_adunit[j]		 - c * c * TSMean[i].ILI_adunit[j]			* TSMean[i].ILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].SARI_adunit[j]		 - c * c * TSMean[i].SARI_adunit[j]			* TSMean[i].SARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].Critical_adunit[j]	 - c * c * TSMean[i].Critical_adunit[j]		* TSMean[i].Critical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].CritRecov_adunit[j] - c * c * TSMean[i].CritRecov_adunit[j]	* TSMean[i].CritRecov_adunit[j]);
				    
				//// incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incI_adunit[j]			- c * c * TSMean[i].incI_adunit[j]				* TSMean[i].incI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incMild_adunit[j]		- c * c * TSMean[i].incMild_adunit[j]			* TSMean[i].incMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incILI_adunit[j]		- c * c * TSMean[i].incILI_adunit[j]			* TSMean[i].incILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incSARI_adunit[j]		- c * c * TSMean[i].incSARI_adunit[j]			* TSMean[i].incSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incCritical_adunit[j]	- c * c * TSMean[i].incCritical_adunit[j]		* TSMean[i].incCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incCritRecov_adunit[j]	- c * c * TSMean[i].incCritRecov_adunit[j]		* TSMean[i].incCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].incD_adunit[j]			- c * c * TSMean[i].incD_adunit[j]				* TSMean[i].incD_adunit[j]);
				    
				//// cumulative incidence
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].cumMild_adunit[j]		- c * c * TSMean[i].cumMild_adunit[j]			* TSMean[i].cumMild_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].cumILI_adunit[j]		- c * c * TSMean[i].cumILI_adunit[j]			* TSMean[i].cumILI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].cumSARI_adunit[j]		- c * c * TSMean[i].cumSARI_adunit[j]			* TSMean[i].cumSARI_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].cumCritical_adunit[j]	- c * c * TSMean[i].cumCritical_adunit[j]		* TSMean[i].cumCritical_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].cumCritRecov_adunit[j]	- c * c * TSMean[i].cumCritRecov_adunit[j]		* TSMean[i].cumCritRecov_adunit[j]);
				for (j = 0; j < P.NumAdunits; j++)		fprintf(dat, "\t%.10f", c * TSVar[i].cumD_adunit[j]			- c * c * TSMean[i].cumD_adunit[j]				* TSMean[i].cumD_adunit[j]);

				if(i != P.NumSamples - 1) fprintf(dat, "\n");
			}
			fclose(dat);
		}
	}
}

void SaveRandomSeeds(void)
{
	/* function: SaveRandomSeeds(void)
	 *
	 * Purpose: outputs the random seeds used for each run to a file
	 * Parameter: none
	 * Returns: none
	 *
	 * Author: ggilani, 09/03/17
	 */
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.seeds.xls", OutFile);
	if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "%li\t%li\n", P.newseed1, P.newseed2);
	fclose(dat);
}
void SaveEvents(void)
{
	/* function: SaveEvents(void)
	 *
	 * Purpose: outputs event log to a csv file if required
	 * Parameters: none
	 * Returns: none
	 *
	 * Author: ggilani, 15/10/2014
	 */
	int i;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.infevents.xls", OutFile);
	if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
	fprintf(dat, "type,t,thread,ind_infectee,cell_infectee,listpos_infectee,adunit_infectee,x_infectee,y_infectee,t_infector,ind_infector,cell_infector\n");
	for (i = 0; i < *nEvents; i++)
	{
		fprintf(dat, "%i\t%.10f\t%i\t%i\t%i\t%i\t%i\t%.10f\t%.10f\t%.10f\t%i\t%i\n",
			InfEventLog[i].type, InfEventLog[i].t, InfEventLog[i].thread, InfEventLog[i].infectee_ind, InfEventLog[i].infectee_cell, InfEventLog[i].listpos, InfEventLog[i].infectee_adunit, InfEventLog[i].infectee_x, InfEventLog[i].infectee_y, InfEventLog[i].t_infector, InfEventLog[i].infector_ind, InfEventLog[i].infector_cell);
	}
	fclose(dat);
}

void LoadSnapshot(void)
{
	FILE* dat;
	int i, j, * CellMemberArray, * CellSuscMemberArray;
	long l;
	long long CM_offset, CSM_offset;
	double t;
	int** Array_InvCDF;
	float* Array_tot_prob, ** Array_cum_trans, ** Array_max_trans;

	if (!(dat = fopen(SnapshotLoadFile, "rb"))) ERR_CRITICAL("Unable to open snapshot file\n");
	fprintf(stderr, "Loading snapshot.");
	if (!(Array_InvCDF = (int**)malloc(P.NCP * sizeof(int*)))) ERR_CRITICAL("Unable to allocate temp cell storage\n");
	if (!(Array_max_trans = (float**)malloc(P.NCP * sizeof(float*)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	if (!(Array_cum_trans = (float**)malloc(P.NCP * sizeof(float*)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	if (!(Array_tot_prob = (float*)malloc(P.NCP * sizeof(float)))) ERR_CRITICAL("Unable to temp allocate cell storage\n");
	for (i = 0; i < P.NCP; i++)
	{
		Array_InvCDF[i] = Cells[i].InvCDF;
		Array_max_trans[i] = Cells[i].max_trans;
		Array_cum_trans[i] = Cells[i].cum_trans;
		Array_tot_prob[i] = Cells[i].tot_prob;
	}

	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.N) ERR_CRITICAL_FMT("Incorrect N (%i %i) in snapshot file.\n", P.N, i);
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.NH) ERR_CRITICAL("Incorrect NH in snapshot file.\n");
	fread_big((void*)&i, sizeof(int), 1, dat); if (i != P.NC) ERR_CRITICAL_FMT("## %i neq %i\nIncorrect NC in snapshot file.", i, P.NC);
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.NCP) ERR_CRITICAL("Incorrect NCP in snapshot file.\n");
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.ncw) ERR_CRITICAL("Incorrect ncw in snapshot file.\n");
	fread_big((void*)& i, sizeof(int), 1, dat); if (i != P.nch) ERR_CRITICAL("Incorrect nch in snapshot file.\n");
	fread_big((void*)& l, sizeof(long), 1, dat); if (l != P.seed1) ERR_CRITICAL("Incorrect seed1 in snapshot file.\n");
	fread_big((void*)& l, sizeof(long), 1, dat); if (l != P.seed2) ERR_CRITICAL("Incorrect seed2 in snapshot file.\n");
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

	zfread_big((void*)Hosts, sizeof(person), (size_t)P.N, dat);
	fprintf(stderr, ".");
	zfread_big((void*)Households, sizeof(household), (size_t)P.NH, dat);
	fprintf(stderr, ".");
	zfread_big((void*)Cells, sizeof(cell), (size_t)P.NC, dat);
	fprintf(stderr, ".");
	zfread_big((void*)Mcells, sizeof(microcell), (size_t)P.NMC, dat);
	fprintf(stderr, ".");
	zfread_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, ".");
	zfread_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.N, dat);
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
	free(Array_tot_prob);
	free(Array_cum_trans);
	free(Array_max_trans);
	free(Array_InvCDF);
	fprintf(stderr, "\n");
	fclose(dat);
}

void SaveSnapshot(void)
{
	FILE* dat;
	int i = 1;

	if (!(dat = fopen(SnapshotSaveFile, "wb"))) ERR_CRITICAL("Unable to open snapshot file\n");

	fwrite_big((void*) & (P.N), sizeof(int), 1, dat);
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
	fwrite_big((void*) & (P.seed1), sizeof(long), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.seed2), sizeof(long), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.TimeStep), sizeof(double), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (P.SnapshotSaveTime), sizeof(double), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (State.CellMemberArray), sizeof(int*), 1, dat);
	fprintf(stderr, "## %i\n", i++);
	fwrite_big((void*) & (State.CellSuscMemberArray), sizeof(int*), 1, dat);
	fprintf(stderr, "## %i\n", i++);

	zfwrite_big((void*)Hosts, sizeof(person), (size_t)P.N, dat);

	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)Households, sizeof(household), (size_t)P.NH, dat);
	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)Cells, sizeof(cell), (size_t)P.NC, dat);
	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)Mcells, sizeof(microcell), (size_t)P.NMC, dat);
	fprintf(stderr, "## %i\n", i++);

	zfwrite_big((void*)State.CellMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, "## %i\n", i++);
	zfwrite_big((void*)State.CellSuscMemberArray, sizeof(int), (size_t)P.N, dat);
	fprintf(stderr, "## %i\n", i++);

	fclose(dat);
}

void UpdateProbs(int DoPlace)
{
	int j;

	if (!DoPlace)
	{
#pragma omp parallel for private(j) schedule(static,500)
		for (j = 0; j < P.NCP; j++)
		{
			CellLookup[j]->tot_prob = 0;
			CellLookup[j]->S0 = CellLookup[j]->S + CellLookup[j]->L + CellLookup[j]->I;
			if ((P.DoDeath) && (!P.DoSIS))
			{
				CellLookup[j]->S0 += CellLookup[j]->n / 5;
				if ((CellLookup[j]->n < 100) || (CellLookup[j]->S0 > CellLookup[j]->n)) CellLookup[j]->S0 = CellLookup[j]->n;
			}
		}
	}
	else
	{
#pragma omp parallel for private(j) schedule(static,500) 
		for (j = 0; j < P.NCP; j++)
		{
			CellLookup[j]->S0 = CellLookup[j]->S;
			CellLookup[j]->tot_prob = 0;
		}
	}
#pragma omp parallel for private(j) schedule(static,500)
	for (j = 0; j < P.NCP; j++)
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
			VariableAndValue = (int)floor(((double)State.trigDC) * P.GlobalIncThreshPop / ((double)P.N));
		else
			VariableAndValue = State.trigDC;
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
void DoOrDontAmendStartTime (double *StartTimeToAmend, double StartTime)
{
	if (*StartTimeToAmend >= 1e10) *StartTimeToAmend = StartTime;
}


void RecordSample(double t, int n)
{
	int i, j, k, S, L, I, R, D, N, cumC, cumTC, cumI, cumR, cumD, cumDC, cumFC; 
	int cumH; //add number of hospitalised, cumulative hospitalisation: ggilani 28/10/14
	int cumCT; //added cumulative number of contact traced: ggilani 15/06/17
	int cumCC; //added cumulative number of cases who are contacts: ggilani 28/05/2019
	int cumDCT; //added cumulative number of cases who are digitally contact traced: ggilani 11/03/20
	int cumHQ, cumAC, cumAH, cumAA, cumACS, cumAPC, cumAPA, cumAPCS, numPC, trigDC,trigAlert;
	int cumC_country[MAX_COUNTRIES]; //add cumulative cases per country
	cell* ct;
	unsigned short int ts;
	float nsy;
	
	//// Severity quantities
	int Mild, ILI, SARI, Critical, CritRecov, cumMild, cumILI, cumSARI, cumCritical, cumCritRecov;
	//int Mild_adunit[MAX_ADUNITS], ILI_adunit[MAX_ADUNITS], SARI_adunit[MAX_ADUNITS], Critical_adunit[MAX_ADUNITS], CritRecov_adunit[MAX_ADUNITS];
	//int cumMild_adunit[MAX_ADUNITS], cumILI_adunit[MAX_ADUNITS], cumSARI_adunit[MAX_ADUNITS], cumCritical_adunit[MAX_ADUNITS], cumCritRecov_adunit[MAX_ADUNITS]; /// ...



	nsy = (float)(DAYS_PER_YEAR / P.SampleStep);
	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	//// initialize to zero
	S = L = I = R = D = cumI = cumC = cumDC = cumTC = cumFC = cumHQ = cumAC = cumAA = cumAH = cumACS = cumAPC = cumAPA = cumAPCS = cumD = cumH = cumCT = cumCC = cumDCT = 0;
	for (i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] = 0;
	if (P.DoSeverity)
		Mild = ILI = SARI = Critical = CritRecov = cumMild = cumILI = cumSARI = cumCritical = cumCritRecov = 0;

#pragma omp parallel for private(i,ct) schedule(static,10000) reduction(+:S,L,I,R,D,cumTC) //added i to private
	for (i = 0; i < P.NCP; i++)
	{
		ct = CellLookup[i];
		S += (int)ct->S;
		L += (int)ct->L;
		I += (int)ct->I;
		R += (int)ct->R;
		D += (int)ct->D;
		cumTC += (int)ct->cumTC;
	}
	cumR = R;
	cumD = D;
	N = S + L + I + R + D;
	if (N != P.N) fprintf(stderr, "## %i #\n", P.N - N);
	State.sumRad2 = 0;
	for (j = 0; j < P.NumThreads; j++)
	{
		cumI += StateT[j].cumI;
		cumC += StateT[j].cumC;
		cumDC += StateT[j].cumDC;
		cumFC += StateT[j].cumFC;
		cumH += StateT[j].cumH; //added cumulative hospitalisation
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
			Mild			+= StateT[j].Mild			;
			ILI				+= StateT[j].ILI			;
			SARI			+= StateT[j].SARI			;
			Critical		+= StateT[j].Critical		;
			CritRecov		+= StateT[j].CritRecov		;
	
			///// cumulative severity states by thread
			cumMild			+= StateT[j].cumMild		;
			cumILI			+= StateT[j].cumILI			;
			cumSARI			+= StateT[j].cumSARI		;
			cumCritical		+= StateT[j].cumCritical	;
			cumCritRecov	+= StateT[j].cumCritRecov	;
		}

		//add up cumulative country counts: ggilani - 12/11/14
		for (i = 0; i < MAX_COUNTRIES; i++) cumC_country[i] += StateT[j].cumC_country[i];
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
	TimeSeries[n].incH = (double)(cumH - State.cumH); //added incidence of hospitalisation
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
	//incidence per country
	for (i = 0; i < MAX_COUNTRIES; i++) TimeSeries[n].incC_country[i] = (double)(cumC_country[i] - State.cumC_country[i]);
	if (P.DoICUTriggers)
	{
		trigDC = cumCritical;
		if (n >= P.TriggersSamplingInterval) trigDC -= (int)TimeSeries[n - P.TriggersSamplingInterval].cumCritical;
	}
	else
	{
		trigDC = cumDC;
		if (n >= P.TriggersSamplingInterval) trigDC -= (int)TimeSeries[n - P.TriggersSamplingInterval].cumDC;
	}
	State.trigDC = trigDC;

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
	State.cumH = cumH; //added cumulative hospitalisation
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

	if (P.DoSeverity)
	{
		//// Record incidence. (Must be done with old State totals)
		TimeSeries[n].incMild			= (double)(cumMild			- State.cumMild			);
		TimeSeries[n].incILI			= (double)(cumILI			- State.cumILI			);
		TimeSeries[n].incSARI			= (double)(cumSARI			- State.cumSARI			);
		TimeSeries[n].incCritical		= (double)(cumCritical		- State.cumCritical		);
		TimeSeries[n].incCritRecov		= (double)(cumCritRecov		- State.cumCritRecov	);

		/////// update state with totals
		State.Mild			= Mild			;
		State.ILI			= ILI			;
		State.SARI			= SARI			;
		State.Critical		= Critical		;
		State.CritRecov		= CritRecov		;
		State.cumMild		= cumMild		;
		State.cumILI		= cumILI		;
		State.cumSARI		= cumSARI		;
		State.cumCritical	= cumCritical	;
		State.cumCritRecov	= cumCritRecov	;

		//// Record new totals for time series. (Must be done with old State totals)
		TimeSeries[n].Mild			= Mild			;
		TimeSeries[n].ILI			= ILI			;
		TimeSeries[n].SARI			= SARI			;
		TimeSeries[n].Critical		= Critical		;
		TimeSeries[n].CritRecov		= CritRecov		;
		TimeSeries[n].cumMild		= cumMild		;
		TimeSeries[n].cumILI		= cumILI		;
		TimeSeries[n].cumSARI		= cumSARI		;
		TimeSeries[n].cumCritical	= cumCritical	;
		TimeSeries[n].cumCritRecov	= cumCritRecov	;

		if (P.DoAdUnits)
			for (i = 0; i <= P.NumAdunits; i++)
			{
				//// Record incidence. Need new total minus old total (same as minus old total plus new total).  
				//// First subtract old total while unchanged. 
				TimeSeries[n].incMild_adunit		[i] = (double)(-State.cumMild_adunit		[i]);
				TimeSeries[n].incILI_adunit			[i] = (double)(-State.cumILI_adunit			[i]);
				TimeSeries[n].incSARI_adunit		[i] = (double)(-State.cumSARI_adunit		[i]);
				TimeSeries[n].incCritical_adunit	[i] = (double)(-State.cumCritical_adunit	[i]);
				TimeSeries[n].incCritRecov_adunit	[i] = (double)(-State.cumCritRecov_adunit	[i]);
				TimeSeries[n].incD_adunit			[i] = (double)(-State.cumD_adunit			[i]);

				//// reset State (don't think StateT) to zero. Don't need to do this with non-admin unit as local variables Mild, cumSARI etc. initialized to zero at beginning of function. Check with Gemma
				State.Mild_adunit			[i] = 0;
				State.ILI_adunit			[i] = 0;
				State.SARI_adunit			[i] = 0;
				State.Critical_adunit		[i] = 0;
				State.CritRecov_adunit		[i] = 0;
				State.cumMild_adunit		[i] = 0;
				State.cumILI_adunit			[i] = 0;
				State.cumSARI_adunit		[i] = 0;
				State.cumCritical_adunit	[i] = 0;
				State.cumCritRecov_adunit	[i] = 0;
				State.cumD_adunit			[i] = 0;

				for (j = 0; j < P.NumThreads; j++)
				{
					//// collate from threads
					State.Mild_adunit			[i] += StateT[j].Mild_adunit			[i];
					State.ILI_adunit			[i] += StateT[j].ILI_adunit				[i];
					State.SARI_adunit			[i] += StateT[j].SARI_adunit			[i];
					State.Critical_adunit		[i] += StateT[j].Critical_adunit		[i];
					State.CritRecov_adunit		[i] += StateT[j].CritRecov_adunit		[i];
					State.cumMild_adunit		[i] += StateT[j].cumMild_adunit			[i];
					State.cumILI_adunit			[i] += StateT[j].cumILI_adunit			[i];
					State.cumSARI_adunit		[i] += StateT[j].cumSARI_adunit			[i];
					State.cumCritical_adunit	[i] += StateT[j].cumCritical_adunit		[i];
					State.cumCritRecov_adunit	[i] += StateT[j].cumCritRecov_adunit	[i];
					State.cumD_adunit			[i] += StateT[j].cumD_adunit			[i];
				}

				//// Record incidence. Need new total minus old total. Add new total 
				TimeSeries[n].incMild_adunit		[i] += (double)(State.cumMild_adunit		[i]);
				TimeSeries[n].incILI_adunit			[i] += (double)(State.cumILI_adunit			[i]);
				TimeSeries[n].incSARI_adunit		[i] += (double)(State.cumSARI_adunit		[i]);
				TimeSeries[n].incCritical_adunit	[i] += (double)(State.cumCritical_adunit	[i]);
				TimeSeries[n].incCritRecov_adunit	[i] += (double)(State.cumCritRecov_adunit	[i]);
				TimeSeries[n].incD_adunit			[i] += (double)(State.cumD_adunit			[i]);

				//// Record new totals for time series. (Must be done with old State totals)
				TimeSeries[n].Mild_adunit			[i] = State.Mild_adunit			[i];
				TimeSeries[n].ILI_adunit			[i] = State.ILI_adunit			[i];
				TimeSeries[n].SARI_adunit			[i] = State.SARI_adunit			[i];
				TimeSeries[n].Critical_adunit		[i] = State.Critical_adunit		[i];
				TimeSeries[n].CritRecov_adunit		[i] = State.CritRecov_adunit	[i];
				TimeSeries[n].cumMild_adunit		[i] = State.cumMild_adunit		[i];
				TimeSeries[n].cumILI_adunit			[i] = State.cumILI_adunit		[i];
				TimeSeries[n].cumSARI_adunit		[i] = State.cumSARI_adunit		[i];
				TimeSeries[n].cumCritical_adunit	[i] = State.cumCritical_adunit	[i];
				TimeSeries[n].cumCritRecov_adunit	[i] = State.cumCritRecov_adunit	[i];
				TimeSeries[n].cumD_adunit			[i] = State.cumD_adunit			[i];
			}
	}

	//update cumulative cases per country
	for (i = 0; i < MAX_COUNTRIES; i++) State.cumC_country[i] = cumC_country[i];
	//update overall state variable for cumulative cases per adunit

	TimeSeries[n].rmsRad = (State.cumI > 0) ? sqrt(State.sumRad2 / ((double)State.cumI)) : 0;
	TimeSeries[n].maxRad = sqrt(State.maxRad2);
	TimeSeries[n].extinct = ((((P.SmallEpidemicCases >= 0) && (State.R <= P.SmallEpidemicCases)) || (P.SmallEpidemicCases < 0)) && (State.I + State.L == 0)) ? 1 : 0;
	for (i = 0; i < NUM_AGE_GROUPS; i++)
	{
		TimeSeries[n].incCa[i] = TimeSeries[n].incIa[i] = TimeSeries[n].incDa[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incCa[i] += (double)StateT[j].cumCa[i];
			TimeSeries[n].incIa[i] += (double)StateT[j].cumIa[i];
			TimeSeries[n].incDa[i] += (double)StateT[j].cumDa[i];
		}
	}

	for (i = 0; i < 2; i++)
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

	for (i = 0; i < INFECT_TYPE_MASK; i++)
	{
		TimeSeries[n].incItype[i] = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			TimeSeries[n].incItype[i] += (double)StateT[j].cumItype[i];
			StateT[j].cumItype[i] = 0;
		}
	}
	if (P.DoAdUnits)
		for (i = 0; i <= P.NumAdunits; i++)
		{
			TimeSeries[n].incI_adunit[i] = TimeSeries[n].incC_adunit[i] = TimeSeries[n].cumT_adunit[i] = TimeSeries[n].incH_adunit[i] = TimeSeries[n].incDC_adunit[i] = TimeSeries[n].incCT_adunit[i] = 0; //added detected cases: ggilani 03/02/15
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
		for (i = 0; i < P.NumAdunits; i++)
			TimeSeries[n].DCT_adunit[i] = (double)AdUnits[i].ndct; //added total numbers of contacts currently isolated due to digital contact tracing: ggilani 11/03/20
	if ((P.DoPlaces) && (t >= P.PlaceCloseTimeStart))
		for (i = 0; i < NUM_PLACE_TYPES; i++)
		{
			numPC = 0;
			for (j = 0; j < P.Nplace[i]; j++)
				if (PLACE_CLOSED(i, j)) numPC++;
			State.NumPlacesClosed[i] = numPC;
			TimeSeries[n].PropPlacesClosed[i] = ((double)numPC) / ((double)P.Nplace[i]);
		}
	for (i = k = 0; i < P.NMC; i++) if (Mcells[i].socdist == 2) k++;
	TimeSeries[n].PropSocDist=((double)k)/((double)P.NMC);
	

	


	if (P.PreControlClusterIdUseDeaths)
	{
		trigAlert = (int)TimeSeries[n].D;
		if (n >= P.PreControlClusterIdDuration) trigAlert -= (int) TimeSeries[n - P.PreControlClusterIdDuration].D;
	}
	else
	{
		trigAlert = State.cumDC;
		if (n >= P.PreControlClusterIdDuration) trigAlert -= (int)TimeSeries[n - P.PreControlClusterIdDuration].cumDC;
	}
	if (trigAlert >= P.PreControlClusterIdCaseThreshold)
	{
		if (P.PreControlClusterIdTime == 0)
		{
			P.PreControlClusterIdTime = t;
			if (P.PreControlClusterIdCalTime >= 0)
			{
				P.PreControlClusterIdHolOffset = P.PreControlClusterIdTime - P.PreControlClusterIdCalTime;
				InterruptRun = 0;
			}
		}
		P.ControlPropCasesId = P.PostAlertControlPropCasesId;

		//// Set Case isolation start time (by admin unit)
		for (i = 0; i < P.NumAdunits; i++)
			if (ChooseTriggerVariableAndValue(i) > ChooseThreshold(i, P.CaseIsolation_CellIncThresh)) //// a little wasteful if doing Global trigs as function called more times than necessary, but worth it for much simpler code. Also this funcion is small portion of runtime.
			{
				if (P.DoInterventionDelaysByAdUnit)
					DoOrDontAmendStartTime(&AdUnits[i].CaseIsolationTimeStart, t + AdUnits[i].CaseIsolationDelay);
				else
					DoOrDontAmendStartTime(&AdUnits[i].CaseIsolationTimeStart, t + P.CaseIsolationTimeStartBase);
			}

		//// Set Household Quarantine start time (by admin unit)
		for (i = 0; i < P.NumAdunits; i++)
			if (ChooseTriggerVariableAndValue(i) > ChooseThreshold(i, P.HHQuar_CellIncThresh)) //// a little wasteful if doing Global trigs as function called more times than necessary, but worth it for much simpler code. Also this funcion is small portion of runtime.
			{
				if (P.DoInterventionDelaysByAdUnit)
					DoOrDontAmendStartTime(&AdUnits[i].HQuarantineTimeStart, t + AdUnits[i].HQuarantineDelay);
				else
					DoOrDontAmendStartTime(&AdUnits[i].HQuarantineTimeStart, t + P.HQuarantineTimeStartBase);
			}

		//// Set DigitalContactTracingTimeStart
		for (i = 0; i < P.NumAdunits; i++)
			if (ChooseTriggerVariableAndValue(i) > ChooseThreshold(i, P.DigitalContactTracing_CellIncThresh)) //// a little wasteful if doing Global trigs as function called more times than necessary, but worth it for much simpler code. Also this funcion is small portion of runtime.
				DoOrDontAmendStartTime(&AdUnits[i].DigitalContactTracingTimeStart	, t + P.DigitalContactTracingTimeStartBase);

		if (P.DoGlobalTriggers)
		{
			int TriggerValue = ChooseTriggerVariableAndValue(0); 
			if (TriggerValue >= ChooseThreshold(0, P.TreatCellIncThresh))
				DoOrDontAmendStartTime(&(P.TreatTimeStart), t + P.TreatTimeStartBase);
			if (TriggerValue >= P.VaccCellIncThresh) DoOrDontAmendStartTime(&P.VaccTimeStart, t + P.VaccTimeStartBase);
			if (TriggerValue >= P.SocDistCellIncThresh)
			{
				DoOrDontAmendStartTime(&P.SocDistTimeStart, t + P.SocDistTimeStartBase);
				//added this for admin unit based intervention delays based on a global trigger: ggilani 17/03/20
				if (P.DoInterventionDelaysByAdUnit)
					for (i = 0; i < P.NumAdunits; i++) 
						DoOrDontAmendStartTime(&AdUnits[i].SocialDistanceTimeStart, t + AdUnits[i].SocialDistanceDelay);
			}
			if (TriggerValue >= P.PlaceCloseCellIncThresh)
			{
				DoOrDontAmendStartTime(&P.PlaceCloseTimeStart, t + P.PlaceCloseTimeStartBase);
				if ((P.PlaceCloseTimeStart2 >= 1e10) && (t >= P.PlaceCloseDuration + P.PlaceCloseTimeStart))
				{
					P.PlaceCloseTimeStart	= t + P.PlaceCloseTimeStartBase2 - P.PlaceCloseTimeStartBase;
					P.PlaceCloseDuration	= P.PlaceCloseDuration2;
				}
				if (P.DoInterventionDelaysByAdUnit)
					for (i = 0; i < P.NumAdunits; i++)
						DoOrDontAmendStartTime(&AdUnits[i].PlaceCloseTimeStart, t + AdUnits[i].PlaceCloseDelay);
			}
			if (TriggerValue >= P.MoveRestrCellIncThresh)
				DoOrDontAmendStartTime(&P.MoveRestrTimeStart, t + P.MoveRestrTimeStartBase);
			if (TriggerValue >= P.KeyWorkerProphCellIncThresh)
				DoOrDontAmendStartTime(&P.KeyWorkerProphTimeStart, t + P.KeyWorkerProphTimeStartBase);
		}
		else
		{
			DoOrDontAmendStartTime(&P.TreatTimeStart			, t + P.TreatTimeStartBase			);
			DoOrDontAmendStartTime(&P.VaccTimeStart				, t + P.VaccTimeStartBase			);
			DoOrDontAmendStartTime(&P.SocDistTimeStart			, t + P.SocDistTimeStartBase		);
			DoOrDontAmendStartTime(&P.PlaceCloseTimeStart		, t + P.PlaceCloseTimeStartBase		);
			DoOrDontAmendStartTime(&P.MoveRestrTimeStart		, t + P.MoveRestrTimeStartBase		);
			DoOrDontAmendStartTime(&P.KeyWorkerProphTimeStart	, t + P.KeyWorkerProphTimeStartBase	);
		}
		DoOrDontAmendStartTime(&P.AirportCloseTimeStart, t + P.AirportCloseTimeStartBase);
	}
	if ((P.PlaceCloseIndepThresh > 0) && (((double)State.cumDC) >= P.PlaceCloseIndepThresh))
		DoOrDontAmendStartTime(&P.PlaceCloseTimeStart, t + P.PlaceCloseTimeStartBase);
	if (t > P.SocDistTimeStart + P.SocDistChangeDelay)
	{
		P.SocDistDurationC = P.SocDistDuration2;
		P.SocDistHouseholdEffectC = P.SocDistHouseholdEffect2;
		P.SocDistSpatialEffectC=P.SocDistSpatialEffect2;
		P.ESocDistHouseholdEffectC = P.ESocDistHouseholdEffect2;
		P.ESocDistSpatialEffectC = P.ESocDistSpatialEffect2;
		for (i = 0; i < P.PlaceTypeNum; i++)
		{
			P.SocDistPlaceEffectC[i] = P.SocDistPlaceEffect2[i];
			P.ESocDistPlaceEffectC[i] = P.ESocDistPlaceEffect2[i];
		}
	}
	




	if (P.OutputBitmap >= 1)
	{
		TSMean = TSMeanNE; TSVar = TSVarNE;
		CaptureBitmap	(n, 0);
		OutputBitmap	(t, 0);
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
		for (i = 0; i < P.N; i++)
		{
			//				i=Cells[b].members[c];
			if (j == 0) j = k = Households[Hosts[i].hh].nh;
			if ((Hosts[i].inf != InfStat_Susceptible) && (Hosts[i].inf != InfStat_ImmuneAtStart))
			{
				if (Hosts[i].latent_time * P.TimeStep <= P.SampleTime)
					TimeSeries[(int)(Hosts[i].latent_time * P.TimeStep / P.SampleStep)].Rdenom++;
				infcountry[Mcells[Hosts[i].mcell].country]++;
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
	for (n = 0; n < P.NumSamples; n++)
	{
		TimeSeries[n].t -= P.PreControlClusterIdTime;
		if (TimeSeries[n].t == 0) lc = n;
		s = 0;
		if (TimeSeries[n].Rdenom == 0) TimeSeries[n].Rdenom = 1e-10;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			TimeSeries[n].Rage[i] /= TimeSeries[n].Rdenom;
		for (i = 0; i < INFECT_TYPE_MASK; i++)
			s += (TimeSeries[n].Rtype[i] /= TimeSeries[n].Rdenom);
		TimeSeries[n].Rdenom = s;
	}
	nf = (sizeof(results) - 3 * sizeof(float*)) / sizeof(double);
	if (!P.DoAdUnits) nf -= MAX_ADUNITS;
	fprintf(stderr, "extinct=%i (%i)\n", (int) TimeSeries[P.NumSamples - 1].extinct, P.NumSamples - 1);
	if (TimeSeries[P.NumSamples - 1].extinct)
	{
		TSMean = TSMeanE; TSVar = TSVarE; P.NRactE++;
	}
	else
	{
		TSMean = TSMeanNE; TSVar = TSVarNE; P.NRactNE++;
	}
	s = 0;
	lc -= 70;
	for (n = 0; n < P.NumSamples; n++)
	{
		if ((n + lc >= 0) && (n + lc < P.NumSamples))
		{
			if (s < TimeSeries[n + lc].incC) { s = TimeSeries[n + lc].incC; t = P.SampleStep * ((double)(n + lc)); }
			res = (double*)&TimeSeries[n + lc];
			res_av = (double*)&TSMean[n];
			res_var = (double*)&TSVar[n];
			for (i = 1; i < nf; i++)
			{
				res_av[i] += res[i];
				res_var[i] += res[i] * res[i];
			}
			if (TSMean[n].cumTmax < TimeSeries[n + lc].cumT) TSMean[n].cumTmax = TimeSeries[n + lc].cumT;
			if (TSMean[n].cumVmax < TimeSeries[n + lc].cumV) TSMean[n].cumVmax = TimeSeries[n + lc].cumV;
		}
		TSMean[n].t += ((double)(n - 70)) * P.SampleStep;
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
	int tn, i, j, k, l, m, p;
	double total_flow, flow;
	ptrdiff_t cl_from, cl_to, cl_from_mcl, cl_to_mcl, mcl_from, mcl_to;
	double pop_dens_from[MAX_ADUNITS], pop_dens_to[MAX_ADUNITS];

#pragma omp parallel for private(tn,i,j,k,l,m,p,total_flow,mcl_from,mcl_to,cl_from,cl_to,cl_from_mcl,cl_to_mcl,pop_dens_from,pop_dens_to,flow) schedule(static) //reduction(+:s,t2)
	for (tn = 0; tn < P.NumThreads; tn++)
	{
		for (i = tn; i < P.NCP; i += P.NumThreads)
		{
			//reset pop density matrix to zero
			for (k = 0; k < P.NumAdunits; k++)
			{
				pop_dens_from[k] = 0.0;
			}

			//find index of cell from which flow travels
			cl_from = CellLookup[i] - Cells;
			cl_from_mcl = (cl_from / P.nch) * P.NMCL * P.nmch + (cl_from % P.nch) * P.NMCL;

			//loop over microcells in these cells to find populations in each admin unit and so flows
			for (k = 0; k < P.NMCL; k++)
			{
				for (l = 0; l < P.NMCL; l++)
				{
					//get index of microcell
					mcl_from = cl_from_mcl + l + k * P.nmch;
					if (Mcells[mcl_from].n > 0)
					{
						//get proportion of each population of cell that exists in each admin unit
						pop_dens_from[Mcells[mcl_from].adunit] += (((double)Mcells[mcl_from].n) / ((double)Cells[cl_from].n));
					}
				}
			}

			for (j = i; j < P.NCP; j++)
			{
				//reset pop density matrix to zero
				for (m = 0; m < P.NumAdunits; m++)
				{
					pop_dens_to[m] = 0.0;
				}

				//find index of cell which flow travels to
				cl_to = CellLookup[j] - Cells;
				cl_to_mcl = (cl_to / P.nch) * P.NMCL * P.nmch + (cl_to % P.nch) * P.NMCL;
				//calculate distance and kernel between the cells
				//total_flow=Cells[cl_from].max_trans[j]*Cells[cl_from].n*Cells[cl_to].n;
				if (j == 0)
				{
					total_flow = Cells[cl_from].cum_trans[j] * Cells[cl_from].n;
				}
				else
				{
					total_flow = (Cells[cl_from].cum_trans[j] - Cells[cl_from].cum_trans[j - 1]) * Cells[cl_from].n;
				}

				//loop over microcells within destination cell
				for (m = 0; m < P.NMCL; m++)
				{
					for (p = 0; p < P.NMCL; p++)
					{
						//get index of microcell
						mcl_to = cl_to_mcl + p + m * P.nmch;
						if (Mcells[mcl_to].n > 0)
						{
							//get proportion of each population of cell that exists in each admin unit
							pop_dens_to[Mcells[mcl_to].adunit] += (((double)Mcells[mcl_to].n) / ((double)Cells[cl_to].n));
						}
					}
				}

				for (m = 0; m < P.NumAdunits; m++)
				{
					for (p = 0; p < P.NumAdunits; p++)
					{
						if (m != p)
						{
							if (AdUnits[m].id / P.CountryDivisor == AdUnits[p].id / P.CountryDivisor)
							{
								flow = total_flow * pop_dens_from[m] * pop_dens_to[p];
							}
							else
							{
								flow = total_flow * pop_dens_from[m] * pop_dens_to[p];
							}
							StateT[tn].origin_dest[m][p] += flow;
							StateT[tn].origin_dest[p][m] += flow;
						}
					}
				}
			}

			////loop over microcells within cell to find the proportion of the cell population in each admin unit
			//k=(cl_from/P.nch)*P.NMCL*P.nmch+(cl_from%P.nch)*P.NMCL;
			//for(l=0;l<P.NMCL;l++)
			//{
			//	for(m=0;m<P.NMCL;m++)
			//	{
			//		mcl_from=k+m+l*P.nmch;
			//		pop_cell_from[Mcells[mcl_from].adunit]+=Mcells[mcl_from].n;
			//	}
			//}
			////loop over cells
			//for(p=(i+1);p<P.NCP;p++)
			//{
			//	//reset population array
			//	for(j=0;j<P.NumAdunits;j++)
			//	{
			//		pop_cell_to[j]=0.0;
			//	}
			//	cl_to=CellLookup[p]-Cells;
			//	//loop over microcells within cell to find the proportion of the cell population in each admin unit
			//	q=(cl_to/P.nch)*P.NMCL*P.nmch+(cl_to%P.nch)*P.NMCL;
			//	for(l=0;l<P.NMCL;l++)
			//	{
			//		for(m=0;m<P.NMCL;m++)
			//		{
			//			mcl_to=q+m+l*P.nmch;
			//			pop_cell_to[Mcells[mcl_to].adunit]+=Mcells[mcl_to].n;
			//		}
			//	}

			//	//find distance and kernel function between cells
			//	dist=dist2_cc_min(Cells+cl_from,Cells+cl_to);
			//	dist_kernel=numKernel(dist);

			//	//add flow between adunits based on how population is distributed
			//	for(l=0;l<P.NumAdunits;l++)
			//	{
			//		for(m=(l+1);m<P.NumAdunits;m++)
			//		{
			//			AdUnits[l].origin_dest[m]+=pop_cell_from[l]*pop_cell_to[m]*dist_kernel;
			//			AdUnits[m].origin_dest[l]+=pop_cell_from[l]*pop_cell_to[m]*dist_kernel;
			//		}
			//	}

		}
	}

	//Sum up flow between adunits across threads
	for (i = 0; i < P.NumAdunits; i++)
	{
		for (j = 0; j < P.NumAdunits; j++)
		{
			for (k = 0; k < P.NumThreads; k++)
			{
				AdUnits[i].origin_dest[j] += StateT[k].origin_dest[i][j];
			}
		}
	}

}

//// Get parameters code (called by ReadParams function)
int GetInputParameter (FILE* dat, FILE* dat2, const char* SItemName, const char* ItemType, void* ItemPtr, int NumItem, int NumItem2, int Offset)
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
		FindFlag = (!strncmp(match, ItemName, strlen(match)));
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
				if ((match[0] == '#') && (match[1] == '1'))
				{
					FindFlag++;
					if (n == 1)
						* ((double*)ItemPtr) = P.clP1;
					else if (n == 2)
						* ((int*)ItemPtr) = (int)P.clP1;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if ((match[0] == '#') && (match[1] == '2'))
				{
					FindFlag++;
					if (n == 1)
						* ((double*)ItemPtr) = P.clP2;
					else if (n == 2)
						* ((int*)ItemPtr) = (int)P.clP2;
					else if (n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
				}
				else if((match[0] == '#') && (match[1] == '3'))
					{
					FindFlag++;
					if(n == 1)
						* ((double*)ItemPtr) = P.clP3;
					else if(n == 2)
						* ((int*)ItemPtr) = (int)P.clP3;
					else if(n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
					}
				else if((match[0] == '#') && (match[1] == '4'))
					{
					FindFlag++;
					if(n == 1)
						* ((double*)ItemPtr) = P.clP4;
					else if(n == 2)
						* ((int*)ItemPtr) = (int)P.clP4;
					else if(n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
					}
				else if((match[0] == '#') && (match[1] == '5'))
					{
					FindFlag++;
					if(n == 1)
						* ((double*)ItemPtr) = P.clP5;
					else if(n == 2)
						* ((int*)ItemPtr) = (int)P.clP5;
					else if(n == 3)
						sscanf(match, "%s", (char*)ItemPtr);
					}
				else if((match[0] == '#') && (match[1] == '6'))
					{
					FindFlag++;
					if(n == 1)
						* ((double*)ItemPtr) = P.clP6;
					else if(n == 2)
						* ((int*)ItemPtr) = (int)P.clP6;
					else if(n == 3)
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
					if ((match[0] != '[') && (!feof(dat)))
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
					if ((match[0] != '[') && (!feof(dat)))
					{
						FindFlag++;
						if (n == 1)
							sscanf(match, "%lf", ((double**)ItemPtr)[j + Offset] + i + Offset); //changed from [j+Offset]+i+Offset to +j+Offset+i, as ItemPtr isn't an array - 01/10: changed it back
						else
							sscanf(match, "%i", ((int**)ItemPtr)[j + Offset] + i + Offset);
					}
					else
					{
						i = NumItem2;
						j = NumItem;
					}
				}
				//Offset=Offset+(NumItem2-1); //added this line to get the correct offset in address position when incrementing j
			} //added these braces
		}
	}
	//	fprintf(stderr,"%s\n",SItemName);
	return FindFlag;
}

