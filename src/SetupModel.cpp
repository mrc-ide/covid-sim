#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "binio.h"
#include "Error.h"
#include "Rand.h"
#include "Kernels.h"
#include "Dist.h"
#include "MachineDefines.h"
#include "Param.h"
#include "SetupModel.h"
#include "Model.h"
#include "ModelMacros.h"
#include "SharedFuncs.h"
#include "InfStat.h"
#include "Bitmap.h"

void* BinFileBuf;
bin_file* BF;
int netbuf[NUM_PLACE_TYPES_NOAIR * 1000000];


///// INITIALIZE / SET UP FUNCTIONS
void SetupModel(char* DensityFile, char* NetworkFile, char* SchoolFile, char* RegDemogFile)
{
	int i, j, k, l, m, i1, i2, j2, l2, m2, tn; //added tn as variable for multi-threaded loops: 28/11/14
	int age; //added age (group): ggilani 09/03/20
	unsigned int rn;
	double t, s, s2, s3, x, y, t2, t3, d, q;
	char buf[2048];
	FILE* dat;

	double maxDigitalContactCoverage, digitalContactTracingProbScaling;

	if (!(Xcg1 = (long*)malloc(MAX_NUM_THREADS * CACHE_LINE_SIZE * sizeof(long)))) ERR_CRITICAL("Unable to allocate ranf storage\n");
	if (!(Xcg2 = (long*)malloc(MAX_NUM_THREADS * CACHE_LINE_SIZE * sizeof(long)))) ERR_CRITICAL("Unable to allocate ranf storage\n");
	setall(P.seed1, P.seed2);

	P.DoBin = -1;
	if (P.DoHeteroDensity)
	{
		if (P.DoAdunitBoundaries)
		{
			fprintf(stderr, "Scanning population density file\n");
			if (!(dat = fopen(DensityFile, "rb"))) ERR_CRITICAL("Unable to open density file\n");
			fread_big(&(P.BinFileLen), sizeof(unsigned int), 1, dat);
			if (P.BinFileLen == 0xf0f0f0f0) //code for first 4 bytes of binary file ## NOTE - SHOULD BE LONG LONG TO COPE WITH BIGGER POPULATIONS
			{
				P.DoBin = 1;
				fread_big(&(P.BinFileLen), sizeof(unsigned int), 1, dat);
				if (!(BinFileBuf = (void*)malloc(P.BinFileLen * sizeof(bin_file)))) ERR_CRITICAL("Unable to allocate binary file buffer\n");
				fread_big(BinFileBuf, sizeof(bin_file), (size_t)P.BinFileLen, dat);
				BF = (bin_file*)BinFileBuf;
				fclose(dat);
			}
			else
			{
				P.DoBin = 0;
				fclose(dat);
				if (!(dat = fopen(DensityFile, "r"))) ERR_CRITICAL("Unable to open density file\n");
				P.BinFileLen = UINT_MAX - 1;
			}
			P.SpatialBoundingBox[0] = P.SpatialBoundingBox[1] = 1e10;
			P.SpatialBoundingBox[2] = P.SpatialBoundingBox[3] = -1e10;
			s2 = 0;
			if (P.DoBin == 0)
			{
				fgets(buf, 2047, dat);
				if (feof(dat)) rn = P.BinFileLen;
			}
			for (rn = 0; rn < P.BinFileLen; rn++)
			{
				if (P.DoBin == 0)
				{
					sscanf(buf, "%lg %lg %lg %i %i", &x, &y, &t, &i2, &l);
					if (l / P.CountryDivisor != i2) //temporarily changed this to 10000 from 100 to work with new admin codes - ggilani 30/05/2018, now changed to CountryDivisor - ggilani 13/05/2019
					{
						//fprintf(stderr,"# %lg %lg %lg %i %i\n",x,y,t,i2,l);
					}

					fgets(buf, 2047, dat);
					if (feof(dat)) rn = P.BinFileLen;
				}
				else
				{
					x = BF[rn].x;
					y = BF[rn].y;
					t = BF[rn].pop;
					i2 = BF[rn].cnt;
					l = BF[rn].ad;
					//					fprintf(stderr,"# %lg %lg %lg %i\t",x,y,t,l);

				}
				m = (l % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor;
				if (P.AdunitLevel1Lookup[m] >= 0)
					if (AdUnits[P.AdunitLevel1Lookup[m]].id / P.AdunitLevel1Mask == l / P.AdunitLevel1Mask)
					{
						AdUnits[P.AdunitLevel1Lookup[m]].cnt_id = i2;
						s2 += t;
						if (x < P.SpatialBoundingBox[0]) P.SpatialBoundingBox[0] = x;
						if (x >= P.SpatialBoundingBox[2]) P.SpatialBoundingBox[2] = x + 1e-6;
						if (y < P.SpatialBoundingBox[1]) P.SpatialBoundingBox[1] = y;
						if (y >= P.SpatialBoundingBox[3]) P.SpatialBoundingBox[3] = y + 1e-6;
					}
			}
			if (!P.DoSpecifyPop) P.N = (int)s2;
			if (P.DoBin == 0) fclose(dat);
		}

#ifdef COUNTRY_THAILAND
		P.width = P.SpatialBoundingBox[2] - P.SpatialBoundingBox[0];
		P.height = P.SpatialBoundingBox[3] - P.SpatialBoundingBox[1];
		P.ncw = (int)(P.width / P.cwidth);
		P.nch = (int)(P.height / P.cwidth);
		P.cwidth = P.width / ((double)P.ncw);
		P.cheight = P.height / ((double)P.nch);
#else
		P.cheight = P.cwidth;
		P.SpatialBoundingBox[0] = floor(P.SpatialBoundingBox[0] / P.cwidth) * P.cwidth;
		P.SpatialBoundingBox[1] = floor(P.SpatialBoundingBox[1] / P.cheight) * P.cheight;
		P.SpatialBoundingBox[2] = ceil(P.SpatialBoundingBox[2] / P.cwidth) * P.cwidth;
		P.SpatialBoundingBox[3] = ceil(P.SpatialBoundingBox[3] / P.cheight) * P.cheight;
		P.width = P.SpatialBoundingBox[2] - P.SpatialBoundingBox[0];
		P.height = P.SpatialBoundingBox[3] - P.SpatialBoundingBox[1];
		P.ncw = 4 * ((int)ceil(P.width / P.cwidth / 4));
		P.nch = 4 * ((int)ceil(P.height / P.cheight / 4));
		P.width = ((double)P.ncw) * P.cwidth;
		P.height = ((double)P.nch) * P.cheight;
		P.SpatialBoundingBox[2] = P.SpatialBoundingBox[0] + P.width;
		P.SpatialBoundingBox[3] = P.SpatialBoundingBox[1] + P.height;
#endif
		P.NC = P.ncw * P.nch;
		fprintf(stderr, "Adjusted bounding box = (%lg, %lg)- (%lg, %lg)\n", P.SpatialBoundingBox[0], P.SpatialBoundingBox[1], P.SpatialBoundingBox[2], P.SpatialBoundingBox[3]);
		fprintf(stderr, "Number of cells = %i (%i x %i)\n", P.NC, P.ncw, P.nch);
		fprintf(stderr, "Population size = %i \n", P.N);
		s = 1;
		P.DoPeriodicBoundaries = 0;
	}
	else
	{
		P.ncw = P.nch = (int)sqrt((double)P.NC);
		P.NC = P.ncw * P.nch;
		fprintf(stderr, "Number of cells adjusted to be %i (%i^2)\n", P.NC, P.ncw);
		s = floor(sqrt((double)P.N));
		P.SpatialBoundingBox[0] = P.SpatialBoundingBox[1] = 0;
		P.SpatialBoundingBox[2] = P.SpatialBoundingBox[3] = s;
		P.N = (int)(s * s);
		fprintf(stderr, "Population size adjusted to be %i (%lg^2)\n", P.N, s);
		P.width = P.height = s;
		P.cwidth = P.width / ((double)P.ncw);
		P.cheight = P.height / ((double)P.nch);
	}
	P.NMC = P.NMCL * P.NMCL * P.NC;
	P.nmcw = P.ncw * P.NMCL;
	P.nmch = P.nch * P.NMCL;
	fprintf(stderr, "Number of micro-cells = %i\n", P.NMC);
	P.scalex = P.BitmapScale;
	P.scaley = P.BitmapAspectScale * P.BitmapScale;
	P.bwidth = (int)(P.width * (P.BoundingBox[2] - P.BoundingBox[0]) * P.scalex);
	P.bwidth = (P.bwidth + 3) / 4;
	P.bwidth *= 4;
	P.bheight = (int)(P.height * (P.BoundingBox[3] - P.BoundingBox[1]) * P.scaley);
	P.bheight += (4 - P.bheight % 4) % 4;
	P.bheight2 = P.bheight + 20; // space for colour legend
	fprintf(stderr, "Bitmap width = %i\nBitmap height = %i\n", P.bwidth, P.bheight);
	P.bminx = (int)(P.width * P.BoundingBox[0] * P.scalex);
	P.bminy = (int)(P.height * P.BoundingBox[1] * P.scaley);
	P.mcwidth = P.cwidth / ((double)P.NMCL);
	P.mcheight = P.cheight / ((double)P.NMCL);
	for (i = 0; i < P.NumSeedLocations; i++)
	{
		P.LocationInitialInfection[i][0] -= P.SpatialBoundingBox[0];
		P.LocationInitialInfection[i][1] -= P.SpatialBoundingBox[1];
	}
	t = dist2_raw(0, 0, P.width, P.height);
	if (P.DoPeriodicBoundaries) t *= 0.25;
	if (!(nKernel = (double*)calloc(NKR + 1, sizeof(double)))) ERR_CRITICAL("Unable to allocate kernel storage\n");
	if (!(nKernelHR = (double*)calloc(NKR + 1, sizeof(double)))) ERR_CRITICAL("Unable to allocate kernel storage\n");
	P.KernelDelta = t / NKR;
	//	fprintf(stderr,"** %i %lg %lg %lg %lg | %lg %lg %lg %lg \n",P.DoUTM_coords,P.SpatialBoundingBox[0],P.SpatialBoundingBox[1],P.SpatialBoundingBox[2],P.SpatialBoundingBox[3],P.width,P.height,t,P.KernelDelta);
	fprintf(stderr, "Coords xmcell=%lg m   ymcell = %lg m\n", sqrt(dist2_raw(P.width / 2, P.height / 2, P.width / 2 + P.mcwidth, P.height / 2)), sqrt(dist2_raw(P.width / 2, P.height / 2, P.width / 2, P.height / 2 + P.mcheight)));
	P.KernelShape = P.MoveKernelShape;
	P.KernelScale = P.MoveKernelScale;
	P.KernelP3 = P.MoveKernelP3;
	P.KernelP4 = P.MoveKernelP4;
	P.KernelType = P.MoveKernelType;
	t2 = 0.0;

	SetupPopulation(DensityFile, SchoolFile, RegDemogFile);
	if (!(TimeSeries = (results*)calloc(P.NumSamples, sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if (!(TSMeanE = (results*)calloc(P.NumSamples, sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if (!(TSVarE = (results*)calloc(P.NumSamples, sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if (!(TSMeanNE = (results*)calloc(P.NumSamples, sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	if (!(TSVarNE = (results*)calloc(P.NumSamples, sizeof(results)))) ERR_CRITICAL("Unable to allocate results storage\n");
	TSMean = TSMeanE; TSVar = TSVarE;

	///// This loops over index l twice just to reset the pointer TSMean from TSMeanE to TSMeanNE (same for TSVar).
	for (l = 0; l < 2; l++)
	{
		for (i = 0; i < P.NumSamples; i++)
		{
			TSMean[i].S = TSMean[i].I = TSMean[i].R = TSMean[i].D = TSMean[i].L =
				TSMean[i].incL = TSMean[i].incI = TSMean[i].incR = TSMean[i].incC = TSMean[i].incDC = TSMean[i].cumDC =
				TSMean[i].incTC = TSMean[i].cumT = TSMean[i].cumTP = TSMean[i].cumUT = TSMean[i].cumV = TSMean[i].H = TSMean[i].incH =
				TSMean[i].incCT = TSMean[i].CT = TSMean[i].incCC = TSMean[i].incDCT = TSMean[i].DCT = //added contact tracing, cases who are contactsTSMean[i].cumTmax = TSMean[i].cumVmax = TSMean[i].incD = TSMean[i].incHQ = TSMean[i].incAC =
				TSMean[i].incAH = TSMean[i].incAA = TSMean[i].incACS = TSMean[i].incAPC =
				TSMean[i].incAPA = TSMean[i].incAPCS = TSMean[i].Rdenom = 0;
			TSVar[i].S = TSVar[i].I = TSVar[i].R = TSVar[i].D = TSVar[i].L =
				TSVar[i].incL = TSVar[i].incI = TSVar[i].incR = TSVar[i].incC = TSVar[i].incTC = TSVar[i].incD = TSVar[i].H = TSVar[i].incH = TSVar[i].incCT = TSVar[i].CT = TSVar[i].incCC = TSMean[i].incDCT = TSVar[i].DCT = 0;
			for (j = 0; j < NUM_PLACE_TYPES; j++) TSMean[i].PropPlacesClosed[j] = TSVar[i].PropPlacesClosed[j] = 0;
			for (j = 0; j < INFECT_TYPE_MASK; j++) TSMean[i].incItype[j] = TSMean[i].Rtype[j] = 0;
			for (j = 0; j < NUM_AGE_GROUPS; j++) TSMean[i].incCa[j] = TSMean[i].incIa[j] = TSMean[i].incDa[j] = TSMean[i].Rage[j] = 0;
			for (j = 0; j < 2; j++)
				TSMean[i].incI_keyworker[j] = TSVar[i].incI_keyworker[j] =
				TSMean[i].incC_keyworker[j] = TSVar[i].incC_keyworker[j] =
				TSMean[i].cumT_keyworker[j] = TSVar[i].cumT_keyworker[j] = 0;
			if (P.DoAdUnits)
				for (j = 0; j <= P.NumAdunits; j++)
					TSMean[i].incI_adunit[j] = TSVar[i].incI_adunit[j] =
					TSMean[i].incC_adunit[j] = TSVar[i].incC_adunit[j] =
					TSMean[i].incD_adunit[j] = TSVar[i].incD_adunit[j] =
					TSMean[i].incDC_adunit[j] = TSVar[i].incDC_adunit[j] =//added detected cases here: ggilani 03/02/15
					TSMean[i].incH_adunit[j] = TSVar[i].incH_adunit[j] =
					TSMean[i].incCT_adunit[j] = TSVar[i].incCT_adunit[j] = //added contact tracing
					TSMean[i].incCC_adunit[j] = TSVar[i].incCC_adunit[j] = //added cases who are contacts: ggilani 28/05/2019
					TSMean[i].incDCT_adunit[j] = TSVar[i].incDCT_adunit[j] = //added digital contact tracing: ggilani 11/03/20
					TSMean[i].cumT_adunit[j] = TSVar[i].cumT_adunit[j] = 0;

			if (P.DoSeverity)
			{
				//// TSMean (each severity for prevalence, incidence and cumulative incidence)
				TSMean[i].Mild = TSMean[i].ILI = TSMean[i].SARI = TSMean[i].Critical = TSMean[i].CritRecov =
					TSMean[i].incMild = TSMean[i].incILI = TSMean[i].incSARI = TSMean[i].incCritical = TSMean[i].incCritRecov =
					TSMean[i].cumMild = TSMean[i].cumILI = TSMean[i].cumSARI = TSMean[i].cumCritical = TSMean[i].cumCritRecov = 0;

				//// TSVar (each severity for prevalence, incidence and cumulative incidence)
				TSVar[i].Mild = TSVar[i].ILI = TSVar[i].SARI = TSVar[i].Critical = TSVar[i].CritRecov =
					TSVar[i].incMild = TSVar[i].incILI = TSVar[i].incSARI = TSVar[i].incCritical = TSVar[i].incCritRecov =
					TSVar[i].cumMild = TSVar[i].cumILI = TSVar[i].cumSARI = TSVar[i].cumCritical = TSVar[i].cumCritRecov = 0;

				//// TSMean admin unit (each severity for prevalence, incidence and cumulative incidence by admin unit)
				if (P.DoAdUnits)
					for (j = 0; j <= P.NumAdunits; j++)
						TSMean[i].Mild_adunit[j] = TSMean[i].ILI_adunit[j] = TSMean[i].SARI_adunit[j] = TSMean[i].Critical_adunit[j] = TSMean[i].CritRecov_adunit[j] =
						TSMean[i].incMild_adunit[j] = TSMean[i].incILI_adunit[j] = TSMean[i].incSARI_adunit[j] = TSMean[i].incCritical_adunit[j] = TSMean[i].incCritRecov_adunit[j] =
						TSMean[i].cumMild_adunit[j] = TSMean[i].cumILI_adunit[j] = TSMean[i].cumSARI_adunit[j] = TSMean[i].cumCritical_adunit[j] = TSMean[i].cumCritRecov_adunit[j] = 0;
			}
		}
		TSMean = TSMeanNE; TSVar = TSVarNE;
	}

	//added memory allocation and initialisation of infection event log, if DoRecordInfEvents is set to 1: ggilani - 10/10/2014
	if (P.DoRecordInfEvents)
	{
		if (!(InfEventLog = (events*)calloc(P.MaxInfEvents, sizeof(events)))) ERR_CRITICAL("Unable to allocate events storage\n");
		if (!(nEvents = (int*)calloc(1, sizeof(int)))) ERR_CRITICAL("Unable to allocate events storage\n");
	}

	P.CellPop2 = ((double)P.N) * ((double)P.N) / (((double)P.NC) * ((double)P.NC));
	SaveAgeDistrib();

	fprintf(stderr, "Initialising kernel...\n");
	InitKernel(0, 1.0);	fprintf(stderr, "Initialising places...\n");
	if (P.DoPlaces)
	{
		if (P.LoadSaveNetwork == 1)
			LoadPeopleToPlaces(NetworkFile);
		else
			AssignPeopleToPlaces();
	}
	if ((P.DoPlaces) && (P.LoadSaveNetwork == 2))
		SavePeopleToPlaces(NetworkFile);
	//SaveDistribs();

	//setall(P.seed3,P.seed4);

	StratifyPlaces();
	for (i = 0; i < P.NC; i++)
	{
		Cells[i].S = Cells[i].n;
		Cells[i].L = Cells[i].I = Cells[i].R = 0;
		//Cells[i].susceptible=Cells[i].members; //added this line
	}
	for (i = 0; i < P.N; i++) Hosts[i].keyworker = 0;
	P.KeyWorkerNum = P.KeyWorkerIncHouseNum = m = l = 0;
	if (P.DoPlaces)
	{
		while ((m < P.KeyWorkerPopNum) && (l < 1000))
		{
			i = (int)(((double)P.N) * ranf_mt(0));
			if (Hosts[i].keyworker)
				l++;
			else
			{
				Hosts[i].keyworker = 1;
				m++;
				P.KeyWorkerNum++;
				P.KeyWorkerIncHouseNum++;
				l = 0;
				if (ranf_mt(0) < P.KeyWorkerHouseProp)
				{
					l2 = Households[Hosts[i].hh].FirstPerson;
					m2 = l2 + Households[Hosts[i].hh].nh;
					for (j2 = l2; j2 < m2; j2++)
						if (!Hosts[j2].keyworker)
						{
							Hosts[j2].keyworker = 1;
							P.KeyWorkerIncHouseNum++;
						}
				}
			}
		}
		for (j = 0; j < NUM_PLACE_TYPES_NOAIR; j++)
		{
			m = l = 0;
			while ((m < P.KeyWorkerPlaceNum[j]) && (l < 1000))
			{
				k = (int)(((double)P.Nplace[j]) * ranf_mt(0));
				for (i2 = 0; (m < P.KeyWorkerPlaceNum[j]) && (i2 < Places[j][k].n); i2++)
				{
					i = Places[j][k].members[i2];
					if ((i < 0) || (i >= P.N)) fprintf(stderr, "## %i # ", i);
					if ((Hosts[i].keyworker) || (ranf_mt(0) >= P.KeyWorkerPropInKeyPlaces[j]))
						l++;
					else
					{
						Hosts[i].keyworker = 1;
						m++;
						P.KeyWorkerNum++;
						P.KeyWorkerIncHouseNum++;
						l = 0;
						l2 = Households[Hosts[i].hh].FirstPerson;
						m2 = l2 + Households[Hosts[i].hh].nh;
						for (j2 = l2; j2 < m2; j2++)
							if ((!Hosts[j2].keyworker) && (ranf_mt(0) < P.KeyWorkerHouseProp))
							{
								Hosts[j2].keyworker = 1;
								P.KeyWorkerIncHouseNum++;
							}
					}
				}
			}
		}
		if (P.KeyWorkerNum > 0) fprintf(stderr, "%i key workers selected in total\n", P.KeyWorkerNum);
		if (P.DoAdUnits)
		{
			for (i = 0; i < P.NumAdunits; i++) AdUnits[i].NP = 0;
			for (j = 0; j < P.PlaceTypeNum; j++)
				if (P.PlaceCloseAdunitPlaceTypes[j] > 0)
				{
					for (k = 0; k < P.Nplace[j]; k++)
						AdUnits[Mcells[Places[j][k].mcell].adunit].NP++;
				}

		}
	}
	fprintf(stderr, "Places intialised.\n");


	//Set up the population for digital contact tracing here... - ggilani 09/03/20
	if (P.DoDigitalContactTracing)
	{

		//if clustering by Households
		if (P.DoHouseholds && P.ClusterDigitalContactUsers)
		{
			//Loop through households
			for (i = 0; i < P.NH; i++)
			{
				if (ranf() < P.PropPopUsingDigitalContactTracing)
				{
					//select this household for digital contact app use
					//loop through household members and check whether they will be selected for use
					i1 = Households[i].FirstPerson;
					i2 = i1 + Households[i].nh;
					for (j = i1; j < i2; j++)
					{
						//get age of host
						age = HOST_AGE_GROUP(j);
						if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;
						//check to see if host will be a user based on age group
						if (ranf() < P.ProportionSmartphoneUsersByAge[age])
						{
							Hosts[j].digitalContactTracingUser = 1;
							P.NDigitalContactUsers++;
							if (!(Hosts[j].contacts = (int*)malloc(MAX_CONTACTS * sizeof(int)))) ERR_CRITICAL("Unable to allocate contacts storage\n");
						}
					}
					P.NDigitalHouseholdUsers++;
				}
			}
			fprintf(stderr, "Number of digital contact tracing households: %i, out of total number of households: %i\n", P.NDigitalHouseholdUsers, P.NH);
			fprintf(stderr, "Number of digital contact tracing users: %i, out of population size: %i\n", P.NDigitalContactUsers, P.N);
		}
		else // Just go through the population and assign people to the digital contact tracing app based on probability by age.
		{
			//Based on the the age distribution of the population and the proportion of smart phone users per age band, there is a theoretical upper maximum proportion of the total population who can be covered by the digital app;
			//First calculate the theoretical maximum
			maxDigitalContactCoverage = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
			{
				maxDigitalContactCoverage += P.PropAgeGroup[0][i] * P.ProportionSmartphoneUsersByAge[i];
			}
			// if proportion of coverage required by user input exceeds maximum theoretical value, display error
			if (P.PropPopUsingDigitalContactTracing > maxDigitalContactCoverage)
			{
				ERR_CRITICAL("Proportion of app users selected is greater than maximum theoretical value.\n");
			}
			else
			{
				//for use with non-clustered
				digitalContactTracingProbScaling = P.PropPopUsingDigitalContactTracing / maxDigitalContactCoverage;
			}

			for (i = 0; i < P.N; i++)
			{
				age = HOST_AGE_GROUP(i);
				if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;

				if (ranf() < (P.ProportionSmartphoneUsersByAge[age] * digitalContactTracingProbScaling))
				{
					Hosts[i].digitalContactTracingUser = 1;
					P.NDigitalContactUsers++;
					if (!(Hosts[i].contacts = (int*)malloc(MAX_CONTACTS * sizeof(int)))) ERR_CRITICAL("Unable to allocate contacts storage\n");
				}
			}
		}
		fprintf(stderr, "Number of digital contact tracing users: %i, out of population size: %i\n", P.NDigitalContactUsers, P.N);
	}


	UpdateProbs(0);
	if (P.DoAirports) SetupAirports();
	if (P.R0scale != 1.0)
	{
		P.HouseholdTrans *= P.R0scale;
		P.R0 *= P.R0scale;
		for (j = 0; j < P.PlaceTypeNum; j++)
			P.PlaceTypeTrans[j] *= P.R0scale;
		fprintf(stderr, "Rescaled transmission coefficients by factor of %lg\n", P.R0scale);
	}
	t = s = t2 = 0;
	for (i = 0; i < MAX_HOUSEHOLD_SIZE; i++)
	{
		t += ((double)(i + 1)) * (P.HouseholdSizeDistrib[0][i] - t2);
		t2 = P.HouseholdSizeDistrib[0][i];
	}
	t2 = s = 0;
	s3 = 1.0;

#pragma omp parallel for private(i,s2,j,k,q,l,d,y,m,tn) schedule(static,1) reduction(+:s,t2) //schedule(static,1000)
	for (tn = 0; tn < P.NumThreads; tn++) //changed this looping to allow for multi-threaded random numbers
	{
		for (i = tn; i < P.N; i += P.NumThreads)
		{
			if (P.InfectiousnessSD == 0)
				Hosts[i].infectiousness = (float)P.AgeInfectiousness[HOST_AGE_GROUP(i)];
			else
				Hosts[i].infectiousness = (float)(P.AgeInfectiousness[HOST_AGE_GROUP(i)] * gen_gamma_mt(P.InfectiousnessGamA, P.InfectiousnessGamR, tn)); //made this multi-threaded: 28/11/14
			q = P.ProportionSymptomatic[HOST_AGE_GROUP(i)];
			if (ranf_mt(tn) < q) //made this multi-threaded: 28/11/14
				Hosts[i].infectiousness = (float)(-P.SymptInfectiousness * Hosts[i].infectiousness);
			j = (int)floor((q = ranf_mt(tn) * CDF_RES)); //made this multi-threaded: 28/11/14
			q -= ((double)j);
			Hosts[i].recovery_time = (unsigned short int) floor(0.5 - (P.InfectiousPeriod * log(q * P.infectious_icdf[j + 1] + (1.0 - q) * P.infectious_icdf[j]) / P.TimeStep));

			if (P.DoHouseholds)
			{
				s2 = P.TimeStep * P.HouseholdTrans * fabs(Hosts[i].infectiousness) * P.HouseholdDenomLookup[Households[Hosts[i].hh].nhr - 1];
				d = 1.0; l = (int)Hosts[i].recovery_time;
				for (k = 0; k < l; k++) { y = 1.0 - s2 * P.infectiousness[k]; d *= ((y < 0) ? 0 : y); }
				l = Households[Hosts[i].hh].FirstPerson;
				m = l + Households[Hosts[i].hh].nh;
				for (k = l; k < m; k++) if ((Hosts[k].inf == InfStat_Susceptible) && (k != i)) s += (1 - d) * P.AgeSusceptibility[HOST_AGE_GROUP(i)];
			}
			q = (P.LatentToSymptDelay > Hosts[i].recovery_time * P.TimeStep) ? Hosts[i].recovery_time * P.TimeStep : P.LatentToSymptDelay;
			s2 = fabs(Hosts[i].infectiousness) * P.RelativeSpatialContact[HOST_AGE_GROUP(i)] * P.TimeStep;
			l = (int)(q / P.TimeStep);
			for (k = 0; k < l; k++) t2 += s2 * P.infectiousness[k];
			s2 *= ((Hosts[i].infectiousness < 0) ? P.SymptSpatialContactRate : 1);
			l = (int)Hosts[i].recovery_time;
			for (; k < l; k++) t2 += s2 * P.infectiousness[k];

		}
	}
	t2 *= (s3 / ((double)P.N));
	s /= ((double)P.N);
	fprintf(stderr, "Household mean size=%lg\nHousehold R0=%lg\n", t, P.R0household = s);
	t = x = y = 0;
	if (P.DoPlaces)
		for (j = 0; j < P.PlaceTypeNum; j++)
			if (j != HOTEL_PLACE_TYPE)
			{
#pragma omp parallel for private(i,k,d,q,s2,s3,t3,l,m,x,y) schedule(static,1000) reduction(+:t)
				for (i = 0; i < P.N; i++)
				{
					k = Hosts[i].PlaceLinks[j];
					if (k >= 0)
					{
						q = (P.LatentToSymptDelay > Hosts[i].recovery_time * P.TimeStep) ? Hosts[i].recovery_time * P.TimeStep : P.LatentToSymptDelay;
						s2 = fabs(Hosts[i].infectiousness) * P.TimeStep * P.PlaceTypeTrans[j];
						x = s2 / P.PlaceTypeGroupSizeParam1[j];
						d = 1.0; l = (int)(q / P.TimeStep);
						for (m = 0; m < l; m++) { y = 1.0 - x * P.infectiousness[m]; d *= ((y < 0) ? 0 : y); }
						s3 = ((double)(Places[j][k].group_size[Hosts[i].PlaceGroupLinks[j]] - 1));
						x *= ((Hosts[i].infectiousness < 0) ? (P.SymptPlaceTypeContactRate[j] * (1 - P.SymptPlaceTypeWithdrawalProp[j])) : 1);
						l = (int)Hosts[i].recovery_time;
						for (; m < l; m++) { y = 1.0 - x * P.infectiousness[m]; d *= ((y < 0) ? 0 : y); }

						t3 = d;
						x = P.PlaceTypePropBetweenGroupLinks[j] * s2 / ((double)Places[j][k].n);
						d = 1.0; l = (int)(q / P.TimeStep);
						for (m = 0; m < l; m++) { y = 1.0 - x * P.infectiousness[m]; d *= ((y < 0) ? 0 : y); }
						x *= ((Hosts[i].infectiousness < 0) ? (P.SymptPlaceTypeContactRate[j] * (1 - P.SymptPlaceTypeWithdrawalProp[j])) : 1);
						l = (int)Hosts[i].recovery_time;
						for (; m < l; m++) { y = 1.0 - x * P.infectiousness[m]; d *= ((y < 0) ? 0 : y); }
						t += (1 - t3 * d) * s3 + (1 - d) * (((double)(Places[j][k].n - 1)) - s3);
					}
				}
				fprintf(stderr, "%lg  ", t / ((double)P.N));
			}
#pragma omp parallel for private(i) schedule(static,500) reduction(+:x,y)
	for (i = 0; i < P.N; i++)
	{
		x += Hosts[i].recovery_time * P.TimeStep;
		y += Hosts[i].recovery_time;
		Hosts[i].recovery_time = 0;
	}
	t /= ((double)P.N);
	x /= ((double)P.N);
	y /= ((double)P.N);
	fprintf(stderr, "R0 for places = %lg\nR0 for random spatial = %lg\nOverall R0=%lg\n", P.R0places = t, P.R0spatial = P.R0 - s - t, P.R0);
	fprintf(stderr, "Mean infectious period (sampled) = %lg (%lg)\n", x, y);
	if (P.DoSI)
		P.LocalBeta = (P.R0 / t2 - s - t);
	else
		P.LocalBeta = (P.R0 - s - t) / t2;
	if ((P.LocalBeta < 0) || (!P.DoSpatial))
	{
		P.LocalBeta = P.R0spatial = 0;
		fprintf(stderr, "Reset spatial R0 to 0\n");
	}
	fprintf(stderr, "LocalBeta = %lg\n", P.LocalBeta);
	TSMean = TSMeanNE; TSVar = TSVarNE;
	fprintf(stderr, "Calculated approx cell probabilities\n");
	for (i = 0; i < INFECT_TYPE_MASK; i++) inftype_av[i] = 0;
	for (i = 0; i < MAX_COUNTRIES; i++) infcountry_av[i] = infcountry_num[i] = 0;
	for (i = 0; i < MAX_SEC_REC; i++)
		for (j = 0; j < MAX_GEN_REC; j++)
			indivR0_av[i][j] = 0;
	for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++)
		for (j = 0; j <= MAX_HOUSEHOLD_SIZE; j++)
			inf_household_av[i][j] = case_household_av[i][j] = 0;
	fprintf(stderr, "Model configuration complete.\n");
	DoInitUpdateProbs = 1;
	for (i = 0; i < P.NC; i++)	Cells[i].tot_treat = 1;  //This makes sure InitModel intialises the cells.
	P.NRactE = P.NRactNE = 0;
	for (i = 0; i < P.N; i++) Hosts[i].esocdist_comply = (ranf() < P.ESocProportionCompliant[HOST_AGE_GROUP(i)]) ? 1 : 0;
	//	else
	//		{
	//		for(i=0;i<P.N;i++) Hosts[i].esocdist_comply=0;
	//		}
	if (P.OutputBitmap)
	{
		InitBMHead();
	}
	if (P.DoMassVacc)
	{
		if (!(State.mvacc_queue = (int*)calloc(P.N, sizeof(int)))) ERR_CRITICAL("Unable to allocate host storage\n");
		for (i = j = 0; i < P.N; i++)
		{
			if ((HOST_AGE_YEAR(i) >= P.VaccPriorityGroupAge[0]) && (HOST_AGE_YEAR(i) <= P.VaccPriorityGroupAge[1]))
			{
				if (ranf() < P.VaccProp)
					State.mvacc_queue[j++] = i;
			}
		}
		k = j;
		for (i = 0; i < P.N; i++)
		{
			if ((HOST_AGE_YEAR(i) < P.VaccPriorityGroupAge[0]) || (HOST_AGE_YEAR(i) > P.VaccPriorityGroupAge[1]))
			{
				if (ranf() < P.VaccProp)
					State.mvacc_queue[j++] = i;
			}
		}
		State.n_mvacc = j;
		fprintf(stderr, "Number to be vaccinated=%i\n", State.n_mvacc);
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < k; j++)
			{
				l = (int)(ranf() * ((double)k));
				m = State.mvacc_queue[j];
				State.mvacc_queue[j] = State.mvacc_queue[l];
				State.mvacc_queue[l] = m;
			}
			for (j = k; j < State.n_mvacc; j++)
			{
				l = k + ((int)(ranf() * ((double)(State.n_mvacc - k))));
				m = State.mvacc_queue[j];
				State.mvacc_queue[j] = State.mvacc_queue[l];
				State.mvacc_queue[l] = m;
			}
		}
		fprintf(stderr, "Configured mass vaccination queue.\n");
	}
	PeakHeightSum = PeakHeightSS = PeakTimeSum = PeakTimeSS = 0;
	i = (P.ncw / 2) * P.nch + P.nch / 2;
	j = (P.ncw / 2 + 2) * P.nch + P.nch / 2;
	fprintf(stderr, "UTM dist horiz=%lg %lg\n", sqrt(dist2_cc(Cells + i, Cells + j)), sqrt(dist2_cc(Cells + j, Cells + i)));
	j = (P.ncw / 2) * P.nch + P.nch / 2 + 2;
	fprintf(stderr, "UTM dist vert=%lg %lg\n", sqrt(dist2_cc(Cells + i, Cells + j)), sqrt(dist2_cc(Cells + j, Cells + i)));
	j = (P.ncw / 2 + 2) * P.nch + P.nch / 2 + 2;
	fprintf(stderr, "UTM dist diag=%lg %lg\n", sqrt(dist2_cc(Cells + i, Cells + j)), sqrt(dist2_cc(Cells + j, Cells + i)));

	//if(P.OutputBitmap)
	//{
	//	CaptureBitmap(0,0);
	//	OutputBitmap(0.0,0);
	//}

}

void SetupPopulation(char* DensityFile, char* SchoolFile, char* RegDemogFile)
{
	int i, j, k, l, m, i2, j2, last_i, mr, ad, tn, *mcl;
	unsigned int rn, rn2;
	double t, s, x, y, xh, yh, maxd, CumAgeDist[NUM_AGE_GROUPS + 1];
	char buf[4096], *col;
	const char delimiters[] = " \t,";
	FILE* dat = NULL, *dat2;
	bin_file rec;


	if (!(Cells = (cell*)calloc(P.NC, sizeof(cell)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(Mcells = (microcell*)calloc(P.NMC, sizeof(microcell)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(mcell_num = (int*)malloc(P.NMC * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(mcell_dens = (double*)malloc(P.NMC * sizeof(double)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(mcell_country = (int*)malloc(P.NMC * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(mcell_adunits = (int*)malloc(P.NMC * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");

	if (P.DoPlaces)
		if (!(Places = (place * *)malloc(P.PlaceTypeNum * sizeof(place*)))) ERR_CRITICAL("Unable to allocate place storage\n");
	for (j = 0; j < P.NMC; j++)
	{
		Mcells[j].n = 0;
		mcell_adunits[j] = -1;
		mcell_dens[j] = 0;
		mcell_num[j] = mcell_country[j] = 0;
	}
	if (P.DoAdUnits)
		for (i = 0; i < MAX_ADUNITS; i++)
			P.PopByAdunit[i][0] = P.PopByAdunit[i][1] = 0;
	if (P.DoHeteroDensity)
	{
		if (!P.DoAdunitBoundaries) P.NumAdunits = 0;
		if (P.DoBin == -1)
		{
			if (!(dat = fopen(DensityFile, "rb"))) ERR_CRITICAL("Unable to open density file\n");
			fread_big(&(P.BinFileLen), sizeof(unsigned int), 1, dat);
			if (P.BinFileLen == 0xf0f0f0f0) //code for first 4 bytes of binary file
			{
				fprintf(stderr, "Reading binary population density file...\n");
				P.DoBin = 1;
				fread_big(&(P.BinFileLen), sizeof(unsigned int), 1, dat);
				if (!(BinFileBuf = (void*)malloc(((size_t)P.BinFileLen) * sizeof(bin_file)))) ERR_CRITICAL("Unable to allocate binary file buffer\n");
				fread_big(BinFileBuf, sizeof(bin_file), (size_t)P.BinFileLen, dat);
				BF = (bin_file*)BinFileBuf;
				fclose(dat);
			}
			else
			{
				fclose(dat);
				P.DoBin = 0;
			}
		}
		if (P.DoBin == 0)
		{
			P.BinFileLen = UINT_MAX - 1;
			fprintf(stderr, "Reading ASCII population density file...\n");
			if (!(dat = fopen(DensityFile, "r"))) ERR_CRITICAL("Unable to open density file\n");
		}
		//		if(!(dat2=fopen("EnvTest.txt","w"))) ERR_CRITICAL("Unable to open test file\n");
		if (P.DoBin == 1)
			fprintf(stderr, "Binary density file contains %i cells.\n", (int)P.BinFileLen);
		else
		{
			fgets(buf, 2047, dat);
			if (feof(dat)) rn = P.BinFileLen;
		}
		for (rn = rn2 = mr = 0; rn < P.BinFileLen; rn++)
		{
			if (P.DoAdUnits)
			{
				if (P.DoBin == 1)
				{
					x = BF[rn].x; y = BF[rn].y; t = BF[rn].pop; i2 = BF[rn].cnt; j2 = BF[rn].ad; //changed from i to rn to loop over indices properly
					rec = BF[rn];
				}
				else
					sscanf(buf, "%lg %lg %lg %i %i", &x, &y, &t, &i2, &j2);
				m = (j2 % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor;
				if (P.DoAdunitBoundaries)
				{
					if (P.AdunitLevel1Lookup[m] >= 0)
					{
						if (j2 / P.AdunitLevel1Mask == AdUnits[P.AdunitLevel1Lookup[m]].id / P.AdunitLevel1Mask)
						{
							k = 1;
							AdUnits[P.AdunitLevel1Lookup[m]].cnt_id = i2;
						}
						else
							k = 0;
					}
					else
						k = 0;
				}
				else
				{
					k = 1;
					if (P.AdunitLevel1Lookup[m] < 0)
					{
						P.AdunitLevel1Lookup[m] = P.NumAdunits;
						AdUnits[P.NumAdunits].id = j2;
						AdUnits[P.NumAdunits].cnt_id = i2;
						P.NumAdunits++;
						if (P.NumAdunits >= MAX_ADUNITS) ERR_CRITICAL("Total number of administrative units exceeds MAX_ADUNITS\n");
					}
					else
					{
						AdUnits[P.AdunitLevel1Lookup[m]].cnt_id = i2;
					}
				}
			}
			else
			{
				k = 1;
				if (P.DoBin == 1)
				{
					x = BF[i].x; y = BF[i].y; t = BF[i].pop; i2 = BF[i].cnt; j2 = BF[i].ad;
					rec = BF[rn];
				}
				else
				{
					sscanf(buf, "%lg %lg %lg %i", &x, &y, &t, &i2);
					j2 = 0;
					rec.x = x; rec.y = y; rec.pop = t; rec.cnt = i2; rec.ad = j2;
				}
			}
			if (P.DoBin == 0)
			{
				fgets(buf, 2047, dat);
				if (feof(dat)) rn = P.BinFileLen;
			}
			if ((k) && (x >= P.SpatialBoundingBox[0]) && (y >= P.SpatialBoundingBox[1]) && (x < P.SpatialBoundingBox[2]) && (y < P.SpatialBoundingBox[3]))
			{
				j = (int)floor((x - P.SpatialBoundingBox[0]) / P.mcwidth + 0.1);
				k = (int)floor((y - P.SpatialBoundingBox[1]) / P.mcheight + 0.1);
				l = j * P.nmch + k;
				if (l < P.NMC)
				{
					mr++;
					mcell_dens[l] += t;
					mcell_country[l] = i2;
					//fprintf(stderr,"mcell %i, country %i, pop %lg\n",l,i2,t);
					mcell_num[l]++;
					if (P.DoAdUnits)
					{
						mcell_adunits[l] = P.AdunitLevel1Lookup[m];
						if (mcell_adunits[l] < 0) fprintf(stderr, "Cell %i has adunits<0\n", l);
						P.PopByAdunit[P.AdunitLevel1Lookup[m]][0] += t;
					}
					else
						mcell_adunits[l] = 0;
					if ((P.OutputDensFile) && (P.DoBin) && (mcell_adunits[l] >= 0))
					{
						if (rn2 < rn) BF[rn2] = rec;
						rn2++;
					}
				}
			}
		}
		//		fclose(dat2);
		fprintf(stderr, "%i valid cells read from density file.\n", mr);
		if ((P.OutputDensFile) && (P.DoBin)) P.BinFileLen = rn2;
		if (P.DoBin == 0)
		{
			fclose(dat);
			if (P.OutputDensFile)
			{
				P.DoBin = 1;
				P.BinFileLen = 0;
				for (l = 0; l < P.NMC; l++)
					if (mcell_adunits[l] >= 0) P.BinFileLen++;
				if (!(BinFileBuf = (void*)malloc(P.BinFileLen * sizeof(bin_file)))) ERR_CRITICAL("Unable to allocate binary file buffer\n");
				BF = (bin_file*)BinFileBuf;
				fprintf(stderr, "Binary density file should contain %i cells.\n", (int)P.BinFileLen);
				rn = 0;
				for (l = 0; l < P.NMC; l++)
					if (mcell_adunits[l] >= 0)
					{
						BF[rn].x = (double)(P.mcwidth * (((double)(l / P.nmch)) + 0.5)) + P.SpatialBoundingBox[0]; //x
						BF[rn].y = (double)(P.mcheight * (((double)(l % P.nmch)) + 0.5)) + P.SpatialBoundingBox[1]; //y
						BF[rn].ad = (P.DoAdUnits) ? (AdUnits[mcell_adunits[l]].id) : 0;
						BF[rn].pop = mcell_dens[l];
						BF[rn].cnt = mcell_country[l];
						rn++;
					}
			}
		}

		if (P.OutputDensFile)
		{
			if (!(dat2 = fopen(OutDensFile, "wb"))) ERR_CRITICAL("Unable to open output density file\n");
			rn = 0xf0f0f0f0;
			fwrite_big((void*)& rn, sizeof(unsigned int), 1, dat2);
			fprintf(stderr, "Saving population density file with NC=%i...\n", (int)P.BinFileLen);
			fwrite_big((void*) & (P.BinFileLen), sizeof(unsigned int), 1, dat2);
			fwrite_big(BinFileBuf, sizeof(bin_file), (size_t)P.BinFileLen, dat2);
			fclose(dat2);
		}
		if (P.DoBin == 1) free(BinFileBuf);
		fprintf(stderr, "Population files read.\n");
		maxd = 0;
		for (i = 0; i < P.NMC; i++)
		{
			if (mcell_num[i] > 0)
			{
				mcell_dens[i] /= ((double)mcell_num[i]);
				Mcells[i].country = (unsigned short)mcell_country[i];
				if (P.DoAdUnits)
					Mcells[i].adunit = mcell_adunits[i];
				else
					Mcells[i].adunit = 0;
			}
			else
				Mcells[i].adunit = -1;
			maxd += mcell_dens[i];
		}
	}
	else
	{
		for (i = 0; i < P.NMC; i++)
		{
			mcell_dens[i] = 1.0;
			Mcells[i].country = 1;
		}
		maxd = ((double)P.NMC);
	}
	if (!P.DoAdUnits) P.NumAdunits = 1;
	if ((P.DoAdUnits) && (P.DoAdunitDemog))
	{
		if (!(State.InvAgeDist = (int**)malloc(P.NumAdunits * sizeof(int*)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		for (i = 0; i < P.NumAdunits; i++)
			if (!(State.InvAgeDist[i] = (int*)malloc(1000 * sizeof(int)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		if (!(dat = fopen(RegDemogFile, "r"))) ERR_CRITICAL("Unable to open regional demography file\n");
		for (k = 0; k < P.NumAdunits; k++)
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				P.PropAgeGroup[k][i] = 0;
			for (i = 0; i < MAX_HOUSEHOLD_SIZE; i++)
				P.HouseholdSizeDistrib[k][i] = 0;
			P.PopByAdunit[k][1] = 0;
		}
		while (!feof(dat))
		{
			fgets(buf, 2047, dat);
			col = strtok(buf, delimiters);
			sscanf(col, "%i", &l);
			m = (l % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor;
			k = P.AdunitLevel1Lookup[m];
			if (k >= 0)
				if (l / P.AdunitLevel1Mask == AdUnits[k].id / P.AdunitLevel1Mask)
				{
					col = strtok(NULL, delimiters);
					sscanf(col, "%lg", &x);
					P.PopByAdunit[k][1] += x;
					t = 0;
					for (i = 0; i < NUM_AGE_GROUPS; i++)
					{
						col = strtok(NULL, delimiters);
						sscanf(col, "%lg", &s);
						P.PropAgeGroup[k][i] += s;
					}
					col = strtok(NULL, delimiters);
					if (P.DoHouseholds)
					{
						sscanf(col, "%lg", &y);
						for (i = 0; i < MAX_HOUSEHOLD_SIZE; i++)
						{
							col = strtok(NULL, delimiters);
							sscanf(col, "%lg", &s);
							P.HouseholdSizeDistrib[k][i] += y * s;
						}
					}
				}
		}
		fclose(dat);
		for (k = 0; k < P.NumAdunits; k++)
		{
			t = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				t += P.PropAgeGroup[k][i];
			CumAgeDist[0] = 0;
			for (i = 1; i <= NUM_AGE_GROUPS; i++)
			{
				P.PropAgeGroup[k][i - 1] /= t;
				CumAgeDist[i] = CumAgeDist[i - 1] + P.PropAgeGroup[k][i - 1];
			}
			for (i = j = 0; i < 1000; i++)
			{
				t = ((double)i) / 1000;
				while (t >= CumAgeDist[j + 1]) j++;
				t = AGE_GROUP_WIDTH * (((double)j) + (t - CumAgeDist[j]) / (CumAgeDist[j + 1] - CumAgeDist[j]));
				State.InvAgeDist[k][i] = (int)t;
			}
			State.InvAgeDist[k][1000 - 1] = NUM_AGE_GROUPS * AGE_GROUP_WIDTH - 1;
			if (P.DoHouseholds)
			{
				t = 0;
				for (i = 0; i < MAX_HOUSEHOLD_SIZE; i++)
					t += P.HouseholdSizeDistrib[k][i];
				P.HouseholdSizeDistrib[k][0] /= t;
				for (i = 1; i < MAX_HOUSEHOLD_SIZE - 1; i++)
					P.HouseholdSizeDistrib[k][i] = P.HouseholdSizeDistrib[k][i] / t + P.HouseholdSizeDistrib[k][i - 1];
				P.HouseholdSizeDistrib[k][MAX_HOUSEHOLD_SIZE - 1] = 1.0;
			}
			else
			{
				for (i = 0; i < MAX_HOUSEHOLD_SIZE - 1; i++)
					P.HouseholdSizeDistrib[k][i] = 1.0;
			}
		}
	}
	else
	{
		if (!(State.InvAgeDist = (int**)malloc(sizeof(int*)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		if (!(State.InvAgeDist[0] = (int*)malloc(1000 * sizeof(int)))) ERR_CRITICAL("Unable to allocate InvAgeDist storage\n");
		CumAgeDist[0] = 0;
		for (i = 1; i <= NUM_AGE_GROUPS; i++)
			CumAgeDist[i] = CumAgeDist[i - 1] + P.PropAgeGroup[0][i - 1];
		for (i = j = 0; i < 1000; i++)
		{
			t = ((double)i) / 1000;
			if (t >= CumAgeDist[j + 1]) j++;
			t = AGE_GROUP_WIDTH * (((double)j) + (t - CumAgeDist[j]) / (CumAgeDist[j + 1] - CumAgeDist[j]));
			State.InvAgeDist[0][i] = (int)t;
		}
		State.InvAgeDist[0][1000 - 1] = NUM_AGE_GROUPS * AGE_GROUP_WIDTH - 1;
	}
	if ((P.DoAdUnits) && (P.DoAdunitDemog) && (P.DoCorrectAdunitPop))
	{
		for (i = 0; i < P.NumAdunits; i++)
			fprintf(stderr, "%i\t%i\t%lg\t%lg\n", i, (AdUnits[i].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor, P.PropAgeGroup[i][0], P.HouseholdSizeDistrib[i][0]);
		maxd = 0;
		for (i = 0; i < P.NMC; i++)
		{
			if (mcell_num[i] > 0)
				mcell_dens[i] *= P.PopByAdunit[mcell_adunits[i]][1] / (1e-10 + P.PopByAdunit[mcell_adunits[i]][0]);
			maxd += mcell_dens[i];
		}
		t = 0;
		for (i = 0; i < P.NumAdunits; i++)
			t += P.PopByAdunit[i][1];
		i = P.N;
		P.N = (int)t;
		fprintf(stderr, "Population size reset from %i to %i\n", i, P.N);
	}
	t = 1.0; m = 0;
	P.NMCP = 0;
	for (i = 0; i < (P.NMC - 1); i++)
	{
		s = mcell_dens[i] / maxd / t;
		if (s > 1.0) s = 1.0;
		m += (Mcells[i].n = (int)ignbin_mt((long)(P.N - m), s, 0));
		t -= mcell_dens[i] / maxd;
		if (Mcells[i].n > 0) P.NMCP++;
		AdUnits[mcell_adunits[i]].n += Mcells[i].n;
	}
	Mcells[P.NMC - 1].n = P.N - m;
	if (Mcells[P.NMC - 1].n > 0) P.NMCP++;
	free(mcell_num);
	free(mcell_country);
	free(mcell_adunits);
	t = 0.0;

	if (!(McellLookup = (microcell * *)malloc(P.NMCP * sizeof(microcell*)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(mcl = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	State.CellMemberArray = mcl;
	P.NCP = 0;
	for (i = i2 = j2 = 0; i < P.NC; i++)
	{
		Cells[i].n = 0;
		k = (i / P.nch) * P.NMCL * P.nmch + (i % P.nch) * P.NMCL;
		Cells[i].members = mcl + j2;
		for (l = 0; l < P.NMCL; l++)
			for (m = 0; m < P.NMCL; m++)
			{
				j = k + m + l * P.nmch;
				if (Mcells[j].n > 0)
				{
					Mcells[j].members = mcl + j2;
					//if(!(Mcells[j].members=(int *) calloc(Mcells[j].n,sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n"); //replaced line above with this to ensure members don't get mixed across microcells
					McellLookup[i2++] = Mcells + j;
					Cells[i].n += Mcells[j].n;
					j2 += Mcells[j].n;
				}
			}
		if (Cells[i].n > 0) P.NCP++;
	}
	fprintf(stderr, "Number of hosts assigned = %i\n", j2);
	if (!P.DoAdUnits) P.AdunitLevel1Lookup[0] = 0;
	fprintf(stderr, "Number of cells with non-zero population = %i\n", P.NCP);
	fprintf(stderr, "Number of microcells with non-zero population = %i\n", P.NMCP);

	if (!(CellLookup = (cell * *)malloc(P.NCP * sizeof(cell*)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	if (!(mcl = (int*)malloc(P.N * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	State.CellSuscMemberArray = mcl;
	i2 = k = 0;
	for (j = 0; j < P.NC; j++)
		if (Cells[j].n > 0)
		{
			CellLookup[i2++] = Cells + j;
			Cells[j].susceptible = mcl + k;
			k += Cells[j].n;
		}
	if (i2 > P.NCP) fprintf(stderr, "######## Over-run on CellLookup array NCP=%i i2=%i ###########\n", P.NCP, i2);
	if (!(RevCellLookup = (int*)malloc(P.NC * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
	i2 = 0;
	for (j = 0; j < P.NC; j++)
		if (Cells[j].n > 0)
			RevCellLookup[j] = i2++;
		else
			RevCellLookup[j] = -1;

	if (!(Hosts = (person*)calloc(P.N, sizeof(person)))) ERR_CRITICAL("Unable to allocate host storage\n");
	fprintf(stderr, "sizeof(person)=%i\n", (int) sizeof(person));
	for (i = 0; i < P.NCP; i++)
	{
		j = (int)(CellLookup[i] - Cells);
		if (Cells[j].n > 0)
		{
			if (!(Cells[j].InvCDF = (int*)malloc(1025 * sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
			if (!(Cells[j].max_trans = (float*)malloc(P.NCP * sizeof(float)))) ERR_CRITICAL("Unable to allocate cell storage\n");
			if (!(Cells[j].cum_trans = (float*)malloc(P.NCP * sizeof(float)))) ERR_CRITICAL("Unable to allocate cell storage\n");
		}
	}
	for (i = 0; i < P.NC; i++)
	{
		Cells[i].cumTC = 0;
		for (j = 0; j < Cells[i].n; j++) Cells[i].members[j] = -1;
	}
	fprintf(stderr, "Cells assigned\n");
	for (i = 0; i <= MAX_HOUSEHOLD_SIZE; i++) denom_household[i] = 0;
	P.NH = 0;
	for (i = j2 = 0; j2 < P.NMCP; j2++)
	{
		j = (int)(McellLookup[j2] - Mcells);
		l = ((j / P.nmch) / P.NMCL) * P.nch + ((j % P.nmch) / P.NMCL);
		ad = ((P.DoAdunitDemog) && (P.DoAdUnits)) ? Mcells[j].adunit : 0;
		for (k = 0; k < Mcells[j].n;)
		{
			m = 1;
			if (P.DoHouseholds)
			{
				s = ranf_mt(0);
				while ((s > P.HouseholdSizeDistrib[ad][m - 1]) && (k + m < Mcells[j].n) && (m < MAX_HOUSEHOLD_SIZE)) m++;
			}
			denom_household[m]++;
			for (i2 = 0; i2 < m; i2++)
			{
				//				fprintf(stderr,"%i ",i+i2);
				Hosts[i + i2].listpos = m; //used temporarily to store household size
				Mcells[j].members[k + i2] = i + i2;
				Cells[l].susceptible[Cells[l].cumTC] = i + i2;
				Cells[l].members[Cells[l].cumTC++] = i + i2;
				Hosts[i + i2].pcell = l;
				Hosts[i + i2].mcell = j;
				Hosts[i + i2].hh = P.NH;
			}
			P.NH++;
			i += m;
			k += m;
		}
	}
	if (!(Households = (household*)malloc(P.NH * sizeof(household)))) ERR_CRITICAL("Unable to allocate household storage\n");
	for (j = 0; j < NUM_AGE_GROUPS; j++) AgeDist[j] = AgeDist2[j] = 0;
	if (P.DoHouseholds) fprintf(stderr, "Household sizes assigned to %i people\n", i);
#pragma omp parallel for private(tn,j2,j,i,k,x,y,xh,yh,i2,m) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)
		for (j2 = tn; j2 < P.NMCP; j2 += P.NumThreads)
		{
			j = (int)(McellLookup[j2] - Mcells);
			x = (double)(j / P.nmch);
			y = (double)(j % P.nmch);
			i = Mcells[j].members[0];
			if (j % 100 == 0)
				fprintf(stderr, "%i=%i (%i %i)            \r", j, Mcells[j].n, Mcells[j].adunit, (AdUnits[Mcells[j].adunit].id % P.AdunitLevel1Mask) / P.AdunitLevel1Divisor);
			for (k = 0; k < Mcells[j].n;)
			{
				m = Hosts[i].listpos;
				xh = P.mcwidth * (ranf_mt(tn) + x);
				yh = P.mcheight * (ranf_mt(tn) + y);
				AssignHouseholdAges(m, i, tn);
				for (i2 = 0; i2 < m; i2++) Hosts[i + i2].listpos = 0;
				if (P.DoHouseholds)
				{
					for (i2 = 0; i2 < m; i2++) {
						Hosts[i + i2].inf = InfStat_Susceptible; //added this so that infection status is set to zero and household r0 is correctly calculated
					}
				}
				Households[Hosts[i].hh].FirstPerson = i;
				Households[Hosts[i].hh].nh = m;
				Households[Hosts[i].hh].nhr = m;
				Households[Hosts[i].hh].loc_x = (float)xh;
				Households[Hosts[i].hh].loc_y = (float)yh;
				i += m;
				k += m;
			}
		}
	if (P.DoCorrectAgeDist)
	{
		double** AgeDistAd, ** AgeDistCorrF, ** AgeDistCorrB;
		if (!(AgeDistAd = (double**)malloc(MAX_ADUNITS * sizeof(double*)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		if (!(AgeDistCorrF = (double**)malloc(MAX_ADUNITS * sizeof(double*)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		if (!(AgeDistCorrB = (double**)malloc(MAX_ADUNITS * sizeof(double*)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		for (i = 0; i < P.NumAdunits; i++)
		{
			if (!(AgeDistAd[i] = (double*)malloc((NUM_AGE_GROUPS + 1) * sizeof(double)))) ERR_CRITICAL("Unable to allocate temp storage\n");
			if (!(AgeDistCorrF[i] = (double*)malloc((NUM_AGE_GROUPS + 1) * sizeof(double)))) ERR_CRITICAL("Unable to allocate temp storage\n");
			if (!(AgeDistCorrB[i] = (double*)malloc((NUM_AGE_GROUPS + 1) * sizeof(double)))) ERR_CRITICAL("Unable to allocate temp storage\n");
		}

		for (i = 0; i < P.NumAdunits; i++)
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				AgeDistAd[i][j] = 0;
		for (i = 0; i < P.N; i++)
		{
			k = (P.DoAdunitDemog) ? Mcells[Hosts[i].mcell].adunit : 0;
			AgeDistAd[k][HOST_AGE_GROUP(i)]++;
		}
		k = (P.DoAdunitDemog) ? P.NumAdunits : 1;
		for (i = 0; i < k; i++)
		{
			s = 0.0;
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				s += AgeDistAd[i][j];
			for (j = 0; j < NUM_AGE_GROUPS; j++)
				AgeDistAd[i][j] /= s;
		}
		for (i = 0; i < k; i++)
		{
			s = t = 0;
			AgeDistCorrB[i][0] = 0;
			for (j = 0; j < NUM_AGE_GROUPS; j++)
			{
				s = t + AgeDistAd[i][j] - P.PropAgeGroup[i][j] - AgeDistCorrB[i][j];
				if (s > 0)
				{
					t = AgeDistCorrF[i][j] = s;
					AgeDistCorrB[i][j + 1] = 0;
				}
				else
				{
					t = AgeDistCorrF[i][j] = 0;
					AgeDistCorrB[i][j + 1] = fabs(s);
				}
				AgeDistCorrF[i][j] /= AgeDistAd[i][j];
				AgeDistCorrB[i][j] /= AgeDistAd[i][j];
			}
		}


#pragma omp parallel for private(tn,j,i,k,m) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (i = tn; i < P.N; i += P.NumThreads)
			{
				m = (P.DoAdunitDemog) ? Mcells[Hosts[i].mcell].adunit : 0;
				j = HOST_AGE_GROUP(i);
				s = ranf_mt(tn);
				if (s < AgeDistCorrF[m][j])
					Hosts[i].age += 5.0;
				else if (s < AgeDistCorrF[m][j] + AgeDistCorrB[m][j])
					Hosts[i].age -= 5.0;
			}
		for (i = 0; i < P.NumAdunits; i++)
		{
			free(AgeDistAd[i]);
			free(AgeDistCorrF[i]);
			free(AgeDistCorrB[i]);
		}
		free(AgeDistAd);
		free(AgeDistCorrF);
		free(AgeDistCorrB);
	}
	for (i = 0; i < P.N; i++) AgeDist[HOST_AGE_GROUP(i)]++;
	fprintf(stderr, "Ages/households assigned\n");

	if (!P.DoRandomInitialInfectionLoc)
	{
		k = (int)(P.LocationInitialInfection[0][0] / P.mcwidth);
		l = (int)(P.LocationInitialInfection[0][1] / P.mcheight);
		j = k * P.nmch + l;

		double rand_r = 0.0; //added these variables so that if initial infection location is empty we can search the 10km neighbourhood to find a suitable cell
		double rand_theta = 0.0;
		int counter = 0;
		if (Mcells[j].n < P.NumInitialInfections[0])
		{
			while (Mcells[j].n < P.NumInitialInfections[0] && counter < 100)
			{
				rand_r = ranf(); rand_theta = ranf();
				rand_r = 0.083 * sqrt(rand_r); rand_theta = 2 * PI * rand_theta; //rand_r is multiplied by 0.083 as this is roughly equal to 10km in decimal degrees
				k = (int)((P.LocationInitialInfection[0][0] + rand_r * cos(rand_theta)) / P.mcwidth);
				l = (int)((P.LocationInitialInfection[0][1] + rand_r * sin(rand_theta)) / P.mcheight);
				j = k * P.nmch + l;
				counter++;
			}
			if (counter < 100)
			{
				P.LocationInitialInfection[0][0] = P.LocationInitialInfection[0][0] + rand_r * cos(rand_theta); //set LocationInitialInfection to actual one used
				P.LocationInitialInfection[0][1] = P.LocationInitialInfection[0][1] + rand_r * sin(rand_theta);
			}
		}
		if (Mcells[j].n < P.NumInitialInfections[0])
			ERR_CRITICAL("Too few people in seed microcell to start epidemic with required number of initial infectionz.\n");
	}
	fprintf(stderr, "Checking cells...\n");
	maxd = 0; last_i = 0;
	for (i = 0; i < P.NMC; i++)
	{
		l = ((i / P.nmch) / P.NMCL) * P.nch + ((i % P.nmch) / P.NMCL);
		if (Cells[l].n == 0)
			mcell_dens[i] = 0;
		else
			last_i = i;
		maxd += mcell_dens[i];
	}
	fprintf(stderr, "Allocating place/age groups...\n");
	for (k = 0; k < NUM_AGE_GROUPS * AGE_GROUP_WIDTH; k++)
	{
		for (l = 0; l < P.PlaceTypeNum; l++)
		{
			PropPlaces[k][l] = PropPlacesC[k][l] = 0.0;
			if ((k < P.PlaceTypeAgeMax[l]) && (k >= P.PlaceTypeAgeMin[l]))
				PropPlaces[k][l] += P.PlaceTypePropAgeGroup[l];
			if ((k < P.PlaceTypeAgeMax2[l]) && (k >= P.PlaceTypeAgeMin2[l]))
				PropPlaces[k][l] += P.PlaceTypePropAgeGroup2[l];
			if ((k < P.PlaceTypeAgeMax3[l]) && (k >= P.PlaceTypeAgeMin3[l]))
				PropPlaces[k][l] += P.PlaceTypePropAgeGroup3[l];
			if (l == HOTEL_PLACE_TYPE)
				PropPlacesC[k][l] = ((l > 0) ? PropPlacesC[k][l - 1] : 0);
			else
				PropPlacesC[k][l] = PropPlaces[k][l] + ((l > 0) ? PropPlacesC[k][l - 1] : 0);
		}
	}
	/*
		for(l=0;l<P.PlaceTypeNum;l++)
			{
			for(k=0;k<NUM_AGE_GROUPS*AGE_GROUP_WIDTH;k++)
				fprintf(stderr, "%i:%lg ",k,PropPlaces[k][l]);
			fprintf(stderr,"\n");
			}
	*/
	/*	if((P.DoAdUnits)&&(P.DoAdunitDemog))
			{for(i=0;i<P.NumAdunits;i++) free(State.InvAgeDist[i]);}
		else
			free(State.InvAgeDist[0]);
		free(State.InvAgeDist);
	*/	P.nsp = 0;
	if ((P.DoSchoolFile) && (P.DoPlaces))
	{
		fprintf(stderr, "Reading school file\n");
		if (!(dat = fopen(SchoolFile, "r"))) ERR_CRITICAL("Unable to open school file\n");
		fscanf(dat, "%i", &P.nsp);
		for (j = 0; j < P.nsp; j++)
		{
			fscanf(dat, "%i %i", &m, &(P.PlaceTypeMaxAgeRead[j]));
			if (!(Places[j] = (place*)malloc(m * sizeof(place)))) ERR_CRITICAL("Unable to allocate place storage\n");
			for (i = 0; i < m; i++)
				if (!(Places[j][i].AvailByAge = (unsigned short int*) malloc(P.PlaceTypeMaxAgeRead[j] * sizeof(unsigned short int)))) ERR_CRITICAL("Unable to allocate place storage\n");
			P.Nplace[j] = 0;
			for (i = 0; i < P.NMC; i++) Mcells[i].np[j] = 0;
		}
		mr = 0;
		while (!feof(dat))
		{
			fscanf(dat, "%lg %lg %i %i", &x, &y, &j, &m);
			for (i = 0; i < P.PlaceTypeMaxAgeRead[j]; i++) fscanf(dat, "%hu", &(Places[j][P.Nplace[j]].AvailByAge[i]));
			Places[j][P.Nplace[j]].loc_x = (float)(x - P.SpatialBoundingBox[0]);
			Places[j][P.Nplace[j]].loc_y = (float)(y - P.SpatialBoundingBox[1]);
			if ((x >= P.SpatialBoundingBox[0]) && (x < P.SpatialBoundingBox[2]) && (y >= P.SpatialBoundingBox[1]) && (y < P.SpatialBoundingBox[3]))
			{
				i = P.nch * ((int)(Places[j][P.Nplace[j]].loc_x / P.cwidth)) + ((int)(Places[j][P.Nplace[j]].loc_y / P.cheight));
				if (Cells[i].n == 0) mr++;
				Places[j][P.Nplace[j]].n = m;
				i = (int)(Places[j][P.Nplace[j]].loc_x / P.mcwidth);
				k = (int)(Places[j][P.Nplace[j]].loc_y / P.mcheight);
				j2 = i * P.nmch + k;
				Mcells[j2].np[j]++;
				Places[j][P.Nplace[j]].mcell = j2;
				P.Nplace[j]++;
				if (P.Nplace[j] % 1000 == 0) fprintf(stderr, "%i read    \r", P.Nplace[j]);
			}
		}
		fclose(dat);
		fprintf(stderr, "%i schools read (%i in empty cells)      \n", P.Nplace[j], mr);
		for (i = 0; i < P.NMC; i++)
			for (j = 0; j < P.nsp; j++)
				if (Mcells[i].np[j] > 0)
				{
					if (!(Mcells[i].places[j] = (int*)malloc(Mcells[i].np[j] * sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
					Mcells[i].np[j] = 0;
				}
		for (j = 0; j < P.nsp; j++)
		{
			t = s = 0;
			for (i = 0; i < P.N; i++)
				t += PropPlaces[HOST_AGE_YEAR(i)][j];
			for (i = 0; i < P.Nplace[j]; i++)
			{
				k = Places[j][i].mcell;
				Mcells[k].places[j][Mcells[k].np[j]++] = i;
				s += (double)Places[j][i].n;
			}
			fprintf(stderr, "School type %i: capacity=%lg demand=%lg\n", j, s, t);
			t /= s;
			for (i = 0; i < P.Nplace[j]; i++)
				Places[j][i].n = (int)ceil(((double)Places[j][i].n) * t);
		}
	}
	if (P.DoPlaces)
	{
		fprintf(stderr, "Configuring places...\n");
#pragma omp parallel for private(tn,j2,i,j,k,t,m,s,x,y,xh,yh) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (j2 = P.nsp + tn; j2 < P.PlaceTypeNum; j2 += P.NumThreads)
			{
				t = 0;
				P.PlaceTypeMaxAgeRead[j2] = 0;
				for (i = 0; i < P.N; i++)
					t += PropPlaces[HOST_AGE_YEAR(i)][j2];
				P.Nplace[j2] = (int)ceil(t / P.PlaceTypeMeanSize[j2]);
				fprintf(stderr, "[%i:%i %g] ", j2, P.Nplace[j2], t);
				if (!(Places[j2] = (place*)calloc(P.Nplace[j2], sizeof(place)))) ERR_CRITICAL("Unable to allocate place storage\n");
				t = 1.0;
				for (m = i = k = 0; i < P.NMC; i++)
				{
					s = mcell_dens[i] / maxd / t;
					if (s > 1.0) s = 1.0;
					if (i == last_i)
						m += (Mcells[last_i].np[j2] = P.Nplace[j2] - m);
					else
						m += (Mcells[i].np[j2] = (int)ignbin_mt((long)(P.Nplace[j2] - m), s, tn));
					t -= mcell_dens[i] / maxd;
					if (Mcells[i].np[j2] > 0)
					{
						if (!(Mcells[i].places[j2] = (int*)malloc(Mcells[i].np[j2] * sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						x = (double)(i / P.nmch);
						y = (double)(i % P.nmch);
						for (j = 0; j < Mcells[i].np[j2]; j++)
						{
							s = ranf_mt(tn);
							xh = P.mcwidth * (ranf_mt(tn) + x);
							yh = P.mcheight * (ranf_mt(tn) + y);
							Places[j2][k].loc_x = (float)xh;
							Places[j2][k].loc_y = (float)yh;
							Places[j2][k].n = 0;
							Places[j2][k].mcell = i;
							Places[j2][k].country = Mcells[i].country;
							Mcells[i].places[j2][j] = k;
							k++;
						}
					}
				}
			}
		for (k = 0; k < NUM_AGE_GROUPS * AGE_GROUP_WIDTH; k++)
			for (l = 1; l < P.PlaceTypeNum; l++)
				if (l != HOTEL_PLACE_TYPE)
				{
					if (PropPlacesC[k][l - 1] < 1)
						PropPlaces[k][l] /= (1 - PropPlacesC[k][l - 1]);
					else if (PropPlaces[k][l] != 0)
						PropPlaces[k][l] = 1.0;
				}
		fprintf(stderr, "Places assigned\n");
	}
	free(mcell_dens);
	l = 0;
	for (j = 0; j < P.NC; j++)
		if (l < Cells[j].n) l = Cells[j].n;
	if (!(SamplingQueue = (int**)malloc(P.NumThreads * sizeof(int*)))) ERR_CRITICAL("Unable to allocate state storage\n");
	P.InfQueuePeakLength = P.N / P.NumThreads / INF_QUEUE_SCALE;
#pragma omp parallel for private(i,k) schedule(static,1)
	for (i = 0; i < P.NumThreads; i++)
	{
		if (!(SamplingQueue[i] = (int*)malloc(2 * (MAX_PLACE_SIZE + CACHE_LINE_SIZE) * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		for (k = 0; k < P.NumThreads; k++)
			if (!(StateT[i].inf_queue[k] = (int*)malloc(P.InfQueuePeakLength * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
		if (!(StateT[i].cell_inf = (float*)malloc((l + 1) * sizeof(float)))) ERR_CRITICAL("Unable to allocate state storage\n");
		if (!(StateT[i].inv_cell_inf = (int*)malloc(1025 * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
	}

	//set up queues and storage for digital contact tracing
	if ((P.DoAdUnits) && (P.DoDigitalContactTracing))
	{
		for (i = 0; i < P.NumAdunits; i++)
		{
			if (!(AdUnits[i].dct_queue = (int*)malloc(P.InfQueuePeakLength * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			if (!(AdUnits[i].dct = (int*)malloc(P.InfQueuePeakLength * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			for (j = 0; j < P.NumThreads; j++)
			{
				if (!(StateT[j].dct_queue[i] = (int*)malloc(P.InfQueuePeakLength * sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			}
		}
	}

	//If outputting origin-destination matrix, set up storage for flow between admin units
	if ((P.DoAdUnits) && (P.DoOriginDestinationMatrix))
	{
		for (i = 0; i < P.NumAdunits; i++)
		{
			if (!(AdUnits[i].origin_dest = (double*)malloc(MAX_ADUNITS * sizeof(double)))) ERR_CRITICAL("Unable to allocate storage for origin destination matrix\n");
			for (j = 0; j < P.NumThreads; j++)
			{
				if (!(StateT[j].origin_dest[i] = (double*)calloc(MAX_ADUNITS, sizeof(double)))) ERR_CRITICAL("Unable to allocate state origin destination matrix storage\n");
			}
			//initialise to zero
			for (j = 0; j < P.NumAdunits; j++)
			{
				AdUnits[i].origin_dest[j] = 0.0;
			}
		}
	}

	for (i = 0; i < P.NC; i++)
	{
		Cells[i].cumTC = 0;
		Cells[i].S = Cells[i].n;
		Cells[i].L = Cells[i].I = 0;
	}
	fprintf(stderr, "Allocated cell and host memory\n");
	fprintf(stderr, "Assigned hosts to cells\n");

}
void SetupAirports(void)
{
	int i, j, k, l, m;
	double x, y, t, tmin;
	indexlist* base, *cur;

	fprintf(stderr, "Assigning airports to microcells\n");
	P.KernelType = P.AirportKernelType;
	P.KernelScale = P.AirportKernelScale;
	P.KernelShape = P.AirportKernelShape;
	P.KernelP3 = P.AirportKernelP3;
	P.KernelP4 = P.AirportKernelP4;
	InitKernel(1, 1.0);
	if (!(Airports[0].DestMcells = (indexlist*)calloc(P.NMCP * NNA, sizeof(indexlist)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	if (!(base = (indexlist*)calloc(P.NMCP * NNA, sizeof(indexlist)))) ERR_CRITICAL("Unable to allocate airport storage\n");
	for (i = 0; i < P.Nairports; i++) Airports[i].num_mcell = 0;
	cur = base;
	for (i = 0; i < P.NMC; i++)
		if (Mcells[i].n > 0)
		{
			Mcells[i].AirportList = cur;
			cur += NNA;
		}
#pragma omp parallel for private(i,j,k,l,x,y,t,tmin) schedule(static,10000)
	for (i = 0; i < P.NMC; i++)
		if (Mcells[i].n > 0)
		{
			if (i % 10000 == 0) fprintf(stderr, "\n%i           ", i);
			x = (((double)(i / P.nmch)) + 0.5) * P.mcwidth;
			y = (((double)(i % P.nmch)) + 0.5) * P.mcheight;
			k = l = 0;
			tmin = 1e20;
			for (j = 0; j < P.Nairports; j++)
				if (Airports[j].total_traffic > 0)
				{
					t = numKernel(dist2_raw(x, y, Airports[j].loc_x, Airports[j].loc_y)) * Airports[j].total_traffic;
					if (k < NNA)
					{
						Mcells[i].AirportList[k].id = j;
						Mcells[i].AirportList[k].prob = (float)t;
						if (t < tmin) { tmin = t; l = k; }
						k++;
					}
					else if (t > tmin)
					{
						Mcells[i].AirportList[l].id = j;
						Mcells[i].AirportList[l].prob = (float)t;
						tmin = 1e20;
						for (k = 0; k < NNA; k++)
							if (Mcells[i].AirportList[k].prob < tmin)
							{
								tmin = Mcells[i].AirportList[k].prob;
								l = k;
							}
					}
				}
			for (j = 0; j < NNA; j++)
				Airports[Mcells[i].AirportList[j].id].num_mcell++;
		}
	cur = Airports[0].DestMcells;
	fprintf(stderr, "Microcell airport lists collated.\n");
	for (i = 0; i < P.Nairports; i++)
	{
		Airports[i].DestMcells = cur;
		cur += Airports[i].num_mcell;
		Airports[i].num_mcell = 0;
	}
#pragma omp parallel for private(i,j,k,l,t,tmin) schedule(static,10000)
	for (i = 0; i < P.NMC; i++)
		if (Mcells[i].n > 0)
		{
			if (i % 10000 == 0) fprintf(stderr, "\n%i           ", i);
			t = 0;
			for (j = 0; j < NNA; j++)
			{
				t += Mcells[i].AirportList[j].prob;
				k = Mcells[i].AirportList[j].id;
#pragma omp critical (airport)
				l = (Airports[k].num_mcell++);
				Airports[k].DestMcells[l].id = i;
				Airports[k].DestMcells[l].prob = Mcells[i].AirportList[j].prob * ((float)Mcells[i].n);
			}
			tmin = 0;
			for (j = 0; j < NNA; j++)
			{
				Mcells[i].AirportList[j].prob = (float)(tmin + Mcells[i].AirportList[j].prob / t);
				tmin = Mcells[i].AirportList[j].prob;
			}
		}
	fprintf(stderr, "Airport microcell lists collated.\n");
	for (i = 0; i < P.Nairports; i++)
		if (Airports[i].total_traffic > 0)
		{
			for (j = 1; j < Airports[i].num_mcell; j++)
				Airports[i].DestMcells[j].prob += Airports[i].DestMcells[j - 1].prob;
			t = Airports[i].DestMcells[Airports[i].num_mcell - 1].prob;
			if (t == 0) t = 1.0;
			for (j = 0; j < Airports[i].num_mcell - 1; j++)
				Airports[i].DestMcells[j].prob = (float)(Airports[i].DestMcells[j].prob / t);
			if (Airports[i].num_mcell > 0) Airports[i].DestMcells[Airports[i].num_mcell - 1].prob = 1.0;
			for (j = l = 0; l <= 1024; l++)
			{
				t = ((double)l) / 1024.0;
				while (Airports[i].DestMcells[j].prob < t) j++;
				Airports[i].Inv_DestMcells[l] = j;
			}
			l = 0;
			for (j = 0; j < Airports[i].num_mcell; j++)
				l += Mcells[Airports[i].DestMcells[j].id].np[HOTEL_PLACE_TYPE];
			if (l < 10)
			{
				fprintf(stderr, "(%i ", l);
				l = 0;
				for (j = 0; j < Airports[i].num_mcell; j++)
					l += Mcells[Airports[i].DestMcells[j].id].n;
				fprintf(stderr, "%i %i) ", Airports[i].num_mcell, l);
			}
		}
	fprintf(stderr, "\nInitialising hotel to airport lookup tables\n");
	free(base);
#pragma omp parallel for private(i,j,l,m,t,tmin) schedule(static,1)
	for (i = 0; i < P.Nairports; i++)
		if (Airports[i].total_traffic > 0)
		{
			m = (int)(Airports[i].total_traffic / HOTELS_PER_1000PASSENGER / 1000);
			if (m < MIN_HOTELS_PER_AIRPORT) m = MIN_HOTELS_PER_AIRPORT;
			fprintf(stderr, "\n%i    ", i);
			tmin = MAX_DIST_AIRPORT_TO_HOTEL * MAX_DIST_AIRPORT_TO_HOTEL * 0.75;
			do
			{
				tmin += 0.25 * MAX_DIST_AIRPORT_TO_HOTEL * MAX_DIST_AIRPORT_TO_HOTEL;
				Airports[i].num_place = 0;
				for (j = 0; j < P.Nplace[HOTEL_PLACE_TYPE]; j++)
					if (dist2_raw(Airports[i].loc_x, Airports[i].loc_y,
						Places[HOTEL_PLACE_TYPE][j].loc_x, Places[HOTEL_PLACE_TYPE][j].loc_y) < tmin)
						Airports[i].num_place++;
			} while (Airports[i].num_place < m);
			if (tmin > MAX_DIST_AIRPORT_TO_HOTEL * MAX_DIST_AIRPORT_TO_HOTEL) fprintf(stderr, "*** %i : %lg %i ***\n", i, sqrt(tmin), Airports[i].num_place);
			if (!(Airports[i].DestPlaces = (indexlist*)calloc(Airports[i].num_place, sizeof(indexlist)))) ERR_CRITICAL("Unable to allocate airport storage\n");
			Airports[i].num_place = 0;
			for (j = 0; j < P.Nplace[HOTEL_PLACE_TYPE]; j++)
				if ((t = dist2_raw(Airports[i].loc_x, Airports[i].loc_y,
					Places[HOTEL_PLACE_TYPE][j].loc_x, Places[HOTEL_PLACE_TYPE][j].loc_y)) < tmin)
				{
					Airports[i].DestPlaces[Airports[i].num_place].prob = (float)numKernel(t);
					Airports[i].DestPlaces[Airports[i].num_place].id = j;
					Airports[i].num_place++;
				}
			t = 0;
			for (j = 0; j < Airports[i].num_place; j++)
			{
				Airports[i].DestPlaces[j].prob = (float)(t + Airports[i].DestPlaces[j].prob);
				t = Airports[i].DestPlaces[j].prob;
			}
			for (j = 0; j < Airports[i].num_place - 1; j++)
				Airports[i].DestPlaces[j].prob = (float)(Airports[i].DestPlaces[j].prob / t);
			if (Airports[i].num_place > 0) Airports[i].DestPlaces[Airports[i].num_place - 1].prob = 1.0;
			for (j = l = 0; l <= 1024; l++)
			{
				t = ((double)l) / 1024.0;
				while (Airports[i].DestPlaces[j].prob < t) j++;
				Airports[i].Inv_DestPlaces[l] = j;
			}
		}
	P.KernelType = P.MoveKernelType;
	P.KernelScale = P.MoveKernelScale;
	P.KernelShape = P.MoveKernelShape;
	P.KernelP3 = P.MoveKernelP3;
	P.KernelP4 = P.MoveKernelP4;
	for (i = 0; i < P.Nplace[HOTEL_PLACE_TYPE]; i++) Places[HOTEL_PLACE_TYPE][i].n = 0;
	InitKernel(0, 1.0);
	fprintf(stderr, "\nAirport initialisation completed successfully\n");
}

#define PROP_OTHER_PARENT_AWAY 0.0


void AssignHouseholdAges(int n, int pers, int tn)
{
	/* Complex household age distribution model
		- picks number of children (nc)
		- tries to space them reasonably
		- picks parental ages to be consistent with childrens' and each other
		- other adults in large housholds are assumed to be grandparents
		- for Thailand, 2 person households are 95% couples without children, 5% 1 parent families
	*/
	int i, j, k, l, nc, ad;
	int a[MAX_HOUSEHOLD_SIZE + 2];

	ad = ((P.DoAdunitDemog) && (P.DoAdUnits)) ? Mcells[Hosts[pers].mcell].adunit : 0;
	if (!P.DoHouseholds)
	{
		for (i = 0; i < n; i++)
			a[i] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
	}
	else
	{
		if (n == 1)
		{
			if (ranf_mt(tn) < ONE_PERS_HOUSE_PROB_OLD)
			{
				do
				{
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while (a[0] < NOCHILD_PERS_AGE);
#else
				while ((a[0] < NOCHILD_PERS_AGE)
					|| (ranf_mt(tn) > (((double)a[0]) - NOCHILD_PERS_AGE + 1) / (OLD_PERS_AGE - NOCHILD_PERS_AGE + 1)));
#endif
			}
			else if ((ONE_PERS_HOUSE_PROB_YOUNG > 0) && (ranf_mt(tn) < ONE_PERS_HOUSE_PROB_YOUNG / (1 - ONE_PERS_HOUSE_PROB_OLD)))
			{
				do
				{
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				} while ((a[0] > YOUNG_AND_SINGLE) || (a[0] < MIN_ADULT_AGE)
					|| (ranf_mt(tn) > 1 - YOUNG_AND_SINGLE_SLOPE * (((double)a[0]) - MIN_ADULT_AGE) / (YOUNG_AND_SINGLE - MIN_ADULT_AGE)));
			}
			else
				while ((a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))]) < MIN_ADULT_AGE);
		}
		else if (n == 2)
		{
			if (ranf_mt(tn) < TWO_PERS_HOUSE_PROB_OLD)
			{
				do
				{
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while (a[0] < NOCHILD_PERS_AGE);
#else
				while ((a[0] < NOCHILD_PERS_AGE)
					|| (ranf_mt(tn) > (((double)a[0]) - NOCHILD_PERS_AGE + 1) / (OLD_PERS_AGE - NOCHILD_PERS_AGE + 1)));
#endif
				do
				{
					a[1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while ((a[1] >= a[0] + MAX_MF_PARTNER_AGE_GAP) || (a[1] < a[0] - MAX_FM_PARTNER_AGE_GAP) || (a[1] < MIN_ADULT_AGE));
#else
				while ((a[1] > a[0] + MAX_MF_PARTNER_AGE_GAP) || (a[1] < a[0] - MAX_FM_PARTNER_AGE_GAP) || (a[1] < NOCHILD_PERS_AGE)
					|| (ranf_mt(tn) > (((double)a[1]) - NOCHILD_PERS_AGE + 1) / (OLD_PERS_AGE - NOCHILD_PERS_AGE + 1)));
#endif
			}
			else if (ranf_mt(tn) < ONE_CHILD_TWO_PERS_PROB / (1 - TWO_PERS_HOUSE_PROB_OLD))
			{
				while ((a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))]) > MAX_CHILD_AGE);
				do
				{
					a[1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while ((a[1] >= a[0] + MAX_PARENT_AGE_GAP) || (a[1] < a[0] + MIN_PARENT_AGE_GAP));
#else
				while ((a[1] > a[0] + MAX_PARENT_AGE_GAP) || (a[1] < a[0] + MIN_PARENT_AGE_GAP) || (a[1] < MIN_ADULT_AGE));
#endif
			}
			else if ((TWO_PERS_HOUSE_PROB_YOUNG > 0) && (ranf_mt(tn) < TWO_PERS_HOUSE_PROB_YOUNG / (1 - TWO_PERS_HOUSE_PROB_OLD - ONE_CHILD_TWO_PERS_PROB)))
			{
				do
				{
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				} while ((a[0] < MIN_ADULT_AGE) || (a[0] > YOUNG_AND_SINGLE)
					|| (ranf_mt(tn) > 1 - YOUNG_AND_SINGLE_SLOPE * (((double)a[0]) - MIN_ADULT_AGE) / (YOUNG_AND_SINGLE - MIN_ADULT_AGE)));
				do
				{
					a[1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while ((a[1] >= a[0] + MAX_MF_PARTNER_AGE_GAP) || (a[1] < a[0] - MAX_FM_PARTNER_AGE_GAP) || (a[1] < MIN_ADULT_AGE));
#else
				while ((a[1] > a[0] + MAX_MF_PARTNER_AGE_GAP) || (a[1] < a[0] - MAX_FM_PARTNER_AGE_GAP) || (a[1] < MIN_ADULT_AGE));
#endif
			}
			else
			{
				do
				{
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				} while (a[0] < MIN_ADULT_AGE);
				do
				{
					a[1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while ((a[1] > a[0] + MAX_MF_PARTNER_AGE_GAP) || (a[1] < a[0] - MAX_FM_PARTNER_AGE_GAP) || (a[1] < MIN_ADULT_AGE));
#else
				while ((a[1] > a[0] + MAX_MF_PARTNER_AGE_GAP) || (a[1] < a[0] - MAX_FM_PARTNER_AGE_GAP) || (a[1] < MIN_ADULT_AGE));
#endif
			}

		}
		else
		{
			if (n == 3)
			{
				if ((ZERO_CHILD_THREE_PERS_PROB > 0) || (TWO_CHILD_THREE_PERS_PROB > 0))
					nc = (ranf_mt(tn) < ZERO_CHILD_THREE_PERS_PROB) ? 0 : ((ranf_mt(tn) < TWO_CHILD_THREE_PERS_PROB) ? 2 : 1);
				else
					nc = 1;
			}
			else if (n == 4)
				nc = (ranf_mt(tn) < ONE_CHILD_FOUR_PERS_PROB) ? 1 : 2;
			else if (n == 5)
				nc = (ranf_mt(tn) < THREE_CHILD_FIVE_PERS_PROB) ? 3 : 2;
			else
#ifdef COUNTRY_INDONESIA
				do
				{
					nc = ignpoi_mt(FRAC_CHILDREN_BIG_HOUSEHOLDS * ((double)n), tn);
				} while ((nc > n - 2) || (nc < 2));
#else
				nc = n - 2 - (int)(3 * ranf_mt(tn));
#endif
			if (nc == 0)
			{
				do
				{
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
					a[1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while ((a[0] < MIN_ADULT_AGE) || (a[1] >= a[0] + MAX_PARENT_AGE_GAP)
					|| (a[1] < a[0] + MIN_PARENT_AGE_GAP) || (a[1] > NOCHILD_PERS_AGE));
#else
				while ((a[1] < MIN_ADULT_AGE) || (a[0] < MIN_ADULT_AGE));
#endif
				do
				{
					a[2] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
				}
#ifdef COUNTRY_THAILAND
				while ((a[2] >= a[1] + MAX_MF_PARTNER_AGE_GAP) || (a[2] < a[1] - MAX_FM_PARTNER_AGE_GAP));
#else
				while ((a[2] >= a[1] + MAX_MF_PARTNER_AGE_GAP) || (a[2] < a[1] - MAX_FM_PARTNER_AGE_GAP));
#endif
			}
			else
			{
#ifdef COUNTRY_THAILAND
				a[0] = 0;
				for (i = 1; i < nc; i++)
					a[i] = a[i - 1] + 1 + ((int)ignpoi_mt(MEAN_CHILD_AGE_GAP - 1, tn));
				j = a[nc - 1] - (MAX_PARENT_AGE_GAP - MIN_PARENT_AGE_GAP);
				if (j > 0)
					j += MAX_PARENT_AGE_GAP;
				else
					j = MAX_PARENT_AGE_GAP;
				do
				{
					a[nc] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
					k = ((nc > 1) ? a[nc - 1] : 0) + a[0];
					l = k - MAX_CHILD_AGE;
					if (ranf_mt(tn) < ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE) l = k - 5;
				} while ((l > 0) || (a[nc] > a[0] + j) || (a[nc] < k + MIN_PARENT_AGE_GAP));
				for (i = 1; i < nc; i++) a[i] += a[0];
				if ((n > nc + 1) && (ranf_mt(tn) > PROP_OTHER_PARENT_AWAY))
				{
					do
					{
						a[nc + 1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
					} while ((a[nc + 1] >= a[nc] + MAX_MF_PARTNER_AGE_GAP)
						|| (a[nc + 1] < a[nc] - MAX_FM_PARTNER_AGE_GAP));
					l = nc + 2;
				}
				else
					l = nc + 1;
#else
				do
				{
					a[0] = 0;
					for (i = 1; i < nc; i++)
						a[i] = a[i - 1] + 1 + ((int)ignpoi_mt(MEAN_CHILD_AGE_GAP - 1, tn));
					a[0] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))] - a[(int)(ranf_mt(tn) * ((double)nc))];
					for (i = 1; i < nc; i++) a[i] += a[0];
					k = (((nc == 1) && (ranf_mt(tn) < ONE_CHILD_PROB_YOUNGEST_CHILD_UNDER_FIVE)) || ((nc == 2) && (ranf_mt(tn) < TWO_CHILDREN_PROB_YOUNGEST_UNDER_FIVE))
						|| ((nc > 2) && (ranf_mt(tn) < PROB_YOUNGEST_CHILD_UNDER_FIVE))) ? 5 : MAX_CHILD_AGE;
				} while ((a[0] < 0) || (a[0] > k) || (a[nc - 1] > MAX_CHILD_AGE));
				j = a[nc - 1] - a[0] - (MAX_PARENT_AGE_GAP - MIN_PARENT_AGE_GAP);
				if (j > 0)
					j += MAX_PARENT_AGE_GAP;
				else
					j = MAX_PARENT_AGE_GAP;
				do
				{
					a[nc] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
					k = a[nc - 1];
					l = k - MAX_CHILD_AGE;
				} while ((a[nc] > a[0] + j) || (a[nc] < k + MIN_PARENT_AGE_GAP) || (a[nc] < MIN_ADULT_AGE));
				if ((n > nc + 1) && (ranf_mt(tn) > PROP_OTHER_PARENT_AWAY))
				{
					do
					{
						a[nc + 1] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))];
					} while ((a[nc + 1] > a[nc] + MAX_MF_PARTNER_AGE_GAP) || (a[nc + 1] < a[nc] - MAX_FM_PARTNER_AGE_GAP)
						|| (a[nc + 1] > a[0] + j) || (a[nc + 1] < k + MIN_PARENT_AGE_GAP) || (a[nc + 1] < MIN_ADULT_AGE));
					l = nc + 2;
				}
				else
					l = nc + 1;
#endif
				if (n > nc + 2)
				{
#ifdef COUNTRY_THAILAND
					j = a[nc] + MIN_PARENT_AGE_GAP;
#else
					j = ((a[nc + 1] > a[nc]) ? a[nc + 1] : a[nc]) + OLDER_GEN_GAP;
#endif
					if (j >= NUM_AGE_GROUPS * AGE_GROUP_WIDTH) j = NUM_AGE_GROUPS * AGE_GROUP_WIDTH - 1;
					if (j < NOCHILD_PERS_AGE) j = NOCHILD_PERS_AGE;
					for (i = nc + 2; i < n; i++)
						while ((a[i] = State.InvAgeDist[ad][(int)(1000.0 * ranf_mt(tn))]) < j);
				}
			}
		}
	}
	for (i = 0; i < n; i++) Hosts[pers + i].age = (unsigned short int) a[i];
}

void AssignPeopleToPlaces(void)
{
	int i, i2, j, j2, k, k2, l, m, m2, tp, f, f2, f3, f4, ic, mx, my, a, cnt, tn, ca, nt, nn;
	int* PeopleArray;
	int* NearestPlaces[MAX_NUM_THREADS];
	double s, t, s2, *NearestPlacesProb[MAX_NUM_THREADS];
	cell* ct;
	int npt;

	npt = NUM_PLACE_TYPES;

	if (P.DoPlaces)
	{
		fprintf(stderr, "Assigning people to places....\n");
		for (i = 0; i < P.NC; i++)
		{
			Cells[i].infected = Cells[i].susceptible;
			if (!(Cells[i].susceptible = (int*)calloc(Cells[i].n, sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
			Cells[i].cumTC = Cells[i].n;
		}

		//PropPlaces initialisation is only valid for non-overlapping places.

		for (i = 0; i < P.N; i++)
		{
			for (tp = 0; tp < npt; tp++) //Changed from 'for(tp=0;tp<P.PlaceTypeNum;tp++)' to try and assign -1 early and avoid problems when using less than the default number of placetypes later 
			{
				Hosts[i].PlaceLinks[tp] = -1;
			}
		}

		for (tp = 0; tp < P.PlaceTypeNum; tp++)
		{
			if (tp != HOTEL_PLACE_TYPE)
			{
				cnt = 0;
				for (a = 0; a < P.NCP; a++)
				{
					i = (int)(CellLookup[a] - Cells);
					Cells[i].n = 0;
					for (j = 0; j < Cells[i].cumTC; j++)
					{
						k = HOST_AGE_YEAR(Cells[i].members[j]);
						f = ((PropPlaces[k][tp] > 0) && (ranf() < PropPlaces[k][tp]));
						if (f)
							for (k = 0; (k < tp) && (f); k++)
								if (Hosts[Cells[i].members[j]].PlaceLinks[k] >= 0) f = 0; //(ranf()<P.PlaceExclusivityMatrix[tp][k]); 
						// Am assuming people can only belong to 1 place (and a hotel) at present
						if (f)
						{
							Cells[i].susceptible[Cells[i].n] = Cells[i].members[j];
							Cells[i].n++;
							cnt++;
						}
					}
					Cells[i].S = Cells[i].n;
					Cells[i].I = 0;
				}
				if (!(PeopleArray = (int*)calloc(cnt, sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
				j2 = 0;
				for (a = 0; a < P.NCP; a++)
				{
					i = (int)(CellLookup[a] - Cells);
					for (j = 0; j < Cells[i].n; j++)
					{
						PeopleArray[j2] = Cells[i].susceptible[j];
						j2++;
					}
				}
				for (i2 = 0; i2 < 2; i2++)
					for (k = 0; k < cnt; k++)
					{
						while ((l = (int)(((double)cnt) * ranf())) == k);
						j2 = PeopleArray[k];
						PeopleArray[k] = PeopleArray[l];
						PeopleArray[l] = j2;
					}
				m = 0;
				if (tp < P.nsp)
				{
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						m += (int)(Places[tp][i].treat_end_time = (unsigned short)Places[tp][i].n);
						Places[tp][i].n = 0;
					}
				}
				else if (P.PlaceTypeSizePower[tp] == 0 && P.PlaceTypeSizeSD[tp] == 0)
				{
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						Places[tp][i].n = 0;
						j = 1 + ((int)ignpoi(P.PlaceTypeMeanSize[tp] - 1));
						if (j > USHRT_MAX - 1) j = USHRT_MAX - 1;
						m += (int)(Places[tp][i].treat_end_time = (unsigned short)j);
					}
				}
				//added this code to allow a place size to be specified according to a lognormal distribution - ggilani 09/02/17
				else if (P.PlaceTypeSizePower[tp] == 0 && P.PlaceTypeSizeSD[tp] > 0)
				{
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						Places[tp][i].n = 0;
						j = (int)gen_lognormal(P.PlaceTypeMeanSize[tp], P.PlaceTypeSizeSD[tp]);
						if (j > USHRT_MAX - 1) j = USHRT_MAX - 1;
						m += (int)(Places[tp][i].treat_end_time = (unsigned short)j);
					}
				}
				else
				{
					s = pow(P.PlaceTypeSizeOffset[tp] / (P.PlaceTypeSizeOffset[tp] + P.PlaceTypeSizeMax[tp] - 1), P.PlaceTypeSizePower[tp]);
					for (i = 0; i < P.Nplace[tp]; i++)
					{
						j = (int)floor(P.PlaceTypeSizeOffset[tp] * pow((1 - s) * ranf() + s, -1 / P.PlaceTypeSizePower[tp]) + 1 - P.PlaceTypeSizeOffset[tp]);
						if (j > USHRT_MAX - 1) j = USHRT_MAX - 1;
						m += (int)(Places[tp][i].treat_end_time = (unsigned short)j);
						Places[tp][i].n = 0;
					}
				}
				if (tp < P.nsp)
				{
					t = ((double)m) / ((double)P.Nplace[tp]);
					fprintf(stderr, "Adjusting place weights by cell (Capacity=%i Demand=%i  Av place size=%lg)\n", m, cnt, t);
					for (i = 0; i < P.Nplace[tp]; i++)
						if (Places[tp][i].treat_end_time > 0)
						{
							j = (int)(Places[tp][i].loc_x / P.cwidth);
							k = j * P.nch + ((int)(Places[tp][i].loc_y / P.cheight));
							Cells[k].I += (int)Places[tp][i].treat_end_time;
						}
					for (k = 0; k < P.NC; k++)
					{
						i = k % P.nch;
						j = k / P.nch;
						f2 = Cells[k].I; f3 = Cells[k].S;
						if ((i > 0) && (j > 0))
						{
							f2 += Cells[(j - 1) * P.nch + (i - 1)].I; f3 += Cells[(j - 1) * P.nch + (i - 1)].S;
						}
						if (i > 0)
						{
							f2 += Cells[j * P.nch + (i - 1)].I; f3 += Cells[j * P.nch + (i - 1)].S;
						}
						if ((i > 0) && (j < P.ncw - 1))
						{
							f2 += Cells[(j + 1) * P.nch + (i - 1)].I; f3 += Cells[(j + 1) * P.nch + (i - 1)].S;
						}
						if (j > 0)
						{
							f2 += Cells[(j - 1) * P.nch + i].I; f3 += Cells[(j - 1) * P.nch + i].S;
						}
						if (j < P.ncw - 1)
						{
							f2 += Cells[(j + 1) * P.nch + i].I; f3 += Cells[(j + 1) * P.nch + i].S;
						}
						if ((i < P.nch - 1) && (j > 0))
						{
							f2 += Cells[(j - 1) * P.nch + (i + 1)].I; f3 += Cells[(j - 1) * P.nch + (i + 1)].S;
						}
						if (i < P.nch - 1)
						{
							f2 += Cells[j * P.nch + (i + 1)].I; f3 += Cells[j * P.nch + (i + 1)].S;
						}
						if ((i < P.nch - 1) && (j < P.ncw - 1))
						{
							f2 += Cells[(j + 1) * P.nch + (i + 1)].I; f3 += Cells[(j + 1) * P.nch + (i + 1)].S;
						}
						Cells[k].L = f3; Cells[k].R = f2;
					}
					m = f2 = f3 = f4 = 0;
					for (k = 0; k < P.NC; k++)
						if ((Cells[k].S > 0) && (Cells[k].I == 0))
						{
							f2 += Cells[k].S; f3++;
							if (Cells[k].R == 0) f4 += Cells[k].S;
						}
					fprintf(stderr, "Demand in cells with no places=%i in %i cells\nDemand in cells with no places <=1 cell away=%i\n", f2, f3, f4);
					for (i = 0; i < P.Nplace[tp]; i++)
						if (Places[tp][i].treat_end_time > 0)
						{
							j = (int)(Places[tp][i].loc_x / P.cwidth);
							k = j * P.nch + ((int)(Places[tp][i].loc_y / P.cheight));
							if ((Cells[k].L > 0) && (Cells[k].R > 0))
							{
								s = ((double)Cells[k].L) / ((double)Cells[k].R);
								Places[tp][i].treat_end_time = (unsigned short)ceil(Places[tp][i].treat_end_time * s);
							}
							m += ((int)Places[tp][i].treat_end_time);
						}
					for (i = 0; i < P.NC; i++) Cells[i].L = Cells[i].I = Cells[i].R = 0;
				}
				t = ((double)m) / ((double)P.Nplace[tp]);
				fprintf(stderr, "Adjusting place weights (Capacity=%i Demand=%i  Av place size=%lg)\n", m, cnt, t);
				for (i = m = 0; i < P.Nplace[tp]; i++)
				{
					s = ((double)Places[tp][i].treat_end_time) * 43 / 40 - 1;
					m += (int)(Places[tp][i].treat_end_time = (unsigned short)(1.0 + ignpoi(s)));
				}
				if (tp < P.nsp)
					s = ((double)cnt) * 1.075;
				else
					s = ((double)cnt) * 1.125;
				j2 = ((int)s) - m;
				for (i = 0; i < j2; i++)
				{
					Places[tp][(int)(((double)P.Nplace[tp]) * ranf())].treat_end_time++;
				}
				j2 = -j2;
				for (i = 0; i < j2; i++)
				{
					while (Places[tp][j = (int)(((double)P.Nplace[tp]) * ranf())].treat_end_time < 2);
					Places[tp][j].treat_end_time--;
				}
				if (P.PlaceTypeNearestNeighb[tp] == 0)
				{
					for (i = 0; i < P.NC; i++) Cells[i].S = 0;
					for (j = 0; j < P.Nplace[tp]; j++)
					{
						i = P.nch * ((int)(Places[tp][j].loc_x / P.cwidth)) + ((int)(Places[tp][j].loc_y / P.cheight));
						Cells[i].S += (int)Places[tp][j].treat_end_time;
					}
					for (i = 0; i < P.NC; i++)
					{
						if (Cells[i].S > Cells[i].cumTC)
						{
							free(Cells[i].susceptible);
							if (!(Cells[i].susceptible = (int*)calloc(Cells[i].S, sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
						}
						Cells[i].S = 0;
					}
					for (j = 0; j < P.Nplace[tp]; j++)
					{
						i = P.nch * ((int)(Places[tp][j].loc_x / P.cwidth)) + ((int)(Places[tp][j].loc_y / P.cheight));
						k = (int)Places[tp][j].treat_end_time;
						for (j2 = 0; j2 < k; j2++)
						{
							Cells[i].susceptible[Cells[i].S] = j;
							Cells[i].S++;
						}
					}
				}
				for (i = 0; i < P.NumThreads; i++)
				{
					if (!(NearestPlaces[i] = (int*)calloc(P.PlaceTypeNearestNeighb[tp] + CACHE_LINE_SIZE, sizeof(int)))) ERR_CRITICAL("Unable to allocate cell storage\n");
					if (!(NearestPlacesProb[i] = (double*)calloc(P.PlaceTypeNearestNeighb[tp] + CACHE_LINE_SIZE, sizeof(double)))) ERR_CRITICAL("Unable to allocate cell storage\n");
				}
				P.KernelType = P.PlaceTypeKernelType[tp];
				P.KernelScale = P.PlaceTypeKernelScale[tp];
				P.KernelShape = P.PlaceTypeKernelShape[tp];
				P.KernelP3 = P.PlaceTypeKernelP3[tp];
				P.KernelP4 = P.PlaceTypeKernelP4[tp];
				InitKernel(1, 1.0);
				UpdateProbs(1);
				ca = 0;
				fprintf(stderr, "Allocating people to place type %i\n", tp);
				a = cnt;
				nt = P.NumThreads;
				nn = P.PlaceTypeNearestNeighb[tp];
				if (P.PlaceTypeNearestNeighb[tp] > 0)
				{
#pragma omp parallel for private(i,i2,j,j2,k,k2,l,m,m2,f,f2,ic,cnt,tn,s,t,mx,my) firstprivate(a,nt,nn) reduction(+:ca) schedule(static,1) //added i3, nh_assigned to private
					for (tn = 0; tn < P.NumThreads; tn++)
					{
						for (j = tn; j < a; j += nt)
						{
							if (j % 1000 == 0) fprintf(stderr, "(%i) %i      \r", tp, j);
							for (i2 = 0; i2 < nn; i2++)
								NearestPlacesProb[tn][i2] = 0;
							l = 1; k = m = m2 = f2 = 0;
							i = PeopleArray[j];
							ic = Hosts[i].mcell;
							mx = ic / P.nmch;
							my = ic % P.nmch;
							if (Hosts[i].PlaceLinks[tp] < 0) //added this so that if any hosts have already be assigned due to their household membership, they will not be reassigned
							{
								while (((k < nn) || (l < 4)) && (l < P.nmcw))
								{
									if ((mx >= 0) && (my >= 0) && (mx < P.nmcw) && (my < P.nmch))
									{
										ic = mx * P.nmch + my;
										if (Mcells[ic].country == Mcells[Hosts[i].mcell].country)
										{
											for (cnt = 0; cnt < Mcells[ic].np[tp]; cnt++)
											{
												if (Mcells[ic].places[tp][cnt] >= P.Nplace[tp]) fprintf(stderr, "#%i %i %i  ", tp, ic, cnt);
												t = dist2_raw(Households[Hosts[i].hh].loc_x, Households[Hosts[i].hh].loc_y,
													Places[tp][Mcells[ic].places[tp][cnt]].loc_x, Places[tp][Mcells[ic].places[tp][cnt]].loc_y);
												s = numKernel(t);
												if (tp < P.nsp)
												{
													t = ((double)Places[tp][Mcells[ic].places[tp][cnt]].treat_end_time);
													if (HOST_AGE_YEAR(i) < P.PlaceTypeMaxAgeRead[tp])
													{
														if ((t > 0) && (Places[tp][Mcells[ic].places[tp][cnt]].AvailByAge[HOST_AGE_YEAR(i)] > 0))
															s *= t;
														else
															s = 0;
													}
													else if (t > 0)
														s *= t;
												}
												k2 = 0; j2 = 0; t = 1e10;
												if (s > 0)
												{
													if (k < nn)
													{
														NearestPlaces[tn][k] = Mcells[ic].places[tp][cnt];
														NearestPlacesProb[tn][k] = s;
														k++;
													}
													else
													{
														for (i2 = 0; i2 < nn; i2++)
														{
															if (NearestPlacesProb[tn][i2] < t)
															{
																t = NearestPlacesProb[tn][i2]; j2 = i2;
															}
														}
														if (s > t)
														{
															NearestPlacesProb[tn][j2] = s;
															NearestPlaces[tn][j2] = Mcells[ic].places[tp][cnt];
														}
													}
												}
											}
										}
									}
									if (m2 == 0)
										mx = mx + 1;
									else if (m2 == 1)
										my = my - 1;
									else if (m2 == 2)
										mx = mx - 1;
									else if (m2 == 3)
										my = my + 1;
									f2 = (f2 + 1) % l;
									if (f2 == 0)
									{
										m2 = (m2 + 1) % 4;
										m = (m + 1) % 2;
										if (m == 0) l++;
									}
								}

								s = 0;
								if (k > nn) fprintf(stderr, "*** k>P.PlaceTypeNearestNeighb[tp] ***\n");
								if (k == 0)
								{
									fprintf(stderr, "# %i %i     \r", i, j);
									Hosts[i].PlaceLinks[tp] = -1;
								}
								else
								{
									for (i2 = 1; i2 < k; i2++)
										NearestPlacesProb[tn][i2] += NearestPlacesProb[tn][i2 - 1];
									s = NearestPlacesProb[tn][k - 1];
									t = ranf_mt(tn);
									f = 0;
									for (i2 = 0; (i2 < k) && (!f); i2++)
									{
										if ((f = (t < NearestPlacesProb[tn][i2] / s)))
										{
											Hosts[i].PlaceLinks[tp] = NearestPlaces[tn][i2];
											ca++;
											if (tp < P.nsp)
											{
#pragma omp critical (places_treat_time)
												Places[tp][Hosts[i].PlaceLinks[tp]].treat_end_time--;
											}
										}
										if (!f) Hosts[i].PlaceLinks[tp] = -1;
										if (NearestPlaces[tn][i2] >= P.Nplace[tp]) fprintf(stderr, "@%i %i %i  ", tp, i, j);
									}
								}
							}
						}

					}
				}
				else
				{
					k2 = cnt - ca;
					m2 = cnt;
					s2 = 0;
					a = k2 / 1000;
					f = k2;
					for (ic = 0; ic <= 30; ic++)
					{
						UpdateProbs(1);
						m2 = f - 1;
						if (ic < 9)
							f = 100 * (9 - ic) * a;
						else if (ic < 18)
							f = 10 * (18 - ic) * a;
						else if (ic < 27)
							f = (27 - ic) * a;
						else
						{
							m2 = k2 - 1; f = 0;
						}
#pragma omp parallel for private(i,i2,j,k,l,m,f2,f3,t,ct,s) reduction(+:ca) /* schedule(dynamic,500)*/ //add s to private variables, added g,g1,g2,i3 and nh_assigned to private variables
						for (i2 = m2; i2 >= f; i2--)
						{
							if (i2 % 1000 == 0)
								fprintf(stderr, "(%i) %i            \r", tp, i2);
							k = PeopleArray[i2];
							i = Hosts[k].pcell;
							f2 = 1;
							f3 = (HOST_AGE_YEAR(k) >= P.PlaceTypeMaxAgeRead[tp]);
							if (Hosts[k].PlaceLinks[tp] < 0)
								while ((f2 > 0) && (f2 < 1000))
								{
									do
									{
										s = ranf();
										l = Cells[i].InvCDF[(int)floor(s * 1024)];
										while (Cells[i].cum_trans[l] < s) l++;
										ct = CellLookup[l];
										m = (int)(ranf() * ((double)ct->S));
										j = -1;
#pragma omp critical
										{
											if (ct->susceptible[m] >= 0)
												if ((f3) || (Places[tp][ct->susceptible[m]].AvailByAge[HOST_AGE_YEAR(k)] > 0))
												{
													j = ct->susceptible[m];
													ct->susceptible[m] = -1;
												}
										}
									} while (j < 0);
									if (j >= P.Nplace[tp])
									{
										fprintf(stderr, "*%i %i: %i %i\n", k, tp, j, P.Nplace[tp]);
										ERR_CRITICAL("Out of bounds place link\n");
									}
									t = dist2_raw(Households[Hosts[k].hh].loc_x, Households[Hosts[k].hh].loc_y, Places[tp][j].loc_x, Places[tp][j].loc_y);
									s = ((double)ct->S) / ((double)ct->S0) * numKernel(t) / Cells[i].max_trans[l];
									if ((P.DoAdUnits) && (P.InhibitInterAdunitPlaceAssignment[tp] > 0))
									{
										if (Mcells[Hosts[k].mcell].adunit != Mcells[Places[tp][j].mcell].adunit) s *= (1 - P.InhibitInterAdunitPlaceAssignment[tp]);
									}
									if (ranf() < s)
									{
#pragma omp critical
										l = (--ct->S);
										if (m < l) ct->susceptible[m] = ct->susceptible[l];
#pragma omp critical (places_treat_time)
										Places[tp][j].treat_end_time--;
										ca++;
										Hosts[k].PlaceLinks[tp] = j;
										f2 = 0;
									}
									else
									{
										ct->susceptible[m] = j;
										f2++;
									}
								}
						}
					}
				}
				fprintf(stderr, "%i hosts assigned to placetype %i\n", ca, tp);
				free(PeopleArray);
				for (i = 0; i < P.Nplace[tp]; i++)
				{
					Places[tp][i].treat_end_time = 0;
					Places[tp][i].n = 0;
				}
				for (i = 0; i < P.NumThreads; i++)
				{
					free(NearestPlacesProb[i]);
					free(NearestPlaces[i]);
				}
			}
		}
		for (i = 0; i < P.NC; i++)
		{
			Cells[i].n = Cells[i].cumTC;
			Cells[i].cumTC = 0;
			Cells[i].S = Cells[i].I = Cells[i].L = Cells[i].R = 0;
			free(Cells[i].susceptible);
			Cells[i].susceptible = Cells[i].infected;
		}
		P.KernelScale = P.MoveKernelScale;
		P.KernelShape = P.MoveKernelShape;
		P.KernelType = P.MoveKernelType;
		P.KernelP3 = P.MoveKernelP3;
		P.KernelP4 = P.MoveKernelP4;
	}

}
void StratifyPlaces(void)
{
	int i, j, k, l, m, n, tn;
	double t;

	if (P.DoPlaces)
	{
		fprintf(stderr, "Initialising groups in places\n");
#pragma omp parallel for private(i,j) schedule(static,500)
		for (i = 0; i < P.N; i++)
			for (j = 0; j < NUM_PLACE_TYPES; j++)
				Hosts[i].PlaceGroupLinks[j] = 0;
		for (j = 0; j < P.PlaceTypeNum; j++)
			for (i = 0; i < P.Nplace[j]; i++)
				Places[j][i].n = 0;
#pragma omp parallel for private(i,j,k,l,m,n,t,tn) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (j = tn; j < P.PlaceTypeNum; j += P.NumThreads)
			{
				if (j == HOTEL_PLACE_TYPE)
				{
					l = 2 * ((int)P.PlaceTypeMeanSize[j]);
					for (i = 0; i < P.Nplace[j]; i++)
					{
						if (!(Places[j][i].members = (int*)calloc(l, sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						Places[j][i].n = 0;
					}
				}
				else
				{
					for (i = 0; i < P.N; i++)
					{
						if (Hosts[i].PlaceLinks[j] >= 0)
							Places[j][Hosts[i].PlaceLinks[j]].n++;
					}
					for (i = 0; i < P.Nplace[j]; i++)
					{
						if (Places[j][i].n > 0)
						{
							if (!(Places[j][i].members = (int*)calloc(Places[j][i].n, sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
						}
						Places[j][i].n = 0;
					}
					for (i = 0; i < P.N; i++)
					{
						k = Hosts[i].PlaceLinks[j];
						if (k >= 0)
						{
							Places[j][k].members[Places[j][k].n] = i;
							Places[j][k].n++;
						}
					}
					for (i = 0; i < P.Nplace[j]; i++)
						if (Places[j][i].n > 0)
						{
							t = ((double)Places[j][i].n) / P.PlaceTypeGroupSizeParam1[j] - 1.0;
							if (t < 0)
								Places[j][i].ng = 1;
							else
								Places[j][i].ng = 1 + (int)ignpoi_mt(t, tn);
							if (!(Places[j][i].group_start = (int*)calloc(Places[j][i].ng, sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
							if (!(Places[j][i].group_size = (int*)calloc(Places[j][i].ng, sizeof(int)))) ERR_CRITICAL("Unable to allocate place storage\n");
							m = Places[j][i].n - Places[j][i].ng;
							for (k = l = 0; k < Places[j][i].ng; k++)
							{
								t = 1 / ((double)(Places[j][i].ng - k));
								Places[j][i].group_start[k] = l;
								Places[j][i].group_size[k] = 1 + ignbin_mt((long)m, t, tn);
								m -= (Places[j][i].group_size[k] - 1);
								l += Places[j][i].group_size[k];
							}
							for (k = 0; k < Places[j][i].n; k++)
							{
								l = (int)(((double)Places[j][i].n) * ranf_mt(tn));
								n = Places[j][i].members[l];
								Places[j][i].members[l] = Places[j][i].members[k];
								Places[j][i].members[k] = n;
							}
							for (k = l = 0; k < Places[j][i].ng; k++)
								for (m = 0; m < Places[j][i].group_size[k]; m++)
								{
									Hosts[Places[j][i].members[l]].PlaceGroupLinks[j] = k;
									l++;
								}
						}
				}
			}

#pragma omp parallel for private(i,j,k,l) schedule(static,1)
		for (i = 0; i < P.NumThreads; i++)
		{
			for (k = 0; k < P.PlaceTypeNum; k++)
			{
				if (P.DoPlaceGroupTreat)
				{
					l = 0;
					for (j = 0; j < P.Nplace[k]; j++)
						l += (int)Places[k][j].ng;
					m = 2 * m / P.NumThreads;
					if (!(StateT[i].p_queue[k] = (int*)calloc(l, sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
					if (!(StateT[i].pg_queue[k] = (int*)calloc(l, sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
				}
				else
				{
					if (!(StateT[i].p_queue[k] = (int*)calloc(P.Nplace[k], sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
					if (!(StateT[i].pg_queue[k] = (int*)calloc(P.Nplace[k], sizeof(int)))) ERR_CRITICAL("Unable to allocate state storage\n");
				}
			}
		}
		fprintf(stderr, "Groups initialised\n");
		/*		s2=t2=0;
				for(j=0;j<P.PlaceTypeNum;j++)
					{
					t=s=0;
					for(i=0;i<P.Nplace[j];i++)
						if(Places[j][i].ng>0)
							{
							for(k=0;k<Places[j][i].ng;k++)
								t+=(double) Places[j][i].group_size[k];
							s+=(double) Places[j][i].ng;
							}
					s2+=s;
					t2+=t;
					fprintf(stderr,"Mean group size for place type %i = %lg\n",j,t/s);
					}
				t=0;
				for(i=0;i<P.N;i++)
					for(j=0;j<P.PlaceTypeNum;j++)
						if(Hosts[i].PlaceLinks[j]>=0)
							t+=(double) Places[j][Hosts[i].PlaceLinks[j]].group_size[Hosts[i].PlaceGroupLinks[j]];
				fprintf(stderr,"Overall mean group size = %lg (%lg)\n",t/((double) P.N),t2/s2);
		*/
	}
}
void LoadPeopleToPlaces(char* NetworkFile)
{
	int i, j, k, l, m, n, npt, i2;
	long s1, s2;
	FILE* dat;
	int fileversion;

	if (!(dat = fopen(NetworkFile, "rb"))) ERR_CRITICAL("Unable to open network file for loading\n");
	fread_big(&fileversion, sizeof(fileversion), 1, dat);
	if (fileversion != NETWORK_FILE_VERSION)
	{
		ERR_CRITICAL("Incompatible network file - please rebuild using '/S:'.\n");
	}

	npt = NUM_PLACE_TYPES_NOAIR;
	fread_big(&i, sizeof(int), 1, dat);
	fread_big(&j, sizeof(int), 1, dat);
	fread_big(&s1, sizeof(long), 1, dat);
	fread_big(&s2, sizeof(long), 1, dat);
	if (i != npt) ERR_CRITICAL("Number of place types does not match saved value\n");
	if (j != P.N) ERR_CRITICAL("Population size does not match saved value\n");
	if ((s1 != P.seed1) || (s2 != P.seed2)) ERR_CRITICAL("Random number seeds do not match saved values\n");
	k = (P.N + 999999) / 1000000;
	for (i = 0; i < P.N; i++)
		for (j = 0; j < P.PlaceTypeNum; j++)
			Hosts[i].PlaceLinks[j] = -1;
	for (i = i2 = 0; i < k; i++)
	{
		l = (i < k - 1) ? 1000000 : (P.N - 1000000 * (k - 1));
		fread_big(&netbuf, sizeof(int), npt * l, dat);
		for (j = 0; j < l; j++)
		{
			n = j * npt;
			for (m = 0; m < npt; m++)
			{
				Hosts[i2].PlaceLinks[m] = netbuf[n + m];
				if (Hosts[i2].PlaceLinks[m] >= P.Nplace[m])
				{
					fprintf(stderr, "*%i %i: %i %i\n", i2, m, Hosts[i2].PlaceLinks[m], P.Nplace[m]);
					ERR_CRITICAL("Out of bounds place link\n");
				}
			}
			i2++;
		}
		fprintf(stderr, "%i loaded            \r", i * 1000000 + l);
	}

	/*	for(i=0;i<P.N;i++)
			{
			if((i+1)%100000==0) fprintf(stderr,"%i loaded            \r",i+1);
			fread_big(&(Hosts[i].PlaceLinks[0]),sizeof(int),P.PlaceTypeNum,dat);
			}
	*/	fprintf(stderr, "\n");
	fclose(dat);
	P.KernelScale = P.MoveKernelScale;
	P.KernelShape = P.MoveKernelShape;
	P.KernelP3 = P.MoveKernelP3;
	P.KernelP4 = P.MoveKernelP4;
}
void SavePeopleToPlaces(char* NetworkFile)
{
	int i, j, npt;
	FILE* dat;
	int fileversion = NETWORK_FILE_VERSION;

	npt = NUM_PLACE_TYPES_NOAIR;
	if (!(dat = fopen(NetworkFile, "wb"))) ERR_CRITICAL("Unable to open network file for saving\n");
	fwrite_big(&fileversion, sizeof(fileversion), 1, dat);

	if (P.PlaceTypeNum > 0)
	{
		fwrite_big(&npt, sizeof(int), 1, dat);
		fwrite_big(&(P.N), sizeof(int), 1, dat);
		fwrite_big(&P.seed1, sizeof(long), 1, dat);
		fwrite_big(&P.seed2, sizeof(long), 1, dat);
		for (i = 0; i < P.N; i++)
		{
			if ((i + 1) % 100000 == 0) fprintf(stderr, "%i saved            \r", i + 1);
			/*			fwrite_big(&(Hosts[i].spatial_norm),sizeof(float),1,dat);
			*/			fwrite_big(&(Hosts[i].PlaceLinks[0]), sizeof(int), npt, dat);
			for (j = 0; j < npt; j++)
				if (Hosts[i].PlaceLinks[j] >= P.Nplace[j])
				{
					fprintf(stderr, "*%i %i: %i %i\n", i, j, Hosts[i].PlaceLinks[j], P.Nplace[j]);
					ERR_CRITICAL("Out of bounds place link\n");
				}
		}
	}

	fprintf(stderr, "\n");
	fflush(dat);
	fclose(dat);
}

void SaveAgeDistrib(void)
{
	int i;
	FILE* dat;
	char outname[1024];

	sprintf(outname, "%s.agedist.xls", OutFile);
	if (!(dat = fopen(outname, "wb"))) ERR_CRITICAL("Unable to open output file\n");
	if (P.DoDeath)
	{
		fprintf(dat, "age\tfreq\tlifeexpect\n");
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "%i\ta%.10f\t%.10f\n", i, AgeDist[i], AgeDist2[i]);
		fprintf(dat, "\np\tlife_expec\tage\n");
		for (i = 0; i <= 1000; i++)
			fprintf(dat, "%.10f\t%.10f\t%i\n", ((double)i) / 1000, P.InvLifeExpecDist[0][i], State.InvAgeDist[0][i]);
	}
	else
	{
		fprintf(dat, "age\tfreq\n");
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			fprintf(dat, "%i\t%.10f\n", i, AgeDist[i]);
	}

	fclose(dat);
}
