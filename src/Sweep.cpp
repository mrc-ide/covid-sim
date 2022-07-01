	#include <climits>
#include <cmath>
#include <cstdlib>

#include "CalcInfSusc.h"
#include "Dist.h"
#include "Error.h"
#include "Files.h"
#include "InfStat.h"
#include "Kernels.h"
#include "Rand.h"
#include "Model.h"
#include "ModelMacros.h"
#include "Param.h"
#include "Sweep.h"
#include "Update.h"
#include <cassert>

// helper functions

void AddToInfectionQueue(const int tn, const int infectee_cell_number, const int infector_index, const int infectee_index, const short int infect_type);
bool AddInfections(const int tn, const int infectee_cell_index, const int infector_index, const int infectee_index,	const short int infect_type);

void TravelReturnSweep(double t)
{
	int l, nr, ner;

	// Convince static analysers that values are set correctly:
	if (!(P.DoAirports && P.HotelPlaceType < P.NumPlaceTypes)) ERR_CRITICAL("DoAirports || HotelPlaceType not set\n");

	if (floor(1 + t + P.ModelTimeStep) != floor(1 + t))
	{
		nr = ner = 0;
		int floorOfTime = (int)floor(t);
		l = 1 + floorOfTime % MAX_TRAVEL_TIME;
		FILE* stderr_shared = stderr;
#pragma omp parallel for reduction(+:nr, ner) schedule(static, 1) default(none) \
			shared(P, Places, Hosts, l, stderr_shared)
		for (int tn = 0; tn < P.NumThreads; tn++)
		{
			for (int j = tn; j < P.Nplace[P.HotelPlaceType]; j += P.NumThreads)
			{
				int n = Places[P.HotelPlaceType][j].n;
				for (int k = n - 1; k >= 0; k--)
				{
					int i = Places[P.HotelPlaceType][j].members[k];
					if (Hosts[i].Travelling == l)
					{
						n--;
						/*						if((n<0)||(Places[P.HotelPlaceType][j].members[n]<0)||(Places[P.HotelPlaceType][j].members[n]>=P.PopSize))
													{Files::xfprintf_stderr("### %i %i %i %i\n",j,k,n,Places[P.HotelPlaceType][j].members[n]);ner++;}
												else if((k<0)||(k>n))
													{Files::xfprintf_stderr("@ %i %i %i %i\n",j,k,n,Places[P.HotelPlaceType][j].members[n]);ner++;}
												else
						*/
						if (k != n)
						{
							Places[P.HotelPlaceType][j].members[k] = Places[P.HotelPlaceType][j].members[n];
						}
						nr++;
						if (Hosts[i].PlaceLinks[P.HotelPlaceType] != j)
						{
							ner++;
							Files::xfprintf(stderr_shared, "(%i %i) ", j, Hosts[i].PlaceLinks[P.HotelPlaceType]);
						}
						Hosts[i].PlaceLinks[P.HotelPlaceType] = -1;
						Hosts[i].Travelling = 0;
					}
				}
				Places[P.HotelPlaceType][j].n = n;
			}
		}
		Files::xfprintf_stderr(" d=%i e=%i>", nr, ner);
	}
}

void TravelDepartSweep(double t)
{
	int d, mps, nld, nad, nsk, bm;
	double nl;

	// Convince static analysers that values are set correctly:
	if (!(P.DoAirports && P.HotelPlaceType < P.NumPlaceTypes)) ERR_CRITICAL("DoAirports || HotelPlaceType not set\n");

	if (floor(1 + t - P.ModelTimeStep) != floor(1 + t))
	{
		bm = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));
		mps = 2 * ((int)P.PlaceTypeMeanSize[P.HotelPlaceType]) - P.NumThreads - 1;
		int floorOfTime = (int)floor(t);
		d = floorOfTime % MAX_TRAVEL_TIME;
		nad = nld = nsk = 0;
#pragma omp parallel for reduction(+:nad,nsk) schedule(static,1) default(none) \
			shared(t, P, Airports, Mcells, Hosts, Places, bm, mps, d)
		for (int tn = 0; tn < P.NumThreads; tn++)
			for (int i = tn; i < P.Nairports; i += P.NumThreads)
				if ((Airports[i].total_traffic > 0) && (Airports[i].num_mcell > 0))
				{
					double s = Airports[i].total_traffic;
					if ((t > P.AirportCloseTimeStart) && (t < P.AirportCloseTimeStart + P.AirportCloseTimeStartBase))
						s *= P.AirportCloseEffectiveness;
					int n = (s > 0) ? ((int)ignpoi_mt((double)s, tn)) : 0;
					int f3 = 0;
					int j = 0;
					while (j < n)
					{
						s = ranf_mt(tn);
						int l = Airports[i].Inv_DestMcells[(int)floor(s * 1024)];
						while (Airports[i].DestMcells[l].prob < s) l++;
						l = Airports[i].DestMcells[l].id;
						int k = (int)(ranf_mt(tn) * ((double)Mcells[l].n));
						int i2 = Mcells[l].members[k];

						// Original:
						// if ((abs(Hosts[i2].inf) < InfStat::InfectiousAsymptomaticNotCase) && (Hosts[i2].inf != InfStat::Case))
						// but also note: above is equivalent to if ((abs(inf) < 2) && (inf != -2)),
						// so if h were -2, it would fail the first case, and the second is redundant.

						if (Hosts[i2].is_not_yet_symptomatic())
						{
							int d2 = HOST_AGE_GROUP(i2);
							if ((P.RelativeTravelRate[d2] == 1) || (ranf_mt(tn) < P.RelativeTravelRate[d2]))
							{
								int f2 = 1;
#pragma omp critical
								{
									if (Hosts[i2].PlaceLinks[P.HotelPlaceType] == -1)
									{
										Hosts[i2].PlaceLinks[P.HotelPlaceType] = -2;
										f2 = 0;
									}
								}
								if (!f2)
								{
									s = ranf_mt(tn);
									l = Airports[i].Inv_prop_traffic[(int)floor(s * 128)];
									while (Airports[i].prop_traffic[l] < s) l++;
									k = Airports[i].conn_airports[l];
									if (bm)
									{
										if (dist2_raw(Airports[i].loc.x, Airports[i].loc.y, Airports[k].loc.x, Airports[k].loc.y) > P.MoveRestrRadius2)
										{
											if (ranf_mt(tn) > P.MoveRestrEffect)
											{
												f2 = 1;
												nsk++;
												j++;
#pragma omp critical
												Hosts[i2].PlaceLinks[P.HotelPlaceType] = -1;
											}
										}
									}
									if (!f2)
									{
										int f = 1;
										do
										{
											s = ranf_mt(tn);
											int m = Airports[k].Inv_DestPlaces[(int)floor(s * 1024)];
											while (Airports[k].DestPlaces[m].prob < s) m++;
											l = Airports[k].DestPlaces[m].id;
											int hp;
#pragma omp critical
											{
												if ((hp = Places[P.HotelPlaceType][l].n) < mps)
												{
													f = 0;
													Places[P.HotelPlaceType][l].n++;
												}
											}
											if (!f)
											{
												f3 = 0;
												Places[P.HotelPlaceType][l].members[hp] = i2;
												d2 = (d + P.InvJourneyDurationDistrib[(int)(ranf_mt(tn) * 1024.0)]) % MAX_TRAVEL_TIME;
												Hosts[i2].PlaceLinks[P.HotelPlaceType] = l;
												Hosts[i2].Travelling = 1 + d2;
												nad++;
												j++;
											}
											f2++;
										} while ((f) && (f2 < 300));
										if (f)
										{
#pragma omp critical
											Hosts[i2].PlaceLinks[P.HotelPlaceType] = -1;
											if (++f3 > 100)
											{
												j++; nsk++;
											}
										}
									}
								}
							}
						}
						else
							j++;
					}
				}
		Files::xfprintf_stderr("<ar=%i as=%i", nad, nsk);
		nl = ((double)P.PlaceTypeMeanSize[P.HotelPlaceType]) * P.HotelPropLocal / P.MeanLocalJourneyTime;
		nsk = 0;
#pragma omp parallel for reduction(+:nld,nsk) schedule(static,1) default(none) \
			shared(P, Places, Cells, CellLookup, Hosts, Households, nl, bm, mps, d)
		for (int tn = 0; tn < P.NumThreads; tn++)
			for (int i = tn; i < P.Nplace[P.HotelPlaceType]; i += P.NumThreads)
			{
				int c = ((int)(Places[P.HotelPlaceType][i].loc.x / P.in_cells_.width)) * P.nch + ((int)(Places[P.HotelPlaceType][i].loc.y / P.in_cells_.height));
				int n = (int)ignpoi_mt(nl * Cells[c].tot_prob, tn);
				if (Places[P.HotelPlaceType][i].n + n > mps)
				{
					nsk += (Places[P.HotelPlaceType][i].n + n - mps);
					n = mps - Places[P.HotelPlaceType][i].n;
				}
				for (int j = 0; j < n; j++)
				{
					int f;
					do
					{
						f = 0;
						double s = ranf_mt(tn);
						int l = Cells[c].InvCDF[(int)floor(s * 1024)];
						while (Cells[c].cum_trans[l] < s) l++;
						Cell* ct = CellLookup[l];
						int m = (int)(ranf_mt(tn) * ((double)ct->S0));
						if (m < (ct->S + ct->L))
						{
							int i2 = ct->susceptible[m];
							int d2 = HOST_AGE_GROUP(i2);
							int f3 = 0;
							if ((Hosts[i2].Travelling == 0) && ((P.RelativeTravelRate[d2] == 1) || (ranf_mt(tn) < P.RelativeTravelRate[d2])))
							{
#pragma omp critical
								{if (Hosts[i2].PlaceLinks[P.HotelPlaceType] == -1) { Hosts[i2].PlaceLinks[P.HotelPlaceType] = -2; f3 = 1; }}
							}
							if (f3)
							{
								double s2 = dist2_raw(Households[Hosts[i2].hh].loc.x, Households[Hosts[i2].hh].loc.y, Places[P.HotelPlaceType][i].loc.x, Places[P.HotelPlaceType][i].loc.y);
								int f2 = 1;
								if ((bm) && (s2 > P.MoveRestrRadius2))
								{
									if (ranf_mt(tn) >= P.MoveRestrEffect)
									{
#pragma omp critical
										Hosts[i2].PlaceLinks[P.HotelPlaceType] = -1;
										nsk++;
										f2 = 0;
									}
								}
								if (f2)
								{
									s = P.KernelLookup.num(s2) / Cells[c].max_trans[l];
									if (ranf_mt(tn) >= s)
									{
#pragma omp critical
										Hosts[i2].PlaceLinks[P.HotelPlaceType] = -1;
										f = 1;
									}
									else
									{
										d2 = (d + P.InvLocalJourneyDurationDistrib[(int)(ranf_mt(tn) * 1024.0)]) % MAX_TRAVEL_TIME;
										int hp = Places[P.HotelPlaceType][i].n;
										Places[P.HotelPlaceType][i].n++;
										Places[P.HotelPlaceType][i].members[hp] = i2;
										Hosts[i2].Travelling = 1 + d2;
										nld++;
#pragma omp critical
										Hosts[i2].PlaceLinks[P.HotelPlaceType] = i;
									}
								}
							}
							else
								f = 1;
						}
						else
							nsk++;
					} while (f);
				}
			}
		Files::xfprintf_stderr(" l=%i ls=%i ", nld, nsk);
	}
}

void InfectSweep(double t, int run) // added run number as argument in order to record it in event log
{
	//// This function takes the day number (t) and run number (run) as inputs. It loops over infected people, and decides whom to infect. Structure is 1) #pragma loop over all cells then 1a) infectious people, which chooses who they will infect, adds them to a queue
	//// Next 2) #pragma loop infects those people from queue (using DoInfect function). This is to avoid race conditions.
	//// Loop 1a) calculates the force of infection exerted by each infected person on (and therefore number of new infections to) i) their house; ii) their place(s); iii) other spatial cells.
	//// Each force of infection includes infectiousness and susceptibility components.
	//// Infectiousness is (broadly) a function of 1 person (their age, treatment status, places, no. people in their household etc.)
	//// Susceptibility is (broadly) a function of 2 people (a person's susceptibility TO ANOTHER PERSON / potential infector)
	//// After loop 1a) over infectious people, spatial infections are doled out.

	int CellQueue;
	unsigned short int TimeStepNow = (unsigned short int) (P.TimeStepsPerDay * t); // TimeStepNow = the timestep number of the start of the current day
	int Day						= (int)floor(t); // integer value of current day (used for indexing).
	double fp					= P.ModelTimeStep / (1 - P.FalsePositiveRate); // fp = false positive
	double seasonality			= (P.DoSeasonality) ? (P.Seasonality[((int)t) % DAYS_PER_YEAR]) : 1.0; // if doing seasonality, pick seasonality from P.Seasonality array using day number in year. Otherwise set to 1.
	
	// Establish if movement restrictions are in place on current day - store in BlanketMoveRestrInPlace, 0:false, 1:true 
	int BlanketMoveRestrInPlace = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));
	// File for storing error reports
	FILE* stderr_shared = stderr;
	
#pragma omp parallel for private(CellQueue) schedule(static,1) default(none) \
		shared(t, P, CellLookup, Hosts, AdUnits, Households, Places, SamplingQueue, Cells, Mcells, StateT, seasonality, TimeStepNow, fp, BlanketMoveRestrInPlace, stderr_shared)
	for (int ThreadNum = 0; ThreadNum < P.NumThreads; ThreadNum++)
	{
		Cell* ThisCell; 
		double Household_Beta; // if doing households, Household_Beta = seasonality * fp * P.HouseholdTrans, else Household_Beta = 0
		double SpatialSeasonal_Beta;
		double SpatialInf_AllPeopleThisCell; ///// spatial infectiousness summed over all infectious people in loop below

		for (int CellIndex = ThreadNum; CellIndex < P.NumPopulatedCells; CellIndex += P.NumThreads) //// loop over (in parallel) all populated cells. Loop 1)
		{
			ThisCell = CellLookup[CellIndex]; // select Cell given by index CellIndex
			SpatialInf_AllPeopleThisCell = 0; ///// spatial infectiousness summed over all infectious people in loop below

			//// quantities that will be used in loop below.
			Person* InfectiousPerson;
			int InfectiousPersonIndex;
			bool DigiContactTrace_ThisPersonNow;

			//// Loop over array of indices of infectious people ThisCell->I in cell ThisCell. Loop 1a)
			for (int InfectiousPersonIndex_ThisCell = 0; InfectiousPersonIndex_ThisCell < ThisCell->I; InfectiousPersonIndex_ThisCell++)
			{
				//// get InfectiousPersonIndex from InfectiousPersonIndex_ThisCell
				InfectiousPersonIndex = ThisCell->infected[InfectiousPersonIndex_ThisCell];
				//// get InfectiousPerson from Hosts (array of people) corresponding to InfectiousPersonIndex, using pointer arithmetic.
				InfectiousPerson = Hosts + InfectiousPersonIndex;
				int AdUnit_ThisPerson = Mcells[InfectiousPerson->mcell].adunit;


				//evaluate flag for digital contact tracing (DigiContactTrace_ThisPersonNow) here at the beginning for each individual
				// DigiContactTrace_ThisPersonNow = 1 if:
				// P.DoDigitalContactTracing = 1 (ie. digital contact tracing functionlity is switched on)
				// AND Day number (t) is greater than the start day for contact tracing in this administrative unit (ie. contact tracing has started)
				// AND Day number (t) is less than the end day for contact tracing in this administrative unit (ie. contact tracing has not ended)
				// AND the selected host is a digital contact tracing user
				// otherwise DigiContactTrace_ThisPersonNow = 0
				DigiContactTrace_ThisPersonNow = ((P.DoDigitalContactTracing) &&
					(t >= AdUnits[AdUnit_ThisPerson].DigitalContactTracingTimeStart) &&
					(t <  AdUnits[AdUnit_ThisPerson].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) &&
					(Hosts[InfectiousPersonIndex].digitalContactTracingUser == 1)); // && (TimeStepNow <= (Hosts[InfectiousPersonIndex].detected_time + P.usCaseIsolationDelay)));

				// BEGIN HOUSEHOLD INFECTIONS


				Household_Beta = (P.DoHouseholds) ? P.Betas[Day][AdUnit_ThisPerson][House] * seasonality * fp : 0;
				if (Household_Beta > 0)
				{
					// For InfectiousPerson's household (InfectiousPerson->hh), 
					// if the number of hosts (nh) in that Household is greater than 1
					// AND the selected host is not travelling
					if ((Households[InfectiousPerson->hh].nh > 1) && (!InfectiousPerson->Travelling))
					{
						int FirstHouseholdMember = Households[InfectiousPerson->hh].FirstPerson;
						int LastHouseholdMember = FirstHouseholdMember + Households[InfectiousPerson->hh].nh;
						// calculate person's infectiousness at household-level
						// using the CalcHouseInf function on the selected cell and timestamp at start of current day
						// then scaling by Household_Beta
						double Household_Infectiousness = Household_Beta * CalcHouseInf(InfectiousPersonIndex, TimeStepNow);

						// Test if any of the individuals in the selected persons household are absent from places
						bool AtLeastOnePersonAbsent = false;
						for (int HouseholdMember = FirstHouseholdMember; (HouseholdMember < LastHouseholdMember) && (!AtLeastOnePersonAbsent); HouseholdMember++) //// loop over household memberts
							for (int PlaceType = 0; (PlaceType < P.NumPlaceTypes) && (!AtLeastOnePersonAbsent); PlaceType++) //// loop over place types
								if (Hosts[HouseholdMember].PlaceLinks[PlaceType] >= 0) //// if person in household has any sort of link to place type
									AtLeastOnePersonAbsent = ((PLACE_CLOSED(PlaceType, Hosts[HouseholdMember].PlaceLinks[PlaceType])) && (HOST_ABSENT(HouseholdMember)));

						// if individuals in the household are absent from places (ie. AtLeastOnePersonAbsent from test immediately above), scale up the infectiousness (Household_Infectiousness) of the household
						if (AtLeastOnePersonAbsent) { Household_Infectiousness *= P.Efficacies[PlaceClosure][House]; }/* NumPCD++;}*/ //// if people in your household are absent from places, person InfectiousPerson/InfectiousPersonIndex is more infectious to them, as they spend more time at home.

						// Loop over household members
						for (int HouseholdMember = FirstHouseholdMember; HouseholdMember < LastHouseholdMember; HouseholdMember++) //// loop over all people in household 
						{
							if (Hosts[HouseholdMember].is_susceptible() && (!Hosts[HouseholdMember].Travelling)) //// if people in household uninfected/susceptible and not travelling
							{
								double Household_FOI = Household_Infectiousness * CalcHouseSusc(HouseholdMember, TimeStepNow, InfectiousPersonIndex);		//// Household force of infection (FOI = infectiousness x susceptibility) from person InfectiousPersonIndex/InfectiousPerson on fellow household member

								// Force of Infection (Household_FOI) > random value between 0 and 1
								if (ranf_mt(ThreadNum) < Household_FOI)
								{
									// explicitly cast to short to resolve level 4 warning
									const short int infect_type = static_cast<short int>(1 + INFECT_TYPE_MASK * (1 + InfectiousPerson->infect_type / INFECT_TYPE_MASK));

									if (AddInfections(ThreadNum, Hosts[HouseholdMember].pcell % P.NumThreads, InfectiousPersonIndex, HouseholdMember, infect_type))
										Hosts[HouseholdMember].infector = InfectiousPersonIndex; //// assign InfectiousPersonIndex as infector of person HouseholdMember
								} // if FOI > random value between 0 and 1
							} // if person in household uninfected/susceptible and not travelling
						} // loop over people in household
					} // if more than one person in household 
				}// if Household_Beta > 0
				// END HOUSHOLD INFECTIONS

				// BEGIN PLACE INFECTIONS
				if (P.DoPlaces) // if places functionality is enabled
				{
					if (!HOST_ABSENT(InfectiousPersonIndex))
					{
						// select microcell (Microcell_ThisPerson) corresponding to selected InfectiousPerson
						Microcell* Microcell_ThisPerson = Mcells + InfectiousPerson->mcell;
						for (int PlaceType = 0; PlaceType < P.NumPlaceTypes; PlaceType++) //// loop over all place types
						{
							// select PlaceLink between selected InfectiousPerson and place from InfectiousPerson's placelinks to place type PlaceType
							int PlaceLink = InfectiousPerson->PlaceLinks[PlaceType];
							if (PlaceLink >= 0)  //// PlaceLink >= 0 means if place type PlaceType is relevant to InfectiousPerson. (Now allowing for partial attendance).
							{
								// Place_Infectiousness
								double Place_Infectiousness = CalcPlaceInf(InfectiousPersonIndex, PlaceType, TimeStepNow);
								Place_Infectiousness		*= P.Betas[Day][AdUnit_ThisPerson][PlaceType] * fp * seasonality / P.PlaceTypeGroupSizeParam1[PlaceType];
								// select microcell of the place linked to InfectiousPerson with link PlaceLink
								Microcell* Microcell_ThisPersonsPlaceLink = Mcells + Places[PlaceType][PlaceLink].mcell;
								// if blanket movement restrictions are in place on current day
								if (BlanketMoveRestrInPlace)
								{
									// if distance between InfectiousPerson's household and linked place
									// is greater than movement restriction radius
									if ((dist2_raw(Households[InfectiousPerson->hh].loc.x, Households[InfectiousPerson->hh].loc.y,
										Places[PlaceType][PlaceLink].loc.x, Places[PlaceType][PlaceLink].loc.y) > P.MoveRestrRadius2))
										Place_Infectiousness *= P.MoveRestrEffect; // multiply infectiousness of place by movement restriction effect
								}
								// else if movement restrictions in effect in either household microcell or place microcell
								else if ((Microcell_ThisPerson->moverest != Microcell_ThisPersonsPlaceLink->moverest) && ((Microcell_ThisPerson->moverest == TreatStat::Treated) || (Microcell_ThisPersonsPlaceLink->moverest == TreatStat::Treated)))
								{
									Place_Infectiousness *= P.MoveRestrEffect; // multiply infectiousness of place by movement restriction effect
								}

								// BEGIN NON-HOTEL INFECTIONS

								// if linked place isn't a hotel and selected host isn't travelling
								if ((PlaceType != P.HotelPlaceType) && (!InfectiousPerson->Travelling))
								{
									// PlaceGroupLink_index is index of group (of place type PlaceType) that selected host is linked to 
									int PlaceGroupLink_index = (InfectiousPerson->PlaceGroupLinks[PlaceType]);

									// calculate infectiousness (PlaceInfectiousness_Scaled_DCT_copy) 
									// which varies if contact tracing is in place
									// if contact tracing isn't in place PlaceInfectiousness_Scaled_DCT_copy is a copy of Place_Infectiousness 
									// if contact tracing is in place, PlaceInfectiousness_Scaled_DCT_copy is Place_Infectiousness  * P.ScalingFactorPlaceDigitalContacts
									// in either case PlaceInfectiousness_Scaled_DCT_copy is capped at 1

									// if contact tracing
									double Place_Infectiousness_DCT_copy, PlaceInfectiousness_Scaled_DCT_copy;
									if (DigiContactTrace_ThisPersonNow)
									{
										// copy Place_Infectiousness
										Place_Infectiousness_DCT_copy = Place_Infectiousness;
										// multiply Place_Infectiousness_DCT_copy by P.ScalingFactorPlaceDigitalContacts 
										PlaceInfectiousness_Scaled_DCT_copy = Place_Infectiousness_DCT_copy * P.ScalingFactorPlaceDigitalContacts;
										// cap Place_Infectiousness_DCT_copy at 1
										if (Place_Infectiousness_DCT_copy > 1) Place_Infectiousness_DCT_copy = 1;
										// cap at 1
										if (PlaceInfectiousness_Scaled_DCT_copy > 1) PlaceInfectiousness_Scaled_DCT_copy = 1;
									}
									else
									{
										// copy Place_Infectiousness to Place_Infectiousness_DCT_copy
										Place_Infectiousness_DCT_copy = Place_Infectiousness;
										// cap Place_Infectiousness_DCT_copy at 1
										if (Place_Infectiousness_DCT_copy > 1) Place_Infectiousness_DCT_copy = 1;
										PlaceInfectiousness_Scaled_DCT_copy = Place_Infectiousness_DCT_copy;
									}

									// if infectiousness is < 0, we have an error - end the program
									int NumPotentialInfecteesPlaceGroup;
									if (PlaceInfectiousness_Scaled_DCT_copy < 0)
									{
										Files::xfprintf(stderr_shared, "@@@ %lg\n", PlaceInfectiousness_Scaled_DCT_copy);
										exit(1);
									}
									// else if infectiousness == 1 (should never be more than 1 due to capping above)
									else if (PlaceInfectiousness_Scaled_DCT_copy >= 1)	//// if place infectiousness above threshold, consider everyone in group a potential infectee...
									{
										// set NumPotentialInfecteesPlaceGroup to be number of people in group in place PlaceType,PlaceLink
										NumPotentialInfecteesPlaceGroup = Places[PlaceType][PlaceLink].group_size[PlaceGroupLink_index];
									}
									else				//// ... otherwise randomly sample (from binomial distribution) number of potential infectees in this place.
									{
										NumPotentialInfecteesPlaceGroup = (int)ignbin_mt((int32_t)Places[PlaceType][PlaceLink].group_size[PlaceGroupLink_index], PlaceInfectiousness_Scaled_DCT_copy, ThreadNum);
									}

									// if potential infectees > 0	
									if (NumPotentialInfecteesPlaceGroup > 0)
									{
										// pick NumPotentialInfecteesPlaceGroup members of place PlaceType,PlaceLink and add them to sampling queue for thread ThreadNum
										SampleWithoutReplacement(ThreadNum, NumPotentialInfecteesPlaceGroup, Places[PlaceType][PlaceLink].group_size[PlaceGroupLink_index]); //// changes thread-specific SamplingQueue.
									}

									// loop over sampling queue of potential infectees
									for (int PotentialInfectee_ThisPlaceGroupIndex = 0; PotentialInfectee_ThisPlaceGroupIndex < NumPotentialInfecteesPlaceGroup; PotentialInfectee_ThisPlaceGroupIndex++)
									{
										// pick potential infectee index PotentialInfectee_PlaceGroup
										int PotentialInfectee_PlaceGroup = Places[PlaceType][PlaceLink].members[Places[PlaceType][PlaceLink].group_start[PlaceGroupLink_index] + SamplingQueue[ThreadNum][PotentialInfectee_ThisPlaceGroupIndex]];
										// calculate place susceptbility based on infectee (PotentialInfectee_PlaceGroup), place type (PlaceType), timestep (TimeStepNow)
										// thread number (ThreadNum)
										double PlaceSusceptibility = CalcPlaceSusc(PotentialInfectee_PlaceGroup, PlaceType, TimeStepNow);

										// allow care home residents to mix more intensely in "groups" (i.e. individual homes) than staff do - to allow for PPE/environmental contamination.
										if ((PlaceType == P.CareHomePlaceType) && ((!Hosts[InfectiousPersonIndex].care_home_resident) || (!Hosts[PotentialInfectee_PlaceGroup].care_home_resident))) PlaceSusceptibility *= P.CareHomeWorkerGroupScaling;
										//these are all place group contacts to be tracked for digital contact tracing - add to StateT queue for contact tracing
										//if infectee is also a user, add them as a contact

										if ((DigiContactTrace_ThisPersonNow) && (Hosts[PotentialInfectee_PlaceGroup].digitalContactTracingUser) && (InfectiousPersonIndex != PotentialInfectee_PlaceGroup) && (!HOST_ABSENT(PotentialInfectee_PlaceGroup)))
										{
											// scale place susceptibility by proportion who self isolate and store as PlaceSusceptibility_DCT_scaled
											double PlaceSusceptibility_DCT_scaled = P.ProportionDigitalContactsIsolate * PlaceSusceptibility;
											// if random number < PlaceSusceptibility_DCT_scaled
											// AND number of contacts of InfectiousPersonIndex(!) is less than maximum digital contact to trace
											if ((Hosts[InfectiousPersonIndex].ncontacts < P.MaxDigitalContactsToTrace) && (ranf_mt(ThreadNum) < PlaceSusceptibility_DCT_scaled))
											{
												Hosts[InfectiousPersonIndex].ncontacts++; //add to number of contacts made
												int AdUnit = Mcells[Hosts[PotentialInfectee_PlaceGroup].mcell].adunit;
												if ((StateT[ThreadNum].ndct_queue[AdUnit] < AdUnits[AdUnit].n))
												{
													//find adunit for contact and add both contact and infectious host to lists - storing both so I can set times later.
													StateT[ThreadNum].dct_queue[AdUnit][StateT[ThreadNum].ndct_queue[AdUnit]++] = { PotentialInfectee_PlaceGroup, InfectiousPersonIndex, TimeStepNow };
												}
												else
													Files::xfprintf(stderr_shared, "No more space in queue! Thread: %i, AdUnit: %i\n", ThreadNum, AdUnit);
											}
										}

										if (Hosts[PotentialInfectee_PlaceGroup].is_susceptible() && (!HOST_ABSENT(PotentialInfectee_PlaceGroup))) //// if person PotentialInfectee_PlaceGroup uninfected and not absent.
										{
											Microcell* MicroCell_PotentialInfectee_PlaceGroup = Mcells + Hosts[PotentialInfectee_PlaceGroup].mcell;
											//downscale PlaceSusceptibility if it has been scaled up do to digital contact tracing
											PlaceSusceptibility *= CalcPersonSusc(PotentialInfectee_PlaceGroup, TimeStepNow, InfectiousPersonIndex) * Place_Infectiousness_DCT_copy / PlaceInfectiousness_Scaled_DCT_copy;

											// if blanket movement restrictions are in place
											if (BlanketMoveRestrInPlace)
											{
												// if potential infectee PotentialInfectee_PlaceGroup's household is further from selected place
												if ((dist2_raw(Households[Hosts[PotentialInfectee_PlaceGroup].hh].loc.x, Households[Hosts[PotentialInfectee_PlaceGroup].hh].loc.y,
													Places[PlaceType][PlaceLink].loc.x, Places[PlaceType][PlaceLink].loc.y) > P.MoveRestrRadius2))
													PlaceSusceptibility *= P.MoveRestrEffect; // multiply susceptibility by movement restriction effect
											}
											// else if movement restrictions are in place in either cell
											else if ((MicroCell_PotentialInfectee_PlaceGroup->moverest != Microcell_ThisPersonsPlaceLink->moverest) && ((MicroCell_PotentialInfectee_PlaceGroup->moverest == TreatStat::Treated) || (Microcell_ThisPersonsPlaceLink->moverest == TreatStat::Treated)))
											{
												// multiply susceptibility by movement restriction effect
												PlaceSusceptibility *= P.MoveRestrEffect;
											}

											// if either susceptiblity is 100% or sample probability PlaceSusceptibility
											if ((PlaceSusceptibility == 1) || (ranf_mt(ThreadNum) < PlaceSusceptibility))
											{
												// explicitly cast to short to resolve level 4 warning
												const short int infect_type = static_cast<short int> (2 + PlaceType + INFECT_TYPE_MASK * (1 + InfectiousPerson->infect_type / INFECT_TYPE_MASK));

												AddInfections(ThreadNum, Hosts[PotentialInfectee_PlaceGroup].pcell % P.NumThreads, InfectiousPersonIndex, PotentialInfectee_PlaceGroup, infect_type);
											}
										}
									}
								}
								// END NON-HOTEL INFECTIONS

								// BEGIN HOTEL INFECTIONS
								// if InfectiousPerson is not travelling or selected link is to a hotel
								if ((PlaceType == P.HotelPlaceType) || (!InfectiousPerson->Travelling))
								{
									Place_Infectiousness *= P.PlaceTypePropBetweenGroupLinks[PlaceType] * P.PlaceTypeGroupSizeParam1[PlaceType] / ((double)Places[PlaceType][PlaceLink].n);
									if (Place_Infectiousness > 1) Place_Infectiousness = 1;
									// if contact tracing in place, multiply Place_Infectiousness_scaled = Place_Infectiousness*scalingfactor, otherwise Place_Infectiousness_scaled = Place_Infectiousness
									double Place_Infectiousness_scaled = (DigiContactTrace_ThisPersonNow) ? (Place_Infectiousness * P.ScalingFactorPlaceDigitalContacts) : Place_Infectiousness;
									// Place_Infectiousness_scaled shouldn't be less than 0 so generate error if it is
									int NumPotentialInfecteesHotel;
									if (Place_Infectiousness_scaled < 0)
									{
										ERR_CRITICAL_FMT("@@@ %lg\n", Place_Infectiousness);
									}
									// if Place_Infectiousness_scaled >=1, everyone in the hotel is a potential infectee
									else if (Place_Infectiousness_scaled >= 1)
										NumPotentialInfecteesHotel = Places[PlaceType][PlaceLink].n;
									// if Place_Infectiousness_scaled between 0 and 1, decide number of potential infectees based on
									// using ignbin_mt function
									else
										NumPotentialInfecteesHotel = (int)ignbin_mt((int32_t)Places[PlaceType][PlaceLink].n, Place_Infectiousness_scaled, ThreadNum);
									// if more than 0 potential infectees, pick n hosts from the hotel and add to sampling queue
									if (NumPotentialInfecteesHotel > 0) SampleWithoutReplacement(ThreadNum, NumPotentialInfecteesHotel, Places[PlaceType][PlaceLink].n);
									// loop over the sampling queue
									for (int PotentialInfectee_HotelIndex = 0; PotentialInfectee_HotelIndex < NumPotentialInfecteesHotel; PotentialInfectee_HotelIndex++)
									{
										// select potential infectee from sampling queue
										int PotentialInfectee_Hotel = Places[PlaceType][PlaceLink].members[SamplingQueue[ThreadNum][PotentialInfectee_HotelIndex]];
										// calculate place susceptibility PlaceSusceptibility
										double PlaceSusceptibility = CalcPlaceSusc(PotentialInfectee_Hotel, PlaceType, TimeStepNow);
										// use group structure to model multiple care homes with shared staff - in which case residents of one "group" don't mix with those in another, only staff do.
										if ((Hosts[InfectiousPersonIndex].care_home_resident) && (Hosts[PotentialInfectee_Hotel].care_home_resident) && (Hosts[InfectiousPersonIndex].PlaceGroupLinks[PlaceType] != Hosts[PotentialInfectee_Hotel].PlaceGroupLinks[PlaceType])) PlaceSusceptibility *= P.CareHomeResidentPlaceScaling;
										// allow care home staff to have lowere contacts in care homes - to allow for PPE/environmental contamination.
										if ((PlaceType == P.CareHomePlaceType) && ((!Hosts[InfectiousPersonIndex].care_home_resident) || (!Hosts[PotentialInfectee_Hotel].care_home_resident))) PlaceSusceptibility *= P.CareHomeWorkerGroupScaling;

										//these are all place group contacts to be tracked for digital contact tracing - add to StateT queue for contact tracing

										//if infectee is also a user, add them as a contact

										// if contact tracing in place AND potential infectee PotentialInfectee_Hotel is a contact tracing user AND PotentialInfectee_Hotel isn't absent AND PotentialInfectee_Hotel isn't InfectiousPersonIndex (suspect this should be si)

										if ((DigiContactTrace_ThisPersonNow) && (Hosts[PotentialInfectee_Hotel].digitalContactTracingUser) && (InfectiousPersonIndex != PotentialInfectee_Hotel) && (!HOST_ABSENT(PotentialInfectee_Hotel)))
										{
											// PlaceSusceptibility_DCT_scaled = place susceptibility * proportion of digital contacts who self isolate
											double PlaceSusceptibility_DCT_scaled = P.ProportionDigitalContactsIsolate * PlaceSusceptibility;
											// if number of contacts of infectious person < maximum and random number < PlaceSusceptibility_DCT_scaled
											if ((Hosts[InfectiousPersonIndex].ncontacts < P.MaxDigitalContactsToTrace) && (ranf_mt(ThreadNum) < PlaceSusceptibility_DCT_scaled))
											{
												Hosts[InfectiousPersonIndex].ncontacts++; //add to number of contacts made
												int ad = Mcells[Hosts[PotentialInfectee_Hotel].mcell].adunit;
												// find adunit for contact and add both contact and infectious host to lists - storing both so I can set times later.
												if ((StateT[ThreadNum].ndct_queue[ad] < AdUnits[ad].n))
													StateT[ThreadNum].dct_queue[ad][StateT[ThreadNum].ndct_queue[ad]++] = { PotentialInfectee_Hotel, InfectiousPersonIndex, TimeStepNow };
												else Files::xfprintf(stderr_shared, "No more space in queue! Thread: %i, AdUnit: %i\n", ThreadNum, ad);
											}
										}

										// if potential infectee PotentialInfectee_Hotel uninfected and not absent.
										if (Hosts[PotentialInfectee_Hotel].is_susceptible() && (!HOST_ABSENT(PotentialInfectee_Hotel)))
										{
											// MicroCell_PotentialInfectee_Hotel = microcell of potential infectee
											Microcell* MicroCell_PotentialInfectee_Hotel = Mcells + Hosts[PotentialInfectee_Hotel].mcell;

											//if doing digital contact tracing, scale down susceptibility here
											PlaceSusceptibility *= CalcPersonSusc(PotentialInfectee_Hotel, TimeStepNow, InfectiousPersonIndex) * Place_Infectiousness / Place_Infectiousness_scaled;
											// if blanket movement restrictions are in place
											if (BlanketMoveRestrInPlace)
											{
												// if potential infectees household is farther away from hotel than restriction radius
												if ((dist2_raw(Households[Hosts[PotentialInfectee_Hotel].hh].loc.x, Households[Hosts[PotentialInfectee_Hotel].hh].loc.y,
													Places[PlaceType][PlaceLink].loc.x, Places[PlaceType][PlaceLink].loc.y) > P.MoveRestrRadius2))
												{
													// multiply susceptibility by movement restriction effect
													PlaceSusceptibility *= P.MoveRestrEffect;
												}
											}
											// else if movement restrictions are in place in potential infectee's cell or hotel's cell
											else if ((MicroCell_PotentialInfectee_Hotel->moverest != Microcell_ThisPersonsPlaceLink->moverest) && ((MicroCell_PotentialInfectee_Hotel->moverest == TreatStat::Treated) || (Microcell_ThisPersonsPlaceLink->moverest == TreatStat::Treated)))
											{
												// multiply susceptibility by movement restriction effect
												PlaceSusceptibility *= P.MoveRestrEffect;
											}

											// ** do infections **

											// is susceptibility is 1 (ie infect everyone) or random number is less than susceptibility
											if ((PlaceSusceptibility == 1) || (ranf_mt(ThreadNum) < PlaceSusceptibility))
											{
												// explicitly cast to short to resolve level 4 warning
												const short int infect_type = static_cast<short int> (2 + PlaceType + MAX_NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + InfectiousPerson->infect_type / INFECT_TYPE_MASK));

												AddInfections(ThreadNum, Hosts[PotentialInfectee_Hotel].pcell % P.NumThreads, InfectiousPersonIndex, PotentialInfectee_Hotel, infect_type);
											} // susceptibility test
										} // PotentialInfectee_Hotel uninfected and not absent.
									} // loop over sampling queue
								} // selected host InfectiousPerson is not travelling or selected link is to a hotel

								// ** END HOTEL INFECTIONS ** 

							} // if place link relevant
						} // loop over place types
					} // if host isn't absent
				} // if places functionality enabled

				// END PLACE INFECTIONS

				// BEGIN SPATIAL INFECTIONS

				//// First determine spatial FOI component (SpatialInf_AllPeopleThisCell)

				// if seasonality beta > 0
				// do spatial infections 
				//// ie sum spatial infectiousness over all infected people, the infections from which are allocated after loop over infected people.

				SpatialSeasonal_Beta = seasonality * fp * P.Betas[Day][AdUnit_ThisPerson][Spatial];
				if (SpatialSeasonal_Beta > 0)
				{
					double SpatialInf_ThisPerson;
					if (InfectiousPerson->Travelling) //// if host currently away from their cell, they cannot add to their cell's spatial infectiousness.
						SpatialInf_ThisPerson = 0;
					else
					{
						// calculate SpatialInf_ThisPerson based on host and timestep
						SpatialInf_ThisPerson = CalcSpatialInf(InfectiousPersonIndex, TimeStepNow);
						// if do digital contact tracing, scale up spatial infectiousness of infectives who are using the app and will be detected
						if (DigiContactTrace_ThisPersonNow)
							SpatialInf_ThisPerson *= P.ScalingFactorSpatialDigitalContacts;
					}
					// test if selected person InfectiousPerson is linked to a place that is closed, f=0 means no links to closed places, otherwise f=1
					int PlaceClosedFlag = 0; // initialise f as 0

											 // If place functionality switched on
					if (P.DoPlaces)
						for (int PlaceType = 0; (PlaceType < P.NumPlaceTypes) && (!PlaceClosedFlag); PlaceType++) // loop over place types until closed place is found
							if (InfectiousPerson->PlaceLinks[PlaceType] >= 0) //// if person has a link to place of type PlaceType...
								PlaceClosedFlag = PLACE_CLOSED(PlaceType, InfectiousPerson->PlaceLinks[PlaceType]); //// find out if that place of type PlaceType is closed.

					if ((PlaceClosedFlag) && (HOST_ABSENT(InfectiousPersonIndex))) //// if place is closed and person is absent then adjust the spatial infectiousness (similar logic to household infectiousness: place closure affects spatial infectiousness
					{
						SpatialInf_ThisPerson *= P.Efficacies[PlaceClosure][Spatial];
						/* NumPCD++; */
						SpatialInf_AllPeopleThisCell += SpatialInf_ThisPerson;
						StateT[ThreadNum].cell_inf[InfectiousPersonIndex_ThisCell] = (float)-SpatialInf_AllPeopleThisCell;
					}
					else
					{
						SpatialInf_AllPeopleThisCell += SpatialInf_ThisPerson;
						StateT[ThreadNum].cell_inf[InfectiousPersonIndex_ThisCell] = (float)SpatialInf_AllPeopleThisCell;
					}
				}
			} // loop over infectious people in cell


			//// Now allocate spatial infections using Force Of Infection (SpatialInf_AllPeopleThisCell) calculated above
			if (SpatialInf_AllPeopleThisCell > 0) //// if spatial infectiousness positive
			{
				// decide how many potential cell to cell infections this cell could cause  
				int NumPotentialCelltoCellInfections = (int)ignpoi_mt(SpatialInf_AllPeopleThisCell * SpatialSeasonal_Beta * ((double)ThisCell->tot_prob), ThreadNum); //// number people this cell's population might infect elsewhere. poisson random number based on spatial infectiousness s5, SpatialSeasonal_Beta (seasonality) and this cell's "probability" (guessing this is a function of its population and geographical size).
				// NumInfectiousThisCell = number of infectious people in cell ThisCell
				int NumInfectiousThisCell = ThisCell->I;

				if (NumPotentialCelltoCellInfections > 0) //// this block normalises cumulative infectiousness cell_inf by person. s5 is the total cumulative spatial infectiousness. Reason is so that infector can be chosen using ranf_mt, which returns random number between 0 and 1.
				{
					//// normalise by cumulative spatial infectiousness.
					for (int InfectiousPerson = 0; InfectiousPerson < NumInfectiousThisCell - 1; InfectiousPerson++) StateT[ThreadNum].cell_inf[InfectiousPerson] /= ((float)SpatialInf_AllPeopleThisCell);
					//// does same as the above loop just a slightly faster calculation. i.e. StateT[ThreadNum].cell_inf[NumInfectiousThisCell - 1] / SpatialInf_AllPeopleThisCell would equal 1 or -1 anyway.
					StateT[ThreadNum].cell_inf[NumInfectiousThisCell - 1] = (StateT[ThreadNum].cell_inf[NumInfectiousThisCell - 1] < 0) ? -1.0f : 1.0f;
				}

				//// loop over infections to dole out. roughly speaking, this determines which infectious person in cell ThisCell infects which person elsewhere.
				for (int Infection = 0; Infection < NumPotentialCelltoCellInfections; Infection++)
				{
					//// decide on infector PotentialInfector_Spatial/InfectiousPerson from cell ThisCell.
					int PotentialInfector_CellIndex; // PotentialInfector_CellIndex = index of infector

					if (NumInfectiousThisCell == 1) PotentialInfector_CellIndex = 0; // if only one infectious person in cell, infector index is first in cell (person 0)
					// if more than one infectious person in cell pick an infectious person (given by index PotentialInfector_CellIndex)
					//// roughly speaking, this determines which infectious person in cell ThisCell infects which person elsewhere
					else
					{
						int StepSize;
						double RandomNum = ranf_mt(ThreadNum);	///// choose random number between 0 and 1
						PotentialInfector_CellIndex = StepSize = NumInfectiousThisCell / 2;		///// assign PotentialInfector_CellIndex and StepSize to be halfway between zero and NumInfectiousThisCell.
						int KeepSearching = 1;
						do
						{
							if (StepSize > 1) StepSize /= 2; //// amount StepSize to change PotentialInfector_CellIndex by reduced by half. Basically a binary search saying, keep amending potential infector PotentialInfector_CellIndex until either PotentialInfector_CellIndex less than zero or more than number of infected people until you find PotentialInfector_CellIndex s.t. spatial infectiousness "matches" RandomNum.
							if ((PotentialInfector_CellIndex > 0) && (fabs(StateT[ThreadNum].cell_inf[PotentialInfector_CellIndex - 1]) >= RandomNum))
							{
								PotentialInfector_CellIndex -= StepSize;
								if (PotentialInfector_CellIndex == 0)							KeepSearching = 0;
							}
							else if ((PotentialInfector_CellIndex < NumInfectiousThisCell - 1) && (fabs(StateT[ThreadNum].cell_inf[PotentialInfector_CellIndex]) < RandomNum))
							{
								PotentialInfector_CellIndex += StepSize;
								if (PotentialInfector_CellIndex == NumInfectiousThisCell - 1)	KeepSearching = 0;
							}
							else									KeepSearching = 0;
						} while (KeepSearching);
					}
					int InfectorPlaceClosedFlag = (StateT[ThreadNum].cell_inf[PotentialInfector_CellIndex] < 0); //// flag for whether infector PotentialInfector_CellIndex had their place(s) closed. <0 (true) = place closed / >=0 (false) = place not closed. Set in if (SpatialSeasonal_Beta > 0) part of loop over infectious people.

					int PotentialInfector_Index = ThisCell->infected[PotentialInfector_CellIndex];
					// PotentialInfector_Spatial is the jth infected person in the cell
					Person* PotentialInfector_Spatial = Hosts + PotentialInfector_Index;

					//calculate flag (DigiContactTrace_ThisPersonNow) for digital contact tracing here at the beginning for each individual infector
					bool DigiContactTrace_ThisPersonNow = ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[PotentialInfector_Spatial->mcell].adunit].DigitalContactTracingTimeStart)
						&& (t < AdUnits[Mcells[PotentialInfector_Spatial->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[PotentialInfector_Index].digitalContactTracingUser == 1)); // && (TimeStepNow <= (Hosts[PotentialInfector_Spatial].detected_time + P.usCaseIsolationDelay)));

					//// decide on infectee

					int KeepSearchingForCellToInfect = 1;  // do the following while KeepSearchingForCellToInfect = 1
					do
					{
						//// chooses which cell person will infect
						// pick RandomNum between 0 and 1
						double RandomNum = ranf_mt(ThreadNum);
						// generate l using InvCDF of selected cell and random integer between 0 and 1024
						int l = ThisCell->InvCDF[(int)floor(RandomNum * 1024)];
						// loop over ThisCell->cum_trans array until find a value >= RandomNum
						while (ThisCell->cum_trans[l] < RandomNum) l++;
						// selecte the cell corresponding to l
						Cell* ct = CellLookup[l];

						///// pick random person SusceptiblePerson within susceptibles of cell ct (S0 initial number susceptibles within cell).
						int SusceptiblePerson = (int)(ranf_mt(ThreadNum) * ((double)ct->S0));
						int PotentialInfectee_Spatial = ct->susceptible[SusceptiblePerson];

						double PotentialInfectorInfecteeDistSquared = dist2(Hosts + PotentialInfectee_Spatial, Hosts + PotentialInfector_Index); /// calculate distance squared between PotentialInfectee_Spatial and PotentialInfector_Spatial
						double AcceptProb = P.KernelLookup.num(PotentialInfectorInfecteeDistSquared) / ThisCell->max_trans[l]; //// acceptance probability

						// initialise KeepSearchingForCellToInfect = 0 (KeepSearchingForCellToInfect = 1 is the while condition for this loop)
						KeepSearchingForCellToInfect = 0;
						// if random number greater than acceptance probablility or infectee is dead
						if ((ranf_mt(ThreadNum) >= AcceptProb) || Hosts[PotentialInfectee_Spatial].is_dead()) //// if rejected, or infectee PotentialInfectee_Spatial/SusceptiblePerson already dead, ensure do-while evaluated again (i.e. choose a new infectee).
						{
							// set KeepSearchingForCellToInfect = 1 so loop continues (i.e. another PotentialInfectee_Spatial will be chosen)
							KeepSearchingForCellToInfect = 1;
						}
						else
						{
							//// if potential infectee not travelling, and either is not part of cell ThisCell or doesn't share a household with infector.
							if ((!Hosts[PotentialInfectee_Spatial].Travelling) && ((ThisCell != ct) || (Hosts[PotentialInfectee_Spatial].hh != PotentialInfector_Spatial->hh)))
							{
								// pick microcell of infector (Microcell_PotentialInfector)
								Microcell* Microcell_PotentialInfector = Mcells + PotentialInfector_Spatial->mcell;
								// pick microcell of infectee (mt)
								Microcell* mt = Mcells + Hosts[PotentialInfectee_Spatial].mcell;
								double Spatial_Susc = CalcSpatialSusc(PotentialInfectee_Spatial, TimeStepNow);
								// Care home residents may have fewer contacts
								if ((Hosts[PotentialInfectee_Spatial].care_home_resident) || (Hosts[PotentialInfector_Index].care_home_resident)) Spatial_Susc *= P.CareHomeResidentSpatialScaling;
								//so this person is a contact - but might not be infected. if we are doing digital contact tracing, we want to add the person to the contacts list, if both are users
								if (DigiContactTrace_ThisPersonNow)
								{
									//if infectee is also a user, add them as a contact
									if (Hosts[PotentialInfectee_Spatial].digitalContactTracingUser && (PotentialInfector_Index != PotentialInfectee_Spatial))
									{
										if ((Hosts[PotentialInfector_Index].ncontacts < P.MaxDigitalContactsToTrace) && (ranf_mt(ThreadNum) < Spatial_Susc * P.ProportionDigitalContactsIsolate))
										{
											Hosts[PotentialInfector_Index].ncontacts++; //add to number of contacts made
											int ad = Mcells[Hosts[PotentialInfectee_Spatial].mcell].adunit;
											if ((StateT[ThreadNum].ndct_queue[ad] < AdUnits[ad].n))
											{
												//find adunit for contact and add both contact and infectious host to lists - storing both so I can set times later.
												StateT[ThreadNum].dct_queue[ad][StateT[ThreadNum].ndct_queue[ad]++] = { PotentialInfectee_Spatial, PotentialInfector_Index, TimeStepNow };
											}
											else
											{
												Files::xfprintf(stderr_shared, "No more space in queue! Thread: %i, AdUnit: %i\n", ThreadNum, ad);
											}
										}
									}
									//scale down susceptibility so we don't over accept
									Spatial_Susc /= P.ScalingFactorSpatialDigitalContacts;
								}


								if (SusceptiblePerson < ct->S)  // only consider susceptible people as possible infectees
								{
									Spatial_Susc *= CalcPersonSusc(PotentialInfectee_Spatial, TimeStepNow, PotentialInfector_Index);
									Spatial_Susc *= P.WAIFW_Matrix_SpatialOnly[HOST_AGE_GROUP(PotentialInfectee_Spatial)][HOST_AGE_GROUP(PotentialInfector_Index)]; //// 
									if (BlanketMoveRestrInPlace)
									{
										if ((dist2_raw(Households[PotentialInfector_Spatial->hh].loc.x, Households[PotentialInfector_Spatial->hh].loc.y,
											Households[Hosts[PotentialInfectee_Spatial].hh].loc.x, Households[Hosts[PotentialInfectee_Spatial].hh].loc.y) > P.MoveRestrRadius2))
											Spatial_Susc *= P.MoveRestrEffect;
									}
									else if ((mt->moverest != Microcell_PotentialInfector->moverest) && ((mt->moverest == TreatStat::Treated) || (Microcell_PotentialInfector->moverest == TreatStat::Treated)))
										Spatial_Susc *= P.MoveRestrEffect;
									if ((!InfectorPlaceClosedFlag) && (HOST_ABSENT(PotentialInfectee_Spatial))) //// if infector did not have place closed, loop over place types of PotentialInfectee_Spatial to see if their places had closed. If they had, amend their susceptibility.
									{
										for (int PlaceType = KeepSearchingForCellToInfect = 0; (PlaceType < P.NumPlaceTypes) && (!KeepSearchingForCellToInfect); PlaceType++)
											if (Hosts[PotentialInfectee_Spatial].PlaceLinks[PlaceType] >= 0)
												KeepSearchingForCellToInfect = PLACE_CLOSED(PlaceType, Hosts[PotentialInfectee_Spatial].PlaceLinks[PlaceType]);
										if (KeepSearchingForCellToInfect) { Spatial_Susc *= P.Efficacies[PlaceClosure][Spatial]; } /* NumPCD++;} */
										KeepSearchingForCellToInfect = 0;
									}
									if ((Spatial_Susc == 1) || (ranf_mt(ThreadNum) < Spatial_Susc)) //// accept/reject
									{
										CellQueue = ((int)(ct - Cells)) % P.NumThreads;

										if (Hosts[PotentialInfectee_Spatial].is_susceptible())
										{
											// explicitly cast to short to resolve level 4 warning
											const short int infect_type = static_cast<short int>(2 + 2 * MAX_NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + PotentialInfector_Spatial->infect_type / INFECT_TYPE_MASK));

											AddInfections(ThreadNum, CellQueue, PotentialInfector_Index, PotentialInfectee_Spatial, infect_type);
										}
									}
								} // SusceptiblePerson < susceptible people in target cell
							} // //// if potential infectee not travelling, and either is not part of cell ThisCell or doesn't share a household with infector
						} // infectee isn't dead
					} while (KeepSearchingForCellToInfect);
				} // loop over infections doled out by cell
			} // SpatialInf_AllPeopleThisCell > 0
		}
	}


#pragma omp parallel for schedule(static,1) default(none) \
		shared(t, run, P, StateT, Hosts, TimeStepNow)
	for (int j = 0; j < P.NumThreads; j++)
	{
		for (int k = 0; k < P.NumThreads; k++)
		{
			for (int i = 0; i < StateT[k].n_queue[j]; i++)
			{
				int infector			= StateT[k].inf_queue[j][i].infector;
				int infectee			= StateT[k].inf_queue[j][i].infectee;
				short int infect_type	= StateT[k].inf_queue[j][i].infect_type;
				Hosts[infectee].infector = infector;
				Hosts[infectee].infect_type = infect_type;
				if (infect_type == -1) //// i.e. if host doesn't have an infector
					DoFalseCase(infectee, t, TimeStepNow, j);
				else
					DoInfect(infectee, t, j, run);
			}
			StateT[k].n_queue[j] = 0;
		}
	}
}

void IncubRecoverySweep(double t)
{
	double ht;
	unsigned short int TimeStepNow; //// this timestep
	TimeStepNow = (unsigned short int) (P.TimeStepsPerDay * t);

	if (P.DoPlaces)
		for (int HolidayNumber = 0; HolidayNumber < P.NumHolidays; HolidayNumber++)
		{
			ht = P.HolidayStartTime[HolidayNumber] + P.HolidaysStartDay_SimTime;
			if ((t + P.ModelTimeStep >= ht) && (t < ht))
			{
//				Files::xfprintf_stderr("Holiday %HolidayNumber t=%lg\n", HolidayNumber, t);
				for (int PlaceType = 0; PlaceType < P.NumPlaceTypes; PlaceType++)
				{
#pragma omp parallel for schedule(static,1) default(none) shared(P, Places, Hosts, HolidayNumber, PlaceType, ht)
					for (int ThreadNum = 0; ThreadNum < P.NumThreads; ThreadNum++)
						for (int PlaceNumber = ThreadNum; PlaceNumber < P.Nplace[PlaceType]; PlaceNumber += P.NumThreads)
						{
							if ((P.HolidayEffect[PlaceType] < 1) /*Proportion of Places of this type closed < 1*/ && ((P.HolidayEffect[PlaceType] == 0) || (ranf_mt(ThreadNum) >= P.HolidayEffect[PlaceType])))
							{
								int HolidayStart	= (int)(ht * P.TimeStepsPerDay);
								int HolidayEnd		= (int)((ht + P.HolidayDuration[HolidayNumber]) * P.TimeStepsPerDay);
								if (Places[PlaceType][PlaceNumber].close_start_time > HolidayStart)  Places[PlaceType][PlaceNumber].close_start_time	= (unsigned short) HolidayStart;
								if (Places[PlaceType][PlaceNumber].close_end_time	< HolidayEnd)	 Places[PlaceType][PlaceNumber].close_end_time		= (unsigned short) HolidayEnd;

								for (int PlaceMember = 0; PlaceMember < Places[PlaceType][PlaceNumber].n; PlaceMember++)
								{
									if (Hosts[Places[PlaceType][PlaceNumber].members[PlaceMember]].absent_start_time	> HolidayStart	) Hosts[Places[PlaceType][PlaceNumber].members[PlaceMember]].absent_start_time	= (unsigned short) HolidayStart;
									if (Hosts[Places[PlaceType][PlaceNumber].members[PlaceMember]].absent_stop_time		< HolidayEnd	) Hosts[Places[PlaceType][PlaceNumber].members[PlaceMember]].absent_stop_time	= (unsigned short) HolidayEnd;
								}
							}
						}
				}
			}
		}

#pragma omp parallel for schedule(static,1) default(none) shared(t, P, CellLookup, Hosts, AdUnits, Mcells, StateT, TimeStepNow)
	for (int ThreadNum = 0; ThreadNum < P.NumThreads; ThreadNum++)	//// loop over threads
		for (int CellIndex = ThreadNum; CellIndex < P.NumPopulatedCells; CellIndex += P.NumThreads)	//// loop/step over populated cells
		{
			Cell* ThisCell = CellLookup[CellIndex]; //// find (pointer-to) ThisCell.
			for (int LatentPerson = ((int)ThisCell->L - 1); LatentPerson >= 0; LatentPerson--) //// loop backwards over latently infected people, hence it starts from L - 1 and goes to zero. Runs backwards because of pointer swapping?
				if (TimeStepNow == Hosts[ThisCell->latent[LatentPerson]].latent_time) //// if now after time at which person became infectious (latent_time a slight misnomer).
					DoIncub(ThisCell->latent[LatentPerson], TimeStepNow, ThreadNum); //// move infected person from latently infected (L) to infectious (I), but not symptomatic

			for (int InfeciousPersonIndexWithinCell = ThisCell->I - 1; InfeciousPersonIndexWithinCell >= 0; InfeciousPersonIndexWithinCell--) ///// loop backwards over Infectious people. Runs backwards because of pointer swapping?
			{
				int InfectiousPersonIndex = ThisCell->infected[InfeciousPersonIndexWithinCell];	//// person index
				Person* InfectiousPerson = Hosts + InfectiousPersonIndex;	//// person

				unsigned short int CaseTime; //// time at which person becomes case (i.e. moves from infectious and asymptomatic to infectious and symptomatic).
				CaseTime = InfectiousPerson->latent_time + ((int)(P.LatentToSymptDelay / P.ModelTimeStep)); //// time that person si/ci becomes case (symptomatic)...
				if ((P.DoSymptoms) && (TimeStepNow == CaseTime)) //// ... if now is that time...
					DoCase(InfectiousPersonIndex, t, TimeStepNow, ThreadNum);		  //// ... change infectious (but asymptomatic) person to infectious and symptomatic. If doing severity, this contains DoMild and DoILI.

				if (P.DoSeverity) // Don't be tempted to if/else or switch following code as there are edge cases where e.g. SARI_time = recovery_or_death_time, and so they're not mutually exclusive.
				{
					if (TimeStepNow == InfectiousPerson->SARI_time)		DoSARI(InfectiousPersonIndex, ThreadNum);	//// see if you can dispense with inequalities by initializing SARI_time, Critical_time etc. to USHRT_MAX
					if (TimeStepNow == InfectiousPerson->Critical_time)	DoCritical(InfectiousPersonIndex, ThreadNum);
					if (TimeStepNow == InfectiousPerson->Stepdown_time)	DoRecoveringFromCritical(InfectiousPersonIndex, ThreadNum);
					if (TimeStepNow == InfectiousPerson->recovery_or_death_time)
					{
						if (InfectiousPerson->to_die)
							DoDeath_FromCriticalorSARIorILI	(InfectiousPersonIndex, ThreadNum);
						else
							DoRecover_FromSeverity(InfectiousPersonIndex, ThreadNum);
					}
				}

				//Adding code to assign recovery or death when leaving the infectious class: ggilani - 22/10/14
				if (TimeStepNow == InfectiousPerson->recovery_or_death_time)
				{
					if (!InfectiousPerson->to_die) //// if person si recovers and this timestep is after they've recovered
					{
						DoRecover(InfectiousPersonIndex, ThreadNum);
						//StateT[ThreadNum].inf_queue[0][StateT[ThreadNum].n_queue[0]++] = ci; //// add them to end of 0th thread of inf queue. 
					}
					else /// if they die and this timestep is after they've died.
					{
						if (HOST_TREATED(InfectiousPersonIndex) && (ranf_mt(ThreadNum) < P.TreatDeathDrop))
							DoRecover(InfectiousPersonIndex, ThreadNum);
						else
							DoDeath(InfectiousPersonIndex, ThreadNum);
					}

					//once host recovers, will no longer make contacts for contact tracing - if we are doing contact tracing and case was infectious when contact tracing was active, increment state vector
					if ((P.DoDigitalContactTracing) && (Hosts[InfectiousPersonIndex].latent_time>= AdUnits[Mcells[Hosts[InfectiousPersonIndex].mcell].adunit].DigitalContactTracingTimeStart) && (Hosts[InfectiousPersonIndex].recovery_or_death_time < AdUnits[Mcells[Hosts[InfectiousPersonIndex].mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[InfectiousPersonIndex].digitalContactTracingUser == 1) && (P.OutputDigitalContactDist))
					{
						if (Hosts[InfectiousPersonIndex].ncontacts > MAX_CONTACTS) Hosts[InfectiousPersonIndex].ncontacts = MAX_CONTACTS;
						//increment bin in State corresponding to this number of contacts
						StateT[ThreadNum].contact_dist[Hosts[InfectiousPersonIndex].ncontacts]++;
					}
				}
			}
		}
}

void DigitalContactTracingSweep(double t)
{
	/**
	 * Function: DigitalContactTracingSweep
	 *
	 * Purpose: to update and count the number of people in each admin unit who are being digitally contact traced each day and remove those who no longer need to be traced
	 * @param t is a double representing the actual simulation time (not the integer timestep)
	 * @return void
	 *
	 * Author: ggilani, 10/03/20 - updated 24/03/20, 14/04/2020
	 */
	unsigned short int TimeStepNow;

	//find current time step
	TimeStepNow = (unsigned short int) (P.TimeStepsPerDay * t);

	FILE* stderr_shared = stderr;
#pragma omp parallel for schedule(static,1) default(none) \
		shared(t, P, AdUnits, StateT, Hosts, TimeStepNow, stderr_shared)
	for (int tn = 0; tn < P.NumThreads; tn++)
	{
		for (int i = tn; i < P.NumAdunits; i += P.NumThreads)
		{
			if (t >= AdUnits[i].DigitalContactTracingTimeStart)
			{
				for (int j = 0; j < P.NumThreads; j++)
				{
					for (int k = 0; k < StateT[j].ndct_queue[i];)
					{
						//start by finding theoretical start and end isolation times for each contact;
						//these are calculated here for each time step instead of InfectSweep when contact event is added as trigger times will be updated for asymptomatic cases detected by testing.
						int infector = StateT[j].dct_queue[i][k].index;
						int contact = StateT[j].dct_queue[i][k].contact;
						unsigned short int contact_time = StateT[j].dct_queue[i][k].contact_time;

						unsigned short int dct_start_time, dct_end_time;
						//this condition is only ever met when a symptomatic case is detected in DoDetectedCase and is not already an index case. If they have already
						//been made an index case due to testing, then this won't occur again for them.
						if (infector==-1)
						{
							//i.e. this is an index case that has been detected by becoming symptomatic and added to the digital contact tracing queue
							dct_start_time = Hosts[contact].dct_trigger_time; //trigger time for these cases is set in DoIncub and already accounts for delay between onset and isolation
							dct_end_time = dct_start_time + (unsigned short int)(P.LengthDigitalContactIsolation * P.TimeStepsPerDay);

						}
						else //We are looking at actual contact events between infectious hosts and their contacts.
						{
							//trigger times are either set in DoDetectedCase or in the loop below (for asymptomatic and presymptomatic cases that are picked up via testing
							//If the contact's index case has a trigger time that means that they have been detected, and we can calculate start and end isolation times for the contact.
							if (Hosts[infector].dct_trigger_time < (USHRT_MAX - 1))
							{
								if (contact_time > Hosts[infector].dct_trigger_time)
								{
									//if the contact time was made after host detected, we should use the later time
									dct_start_time = contact_time + (unsigned short int) (P.DigitalContactTracingDelay * P.TimeStepsPerDay);
								}
								else
								{
									//if the contact time was made before or at the same time as detection, use the trigger time instead
									dct_start_time = Hosts[infector].dct_trigger_time + (unsigned short int) (P.DigitalContactTracingDelay * P.TimeStepsPerDay);
								}
								dct_end_time = dct_start_time + (unsigned short int)(P.LengthDigitalContactIsolation * P.TimeStepsPerDay);
							}
							else
							{
								dct_start_time = USHRT_MAX - 1; //for contacts of asymptomatic or presymptomatic cases - they won't get added as their index case won't know that they are infected (unless explicitly tested)
								//but we keep them in the queue in case their index case is detected as the contact of someone else and gets their trigger time set
								//set dct_end_time to recovery time of infector, in order to remove from queue if their infector isn't detected before they recover.
								dct_end_time = Hosts[infector].recovery_or_death_time;
							}
						}

						//if we've reached the start time for isolation
						if (dct_start_time == TimeStepNow)
						{
							//if the host has been detected due to being symptomatic, they are now an index case - set this variable now. For index cases detected by testing, this will be set on testing
							if ((infector==-1) && (Hosts[contact].index_case_dct == 0)) //don't really need the second condition as the first should only be true when the second isn't (due to how this contact is logged in DoDetectedCase)
							{
								Hosts[contact].index_case_dct = 1; //assign them as an index case
							}

							//if contact is not being traced at all
							if (Hosts[contact].digitalContactTraced == 0)
							{
								//move into the contact tracing list for that admin unit, set start and end times, update flag and remove from queue
								if (AdUnits[i].ndct < AdUnits[i].n) //AdUnits[i].n is length of queue
								{
									Hosts[contact].dct_start_time = dct_start_time;
									Hosts[contact].dct_end_time = dct_end_time;
									Hosts[contact].digitalContactTraced = 1;
									// At this point, we do testing on index cases who have been picked up on symptoms alone, in order to figure out whether and when
									// to remove their contacts (if P.RemoveContactsOfNegativeIndexCase). It's much harder to do it in the next loop as we don't have all
									// the information about the contact event there and would need to loop over all contacts again to look for their index case
									// This would cause race conditions due to having a loop over adunits within threaded loop over admin units
									// Only set test times if P.DoDCTTest. If P.DoDCTTest==0, but we are finding contacts of contacts, we check to see if contacts should become index cases every day they are in isolation.
									if (P.DoDCTTest)
									{
										if (Hosts[contact].index_case_dct == 1)
										{
											//set testing time (which has a different delay to contact testing delay), but no need to set index_case link
											Hosts[contact].dct_test_time = dct_start_time + (unsigned short int)(P.DelayToTestIndexCase * P.TimeStepsPerDay);
											//if host is infectious at test time
											if ((Hosts[contact].dct_test_time >= Hosts[contact].latent_time) && (Hosts[contact].dct_test_time < Hosts[contact].recovery_or_death_time))
											{
												//if false negative, remove from queue by setting the end time to the test time
												if ((P.SensitivityDCT == 0) || ((P.SensitivityDCT < 1) && (ranf_mt(tn) >= P.SensitivityDCT)))
												{
													Hosts[contact].dct_end_time = Hosts[contact].dct_test_time;
													//set index_dct_flag to 2 to indicate that contacts should be removed, if we are removing based on negative test result of index case
													if (P.RemoveContactsOfNegativeIndexCase) Hosts[contact].index_case_dct = 2;
												}
											}
											//if host is non-infectious)
											else
											{
												//if true negative, remove from list
												if ((P.SpecificityDCT == 1) || ((P.SpecificityDCT > 0) && (ranf_mt(tn) < P.SpecificityDCT)))
												{
													//again mark them to be removed from list at test time rather than end_time, and change index_case_dct flag
													Hosts[contact].dct_end_time = Hosts[contact].dct_test_time;
													if (P.RemoveContactsOfNegativeIndexCase) Hosts[contact].index_case_dct = 2;
												}
											}
										}
										else if (Hosts[contact].index_case_dct == 0)
										{
											//if their infector is set to be removed from the list at test time, and their contacts will also be removed at this stage
											if ((Hosts[infector].index_case_dct == 2) && (P.RemoveContactsOfNegativeIndexCase))
											{
												//set end time to match end time of infector
												Hosts[contact].dct_end_time = Hosts[infector].dct_end_time;
											}
											else
											{
												//set testing time
												Hosts[contact].dct_test_time = dct_start_time + (unsigned short int)(P.DelayToTestDCTContacts * P.TimeStepsPerDay);
											}
										}
									}

									//actually put them in the queue
									AdUnits[i].dct[AdUnits[i].ndct] = contact;
									//update number of people in queue
									AdUnits[i].ndct++;
									//increment state variables
									StateT[tn].cumDCT_adunit[i]++;
									StateT[tn].cumDCT++;

									//now remove this case from the queues
									StateT[j].dct_queue[i][k] = StateT[j].dct_queue[i][StateT[j].ndct_queue[i] - 1];
									StateT[j].dct_queue[i][StateT[j].ndct_queue[i] - 1] = { contact,infector,contact_time };
									StateT[j].ndct_queue[i]--;
								}
								else
								{
									Files::xfprintf(stderr_shared, "No more space in queue! AdUnit: %i, ndct=%i, max queue length: %i\n", i, AdUnits[i].ndct, AdUnits[i].n);
									Files::xfprintf(stderr_shared, "Error!\n");
									k++;
								}
							}
							//else if contact is already being contact traced
							else if (Hosts[contact].digitalContactTraced == 1)
							{
								if (P.DoDCTTest)
								{
									//if case has been detected due to being symptomatic, then we will update their testing time if they would be tested earlier based on being an index case as opposed to being a contact of another case
									//If they are already being contact traced and testing is on, they should have been set a test_time
									if ((Hosts[contact].index_case_dct == 1) && (Hosts[contact].dct_test_time > (dct_start_time + (unsigned short int)(P.DelayToTestIndexCase * P.TimeStepsPerDay))))
									{
										Hosts[contact].dct_test_time = dct_start_time + (unsigned short int)(P.DelayToTestIndexCase * P.TimeStepsPerDay);
										//update end time (which is always at least equal to, but may be later that the current one)
										Hosts[contact].dct_end_time = dct_end_time;
										//check to see if test will be negative, if so, tag them for early removal and update index_dct_flag
										if ((Hosts[contact].dct_test_time >= Hosts[contact].latent_time) && (Hosts[contact].dct_test_time < Hosts[contact].recovery_or_death_time))
										{
											//if false negative, remove from
											if ((P.SensitivityDCT == 0) || ((P.SensitivityDCT < 1) && (ranf_mt(tn) >= P.SensitivityDCT)))
											{
												Hosts[contact].dct_end_time = Hosts[contact].dct_test_time;
												//set index_dct_flag to 2 to indicate that contacts should be removed
												if (P.RemoveContactsOfNegativeIndexCase) Hosts[contact].index_case_dct = 2;
											}
										}
										//if host is non-infectious
										else
										{
											//if true negative, remove from list
											if ((P.SpecificityDCT == 1) || ((P.SpecificityDCT > 0) && (ranf_mt(tn) < P.SpecificityDCT)))
											{
												//again mark them to be removed from list at test time rather than end_time, and change index_case_dct flag
												Hosts[contact].dct_end_time = Hosts[contact].dct_test_time;
												if (P.RemoveContactsOfNegativeIndexCase) Hosts[contact].index_case_dct = 2;
											}
										}
									}
									else
									{
										//we don't want to remove this contact if they are also linked to another case - their testing time shouldn't change.
										//but we'll only extend their end time if they wouldn't potentially be removed by having a negative contact
										if ((!P.RemoveContactsOfNegativeIndexCase) || ((P.RemoveContactsOfNegativeIndexCase) && (Hosts[infector].index_case_dct == 1)))
										{
											//extend end time
											Hosts[contact].dct_end_time = dct_end_time;
										}
										//otherwise if contact would have been removed if they didn't have another contact, we keep their original end time
									}

								}
								else
								{
									//just extend the isolation end time, but we're not going to update testing time or as we still want the testing time to be dependent on the earlier contact.
									Hosts[contact].dct_end_time = dct_end_time; //we could choose to not extend the time for cases who are index cases. If they are tested and are negative, they'd be removed earlier anyway. If positive, they will stay isolated for a bit longer

								}
								//now remove this case from the queue
								StateT[j].dct_queue[i][k] = StateT[j].dct_queue[i][StateT[j].ndct_queue[i] - 1];
								StateT[j].dct_queue[i][StateT[j].ndct_queue[i] - 1] = { contact,infector,contact_time };
								StateT[j].ndct_queue[i]--;
							}
						}
						//if contact of an asymptomatic host has passed the recovery time of their asymptomatic index, they would no longer be identified by testing of their index case - remove from the queue so they don't stay here forever
						else if ((dct_start_time == (USHRT_MAX - 1)) && (dct_end_time == TimeStepNow))
						{
							//now remove this case from the queue
							StateT[j].dct_queue[i][k] = StateT[j].dct_queue[i][StateT[j].ndct_queue[i] - 1];
							StateT[j].dct_queue[i][StateT[j].ndct_queue[i] - 1] = { contact,infector,contact_time };
							StateT[j].ndct_queue[i]--;
						}
						else
						{
							k++;
						}

					}
				}
			}
		}
	}

#pragma omp parallel for schedule(static,1) default(none) \
		shared(t, P, AdUnits, Hosts, TimeStepNow)
	for (int tn = 0; tn < P.NumThreads; tn++)
	{
		for (int i = tn; i < P.NumAdunits; i += P.NumThreads)
		{
			if (t >= AdUnits[i].DigitalContactTracingTimeStart)
			{
				for (int j = 0; j < AdUnits[i].ndct;)
				{
					int contact = AdUnits[i].dct[j];

					//first do testing of index cases and their contacts
					if (P.DoDCTTest)
					{
						if ((Hosts[contact].dct_test_time == TimeStepNow) && (Hosts[contact].index_case_dct == 0))
						{
							//if host is positive
							if (Hosts[contact].is_infectious_asymptomatic_not_case() ||
								Hosts[contact].is_case() ||
								Hosts[contact].is_infectious_almost_symptomatic())
							{
								//if the test is a false negative
								if ((P.SensitivityDCT == 0) || ((P.SensitivityDCT < 1) && (ranf_mt(tn) >= P.SensitivityDCT)))
								{
									Hosts[contact].dct_end_time = TimeStepNow;
								}
								//else if a true positive
								else if (P.FindContactsOfDCTContacts)
								{
									//set them to be an index case
									Hosts[contact].index_case_dct = 1;
									//set trigger time to pick up their contacts in the next time step
									Hosts[contact].dct_trigger_time = TimeStepNow + 1; //added the +1 here so that if there are no delays, the contacts will still get picked up correctly
									//if they are an infectious, asymptomatic non-case, call DoDetectedCase in order to trigger HQ and PC too.
									if (Hosts[contact].is_infectious_asymptomatic_not_case())
									{
										DoDetectedCase(contact, t, TimeStepNow, tn);
										Hosts[contact].detected = 1; Hosts[contact].detected_time = TimeStepNow;
									}
								}
							}
							//or if host is negative
							else
							{
								//and is a true negative
								if ((P.SpecificityDCT == 1) || ((P.SpecificityDCT > 0) && (ranf_mt(tn) < P.SpecificityDCT)))
								{
									Hosts[contact].dct_end_time = TimeStepNow;
								}
								//can't track contacts of false positives as they don't make any contacts in InfectSweep
							}
						}
					}
					else if (P.FindContactsOfDCTContacts)
					{
						//check every day to see if contacts become index cases - but they have to be infectious. Otherwise we could set the trigger time and cause their contacts to be traced when they are not being traced themselves.
						if ((Hosts[contact].index_case_dct == 0) && (
							Hosts[contact].is_infectious_almost_symptomatic() ||
							Hosts[contact].is_case() ||
							Hosts[contact].is_infectious_almost_symptomatic()))
							//if ((Hosts[contact].dct_test_time == TimeStepNow) && (Hosts[contact].index_case_dct == 0) && ((abs(Hosts[contact].inf) == 2) || (Hosts[contact].inf == -1)))
						{
							//set them to be an index case
							Hosts[contact].index_case_dct = 1;
							//set trigger time to pick up their contacts in the next time step
							Hosts[contact].dct_trigger_time = TimeStepNow + 1; //added the +1 here so that if there are no delays, the contacts will still get picked up correctly
							//if they are asymptomatic, i.e. specifically if they have inf flag 2, call DoDetectedCase in order to trigger HQ and PC too.
							if (Hosts[contact].is_infectious_asymptomatic_not_case())
							{
								DoDetectedCase(contact, t, TimeStepNow, tn);
								Hosts[contact].detected = 1; Hosts[contact].detected_time = TimeStepNow;
							}
						}
					}

					//now remove hosts who have reached the end of their isolation time
					if (Hosts[contact].dct_end_time == TimeStepNow)
					{
						//stop contact tracing this host
						Hosts[contact].digitalContactTraced = 0;
						//remove index_case_dct flag to 0;
						if (Hosts[contact].index_case_dct)
						{
							Hosts[contact].index_case_dct = 0;
							//Hosts[contact].dct_trigger_time = USHRT_MAX - 1;
						}

						//remove from list
						//k = contact;
						AdUnits[i].dct[j] = AdUnits[i].dct[AdUnits[i].ndct - 1];
						AdUnits[i].dct[AdUnits[i].ndct - 1] = contact;
						AdUnits[i].ndct--;
					}
					else
					{
						j++;
					}
				}
			}
		}
	}
}

int TreatSweep(double t)
{
	///// function loops over microcells to decide which cells are treated (either with treatment, vaccine, social distancing, movement restrictions etc.)

	int TreatFlag = 0, TreatFlag1 = 0; //// Function returns TreatFlag. If TreatFlag == 0, function no longer called. Anytime any treatment used, TreatFlag set to 1. 
	int nckwp;

	//// time steps
	unsigned short int TimeStepNow = (unsigned short int) (P.TimeStepsPerDay * t);	////  time-step now.
	unsigned short int t_TreatStart;												////  time-step treatment begin
	unsigned short int t_TreatEnd;													////  time-step treatment finish
	unsigned short int t_VacStart;													////  time-step vaccination begin
	unsigned short int t_PlaceClosure_End;											////  time-step place closure finish
	unsigned short int t_MoveRestrict_Start;										////  time-step movement restriction begin
	unsigned short int t_MoveRestrict_End;											////  time-step movement restriction finish
	unsigned short int t_SocDist_End;												////  time-step social distancing finish
	unsigned short int t_KeyWorkerPlaceClosure_End;									////  time-step key worker place closure finish
	double radius;

	int global_trig = 0;
	if (P.DoGlobalTriggers)
	{
		if (P.DoPerCapitaTriggers)	global_trig = (int)floor(((double)State.trigDetectedCases) * P.GlobalIncThreshPop / ((double)P.PopSize));
		else						global_trig = State.trigDetectedCases;
	}

	///// block loops over places (or place groups if P.DoPlaceGroupTreat == 1) and determines whom to prophylactically treat
	if ((P.DoPlaces) && (t >= P.TreatTimeStart) && (t < P.TreatTimeStart + P.TreatPlaceGeogDuration) && (State.cumT < P.TreatMaxCourses))
	{
		t_TreatEnd = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatDelayMean + P.TreatProphCourseLength));

#pragma omp parallel for private(TreatFlag) reduction(+:TreatFlag1) schedule(static,1) default(none) \
			shared(P, StateT, Places, Hosts, TimeStepNow, t_TreatEnd)
		for (int Thread = 0; Thread < P.NumThreads; Thread++)
			for (int PlaceType = 0; PlaceType < P.NumPlaceTypes; PlaceType++)
			{
				for (int PlaceNumQueueIndex = 0; PlaceNumQueueIndex < StateT[Thread].np_queue[PlaceType]; PlaceNumQueueIndex++) //// loop over all places IN QUEUE, not all a places
				{
					int PlaceNum = StateT[Thread].p_queue[PlaceType][PlaceNumQueueIndex]; //// note PlaceNum is index of place, not index of place in place queue.
					if (P.DoPlaceGroupTreat)
					{
						int PlaceGroupIndex = StateT[Thread].pg_queue[PlaceType][PlaceNumQueueIndex];
						TreatFlag = PlaceGroupIndex; //// keep this as a flag 
						for (int PG_member = ((int)Places[PlaceType][PlaceNum].group_start[PlaceGroupIndex]); PG_member < ((int)(Places[PlaceType][PlaceNum].group_start[PlaceGroupIndex] + Places[PlaceType][PlaceNum].group_size[PlaceGroupIndex])); PG_member++) // loop over people in place group.
						{
							/*if((Places[PlaceType][PlaceNum].members[PG_member]<0)||(Places[PlaceType][PlaceNum].members[PG_member]>P.PopSize-1))
								Files::xfprintf_stderr("\n*** npq=%i gn=%i h=%i PG_member=%i PlaceType=%i PlaceNum=%i PlaceGroupIndex=%i s=%i n=%i ***\n",
									StateT[Thread].np_queue[PlaceType],
									Places[PlaceType][PlaceNum].n,
									Places[PlaceType][PlaceNum].members[PG_member],
									PG_member,PlaceType,PlaceNum,PlaceGroupIndex,
									(int) Places[PlaceType][PlaceNum].group_start[PlaceGroupIndex],
									(int) Places[PlaceType][PlaceNum].group_size[PlaceGroupIndex]);
							else
							*/
							if ((!HOST_TO_BE_TREATED(Places[PlaceType][PlaceNum].members[PG_member])) && ((P.TreatPlaceTotalProp[PlaceType] == 1) || (ranf_mt(Thread) < P.TreatPlaceTotalProp[PlaceType])))
								DoProph(Places[PlaceType][PlaceNum].members[PG_member], TimeStepNow, Thread);
						}
					}
					else
					{
						if ((Places[PlaceType][PlaceNum].treat) && (!PLACE_TREATED(PlaceType, PlaceNum)))
						{
							TreatFlag1 = 1;
							Places[PlaceType][PlaceNum].treat_end_time = t_TreatEnd;
							for (int PG_member = 0; PG_member < Places[PlaceType][PlaceNum].n; PG_member++)
								if (!HOST_TO_BE_TREATED(Places[PlaceType][PlaceNum].members[PG_member]))
								{
									if ((P.TreatPlaceTotalProp[PlaceType] == 1) || (ranf_mt(Thread) < P.TreatPlaceTotalProp[PlaceType]))
										DoProph(Places[PlaceType][PlaceNum].members[PG_member], TimeStepNow, Thread);
								}
						}
						Places[PlaceType][PlaceNum].treat = 0;
					}
				}
				StateT[Thread].np_queue[PlaceType] = 0;
			}
	}

	///// block vaccinates everyone in mass vaccination queue. Don't know why loop is done twice (although State.mvacc_cum is reset at end so will relate to that)
	if ((P.DoMassVacc) && (t >= P.VaccTimeStart))
		for (int j = 0; j < 2; j++)
		{
			int m = (int)P.VaccMaxCourses;
			if (m > State.n_mvacc) m = State.n_mvacc;
#pragma omp parallel for schedule(static,1000) default(none) \
				shared(State, m, TimeStepNow)
			for (int i = State.mvacc_cum; i < m; i++)
				DoVacc(State.mvacc_queue[i], TimeStepNow);
			State.mvacc_cum = m;
		}

	//// Main block of function - assigns various start and end times for various "treatments" (e.g. vaccination, treatment, place closure, social distancing)
	//// Then Loops over microcells seeing which if any of these start and end times are relevant for each microcell. If they are, various "TreatStat" treatment status flags are changed.
	//// When flagging microcells as "ToBeTreated", for some interventions (treatment, vaccination),
	//// surrounding microcells also flagged for as "ToBeTreated".
	//// Function returns TreatFlag. By default, TreatFlag is set to 0. Anytime any of these steps are performed, TreatFlag is set to 1.
	//// Therefore if nothing happens over all microcells then this function TreatSweep is no longer called. 
	if ((t >= P.TreatTimeStart) || (t >= P.VaccTimeStartGeo) || (t >= P.PlaceCloseTimeStart) || (t >= P.MoveRestrTimeStart) || (t >= P.SocDistTimeStart) || (t >= P.KeyWorkerProphTimeStart)) //changed this to start time geo
	{
		t_TreatStart				= (unsigned short int) (P.TimeStepsPerDay		* (t + P.TreatDelayMean));
		t_TreatEnd					= (unsigned short int) (P.TimeStepsPerDay		* (t + P.TreatProphCourseLength) - 1);
		t_VacStart					= (unsigned short int) (P.TimeStepsPerDay		* (t + P.VaccDelayMean));
		t_PlaceClosure_End			= (unsigned short int) ceil(P.TimeStepsPerDay	* (t + P.PlaceCloseDelayMean + P.PlaceCloseDuration));
		t_MoveRestrict_Start		= (unsigned short int) floor(P.TimeStepsPerDay	* (t + P.MoveDelayMean));
		t_MoveRestrict_End			= (unsigned short int) ceil(P.TimeStepsPerDay	* (t + P.MoveRestrDuration));
		t_SocDist_End				= (unsigned short int) ceil(P.TimeStepsPerDay	* (t + P.SocDistDurationCurrent));
		t_KeyWorkerPlaceClosure_End = (unsigned short int) ceil(P.TimeStepsPerDay	* (t + P.KeyWorkerProphRenewalDuration));
		nckwp = (int)ceil(P.KeyWorkerProphDuration / P.TreatProphCourseLength);

#pragma omp parallel for private(radius) reduction(+:TreatFlag) schedule(static,1) default(none) \
			shared(t, P, Hosts, Mcells, McellLookup, AdUnits, State, global_trig, TimeStepNow, t_TreatEnd, t_TreatStart, t_VacStart, t_PlaceClosure_End, t_MoveRestrict_End, t_MoveRestrict_Start, t_SocDist_End, t_KeyWorkerPlaceClosure_End, nckwp)
		for (int ThreadNum = 0; ThreadNum < P.NumThreads; ThreadNum++)
			for (int PopulatedMicroCellNum = ThreadNum; PopulatedMicroCellNum < P.NumPopulatedMicrocells; PopulatedMicroCellNum += P.NumThreads) //// loop over populated microcells
			{
				int mcellnum	= (int)(McellLookup[PopulatedMicroCellNum] - Mcells); //// microcell number
				int adi			= (P.DoAdUnits) ? Mcells[mcellnum].adunit : -1;
				int AdminUnit	= (P.DoAdUnits) ? AdUnits[adi].id : 0;

				//// Code block goes through various types of treatments/interventions (vaccination/movement restrictions etc.),
				//// assesses whether various triggers (counts) are over a certain threshold, (specified in ReadParams)
				//// and then implements those treatments by setting various flags (i.e. .treat/ .vacc etc.) by microcell.
				//// Further, this block assigns all microcells that are within this admin unit (and around this microcell) to be treated, using the flags set to avoid duplication.

				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
				//// **** //// **** //// **** //// **** TREATMENT
				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

				// if Microcell already treated but TimeStepNow after treatment end time.
				if ((Mcells[mcellnum].treat == TreatStat::Treated		) && (TimeStepNow >= Mcells[mcellnum].treat_end_time))
				{
					TreatFlag						= 1;
					Mcells[mcellnum].treat			= TreatStat::Untreated;
				}
				// if Microcell flagged as to be treated TimeStepNow after treatment start time.
				if ((Mcells[mcellnum].treat == TreatStat::ToBeTreated	) && (TimeStepNow >= Mcells[mcellnum].treat_start_time))
				{
					TreatFlag						= 1;
					Mcells[mcellnum].treat			= TreatStat::Treated;
					Mcells[mcellnum].treat_trig		= 0;					// reset trigger
					Mcells[mcellnum].treat_end_time	= t_TreatEnd;

					for (int MicrocellMember = 0; MicrocellMember < Mcells[mcellnum].n; MicrocellMember++) // loop over microcell members
					{
						int PersonInMicrocell = Mcells[mcellnum].members[MicrocellMember];
						if ((!HOST_TO_BE_TREATED(PersonInMicrocell)) && ((P.TreatPropRadial == 1) || (ranf_mt(ThreadNum) < P.TreatPropRadial)))
							DoProphNoDelay(PersonInMicrocell, TimeStepNow, ThreadNum, 1);
					}
				}

				// Has trigger threshold been reached (treatment)?
				bool TrigThreshReached_Treatment = false; 
				if (P.DoGlobalTriggers)
					TrigThreshReached_Treatment = (global_trig >= P.TreatCellIncThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : (int)P.TreatCellIncThresh;
					TrigThreshReached_Treatment = (State.trigDC_adunit[adi] > trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : (int)P.TreatCellIncThresh;
					TrigThreshReached_Treatment = (Mcells[mcellnum].treat_trig >= trig_thresh);
				}

				// Block below spirals around microcell, changing treatment status variables of nearby microcells, mostly flagging them for treatment in the lines above
				if ((t >= P.TreatTimeStart) && (Mcells[mcellnum].treat == TreatStat::Untreated) && (TrigThreshReached_Treatment) && (P.TreatRadius2 > 0))
				{
					MicroCellPosition Nearby_mcellposition	= P.get_micro_cell_position_from_cell_index(mcellnum); // initialize Nearby_mcellposition to be position of microcell.
					Direction CurrentDirection				= Direction::Right; // initialize CurrentDirection to go to the right

					int NearbyMCell	= mcellnum; // initialize NearbyMCell to mcellnum
					int PeopleRequiringRadialProphylaxis = 0;
					int i = 0, MicroCellCounter = 0, ColumnCounter = 1;
					bool AskAgainIfStillTreating = false, StillTreating = true;

					if ((!P.TreatByAdminUnit) || (AdminUnit > 0))
					{
						int ad2 = AdminUnit / P.TreatAdminUnitDivisor;
						do
						{
							// depending on characteristics of the nearby Microcells (starting with this microcell), alter their TreatStat variables (start time, TreatStat etc.)
							if (P.is_in_bounds(Nearby_mcellposition))
							{
								bool TreatThisMicroCell = false; 
								if (P.TreatByAdminUnit)
									TreatThisMicroCell = (AdUnits[Mcells[NearbyMCell].adunit].id / P.TreatAdminUnitDivisor == ad2);
								else
									TreatThisMicroCell = ((radius = dist2_mm(Mcells + mcellnum, Mcells + NearbyMCell)) < P.TreatRadius2);

								if (TreatThisMicroCell)
								{
									TreatFlag = 1;
									AskAgainIfStillTreating = 1;
									if ((Mcells[NearbyMCell].n > 0) && (Mcells[NearbyMCell].treat == TreatStat::Untreated))
									{
										Mcells[NearbyMCell].treat_start_time	= t_TreatStart;
										Mcells[NearbyMCell].treat				= TreatStat::ToBeTreated;
										PeopleRequiringRadialProphylaxis		+= Mcells[NearbyMCell].n; // add total number of people in microcell to PeopleRequiringRadialProphylaxis
									}
								}
							}

							// Change microcell position
							Nearby_mcellposition += CurrentDirection; // (depending on CurrentDirection, increment or decrememnt x or y coordinates of microcell position
							MicroCellCounter = (MicroCellCounter + 1) % ColumnCounter;
							if (MicroCellCounter == 0) // i.e. if MicroCellCounter before modulo arithmetic was some multiple of ColumnCounter (i.e. if you've got to the end of the row or column depending of direction). 
							{
								CurrentDirection = rotate_left(CurrentDirection); // up -> left, left -> down, down -> right, right -> up
								i = (i + 1) % 2;
								if (i == 0) ColumnCounter++;
								if (CurrentDirection == Direction::Up)
								{
									StillTreating = AskAgainIfStillTreating;
									AskAgainIfStillTreating = 0; // reset AskAgainIfStillTreating to 0 (will be set to 1 again in block above if Still treating - i.e. radius and doses not exceeded and not all nearby admin units already treated). 
								}
							}
							// Change microcell 
							NearbyMCell = P.get_micro_cell_index_from_position(Nearby_mcellposition);

						} while ((StillTreating) && (PeopleRequiringRadialProphylaxis < P.TreatMaxCoursesPerCase));
					}
				}


				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
				//// **** //// **** //// **** //// **** VACCINATION
				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****


				// if Microcell flagged as to be treated TimeStepNow after vaccination start time.
				//// vaccinates proportion VaccProp of people in microcell (or at least adds them to geovacc_queue).
				if ((Mcells[mcellnum].vacc == TreatStat::ToBeTreated) && (TimeStepNow >= Mcells[mcellnum].vacc_start_time))
				{
					TreatFlag					= 1;
					Mcells[mcellnum].vacc		= TreatStat::Treated;
					Mcells[mcellnum].vacc_trig	= 0;
					//if(State.cumVG+P.NumThreads*Mcells[mcellnum].n<P.VaccMaxCourses) //changed to VG - commented this out for now, we'll add everyone to queues and deal with the number of doses available in the vaccination function
					{
						for (int MicrocellMember = 0; MicrocellMember < Mcells[mcellnum].n; MicrocellMember++)
						{
							int PersonInMicrocell = Mcells[mcellnum].members[MicrocellMember];
							//#pragma omp critical (state_cumV_daily) //added this
							if (((P.VaccProp == 1) || (ranf_mt(ThreadNum) < P.VaccProp)))
								DoVaccNoDelay(PersonInMicrocell,TimeStepNow); //add to the queue
						}
					}
				}

				// Has trigger threshold been reached (vaccination)?
				bool TrigThreshReached_Vacc = false;
				if (P.DoGlobalTriggers)
					TrigThreshReached_Vacc = (global_trig >= P.VaccCellIncThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : (int)P.VaccCellIncThresh;
					TrigThreshReached_Vacc = (State.trigDC_adunit[adi] > trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : (int)P.VaccCellIncThresh;
					TrigThreshReached_Vacc = (Mcells[mcellnum].treat_trig >= trig_thresh);
				}

				// Block below spirals around microcell, changing treatment status variables of nearby microcells, mostly flagging them for vaccination in the lines above
				if ((!P.DoMassVacc) && (P.VaccRadius2 > 0) && (t >= P.VaccTimeStartGeo) && (Mcells[mcellnum].vacc == TreatStat::Untreated) && (TrigThreshReached_Vacc)) //changed from VaccTimeStart to VaccTimeStarGeo
				{
					MicroCellPosition Nearby_mcellposition	= P.get_micro_cell_position_from_cell_index(mcellnum); // initialize Nearby_mcellposition to be position of microcell.
					Direction CurrentDirection				= Direction::Right;// initialize CurrentDirection to go to the right

					int NearbyMCell = mcellnum;  // initialize NearbyMCell to mcellnum
					int i = 0, MicroCellCounter = 0, ColumnCounter = 1;
					bool AskAgainIfStillVaccinating = false, StillVaccinating = true;

					if ((!P.VaccByAdminUnit) || (AdminUnit > 0))
					{
						int ad2 = AdminUnit / P.VaccAdminUnitDivisor;
						do
						{
							// depending on characteristics of the nearby Microcells (starting with this microcell), alter their TreatStat variables (start time, TreatStat etc.)
							if (P.is_in_bounds(Nearby_mcellposition))
							{
								bool VaccinateThisMicroCell = false;
								if (P.VaccByAdminUnit)
								{
									VaccinateThisMicroCell = (AdUnits[Mcells[NearbyMCell].adunit].id / P.VaccAdminUnitDivisor == ad2);
									radius = 1e20;
								}
								else
									VaccinateThisMicroCell = ((radius = dist2_mm(Mcells + mcellnum, Mcells + NearbyMCell)) < P.VaccRadius2);

								if (VaccinateThisMicroCell)
								{
									TreatFlag = 1;
									AskAgainIfStillVaccinating = 1;
									if (radius < P.VaccMinRadius2)
										Mcells[NearbyMCell].vacc = TreatStat::DontTreatAgain;
									else if ((Mcells[NearbyMCell].n > 0) && (Mcells[NearbyMCell].vacc == TreatStat::Untreated))
									{
										Mcells[NearbyMCell].vacc_start_time = t_VacStart;
										Mcells[NearbyMCell].vacc			= TreatStat::ToBeTreated;
									}
								}
							}

							// Change microcell position
							Nearby_mcellposition += CurrentDirection; // (depending on CurrentDirection, increment or decrememnt x or y coordinates of microcell position
							MicroCellCounter = (MicroCellCounter + 1) % ColumnCounter;
							if (MicroCellCounter == 0) // i.e. if MicroCellCounter before modulo arithmetic was some multiple of ColumnCounter (i.e. if you've got to the end of the row or column depending of direction). 
							{
								CurrentDirection = rotate_left(CurrentDirection); // up -> left, left -> down, down -> right, right -> up
								i = (i + 1) % 2;
								if (i == 0) ColumnCounter++;
								if (CurrentDirection == Direction::Up)
								{
									StillVaccinating = AskAgainIfStillVaccinating;
									AskAgainIfStillVaccinating = 0; // reset AskAgainIfStillTreating to 0 (will be set to 1 again in block above if Still treating - i.e. radius and doses not exceeded and not all nearby admin units already treated). 
								}
							}
							// Change microcell 
							NearbyMCell = P.get_micro_cell_index_from_position(Nearby_mcellposition);

						} while (StillVaccinating);
					}
				}

				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
				//// **** //// **** //// **** //// **** PLACE CLOSURE
				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

				// Has off-trigger threshold been reached? I.e. so places can open again?
				///// note that here question is whether trigger lower than stop threshold. A few blocks down meaning changes to almost the opposite: asking whether trigger has exceeded threshold in order to close places for first time.

				bool BelowTrigThreshold_PlaceOpen = false;
				if (P.DoGlobalTriggers)
					BelowTrigThreshold_PlaceOpen = (global_trig < P.PlaceCloseCellIncStopThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.PlaceCloseCellIncStopThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncStopThresh;
					BelowTrigThreshold_PlaceOpen = (State.trigDC_adunit[adi] < trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.PlaceCloseCellIncStopThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncStopThresh;
					BelowTrigThreshold_PlaceOpen = (Mcells[mcellnum].treat_trig < trig_thresh);
				}

				// if Microcell already treated but TimeStepNow after placeclose end time.
				if ((Mcells[mcellnum].placeclose == TreatStat::Treated) && ((BelowTrigThreshold_PlaceOpen) || (TimeStepNow >= Mcells[mcellnum].place_end_time))) //// if place closure has started, the places in this microcell are closed, and either stop threshold has been reached or place_end_time has passed, go through block
				{
					if (P.DoPlaceCloseOnceOnly)
						Mcells[mcellnum].placeclose = TreatStat::DontTreatAgain;
					else
						Mcells[mcellnum].placeclose = TreatStat::Untreated;

					TreatFlag						= 1;
					Mcells[mcellnum].place_end_time = TimeStepNow;
					Mcells[mcellnum].place_trig		= 0;

					if (BelowTrigThreshold_PlaceOpen)
						for (int PlaceType = 0; PlaceType < P.NumPlaceTypes; PlaceType++)
							if (PlaceType != P.HotelPlaceType)
								for (int PlaceNumber = 0; PlaceNumber < Mcells[mcellnum].NumPlacesByType[PlaceType]; PlaceNumber++)
									DoPlaceOpen(PlaceType, Mcells[mcellnum].places[PlaceType][PlaceNumber], TimeStepNow);
				}

				if ((P.DoPlaces) && (t >= P.PlaceCloseTimeStart) && (Mcells[mcellnum].placeclose == TreatStat::Untreated)) //// if doing places, time now is after policy has begun, but place hasn't closed yet.
				{
					///// note that here TrigThreshReached_PlaceClosure bool asks whether trigger has exceeded threshold in order to close places for first time.A few blocks up meaning was almost the opposite: asking whether trigger lower than stop threshold.

					bool TrigThreshReached_PlaceClosure = false;
					if (P.DoGlobalTriggers)
						TrigThreshReached_PlaceClosure = (global_trig >= P.PlaceCloseCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
						TrigThreshReached_PlaceClosure = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
						TrigThreshReached_PlaceClosure = (Mcells[mcellnum].treat_trig >= trig_thresh);
					}

					if (((P.PlaceCloseByAdminUnit) && (AdUnits[Mcells[mcellnum].adunit].place_close_trig < USHRT_MAX - 1)
						&& (((double)AdUnits[Mcells[mcellnum].adunit].place_close_trig) / ((double)AdUnits[Mcells[mcellnum].adunit].NP) > P.PlaceCloseAdunitPropThresh))
						|| ((!P.PlaceCloseByAdminUnit) && (TrigThreshReached_PlaceClosure)))
					{
						//							if(P.PlaceCloseByAdminUnit) AdUnits[Mcells[mcellnum].adunit].place_close_trig=USHRT_MAX-1; // This means schools only close once
						int interventionFlag; //added this as a way to help filter out when interventions start
						interventionFlag = 1;
						if ((P.DoInterventionDelaysByAdUnit) && ((t <= AdUnits[Mcells[mcellnum].adunit].PlaceCloseTimeStart) || (t >= (AdUnits[Mcells[mcellnum].adunit].PlaceCloseTimeStart + AdUnits[Mcells[mcellnum].adunit].PlaceCloseDuration))))
							interventionFlag = 0;

						if ((interventionFlag == 1) && ((!P.PlaceCloseByAdminUnit) || (AdminUnit > 0)))
						{
							if ((Mcells[mcellnum].n > 0) && (Mcells[mcellnum].placeclose == TreatStat::Untreated))
							{
								//if doing intervention delays and durations by admin unit based on global triggers
								if (P.DoInterventionDelaysByAdUnit)
									Mcells[mcellnum].place_end_time = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.PlaceCloseDelayMean + AdUnits[Mcells[mcellnum].adunit].PlaceCloseDuration));
								else
									Mcells[mcellnum].place_end_time = t_PlaceClosure_End;

								Mcells[mcellnum].place_trig = 0;
								Mcells[mcellnum].placeclose = TreatStat::Treated;
								for (int PlaceType = 0; PlaceType < P.NumPlaceTypes; PlaceType++)
									if (PlaceType != P.HotelPlaceType)
										for (int PlaceNumber = 0; PlaceNumber < Mcells[mcellnum].NumPlacesByType[PlaceType]; PlaceNumber++)
											DoPlaceClose(PlaceType, Mcells[mcellnum].places[PlaceType][PlaceNumber], TimeStepNow, ThreadNum, 1);
							}
						}
					}
				}


				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
				//// **** //// **** //// **** //// **** MOVEMENT RESTRICTIONS
				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

				// if Microcell already treated but TimeStepNow after movement restriction end time.
				if ((Mcells[mcellnum].moverest == TreatStat::Treated) && (TimeStepNow >= Mcells[mcellnum].move_end_time))
				{
					TreatFlag = 1;
					if (P.DoMoveRestrOnceOnly)
						Mcells[mcellnum].moverest = TreatStat::DontTreatAgain;
					else
						Mcells[mcellnum].moverest = TreatStat::Untreated;

				}
				// if Microcell flagged as to be treated TimeStepNow after movement restriction start time.
				if ((Mcells[mcellnum].moverest == TreatStat::ToBeTreated) && (TimeStepNow >= Mcells[mcellnum].move_start_time))
				{
					TreatFlag						= 1;
					Mcells[mcellnum].moverest		= TreatStat::Treated;
					Mcells[mcellnum].move_trig		= 0;					// reset trigger
					Mcells[mcellnum].move_end_time	= t_MoveRestrict_End;
				}

				// Has trigger threshold been reached (movement restriction)?
				bool TrigThreshReached_MoveRestrict = false;
				if (P.DoGlobalTriggers)
					TrigThreshReached_MoveRestrict = (global_trig >= P.MoveRestrCellIncThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
					TrigThreshReached_MoveRestrict = (State.trigDC_adunit[adi] > trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
					TrigThreshReached_MoveRestrict = (Mcells[mcellnum].treat_trig >= trig_thresh);
				}

				// Block below spirals around microcell, changing treatment status variables of nearby microcells, mostly flagging them for treatment in the lines above
				if ((t >= P.MoveRestrTimeStart) && (Mcells[mcellnum].moverest == TreatStat::Untreated) && (TrigThreshReached_MoveRestrict))
				{
					MicroCellPosition Nearby_mcellposition	= P.get_micro_cell_position_from_cell_index(mcellnum); // initialize Nearby_mcellposition to be position of microcell.
					Direction CurrentDirection				= Direction::Right; // initialize CurrentDirection to go to the right

					int NearbyMCell = mcellnum; // initialize NearbyMCell to mcellnum
					int i = 0, MicroCellCounter = 0, ColumnCounter = 1;
					bool AskAgainIfStillTreating = false, StillTreating = true;

					if ((!P.MoveRestrByAdminUnit) || (AdminUnit > 0))
					{
						int ad2 = AdminUnit / P.MoveRestrAdminUnitDivisor;
						do
						{
							// depending on characteristics of the nearby Microcells (starting with this microcell), alter their TreatStat variables (start time, TreatStat etc.)
							if (P.is_in_bounds(Nearby_mcellposition))
							{
								bool TreatThisMicroCell = false;
								if (P.MoveRestrByAdminUnit)
									TreatThisMicroCell = (AdUnits[Mcells[NearbyMCell].adunit].id / P.MoveRestrAdminUnitDivisor == ad2);
								else
									TreatThisMicroCell = ((radius = dist2_mm(Mcells + mcellnum, Mcells + NearbyMCell)) < P.MoveRestrRadius2);
								if (TreatThisMicroCell)
								{
									TreatFlag = 1;
									AskAgainIfStillTreating = 1;
									if ((Mcells[NearbyMCell].n > 0) && (Mcells[NearbyMCell].moverest == TreatStat::Untreated))
									{
										Mcells[NearbyMCell].move_start_time = t_MoveRestrict_Start;
										Mcells[NearbyMCell].moverest		= TreatStat::ToBeTreated;
									}
								}
							}

							// Change microcell position
							Nearby_mcellposition += CurrentDirection; // (depending on CurrentDirection, increment or decrememnt x or y coordinates of microcell position
							MicroCellCounter = (MicroCellCounter + 1) % ColumnCounter;
							if (MicroCellCounter == 0) // i.e. if MicroCellCounter before modulo arithmetic was some multiple of ColumnCounter (i.e. if you've got to the end of the row or column depending of direction). 
							{
								CurrentDirection = rotate_left(CurrentDirection); // up -> left, left -> down, down -> right, right -> up
								i = (i + 1) % 2;
								if (i == 0) ColumnCounter++;
								if (CurrentDirection == Direction::Up)
								{
									StillTreating = AskAgainIfStillTreating;
									AskAgainIfStillTreating = 0; // reset AskAgainIfStillTreating to 0 (will be set to 1 again in block above if Still treating - i.e. radius and doses not exceeded and not all nearby admin units already treated). 
								}
							}
							// Change microcell 
							NearbyMCell = P.get_micro_cell_index_from_position(Nearby_mcellposition);

						} while (StillTreating);
					}
				}

				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
				//// **** //// **** //// **** //// **** SOCIAL DISTANCING
				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

				// Has off-trigger threshold been reached? i.e. so social distancing can stop?
				///// note that here question is whether trigger lower than stop threshold. A few blocks down meaning changes to almost the opposite: asking whether trigger has exceeded threshold in order to start social distancing..

				bool BelowTrigThreshold_SocDist = false;
				if (P.DoGlobalTriggers)
					BelowTrigThreshold_SocDist = (global_trig < P.SocDistCellIncStopThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.SocDistCellIncStopThresh)) / P.IncThreshPop)) : P.SocDistCellIncStopThresh;
					BelowTrigThreshold_SocDist = (State.trigDC_adunit[adi] < trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.SocDistCellIncStopThresh)) / P.IncThreshPop)) : P.SocDistCellIncStopThresh;
					BelowTrigThreshold_SocDist = (Mcells[mcellnum].treat_trig < trig_thresh);
				}

				//// if: policy of social distancing has started AND this microcell cell has been labelled to as undergoing social distancing, AND either trigger not reached (note definition of f2 changes in next few lines) or end time has passed.
				if ((t >= P.SocDistTimeStart) && (Mcells[mcellnum].socdist == TreatStat::Treated) && ((BelowTrigThreshold_SocDist) || (TimeStepNow >= Mcells[mcellnum].socdist_end_time)))
				{
					TreatFlag = 1;
					if (P.DoSocDistOnceOnly)
						Mcells[mcellnum].socdist		= TreatStat::DontTreatAgain;
					else
						Mcells[mcellnum].socdist		= TreatStat::Untreated;


					Mcells[mcellnum].socdist_trig		= 0;	//// reset trigger
					Mcells[mcellnum].socdist_end_time	= TimeStepNow; //// record end time.
				}

				///// note that here TrigThreshReached_SocDist bool asks whether trigger has exceeded threshold in order to social distance for first time. A few blocks up meaning was almost the opposite: asking whether trigger lower than stop threshold.
				bool TrigThreshReached_SocDist = false;
				if (P.DoGlobalTriggers)
					TrigThreshReached_SocDist = (global_trig >= P.SocDistCellIncThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
					TrigThreshReached_SocDist = (State.trigDC_adunit[adi] >= trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
					TrigThreshReached_SocDist = (Mcells[mcellnum].treat_trig >= trig_thresh);
				}

				if ((t >= P.SocDistTimeStart) && (Mcells[mcellnum].socdist == TreatStat::Untreated) && (TrigThreshReached_SocDist))
				{
					// some code to try and deal with intervention delays and durations by admin unit based on global triggers
					int interventionFlag; //added this as a way to help filter out when interventions start
					interventionFlag = 1;

					if (P.DoInterventionDelaysByAdUnit)
						if ((t <= AdUnits[Mcells[mcellnum].adunit].SocialDistanceTimeStart) ||
							(t >= (AdUnits[Mcells[mcellnum].adunit].SocialDistanceTimeStart + AdUnits[Mcells[mcellnum].adunit].SocialDistanceDuration))) //// i.e. if outside window of social distancing for this admin unit.
							interventionFlag = 0;

					if (interventionFlag == 1)
						if ((Mcells[mcellnum].n > 0) && (Mcells[mcellnum].socdist == TreatStat::Untreated)) //// if microcell populated and not currently undergoing social distancing
						{
							Mcells[mcellnum].socdist		= TreatStat::Treated; //// update flag to denote that cell is undergoing social distancing
							Mcells[mcellnum].socdist_trig	= 0; /// reset trigger

							//// set (admin-specific) social distancing end time.
							if (P.DoInterventionDelaysByAdUnit)
								Mcells[mcellnum].socdist_end_time = (unsigned short int) ceil(P.TimeStepsPerDay * (t + AdUnits[Mcells[mcellnum].adunit].SocialDistanceDuration));
							else
								Mcells[mcellnum].socdist_end_time = t_SocDist_End;
						}
				}


				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
				//// **** //// **** //// **** //// **** KEY-WORKER PROPHYLAXIS
				//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

				// if Microcell already treated key-workers but TimeStepNow after treatment end time.
				if ((Mcells[mcellnum].keyworkerproph == 2) && (TimeStepNow >= Mcells[mcellnum].keyworkerproph_end_time))
				{
					TreatFlag = 1;
					Mcells[mcellnum].keyworkerproph = P.DoKeyWorkerProphOnceOnly;
				}

				// Has trigger threshold been reached (KeyWorkerProph)?
				bool TrigThreshReached_KeyWorkerProph = false;
				if (P.DoGlobalTriggers)
					TrigThreshReached_KeyWorkerProph = (global_trig >= P.KeyWorkerProphCellIncThresh);
				else if (P.DoAdminTriggers)
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
					TrigThreshReached_KeyWorkerProph = (State.trigDC_adunit[adi] > trig_thresh);
				}
				else
				{
					int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
					TrigThreshReached_KeyWorkerProph = (Mcells[mcellnum].treat_trig >= trig_thresh);
				}

				// Block below spirals around microcell, changing treatment status variables of nearby microcells, mostly flagging them for treatment in the lines above
				if ((P.DoPlaces) && (t >= P.KeyWorkerProphTimeStart) && (Mcells[mcellnum].keyworkerproph == 0) && (TrigThreshReached_KeyWorkerProph))
				{
					MicroCellPosition Nearby_mcellposition	= P.get_micro_cell_position_from_cell_index(mcellnum); // initialize Nearby_mcellposition to be position of microcell.
					Direction CurrentDirection				= Direction::Right; // initialize CurrentDirection to go to the right

					int NearbyMCell = mcellnum; // initialize NearbyMCell to mcellnum
					int i = 0, MicroCellCounter = 0, ColumnCounter = 1;
					bool AskAgainIfStillTreating = false, StillTreating = true;

					do
					{
						// depending on characteristics of the nearby Microcells (starting with this microcell), alter their TreatStat variables (start time, TreatStat etc.)
						if (P.is_in_bounds(Nearby_mcellposition))
							if (dist2_mm(Mcells + mcellnum, Mcells + NearbyMCell) < P.KeyWorkerProphRadius2)
							{
								TreatFlag = 1;
								AskAgainIfStillTreating = 1;

								if ((Mcells[NearbyMCell].n > 0) && (Mcells[NearbyMCell].keyworkerproph == 0))
								{
									Mcells[NearbyMCell].keyworkerproph			= 2;
									Mcells[NearbyMCell].keyworkerproph_trig		= 0;							// reset trigger
									Mcells[NearbyMCell].keyworkerproph_end_time = t_KeyWorkerPlaceClosure_End;

									for (int MicrocellMember = 0; MicrocellMember < Mcells[NearbyMCell].n; MicrocellMember++)
									{
										int PerosonInMicrocell = Mcells[NearbyMCell].members[MicrocellMember];
										if ((Hosts[PerosonInMicrocell].keyworker) && (!HOST_TO_BE_TREATED(PerosonInMicrocell)))
											DoProphNoDelay(PerosonInMicrocell, TimeStepNow, ThreadNum, nckwp);
									}
								}
							}

						// Change microcell position
						Nearby_mcellposition += CurrentDirection; // (depending on CurrentDirection, increment or decrememnt x or y coordinates of microcell position
						MicroCellCounter = (MicroCellCounter + 1) % ColumnCounter;
						if (MicroCellCounter == 0)
						{
							CurrentDirection = rotate_left(CurrentDirection); // up -> left, left -> down, down -> right, right -> up
							i = (i + 1) % 2;
							if (i == 0) ColumnCounter++;
							if (CurrentDirection == Direction::Up)
							{
								StillTreating = AskAgainIfStillTreating;
								AskAgainIfStillTreating = 0; // reset AskAgainIfStillTreating to 0 (will be set to 1 again in block above if Still treating - i.e. radius and doses not exceeded and not all nearby admin units already treated). 
							}
						}
						// Change microcell 
						NearbyMCell = P.get_micro_cell_index_from_position(Nearby_mcellposition);

					} while (StillTreating);
				}

			} // End of PopulatedMicroCellNum loop, ThreadNum loop, and pragma.

		for (int i = 0; i < P.NumThreads; i++)
		{
			State.cumT += StateT[i].cumT;
			State.cumTP += StateT[i].cumTP;
			State.cumUT += StateT[i].cumUT;
			//State.cumV+=StateT[i].cumV;
			StateT[i].cumT = StateT[i].cumUT = StateT[i].cumTP = 0; //StateT[i].cumV=0;
		}
	}

	TreatFlag += TreatFlag1;
	return (TreatFlag > 0);
}

/**
 * Function: AddToInfectionQueue
 *
 * Purpose: add to the infection queue
 * @param tn - thread number
 * @param infectee_cell_number
 * @param infector_index
 * @param infectee_index
 * @param infect_type
 * @return void
 */
void AddToInfectionQueue(const int tn, const int infectee_cell_number, const int infector_index, const int infectee_index, const short int infect_type)
{
	StateT[tn].inf_queue[infectee_cell_number][StateT[tn].n_queue[infectee_cell_number]++] = { infector_index, infectee_index, infect_type };
}

/**
 * Function: AddInfections
 *
 * Purpose: add to the infection queue
 * @param tn - thread number
 * @param infectee_cell_index -  cell number of potential infectee
 * @param infector_index
 * @param infectee_index
 * @param infect_type
 * @return true if a positive entry is added to the queue
 */
bool AddInfections(const int tn, const int infectee_cell_index, const int infector_index, const int infectee_index,
	const short int infect_type)
{
	bool positiveEntryAdded = false;
	assert(infectee_cell_index >= 0);
	assert(infector_index >= 0);
	assert(infectee_index >= 0);
	assert(infect_type >= 0);
	
	// if there is space in queue for this thread
	if ((StateT[tn].n_queue[infectee_cell_index] < P.InfQueuePeakLength))
	{
		// if random number < false positive rate
		if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
			// add false positive to infection queue
			AddToInfectionQueue(tn, infectee_cell_index, -1, infectee_index, -1);
		else
		{
			// infect infectee_index

			// infect_type: first 4 bits store type of infection
			//				1= household
			//				2..MAX_NUM_PLACE_TYPES+1 = within-class/work-group place based transmission
			//				MAX_NUM_PLACE_TYPES+2..2*MAX_NUM_PLACE_TYPES+1 = between-class/work-group place based transmission
			//				2*MAX_NUM_PLACE_TYPES+2 = "spatial" transmission (spatially local random mixing)
			// bits >4 store the generation of infection

			AddToInfectionQueue(tn, infectee_cell_index, infector_index, infectee_index, infect_type);
			positiveEntryAdded = true;
		}
	}
	return positiveEntryAdded;
}
