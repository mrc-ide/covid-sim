#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "CalcInfSusc.h"
#include "Dist.h"
#include "Error.h"
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
	if (!(P.DoAirports && P.HotelPlaceType < P.PlaceTypeNum)) ERR_CRITICAL("DoAirports || HotelPlaceType not set\n");

	if (floor(1 + t + P.TimeStep) != floor(1 + t))
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
													{fprintf(stderr,"### %i %i %i %i\n",j,k,n,Places[P.HotelPlaceType][j].members[n]);ner++;}
												else if((k<0)||(k>n))
													{fprintf(stderr,"@ %i %i %i %i\n",j,k,n,Places[P.HotelPlaceType][j].members[n]);ner++;}
												else
						*/
						if (k != n)
						{
							Places[P.HotelPlaceType][j].members[k] = Places[P.HotelPlaceType][j].members[n];
						}
						nr++;
						if (Hosts[i].PlaceLinks[P.HotelPlaceType] != j)
						{
							ner++; fprintf(stderr_shared, "(%i %i) ", j, Hosts[i].PlaceLinks[P.HotelPlaceType]);
						}
						Hosts[i].PlaceLinks[P.HotelPlaceType] = -1;
						Hosts[i].Travelling = 0;
					}
				}
				Places[P.HotelPlaceType][j].n = n;
			}
		}
		fprintf(stderr, " d=%i e=%i>", nr, ner);
	}
}

void TravelDepartSweep(double t)
{
	int d, mps, nld, nad, nsk, bm;
	double nl;

	// Convince static analysers that values are set correctly:
	if (!(P.DoAirports && P.HotelPlaceType < P.PlaceTypeNum)) ERR_CRITICAL("DoAirports || HotelPlaceType not set\n");

	if (floor(1 + t - P.TimeStep) != floor(1 + t))
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
										if (P.distance_->distance_squared(Airports[i].loc, Airports[k].loc) > P.MoveRestrRadius2)
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
		fprintf(stderr, "<ar=%i as=%i", nad, nsk);
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
								double s2 = P.distance_->distance_squared(Households[Hosts[i2].hh].loc, Places[P.HotelPlaceType][i].loc);
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
		fprintf(stderr, " l=%i ls=%i ", nld, nsk);
	}
}

void InfectSweep(double t, int run) //added run number as argument in order to record it in event log
{
	//// This function takes the day number (t) and run number (run) as inputs. It loops over infected people, and decides whom to infect. Structure is 1) #pragma loop over all cells then 1a) infectious people, which chooses who they will infect, adds them to a queue
	//// Next 2) #pragma loop infects those people from queue (using DoInfect function). This is to avoid race conditions.
	//// Loop 1a) calculates the force of infection exerted by each infected person on (and therefore number of new infections to) i) their house; ii) their place(s); iii) other spatial cells.
	//// Each force of infection includes infectiousness and susceptibility components.
	//// Infectiousness is (broadly) a function of 1 person (their age, treatment status, places, no. people in their household etc.)
	//// Susceptibility is (broadly) a function of 2 people (a person's susceptibility TO ANOTHER PERSON / potential infector)
	//// After loop 1a) over infectious people, spatial infections are doled out.

	int n; //// number of people you could potentially infect in your place group, then number of potential spatial infections doled out by cell on other cells.
	int f, f2, cq /*cell queue*/, bm/*movement restrictions in place*/, ci /*person index*/;
	double seasonality, sbeta, hbeta;
	//// various quantities of force of infection, including "infectiousness" and "susceptibility" components
	double s; // household Force Of Infection (FOI) on fellow household member, then place susceptibility, then random number for spatial infections allocation* / ;
	double s2; // spatial infectiousness, then distance in spatial infections allocation
	double s3, s3_scaled; // household, then place infectiousness
	double s4, s4_scaled; // place infectiousness (copy of s3 as some code commented out
	double s5; //// total spatial infectiousness summed over all infectious people in cell.
	double s6;
	double fp; //// false positive
	unsigned short int ts;

	// if not doing seasonality
	if (!P.DoSeasonality)
	{	
		// set seasonality to 1
		seasonality = 1.0;
	}
	else{
		// otherwise pick seasonality from P.Seasonality array using day number in year
		seasonality = P.Seasonality[((int)t) % DAYS_PER_YEAR];
	}
	
	// ts = the timestep number of the start of the current day
	ts = (unsigned short int) (P.TimeStepsPerDay * t); 
	// fp = false positive
	fp = P.TimeStep / (1 - P.FalsePositiveRate);
	// sbeta seasonality beta
	sbeta = seasonality * fp * P.LocalBeta;
	
	// hbeta = household beta
	// if doing households, hbeta = seasonality * fp * P.HouseholdTrans, else hbeta = 0
	hbeta = (P.DoHouseholds) ? (seasonality * fp * P.HouseholdTrans) : 0;
	
	// Establish if movement restrictions are in place on current day - store in bm, 0:false, 1:true 
	bm = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));
	// File for storing error reports
	FILE* stderr_shared = stderr;
	
#pragma omp parallel for private(n,f,f2,s,s2,s3,s4,s5,s6,cq,ci,s3_scaled,s4_scaled) schedule(static,1) default(none) \
		shared(t, P, CellLookup, Hosts, AdUnits, Households, Places, SamplingQueue, Cells, Mcells, StateT, hbeta, sbeta, seasonality, ts, fp, bm, stderr_shared)
	for (int tn = 0; tn < P.NumThreads; tn++)
		for (int b = tn; b < P.NCP; b += P.NumThreads) //// loop over (in parallel) all populated cells. Loop 1)
		{
			Cell* c = CellLookup[b]; // select Cell given by index b
			s5 = 0; ///// spatial infectiousness summed over all infectious people in loop below
			
			//// Loop over array of indices of infectious people c->I in cell c. Loop 1a)
			for (int j = 0; j < c->I; j++) 
			{
				//// get person index ci of j'th infected person in cell
				ci = c->infected[j];
				//// get person si from Hosts (array of people) corresponding to ci, using pointer arithmetic.
				Person* si = Hosts + ci;

				//evaluate flag for digital contact tracing (fct) here at the beginning for each individual
				// fct = 1 if:
				// P.DoDigitalContactTracing = 1 (ie. digital contact tracing functionlity is switched on)
				// AND Day number (t) is greater than the start day for contact tracing in this administrative unit (ie. contact tracing has started)
				// AND Day number (t) is less than the end day for contact tracing in this administrative unit (ie. contact tracing has not ended)
				// AND the selected host is a digital contact tracing user
				// otherwise fct = 0
				int fct = ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart)
				&& (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1)); // && (ts <= (Hosts[ci].detected_time + P.usCaseIsolationDelay)));

				// BEGIN HOUSEHOLD INFECTIONS
				
				//// Household Force Of Infection (FOI) component
				
				// hbeta =  seasonality * fp * P.HouseholdTrans or 0 depending on whether households functionality is on or off - see start of function 
				
				// if household beta (hbeta) > 0 
				if (hbeta > 0)
				{
					// For selected host si's household (si->hh), 
					// if the number of hosts (nh) in that Household is greater than 1
					// AND the selected host is not travelling
					if ((Households[si->hh].nh > 1) && (!si->Travelling))
					{
						int l = Households[si->hh].FirstPerson;
						int m = l + Households[si->hh].nh;
						// calculate infectiousness of selected household (s3)
						// using the CalcHouseInf function on the selected cell and timestamp at start of current day
						// then scaling by hbeta
						s3 = hbeta * CalcHouseInf(ci, ts);  

						// Test if any of the individuals in the selected persons household are absent from places
						// f=0 means noone absent, f=1 means at least one absent
						f = 0; // initialise f to be 0
						for (int i3 = l; (i3 < m) && (!f); i3++){ //// loop over people in household
							for (int i2 = 0; (i2 < P.PlaceTypeNum) && (!f); i2++){ //// loop over place types
								if (Hosts[i3].PlaceLinks[i2] >= 0){ //// if person in household has any sort of link to place type
									// if person is absent set f=1
									f = ((PLACE_CLOSED(i2, Hosts[i3].PlaceLinks[i2]))&&(HOST_ABSENT(i3)));
								}
							}
						}

						// if individuals in the household are absent from places (ie. f==1 from test immediately above), scale up the infectiousness (s3) of the household
						if (f) { s3 *= P.PlaceCloseHouseholdRelContact; }/* NumPCD++;}*/ //// if people in your household are absent from places, person si/ci is more infectious to them, as they spend more time at home.
						
						// Loop from l (the index of the first person in the household) to m-1 (the index of the last person in the household)
						// ie. loop over everyone in the household
						for (int i3 = l; i3 < m; i3++) //// loop over all people in household (note goes from l to m - 1)
						{
							if (Hosts[i3].is_susceptible() && (!Hosts[i3].Travelling)) //// if people in household uninfected/susceptible and not travelling
							{
								s = s3 * CalcHouseSusc(i3, ts, ci);		//// FOI ( = infectiousness x susceptibility) from person ci/si on fellow household member i3
								
								// Force of Infection (s) > random value between 0 and 1
								if (ranf_mt(tn) < s) 
								{
									// explicitly cast to short to resolve level 4 warning
									const short int infect_type = static_cast<short int>(1 + INFECT_TYPE_MASK * (1 + si->
										infect_type / INFECT_TYPE_MASK));

									if (AddInfections(tn, Hosts[i3].pcell % P.NumThreads, ci, i3, infect_type))
									{
										// ** infect household member i3 **
										Hosts[i3].infector = ci; //// assign person ci as infector of person i3
									}
								}// if FOI > s
							} // if person in household uninfected/susceptible and not travelling
						}// loop over people in household
					} // if more than one person in household 
				}// if hbeta > 0
				
				// END HOUSHOLD INFECTIONS
				
				// BEGIN PLACE INFECTIONS
				
				// Still with infected person (si) = Hosts[ci]
				// if places functionality is enabled
				if (P.DoPlaces)
				{
					// if host with index ci isn't absent
					if (!HOST_ABSENT(ci))
					{
						// select microcell (mi) corresponding to selected host (si)
						Microcell* mi = Mcells + si->mcell;
						for (int k = 0; k < P.PlaceTypeNum; k++) //// loop over all place types
						{
							// select link (l) between selected host (si) and place from si's placelinks to place type k
							int l = si->PlaceLinks[k];
							if (l >= 0)  //// l>=0 means if place type k is relevant to person si. (Now allowing for partial attendance).
							{
								// infectiousness of place (s3)
								// = false positive rate * seasonality * place infectiousness
								s3 = fp * seasonality * CalcPlaceInf(ci, k, ts);
								// select microcell of the place linked to host si with link l
								Microcell* mp = Mcells + Places[k][l].mcell;
								// if blanket movement restrictions are in place on current day
								if (bm)
								{
									// if distance between si's household and linked place
									// is greater than movement restriction radius
									if ((P.distance_->distance_squared(Households[si->hh].loc, Places[k][l].loc) > P.MoveRestrRadius2))
									{
										// multiply infectiousness of place by movement restriction effect
										s3 *= P.MoveRestrEffect;
									}
								}
								// else if movement restrictions in effect in either household microcell or place microcell
								else if ((mi->moverest != mp->moverest) && ((mi->moverest == TreatStat::Treated) || (mp->moverest == TreatStat::Treated)))
								{
									// multiply infectiousness of place by movement restriction effect
									s3 *= P.MoveRestrEffect;
								}
								
								// BEGIN NON-HOTEL INFECTIONS
								
								// if linked place isn't a hotel and selected host isn't travelling
								if ((k != P.HotelPlaceType) && (!si->Travelling))
								{
									// i2 is index of group (of place type k) that selected host is linked to 
									int i2 = (si->PlaceGroupLinks[k]);
									
									// calculate infectiousness (s4_scaled) 
									// which varies if contact tracing is in place
									// if contact tracing isn't in place s4_scaled is a copy of s3 
									// if contact tracing is in place, s4_scaled is s3  * P.ScalingFactorPlaceDigitalContacts
									// in either case s4_scaled is capped at 1
									
									// if contact tracing
									if (fct)
									{
										// copy s3
										s4 = s3;
										// multiply s4 by P.ScalingFactorPlaceDigitalContacts 
										s4_scaled = s4 *P.ScalingFactorPlaceDigitalContacts;
										// cap s4 at 1
										if (s4 > 1) s4 = 1;
										// cap at 1
										if (s4_scaled > 1) s4_scaled = 1;
									}
									else
									{
										// copy s3 to s4
										s4 = s3;
										// cap s4 at 1
										if (s4 > 1) s4 = 1;
										s4_scaled = s4;
									}

									// if infectiousness is < 0, we have an error - end the program
									if (s4_scaled < 0)
									{
										fprintf(stderr_shared, "@@@ %lg\n", s4_scaled);
										exit(1);
									}
									// else if infectiousness == 1 (should never be more than 1 due to capping above)
									else if (s4_scaled >= 1)	//// if place infectiousness above threshold, consider everyone in group a potential infectee...
									{
										// set n to be number of people in group in place k,l
										n = Places[k][l].group_size[i2];
									}
									else				//// ... otherwise randomly sample (from binomial distribution) number of potential infectees in this place.
									{
										n = (int)ignbin_mt((int32_t)Places[k][l].group_size[i2], s4_scaled, tn);
									}
									
									// if potential infectees > 0	
									if (n > 0) 
									{
										// pick n members of place k,l and add them to sampling queue for thread tn
										SampleWithoutReplacement(tn, n, Places[k][l].group_size[i2]); //// changes thread-specific SamplingQueue.
									}
									
									// loop over sampling queue of potential infectees
									for (int m = 0; m < n; m++)
									{
										// pick potential infectee index i3
										int i3 = Places[k][l].members[Places[k][l].group_start[i2] + SamplingQueue[tn][m]];
										// calculate place susceptbility based on infectee (i3), place type (k), timestep (ts)
										// cell (ci) and thread number (tn)
										s = CalcPlaceSusc(i3, k, ts);

										// allow care home residents to mix more intensely in "groups" (i.e. individual homes) than staff do - to allow for PPE/environmental contamination.
										if ((k==P.CareHomePlaceType)&&((!Hosts[ci].care_home_resident)||(!Hosts[i3].care_home_resident))) s *= P.CareHomeWorkerGroupScaling;
										//these are all place group contacts to be tracked for digital contact tracing - add to StateT queue for contact tracing
										//if infectee is also a user, add them as a contact
										
										if ((fct) && (Hosts[i3].digitalContactTracingUser) && (ci != i3) && (!HOST_ABSENT(i3)))
										{
											// scale place susceptibility by proportion who self isolate and store as s6
											s6 = P.ProportionDigitalContactsIsolate * s;
											// if random number < s6
											// AND number of contacts of ci(!) is less than maximum digital contact to trace
											if ((Hosts[ci].ncontacts < P.MaxDigitalContactsToTrace) && (ranf_mt(tn) <s6))
											{
												Hosts[ci].ncontacts++; //add to number of contacts made
												int ad = Mcells[Hosts[i3].mcell].adunit;
												if ((StateT[tn].ndct_queue[ad] < AdUnits[ad].n))
												{
													//find adunit for contact and add both contact and infectious host to lists - storing both so I can set times later.
													StateT[tn].dct_queue[ad][StateT[tn].ndct_queue[ad]++] = { i3,ci,ts };
												}
												else
												{
													fprintf(stderr_shared, "No more space in queue! Thread: %i, AdUnit: %i\n", tn, ad);
												}
											}
										}

										if (Hosts[i3].is_susceptible() && (!HOST_ABSENT(i3))) //// if person i3 uninfected and not absent.
										{
											Microcell* mt = Mcells + Hosts[i3].mcell;
											//downscale s if it has been scaled up do to digital contact tracing
											s *= CalcPersonSusc(i3, ts, ci)*s4/s4_scaled;

											// if blanket movement restrictions are in place
											if (bm)
											{
												// if potential infectee i3's household is further from selected place
												if ((P.distance_->distance_squared(Households[Hosts[i3].hh].loc, Places[k][l].loc) > P.MoveRestrRadius2))
												{
													// multiply susceptibility by movement restriction effect
													s *= P.MoveRestrEffect;
												}
											}
											// else if movement restrictions are in place in either cell
											else if ((mt->moverest != mp->moverest) && ((mt->moverest == TreatStat::Treated) || (mp->moverest == TreatStat::Treated)))
											{
												// multiply susceptibility by movement restriction effect
												s *= P.MoveRestrEffect;
											}

											// if either susceptiblity is 100% or sample probability s
											if ((s == 1) || (ranf_mt(tn) < s))
											{
												// explicitly cast to short to resolve level 4 warning
												const short int infect_type = static_cast<short int> (2 + k + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK));

												AddInfections(tn, Hosts[i3].pcell % P.NumThreads, ci, i3, infect_type);
											}
										}
									}
								}
								
								// END NON-HOTEL INFECTIONS
								
								// BEGIN HOTEL INFECTIONS
								
								// if selected host si is not travelling or selected link is to a hotel
								if ((k == P.HotelPlaceType) || (!si->Travelling))
								{
									s3 *= P.PlaceTypePropBetweenGroupLinks[k] * P.PlaceTypeGroupSizeParam1[k] / ((double)Places[k][l].n);
									if (s3 > 1) s3 = 1;
									// if contact tracing in place, multiply s3_scaled = s3*scalingfactor, otherwise s3_scaled = s3
									s3_scaled = (fct) ? (s3 * P.ScalingFactorPlaceDigitalContacts) : s3;
									// s3_scales shouldn't be less than 0 so generate error if it is
									if (s3_scaled < 0)
									{
										ERR_CRITICAL_FMT("@@@ %lg\n", s3);
									}
									// if s3_scaled >=1, everyone in the hotel is a potential infectee
									else if (s3_scaled >= 1)
										n = Places[k][l].n;
									// if s3_scaled between 0 and 1, decide number of potential infectees based on
									// using ignbin_mt function
									else
										n = (int)ignbin_mt((int32_t)Places[k][l].n, s3_scaled, tn);
									// if more than 0 potential infectees, pick n hosts from the hotel and add to sampling queue
									if (n > 0) SampleWithoutReplacement(tn, n, Places[k][l].n);
									// loop over the sampling queue
									for (int m = 0; m < n; m++)
									{
										// select potential infectee from sampling queue
										int i3 = Places[k][l].members[SamplingQueue[tn][m]];
										// calculate place susceptibility s
										s = CalcPlaceSusc(i3, k, ts);
										// use group structure to model multiple care homes with shared staff - in which case residents of one "group" don't mix with those in another, only staff do.
										if ((Hosts[ci].care_home_resident) && (Hosts[i3].care_home_resident) && (Hosts[ci].PlaceGroupLinks[k]!= Hosts[i3].PlaceGroupLinks[k])) s *= P.CareHomeResidentPlaceScaling;
										// allow care home staff to have lowere contacts in care homes - to allow for PPE/environmental contamination.
										if ((k == P.CareHomePlaceType) && ((!Hosts[ci].care_home_resident) || (!Hosts[i3].care_home_resident))) s *= P.CareHomeWorkerGroupScaling;
										
										//these are all place group contacts to be tracked for digital contact tracing - add to StateT queue for contact tracing

										//if infectee is also a user, add them as a contact
										
										// if contact tracing in place AND potential infectee i3 is a contact tracing user AND i3 isn't absent AND i3 isn't ci (suspect this should be si)

										if ((fct) && (Hosts[i3].digitalContactTracingUser) && (ci != i3) && (!HOST_ABSENT(i3)))
										{
											// s6 = place susceptibility * proportion of digital contacts who self isolate
											s6 = P.ProportionDigitalContactsIsolate * s;
											// if number of contacts of infectious person < maximum and random number < s6
											if ((Hosts[ci].ncontacts < P.MaxDigitalContactsToTrace) && (ranf_mt(tn) < s6))
											{
												Hosts[ci].ncontacts++; //add to number of contacts made
												int ad = Mcells[Hosts[i3].mcell].adunit;
												if ((StateT[tn].ndct_queue[ad] < AdUnits[ad].n))
												{
													//find adunit for contact and add both contact and infectious host to lists - storing both so I can set times later.
													StateT[tn].dct_queue[ad][StateT[tn].ndct_queue[ad]++] = { i3,ci,ts };
												}
												else
												{
													fprintf(stderr_shared, "No more space in queue! Thread: %i, AdUnit: %i\n", tn, ad);
												}
											}
										}

										// if potential infectee i3 uninfected and not absent.
										if (Hosts[i3].is_susceptible() && (!HOST_ABSENT(i3)))
										{
											// mt = microcell of potential infectee
											Microcell* mt = Mcells + Hosts[i3].mcell;

											//if doing digital contact tracing, scale down susceptibility here
											s*= CalcPersonSusc(i3, ts, ci)*s3/s3_scaled;
											// if blanket movement restrictions are in place
											if (bm)
											{
												// if potential infectees household is farther away from hotel than restriction radius
												if ((P.distance_->distance_squared(Households[Hosts[i3].hh].loc, Places[k][l].loc) > P.MoveRestrRadius2))
												{
													// multiply susceptibility by movement restriction effect
													s *= P.MoveRestrEffect;
												}
											}
											// else if movement restrictions are in place in potential infectee's cell or hotel's cell
											else if ((mt->moverest != mp->moverest) && ((mt->moverest == TreatStat::Treated) || (mp->moverest == TreatStat::Treated)))
											{
												// multiply susceptibility by movement restriction effect
												s *= P.MoveRestrEffect;
											}
											
											// ** do infections **
											
											// is susceptibility is 1 (ie infect everyone) or random number is less than susceptibility
											if ((s == 1) || (ranf_mt(tn) < s))
											{
												// explicitly cast to short to resolve level 4 warning
												const short int infect_type = static_cast<short int> (2 + k + NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK));
												
												AddInfections(tn, Hosts[i3].pcell% P.NumThreads, ci, i3, infect_type);
											}// susceptibility test
										}// potential infectee i3 uninfected and not absent.
									}// loop over sampling queue
								}// selected host si is not travelling or selected link is to a hotel
								
								// ** END HOTEL INFECTIONS ** 
								
							}// if place link relevant
						}// loop over place types
					}// if host isn't absent
				}// if places functionality enabled
				
				// END PLACE INFECTIONS
				
				// BEGIN SPATIAL INFECTIONS
				
				//// First determine spatial FOI component (s5)
				
				// if seasonality beta > 0
				// do spatial infections 
				//// ie sum spatial infectiousness over all infected people, the infections from which are allocated after loop over infected people.
				if (sbeta > 0) 
				{
					if (si->Travelling) //// if host currently away from their cell, they cannot add to their cell's spatial infectiousness.
					{
						s2 = 0; 
						f = 0;
					}
					else
					{
						// calculate spatial infectiousness (s2) based on host and timestep
						s2 = CalcSpatialInf(ci, ts);
						//if do digital contact tracing, scale up spatial infectiousness of infectives who are using the app and will be detected
						if (fct)
						{
							s2 *= P.ScalingFactorSpatialDigitalContacts;
						}
					}
					// test if selected person si is linked to a place that is closed, f=0 means no links to closed places, otherwise f=1
					f = 0; // initialise f as 0
					// If place functionality switched on
					if (P.DoPlaces)
					{
						// loop over place types until closed place is found
						for (int i3 = 0; (i3 < P.PlaceTypeNum) && (!f); i3++)
						{
							if (si->PlaceLinks[i3] >= 0) //// if person has a link to place of type i3...
							{
								// if place is closed set f=1
								f = PLACE_CLOSED(i3, si->PlaceLinks[i3]); //// find out if that place of type i3 is closed.
							}
						}
					}// if doing places

					if((f) && (HOST_ABSENT(ci))) //// if place is closed and person is absent then adjust the spatial infectiousness (similar logic to household infectiousness: place closure affects spatial infectiousness
					{
						s2 *= P.PlaceCloseSpatialRelContact;
						/* NumPCD++; */
						s5 += s2;
						StateT[tn].cell_inf[j] = (float)-s5;
					}
					else
					{
						s5 += s2;
						StateT[tn].cell_inf[j] = (float)s5;
					}
				}
			} // loop over infectious people in cell
			
			
			//// Now allocate spatial infections using Force Of Infection (s5) calculated above
			if (s5 > 0) //// if spatial infectiousness positive
			{
				
				// decide how many potential cell to cell infections this cell could cause  
				n = (int)ignpoi_mt(s5 * sbeta * ((double)c->tot_prob), tn); //// number people this cell's population might infect elsewhere. poisson random number based on spatial infectiousness s5, sbeta (seasonality) and this cell's "probability" (guessing this is a function of its population and geographical size).
				// i2 = number of infectious people in cell c
				int i2 = c->I;
				
				if (n > 0) //// this block normalises cumulative infectiousness cell_inf by person. s5 is the total cumulative spatial infectiousness. Reason is so that infector can be chosen using ranf_mt, which returns random number between 0 and 1.
				{
					//// normalise by cumulative spatial infectiousness.
					for (int j = 0; j < i2 - 1; j++) StateT[tn].cell_inf[j] /= ((float) s5);
					//// does same as the above loop just a slightly faster calculation. i.e. StateT[tn].cell_inf[i2 - 1] / s5 would equal 1 or -1 anyway.
					StateT[tn].cell_inf[i2 - 1] = (StateT[tn].cell_inf[i2 - 1] < 0) ? -1.0f : 1.0f;
				}
				
				//// loop over infections to dole out. roughly speaking, this determines which infectious person in cell c infects which person elsewhere.
				for (int k = 0; k < n; k++)  
				{
					//// decide on infector ci/si from cell c.
					int j; // j = index of infector
					// if only one infectious person in cell
					if (i2 == 1)
					{
						j = 0; // infector index is first in cell (person 0)
					}
					// if more than one infectious person in cell pick an infectious person (given by index j)
					//// roughly speaking, this determines which infectious person in cell c infects which person elsewhere
					else				
					{
						int m;
						s = ranf_mt(tn);	///// choose random number between 0 and 1
						j = m = i2 / 2;		///// assign j and m to be halfway between zero and number of infected people i2 = c->I.
						f = 1;
						do
						{
							if (m > 1) m /= 2; //// amount m to change j by reduced by half. Looks like a binary search. Basically saying, keep amending potential infector j until either j less than zero or more than number of infected people until you find j s.t. spatial infectiousness "matches" s.
							if ((j > 0) && (fabs(StateT[tn].cell_inf[j - 1]) >= s))
							{
								j -= m;
								if (j == 0)			f = 0;
							}
							else if ((j < i2 - 1) && (fabs(StateT[tn].cell_inf[j]) < s))
							{
								j += m;
								if (j == i2 - 1)	f = 0;
							}
							else					f = 0;
						} while (f);
					}
					f = (StateT[tn].cell_inf[j] < 0); //// flag for whether infector j had their place(s) closed. <0 (true) = place closed / >=0 (false) = place not closed. Set in if (sbeta > 0) part of loop over infectious people.
					// ci is the index of the jth infectious person in the cell
					ci = c->infected[j];
					// si is the jth selected person in the cell
					Person* si = Hosts + ci;

					//calculate flag (fct) for digital contact tracing here at the beginning for each individual infector
					int fct = ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart)
						&& (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1)); // && (ts <= (Hosts[ci].detected_time + P.usCaseIsolationDelay)));


					//// decide on infectee
					
					// do the following while f2=0 
					do
					{
						//// chooses which cell person will infect
						// pick random s between 0 and 1
						s = ranf_mt(tn);
						// generate l using InvCDF of selected cell and random integer between 0 and 1024
						int l = c->InvCDF[(int)floor(s * 1024)];
						// loop over c->cum_trans array until find a value >= random number s
						while (c->cum_trans[l] < s) l++;
						// selecte the cell corresponding to l
						Cell* ct = CellLookup[l];

						///// pick random person m within susceptibles of cell ct (S0 initial number susceptibles within cell).
						int m = (int)(ranf_mt(tn) * ((double)ct->S0));
						int i3 = ct->susceptible[m];
						
						s2 = dist2(Hosts + i3, Hosts + ci); /// calculate distance squared between this susceptible person and person ci/si identified earlier
						s = P.KernelLookup.num(s2) / c->max_trans[l]; //// acceptance probability
						
						// initialise f2=0 (f2=1 is the while condition for this loop)
						f2 = 0;
						// if random number greater than acceptance probablility or infectee is dead
						if ((ranf_mt(tn) >= s) || Hosts[i3].is_dead()) //// if rejected, or infectee i3/m already dead, ensure do-while evaluated again (i.e. choose a new infectee).
						{
							// set f2=1 so loop continues 
							f2 = 1;
						}
						else
						{
							//// if potential infectee not travelling, and either is not part of cell c or doesn't share a household with infector.
							if ((!Hosts[i3].Travelling) && ((c != ct) || (Hosts[i3].hh != si->hh))) 
							{
								// pick microcell of infector (mi)
								Microcell* mi = Mcells + si->mcell;
								// pick microcell of infectee (mt)
								Microcell* mt = Mcells + Hosts[i3].mcell;
								s = CalcSpatialSusc(i3, ts);
								// Care home residents may have fewer contacts
								if ((Hosts[i3].care_home_resident) || (Hosts[ci].care_home_resident)) s *= P.CareHomeResidentSpatialScaling;
								//so this person is a contact - but might not be infected. if we are doing digital contact tracing, we want to add the person to the contacts list, if both are users
								if (fct)
								{
									//if infectee is also a user, add them as a contact
									if (Hosts[i3].digitalContactTracingUser && (ci != i3))
									{
										if ((Hosts[ci].ncontacts < P.MaxDigitalContactsToTrace) && (ranf_mt(tn) < s * P.ProportionDigitalContactsIsolate))
										{
											Hosts[ci].ncontacts++; //add to number of contacts made
											int ad = Mcells[Hosts[i3].mcell].adunit;
											if ((StateT[tn].ndct_queue[ad] < AdUnits[ad].n))
											{
												//find adunit for contact and add both contact and infectious host to lists - storing both so I can set times later.
												StateT[tn].dct_queue[ad][StateT[tn].ndct_queue[ad]++] = { i3,ci,ts };
											}
											else
											{
												fprintf(stderr_shared, "No more space in queue! Thread: %i, AdUnit: %i\n", tn, ad);
											}
										}
									}
									//scale down susceptibility so we don't over accept
									s /= P.ScalingFactorSpatialDigitalContacts;
								}
								
							
								if (m < ct->S)  // only bother trying to infect susceptible people
								{
									s *= CalcPersonSusc(i3, ts, ci);
									if (bm)
									{
										if ((P.distance_->distance_squared(Households[si->hh].loc, Households[Hosts[i3].hh].loc) > P.MoveRestrRadius2))
											s *= P.MoveRestrEffect;
									}
									else if ((mt->moverest != mi->moverest) && ((mt->moverest == TreatStat::Treated) || (mi->moverest == TreatStat::Treated)))
										s *= P.MoveRestrEffect;
									if ((!f)&& (HOST_ABSENT(i3))) //// if infector did not have place closed, loop over place types of infectee i3 to see if their places had closed. If they had, amend their susceptibility.
									{
										for (m = f2 = 0; (m < P.PlaceTypeNum) && (!f2); m++)
											if (Hosts[i3].PlaceLinks[m] >= 0)
											{
												f2 = PLACE_CLOSED(m, Hosts[i3].PlaceLinks[m]);
											}
										if (f2) { s *= P.PlaceCloseSpatialRelContact; }/* NumPCD++;} */
										f2 = 0;
									}
									if ((s == 1) || (ranf_mt(tn) < s)) //// accept/reject
									{
										cq = ((int)(ct - Cells)) % P.NumThreads;

										if (Hosts[i3].is_susceptible())
										{
											// explicitly cast to short to resolve level 4 warning
											const short int infect_type = static_cast<short int>(2 + 2 * NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK));
											
											AddInfections(tn, cq, ci, i3, infect_type);
										}
									}
								}// m < susceptible people in target cell
							}// //// if potential infectee not travelling, and either is not part of cell c or doesn't share a household with infector
						}// infectee isn't dead
					} while (f2);
				}// loop over infections doled out by cell
			}// s5 > 0
		}


#pragma omp parallel for schedule(static,1) default(none) \
		shared(t, run, P, StateT, Hosts, ts)
	for (int j = 0; j < P.NumThreads; j++)
	{
		for (int k = 0; k < P.NumThreads; k++)
		{
			for (int i = 0; i < StateT[k].n_queue[j]; i++)
			{
				int infector = StateT[k].inf_queue[j][i].infector;
				int infectee = StateT[k].inf_queue[j][i].infectee;
				short int infect_type = StateT[k].inf_queue[j][i].infect_type;
				Hosts[infectee].infector = infector;
				Hosts[infectee].infect_type = infect_type;
				if (infect_type == -1) //// i.e. if host doesn't have an infector
					DoFalseCase(infectee, t, ts, j);
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
	unsigned short int ts; //// this timestep
	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	if (P.DoPlaces)
		for (int i = 0; i < P.NumHolidays; i++)
		{
			ht = P.HolidayStartTime[i] + P.HolidaysStartDay_SimTime;
			if ((t + P.TimeStep >= ht) && (t < ht))
			{
//				fprintf(stderr, "Holiday %i t=%lg\n", i, t);
				for (int j = 0; j < P.PlaceTypeNum; j++)
				{
#pragma omp parallel for schedule(static,1) default(none) shared(P, Places, Hosts, i, j, ht)
					for(int tn=0;tn<P.NumThreads;tn++)
						for (int k = tn; k < P.Nplace[j]; k+=P.NumThreads)
						{
							if ((P.HolidayEffect[j] < 1) && ((P.HolidayEffect[j] == 0) || (ranf_mt(tn) >= P.HolidayEffect[j])))
							{
								int l = (int)(ht * P.TimeStepsPerDay);
								if (Places[j][k].close_start_time > l)  Places[j][k].close_start_time = (unsigned short) l;
								int b = (int)((ht + P.HolidayDuration[i]) * P.TimeStepsPerDay);
								if (Places[j][k].close_end_time < b)	  Places[j][k].close_end_time = (unsigned short) b;
								for (int ci = 0; ci < Places[j][k].n; ci++)
								{
									if (Hosts[Places[j][k].members[ci]].absent_start_time > l) Hosts[Places[j][k].members[ci]].absent_start_time = (unsigned short)l;
									if (Hosts[Places[j][k].members[ci]].absent_stop_time < b) Hosts[Places[j][k].members[ci]].absent_stop_time = (unsigned short)b;
								}
							}
						}
				}
			}
		}

#pragma omp parallel for schedule(static,1) default(none) shared(t, P, CellLookup, Hosts, AdUnits, Mcells, StateT, ts)
	for (int tn = 0; tn < P.NumThreads; tn++)	//// loop over threads
		for (int b = tn; b < P.NCP; b += P.NumThreads)	//// loop/step over populated cells
		{
			Cell* c = CellLookup[b]; //// find (pointer-to) cell.
			for (int j = ((int)c->L - 1); j >= 0; j--) //// loop backwards over latently infected people, hence it starts from L - 1 and goes to zero. Runs backwards because of pointer swapping?
				if (ts >= Hosts[c->latent[j]].latent_time) //// if now after time at which person became infectious (latent_time a slight misnomer).
					DoIncub(c->latent[j], ts, tn); //// move infected person from latently infected (L) to infectious (I), but not symptomatic
			//StateT[tn].n_queue[0] = StateT[tn].n_queue[1] = 0;
			for (int j = c->I - 1; j >= 0; j--) ///// loop backwards over Infectious people. Runs backwards because of pointer swapping?
			{
				int ci = c->infected[j];	//// person index
				Person* si = Hosts + ci;		//// person

				unsigned short int tc; //// time at which person becomes case (i.e. moves from infectious and asymptomatic to infectious and symptomatic).
				/* Following line not 100% consistent with DoIncub. All severity time points (e.g. SARI time) are added to latent_time, not latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep))*/
				tc = si->latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep)); //// time that person si/ci becomes case (symptomatic)...
				if ((P.DoSymptoms) && (ts == tc)) //// ... if now is that time...
					DoCase(ci, t, ts, tn);		  //// ... change infectious (but asymptomatic) person to infectious and symptomatic. If doing severity, this contains DoMild and DoILI.

				if (P.DoSeverity)
				{
					if (ts >= si->SARI_time)					DoSARI(ci, tn);	//// see if you can dispense with inequalities by initializing SARI_time, Critical_time etc. to USHRT_MAX
					if (ts >= si->Critical_time)				DoCritical(ci, tn);
					if (ts >= si->RecoveringFromCritical_time)	DoRecoveringFromCritical(ci, tn);
					if (ts >= si->recovery_or_death_time)
					{
						if (si->to_die)
							DoDeath_FromCriticalorSARIorILI(ci, tn);
						else
							DoRecover_FromSeverity(ci, tn);
					}
				}

				//Adding code to assign recovery or death when leaving the infectious class: ggilani - 22/10/14
				if (ts >= si->recovery_or_death_time)
				{
					if (!si->to_die) //// if person si recovers and this timestep is after they've recovered
					{
						DoRecover(ci, tn);
						//StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++] = ci; //// add them to end of 0th thread of inf queue. Don't get why 0 here.
					}
					else /// if they die and this timestep is after they've died.
					{
						if (HOST_TREATED(ci) && (ranf_mt(tn) < P.TreatDeathDrop))
							DoRecover(ci, tn);
						else
							DoDeath(ci, tn);
					}

					//once host recovers, will no longer make contacts for contact tracing - if we are doing contact tracing and case was infectious when contact tracing was active, increment state vector
					if ((P.DoDigitalContactTracing) && (Hosts[ci].latent_time>= AdUnits[Mcells[Hosts[ci].mcell].adunit].DigitalContactTracingTimeStart) && (Hosts[ci].recovery_or_death_time < AdUnits[Mcells[Hosts[ci].mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1) && (P.OutputDigitalContactDist))
					{
						if (Hosts[ci].ncontacts > MAX_CONTACTS) Hosts[ci].ncontacts = MAX_CONTACTS;
						//increment bin in State corresponding to this number of contacts
						StateT[tn].contact_dist[Hosts[ci].ncontacts]++;
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
	unsigned short int ts;

	//find current time step
	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	FILE* stderr_shared = stderr;
#pragma omp parallel for schedule(static,1) default(none) \
		shared(t, P, AdUnits, StateT, Hosts, ts, stderr_shared)
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
						if (dct_start_time == ts)
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
									//At this point, we do testing on index cases who have been picked up on symptoms alone, in order to figure out whether and when
									//to remove their contacts (if P.RemoveContactsOfNegativeIndexCase). It's much harder to do it in the next loop as we don't have all
									//the information about the contact event there and would need to loop over all contacts again to look for their index case
									//This would cause race conditions due to having a loop over adunits within threaded loop over admin units
									//Only set test times if P.DoDCTTest. If P.DoDCTTest==0, but we are finding contacts of contacts, we check to see if contacts should become index cases every day they are in isolation.
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
									fprintf(stderr_shared, "No more space in queue! AdUnit: %i, ndct=%i, max queue length: %i\n", i, AdUnits[i].ndct, AdUnits[i].n);
									fprintf(stderr_shared, "Error!\n");
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
						else if ((dct_start_time == (USHRT_MAX - 1)) && (dct_end_time == ts))
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
		shared(t, P, AdUnits, Hosts, ts)
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
						if ((Hosts[contact].dct_test_time == ts) && (Hosts[contact].index_case_dct == 0))
						{
							//if host is positive
							if (Hosts[contact].is_infectious_asymptomatic_not_case() ||
								Hosts[contact].is_case() ||
								Hosts[contact].is_infectious_almost_symptomatic())
							{
								//if the test is a false negative
								if ((P.SensitivityDCT == 0) || ((P.SensitivityDCT < 1) && (ranf_mt(tn) >= P.SensitivityDCT)))
								{
									Hosts[contact].dct_end_time = ts;
								}
								//else if a true positive
								else if (P.FindContactsOfDCTContacts)
								{
									//set them to be an index case
									Hosts[contact].index_case_dct = 1;
									//set trigger time to pick up their contacts in the next time step
									Hosts[contact].dct_trigger_time = ts + 1; //added the +1 here so that if there are no delays, the contacts will still get picked up correctly
									//if they are an infectious, asymptomatic non-case, call DoDetectedCase in order to trigger HQ and PC too.
									if (Hosts[contact].is_infectious_asymptomatic_not_case())
									{
										DoDetectedCase(contact, t, ts, tn);
										Hosts[contact].detected = 1; Hosts[contact].detected_time = ts;
									}
								}
							}
							//or if host is negative
							else
							{
								//and is a true negative
								if ((P.SpecificityDCT == 1) || ((P.SpecificityDCT > 0) && (ranf_mt(tn) < P.SpecificityDCT)))
								{
									Hosts[contact].dct_end_time = ts;
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
							//if ((Hosts[contact].dct_test_time == ts) && (Hosts[contact].index_case_dct == 0) && ((abs(Hosts[contact].inf) == 2) || (Hosts[contact].inf == -1)))
						{
							//set them to be an index case
							Hosts[contact].index_case_dct = 1;
							//set trigger time to pick up their contacts in the next time step
							Hosts[contact].dct_trigger_time = ts + 1; //added the +1 here so that if there are no delays, the contacts will still get picked up correctly
							//if they are asymptomatic, i.e. specifically if they have inf flag 2, call DoDetectedCase in order to trigger HQ and PC too.
							if (Hosts[contact].is_infectious_asymptomatic_not_case())
							{
								DoDetectedCase(contact, t, ts, tn);
								Hosts[contact].detected = 1; Hosts[contact].detected_time = ts;
							}
						}
					}

					//now remove hosts who have reached the end of their isolation time
					if (Hosts[contact].dct_end_time == ts)
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

	int TreatFlag, TreatFlag1; //// Function returns TreatFlag. If TreatFlag == 0, function no longer called. Anytime any treatment used, TreatFlag set to 1. 
	int f2, f3, f4; //// various fail conditions. Used for other things
	int nckwp;

	//// time steps
	unsigned short int ts;								////  time-step now.
	unsigned short int t_TreatStart;					////  time-step treatment begin
	unsigned short int t_TreatEnd;						////  time-step treatment finish
	unsigned short int t_VacStart;						////  time-step vaccination begin
	unsigned short int t_PlaceClosure_End;				////  time-step place closure finish
	unsigned short int t_MoveRestrict_Start;			////  time-step movement restriction begin
	unsigned short int t_MoveRestrict_End;				////  time-step movement restriction finish
	unsigned short int t_SocDist_End;					////  time-step social distancing finish
	unsigned short int t_KeyWorkerPlaceClosure_End;		////  time-step key worker place closure finish
	int global_trig;
	double r;

	ts = (unsigned short int) (P.TimeStepsPerDay * t);
	TreatFlag = TreatFlag1 = 0;
	if (P.DoGlobalTriggers)
	{
		if (P.DoPerCapitaTriggers)
			global_trig = (int)floor(((double)State.trigDetectedCases) * P.GlobalIncThreshPop / ((double)P.PopSize));
		else
			global_trig = State.trigDetectedCases;
	}
	else
		global_trig = 0;

	///// block loops over places (or place groups if P.DoPlaceGroupTreat == 1) and determines whom to prophylactically treat
	if ((P.DoPlaces) && (t >= P.TreatTimeStart) && (t < P.TreatTimeStart + P.TreatPlaceGeogDuration) && (State.cumT < P.TreatMaxCourses))
	{
		t_TreatEnd = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatDelayMean + P.TreatProphCourseLength));

#pragma omp parallel for private(TreatFlag) reduction(+:TreatFlag1) schedule(static,1) default(none) \
			shared(P, StateT, Places, Hosts, ts, t_TreatEnd)
		for (int Thread = 0; Thread < P.NumThreads; Thread++)
			for (int PlaceType = 0; PlaceType < P.PlaceTypeNum; PlaceType++)
			{
				for (int PlaceNumQueueIndex = 0; PlaceNumQueueIndex < StateT[Thread].np_queue[PlaceType]; PlaceNumQueueIndex++) //// loop over all plaes IN QUEUE, not all a places
				{
					int PlaceNum = StateT[Thread].p_queue[PlaceType][PlaceNumQueueIndex]; //// note PlaceNum is index of place, not index of place in place queue.
					if (P.DoPlaceGroupTreat)
					{
						int PlaceGroupIndex = StateT[Thread].pg_queue[PlaceType][PlaceNumQueueIndex];
						TreatFlag = PlaceGroupIndex; //// keep this as a flag 
						for (int PG_member = ((int)Places[PlaceType][PlaceNum].group_start[PlaceGroupIndex]); PG_member < ((int)(Places[PlaceType][PlaceNum].group_start[PlaceGroupIndex] + Places[PlaceType][PlaceNum].group_size[PlaceGroupIndex])); PG_member++) // loop over people in place group.
						{
							/*if((Places[PlaceType][PlaceNum].members[PG_member]<0)||(Places[PlaceType][PlaceNum].members[PG_member]>P.PopSize-1))
								fprintf(stderr,"\n*** npq=%i gn=%i h=%i PG_member=%i PlaceType=%i PlaceNum=%i PlaceGroupIndex=%i s=%i n=%i ***\n",
									StateT[Thread].np_queue[PlaceType],
									Places[PlaceType][PlaceNum].n,
									Places[PlaceType][PlaceNum].members[PG_member],
									PG_member,PlaceType,PlaceNum,PlaceGroupIndex,
									(int) Places[PlaceType][PlaceNum].group_start[PlaceGroupIndex],
									(int) Places[PlaceType][PlaceNum].group_size[PlaceGroupIndex]);
							else
							*/
							if ((!HOST_TO_BE_TREATED(Places[PlaceType][PlaceNum].members[PG_member])) && ((P.TreatPlaceTotalProp[PlaceType] == 1) || (ranf_mt(Thread) < P.TreatPlaceTotalProp[PlaceType])))
								DoProph(Places[PlaceType][PlaceNum].members[PG_member], ts, Thread);
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
										DoProph(Places[PlaceType][PlaceNum].members[PG_member], ts, Thread);
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
				shared(State, m, ts)
			for (int i = State.mvacc_cum; i < m; i++)
				DoVacc(State.mvacc_queue[i], ts);
			State.mvacc_cum = m;
		}
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

#pragma omp parallel for private(f2,f3,f4,r) reduction(+:TreatFlag) schedule(static,1) default(none) \
			shared(t, P, Hosts, Mcells, McellLookup, AdUnits, State, global_trig, ts, t_TreatEnd, t_TreatStart, t_VacStart, t_PlaceClosure_End, t_MoveRestrict_End, t_MoveRestrict_Start, t_SocDist_End, t_KeyWorkerPlaceClosure_End, nckwp)
		for (int tn = 0; tn < P.NumThreads; tn++)
			for (int bs = tn; bs < P.NMCP; bs += P.NumThreads) //// loop over populated microcells
			{
				int mcellnum = (int)(McellLookup[bs] - Mcells); //// microcell number
				int adi = (P.DoAdUnits) ? Mcells[mcellnum].adunit : -1;
				int ad = (P.DoAdUnits) ? AdUnits[adi].id : 0;

					//// Code block goes through various types of treatments/interventions (vaccination/movement restrictions etc.),
					//// assesses whether various triggers (counts) are over a certain threshold, (specified in ReadParams)
					//// and then implements those treatments by setting various flags (i.e. .treat/ .vacc etc.) by microcell.
					//// Further, this block assigns all microcells that are within this admin unit (and around this microcell) to be treated, using the flags set to avoid duplication.

					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
					//// **** //// **** //// **** //// **** TREATMENT
					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

					if ((Mcells[mcellnum].treat == TreatStat::Treated) && (ts >= Mcells[mcellnum].treat_end_time))
					{
						TreatFlag = 1;
						Mcells[mcellnum].treat = TreatStat::Untreated;
					}
					if ((Mcells[mcellnum].treat == TreatStat::ToBeTreated) && (ts >= Mcells[mcellnum].treat_start_time))
					{
						TreatFlag = 1;
						Mcells[mcellnum].treat				= TreatStat::Treated;
						Mcells[mcellnum].treat_trig		= 0;
						Mcells[mcellnum].treat_end_time	= t_TreatEnd;
						for (int i = 0; i < Mcells[mcellnum].n; i++)
						{
							int l = Mcells[mcellnum].members[i];
							if ((!HOST_TO_BE_TREATED(l)) && ((P.TreatPropRadial == 1) || (ranf_mt(tn) < P.TreatPropRadial)))
								DoProphNoDelay(l, ts, tn, 1);
						}
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.TreatCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : (int)P.TreatCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : (int)P.TreatCellIncThresh;
						f2 = (Mcells[mcellnum].treat_trig >= trig_thresh);
					}
					if ((t >= P.TreatTimeStart) && (Mcells[mcellnum].treat == TreatStat::Untreated) && (f2) && (P.TreatRadius2 > 0))
					{
						MicroCellPosition min = P.get_micro_cell_position_from_cell_index(mcellnum);
						Direction j = Direction::Right;
						int k = mcellnum;
						int maxx = 0;
						int i, m, l;
						i = m = f2 = 0;
						l = f3 = 1;
						if ((!P.TreatByAdminUnit) || (ad > 0))
						{
							int ad2 = ad / P.TreatAdminUnitDivisor;
							do
							{
								if (P.is_in_bounds(min))
								{
									if (P.TreatByAdminUnit)
										f4 = (AdUnits[Mcells[k].adunit].id / P.TreatAdminUnitDivisor == ad2);
									else
										f4 = ((r = dist2_mm(Mcells + mcellnum, Mcells + k)) < P.TreatRadius2);
									if (f4)
									{
										TreatFlag = f2 = 1;
										if ((Mcells[k].n > 0) && (Mcells[k].treat == TreatStat::Untreated))
										{
											Mcells[k].treat_start_time = t_TreatStart;
											Mcells[k].treat = TreatStat::ToBeTreated;
											maxx += Mcells[k].n;
										}
									}
								}
								min += j;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = rotate_left(j);
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == Direction::Up)
									{
										f3 = f2;
										f2 = 0;
									}
								}
								k = P.get_micro_cell_index_from_position(min);
							} while ((f3) && (maxx < P.TreatMaxCoursesPerCase));
						}
					}


					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
					//// **** //// **** //// **** //// **** VACCINATION
					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****


					//// vaccinates proportion VaccProp of people in microcell (or at least adds them to geovacc_queue).
					if ((Mcells[mcellnum].vacc == TreatStat::ToBeTreated) && (ts >= Mcells[mcellnum].vacc_start_time))
					{
						TreatFlag = 1;
						Mcells[mcellnum].vacc_trig = 0;
						//if(State.cumVG+P.NumThreads*Mcells[mcellnum].n<P.VaccMaxCourses) //changed to VG - commented this out for now, we'll add everyone to queues and deal with the number of doses available in the vaccination function
						{
							for (int i = 0; i < Mcells[mcellnum].n; i++)
							{
								int l = Mcells[mcellnum].members[i];
								//#pragma omp critical (state_cumV_daily) //added this
								if (((P.VaccProp == 1) || (ranf_mt(tn) < P.VaccProp)))
								{
									//add to the queue
									DoVaccNoDelay(l,ts);
								}
							}
							Mcells[mcellnum].vacc = TreatStat::Treated;
						}
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.VaccCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : (int)P.VaccCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : (int)P.VaccCellIncThresh;
						f2 = (Mcells[mcellnum].treat_trig >= trig_thresh);
					}
					if ((!P.DoMassVacc) && (P.VaccRadius2 > 0) && (t >= P.VaccTimeStartGeo) && (Mcells[mcellnum].vacc == TreatStat::Untreated) && (f2)) //changed from VaccTimeStart to VaccTimeStarGeo
					{
						MicroCellPosition min = P.get_micro_cell_position_from_cell_index(mcellnum);
						Direction j = Direction::Right;
						int k = mcellnum;
						int i, l, m;
						i = m = f2 = 0;
						l = f3 = 1;
						if ((!P.VaccByAdminUnit) || (ad > 0))
						{
							int ad2 = ad / P.VaccAdminUnitDivisor;
							do
							{
								if (P.is_in_bounds(min))
								{
									if (P.VaccByAdminUnit)
									{
										f4 = (AdUnits[Mcells[k].adunit].id / P.VaccAdminUnitDivisor == ad2);
										r = 1e20;
									}
									else
										f4 = ((r = dist2_mm(Mcells + mcellnum, Mcells + k)) < P.VaccRadius2);
									if (f4)
									{
										TreatFlag = f2 = 1;
										if (r < P.VaccMinRadius2)
											Mcells[k].vacc = TreatStat::DontTreatAgain;
										else if ((Mcells[k].n > 0) && (Mcells[k].vacc == TreatStat::Untreated))
										{
											Mcells[k].vacc_start_time = t_VacStart;
											Mcells[k].vacc = TreatStat::ToBeTreated;
										}
									}
								}
								min += j;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = rotate_left(j);
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == Direction::Up)
									{
										f3 = f2;
										f2 = 0;
									}
								}
								k = P.get_micro_cell_index_from_position(min);
							} while (f3);
						}
					}

					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
					//// **** //// **** //// **** //// **** PLACE CLOSURE
					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

					///// note that here f2 bool asks whether trigger lower than stop threshold. A few blocks down meaning changes to almost the opposite: asking whether trigger has exceeded threshold in order to close places for first time.
					if (P.DoGlobalTriggers)
						f2 = (global_trig < P.PlaceCloseCellIncStopThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.PlaceCloseCellIncStopThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncStopThresh;
						f2 = (State.trigDC_adunit[adi] < trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.PlaceCloseCellIncStopThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncStopThresh;
						f2 = (Mcells[mcellnum].treat_trig < trig_thresh);
					}
					if ((Mcells[mcellnum].placeclose == TreatStat::Treated) && ((f2) || (ts >= Mcells[mcellnum].place_end_time))) //// if place closure has started, the places in this microcell are closed, and either stop threshold has been reached or place_end_time has passed, go through block
					{
						TreatFlag = 1;
						if (P.DoPlaceCloseOnceOnly)
							Mcells[mcellnum].placeclose = TreatStat::DontTreatAgain;
						else
							Mcells[mcellnum].placeclose = TreatStat::Untreated;

						Mcells[mcellnum].place_end_time = ts;
						Mcells[mcellnum].place_trig = 0;
						if (f2)
						{
							for (int j2 = 0; j2 < P.PlaceTypeNum; j2++)
								if (j2 != P.HotelPlaceType)
									for (int i2 = 0; i2 < Mcells[mcellnum].np[j2]; i2++)
										DoPlaceOpen(j2, Mcells[mcellnum].places[j2][i2], ts);
						}
					}

					if ((P.DoPlaces) && (t >= P.PlaceCloseTimeStart) && (Mcells[mcellnum].placeclose == TreatStat::Untreated)) //// if doing places, time now is after policy has begun, but place hasn't closed yet.
					{
						///// note that here f2 bool asks whether trigger has exceeded threshold in order to close places for first time.A few blocks up meaning was almost the opposite: asking whether trigger lower than stop threshold.

						if (P.DoGlobalTriggers)
						{
							f2 = (global_trig >= P.PlaceCloseCellIncThresh);
						}
						else if (P.DoAdminTriggers)
						{
							int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
							f2 = (State.trigDC_adunit[adi] > trig_thresh);
						}
						else
						{
							int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
							f2 = (Mcells[mcellnum].treat_trig >= trig_thresh);
						}
						if (((P.PlaceCloseByAdminUnit) && (AdUnits[Mcells[mcellnum].adunit].place_close_trig < USHRT_MAX - 1)
							&& (((double)AdUnits[Mcells[mcellnum].adunit].place_close_trig) / ((double)AdUnits[Mcells[mcellnum].adunit].NP) > P.PlaceCloseAdunitPropThresh))
							|| ((!P.PlaceCloseByAdminUnit) && (f2)))
						{
							//							if(P.PlaceCloseByAdminUnit) AdUnits[Mcells[mcellnum].adunit].place_close_trig=USHRT_MAX-1; // This means schools only close once
							int interventionFlag; //added this as a way to help filter out when interventions start
							interventionFlag = 1;
							if ((P.DoInterventionDelaysByAdUnit)&&((t <= AdUnits[Mcells[mcellnum].adunit].PlaceCloseTimeStart) || (t >= (AdUnits[Mcells[mcellnum].adunit].PlaceCloseTimeStart + AdUnits[Mcells[mcellnum].adunit].PlaceCloseDuration))))
									interventionFlag = 0;

							if ((interventionFlag == 1) && ((!P.PlaceCloseByAdminUnit) || (ad > 0)))
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
									for (int j2 = 0; j2 < P.PlaceTypeNum; j2++)
										if (j2 != P.HotelPlaceType)
											for (int i2 = 0; i2 < Mcells[mcellnum].np[j2]; i2++)
												DoPlaceClose(j2, Mcells[mcellnum].places[j2][i2], ts, tn, 1);
								}
							}
						}
					}


					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
					//// **** //// **** //// **** //// **** MOVEMENT RESTRICTIONS
					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

					if ((Mcells[mcellnum].moverest == TreatStat::Treated) && (ts >= Mcells[mcellnum].move_end_time))
					{
						TreatFlag = 1;
						if (P.DoMoveRestrOnceOnly)
							Mcells[mcellnum].moverest = TreatStat::DontTreatAgain;
						else
							Mcells[mcellnum].moverest = TreatStat::Untreated;

					}
					if ((Mcells[mcellnum].moverest == TreatStat::ToBeTreated) && (ts >= Mcells[mcellnum].move_start_time))
					{
						TreatFlag = 1;
						Mcells[mcellnum].moverest = TreatStat::Treated;
						Mcells[mcellnum].move_trig = 0;
						Mcells[mcellnum].move_end_time = t_MoveRestrict_End;
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.MoveRestrCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
						f2 = (Mcells[mcellnum].treat_trig >= trig_thresh);
					}

					if ((t >= P.MoveRestrTimeStart) && (Mcells[mcellnum].moverest == TreatStat::Untreated) && (f2))
					{
						MicroCellPosition min = P.get_micro_cell_position_from_cell_index(mcellnum);
						Direction j = Direction::Right;
						int k = mcellnum;
						int i, l, m;
						i = m = f2 = 0;
						l = f3 = 1;
						if ((!P.MoveRestrByAdminUnit) || (ad > 0))
						{
							int ad2 = ad / P.MoveRestrAdminUnitDivisor;
							do
							{
								if (P.is_in_bounds(min))
								{
									if (P.MoveRestrByAdminUnit)
										f4 = (AdUnits[Mcells[k].adunit].id / P.MoveRestrAdminUnitDivisor == ad2);
									else
										f4 = ((r = dist2_mm(Mcells + mcellnum, Mcells + k)) < P.MoveRestrRadius2);
									if (f4)
									{
										TreatFlag = f2 = 1;
										if ((Mcells[k].n > 0) && (Mcells[k].moverest == TreatStat::Untreated))
										{
											Mcells[k].move_start_time = t_MoveRestrict_Start;
											Mcells[k].moverest = TreatStat::ToBeTreated;
										}
									}
								}
								min += j;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = rotate_left(j);
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == Direction::Up) { f3 = f2; f2 = 0; }
								}
								k = P.get_micro_cell_index_from_position(min);
							} while (f3);
						}
					}

					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
					//// **** //// **** //// **** //// **** SOCIAL DISTANCING
					//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****


					if (P.DoGlobalTriggers)
						f2 = (global_trig < P.SocDistCellIncStopThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.SocDistCellIncStopThresh)) / P.IncThreshPop)) : P.SocDistCellIncStopThresh;
						f2 = (State.trigDC_adunit[adi] < trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.SocDistCellIncStopThresh)) / P.IncThreshPop)) : P.SocDistCellIncStopThresh;
						f2 = (Mcells[mcellnum].treat_trig < trig_thresh);
					}

					//// if: policy of social distancing has started AND this microcell cell has been labelled to as undergoing social distancing, AND either trigger not reached (note definition of f2 changes in next few lines) or end time has passed.
					if ((t >= P.SocDistTimeStart) && (Mcells[mcellnum].socdist == TreatStat::Treated) && ((f2) || (ts >= Mcells[mcellnum].socdist_end_time)))
					{
						TreatFlag = 1;
						if (P.DoSocDistOnceOnly)
							Mcells[mcellnum].socdist = TreatStat::DontTreatAgain;
						else
							Mcells[mcellnum].socdist = TreatStat::Untreated;


						Mcells[mcellnum].socdist_trig = 0;	//// reset trigger
						Mcells[mcellnum].socdist_end_time = ts; //// record end time.
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.SocDistCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
						f2 = (State.trigDC_adunit[adi] >= trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
						f2 = (Mcells[mcellnum].treat_trig >= trig_thresh);
					}
					if ((t >= P.SocDistTimeStart) && (Mcells[mcellnum].socdist == TreatStat::Untreated) && (f2))
					{
						//some code to try and deal with intervention delays and durations by admin unit based on global triggers
						int interventionFlag; //added this as a way to help filter out when interventions start
						interventionFlag = 1;

						if (P.DoInterventionDelaysByAdUnit)
							if ((t <= AdUnits[Mcells[mcellnum].adunit].SocialDistanceTimeStart) ||
								(t >= (AdUnits[Mcells[mcellnum].adunit].SocialDistanceTimeStart + AdUnits[Mcells[mcellnum].adunit].SocialDistanceDuration))) //// i.e. if outside window of social distancing for this admin unit.
								interventionFlag = 0;

						if (interventionFlag == 1)
							if ((Mcells[mcellnum].n > 0) && (Mcells[mcellnum].socdist == TreatStat::Untreated)) //// if microcell populated and not currently undergoing social distancing
							{
								Mcells[mcellnum].socdist = TreatStat::Treated; //// update flag to denote that cell is undergoing social distancing
								Mcells[mcellnum].socdist_trig = 0; /// reset trigger
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

					if ((Mcells[mcellnum].keyworkerproph == 2) && (ts >= Mcells[mcellnum].keyworkerproph_end_time))
					{
						TreatFlag = 1;
						Mcells[mcellnum].keyworkerproph = P.DoKeyWorkerProphOnceOnly;
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.KeyWorkerProphCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						int trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[mcellnum].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
						f2 = (Mcells[mcellnum].treat_trig >= trig_thresh);
					}
					if ((P.DoPlaces) && (t >= P.KeyWorkerProphTimeStart) && (Mcells[mcellnum].keyworkerproph == 0) && (f2))
					{
						MicroCellPosition min = P.get_micro_cell_position_from_cell_index(mcellnum);
						Direction j = Direction::Right;
						int k = mcellnum;
						int i, l, m;
						i = m = f2 = 0;
						l = f3 = 1;
						do
						{
							if (P.is_in_bounds(min))
								if (dist2_mm(Mcells + mcellnum, Mcells + k) < P.KeyWorkerProphRadius2)
								{
									TreatFlag = f2 = 1;
									if ((Mcells[k].n > 0) && (Mcells[k].keyworkerproph == 0))
									{
										Mcells[k].keyworkerproph = 2;
										Mcells[k].keyworkerproph_trig = 0;
										Mcells[k].keyworkerproph_end_time = t_KeyWorkerPlaceClosure_End;
										for (int i2 = 0; i2 < Mcells[k].n; i2++)
										{
											int j2 = Mcells[k].members[i2];
											if ((Hosts[j2].keyworker) && (!HOST_TO_BE_TREATED(j2)))
												DoProphNoDelay(j2, ts, tn, nckwp);
										}
									}
								}
							min += j;
							m = (m + 1) % l;
							if (m == 0)
							{
								j = rotate_left(j);
								i = (i + 1) % 2;
								if (i == 0) l++;
								if (j == Direction::Up)
								{
									f3 = f2;
									f2 = 0;
								}
							}
							k = P.get_micro_cell_index_from_position(min);
						} while (f3);
					}

			} // End of bs loop, tn loop, and pragma.
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
			//				2..NUM_PLACE_TYPES+1 = within-class/work-group place based transmission
			//				NUM_PLACE_TYPES+2..2*NUM_PLACE_TYPES+1 = between-class/work-group place based transmission
			//				2*NUM_PLACE_TYPES+2 = "spatial" transmission (spatially local random mixing)
			// bits >4 store the generation of infection

			AddToInfectionQueue(tn, infectee_cell_index, infector_index, infectee_index, infect_type);
			positiveEntryAdded = true;
		}
	}
	return positiveEntryAdded;
}
