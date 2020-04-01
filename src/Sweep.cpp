#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
#include "SharedFuncs.h"
#include "Update.h"


void TravelReturnSweep(double t)
{
	int i, j, k, l, n, nr, ner, tn;

	if (floor(1 + t + P.TimeStep) != floor(1 + t))
	{
		nr = ner = 0;
		k = (int)floor(t);
		l = 1 + k % MAX_TRAVEL_TIME;
#pragma omp parallel for private(i,j,k,n,tn) reduction(+:nr,ner) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
		{
			for (j = tn; j < P.Nplace[HOTEL_PLACE_TYPE]; j += P.NumThreads)
			{
				n = Places[HOTEL_PLACE_TYPE][j].n;
				for (k = n - 1; k >= 0; k--)
				{
					i = Places[HOTEL_PLACE_TYPE][j].members[k];
					if (Hosts[i].Travelling == l)
					{
						n--;
						/*						if((n<0)||(Places[HOTEL_PLACE_TYPE][j].members[n]<0)||(Places[HOTEL_PLACE_TYPE][j].members[n]>=P.N))
													{fprintf(stderr,"### %i %i %i %i\n",j,k,n,Places[HOTEL_PLACE_TYPE][j].members[n]);ner++;}
												else if((k<0)||(k>n))
													{fprintf(stderr,"@ %i %i %i %i\n",j,k,n,Places[HOTEL_PLACE_TYPE][j].members[n]);ner++;}
												else
						*/
						if (k != n)
						{
							Places[HOTEL_PLACE_TYPE][j].members[k] = Places[HOTEL_PLACE_TYPE][j].members[n];
						}
						nr++;
						if (Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE] != j)
						{
							ner++; fprintf(stderr, "(%i %i) ", j, Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE]);
						}
						Hosts[i].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
						Hosts[i].Travelling = 0;
					}
				}
				Places[HOTEL_PLACE_TYPE][j].n = n;
			}
		}
		fprintf(stderr, " d=%i e=%i>", nr, ner);
	}
}

void TravelDepartSweep(double t)
{
	int c, i, i2, j, k, l, d, d2, m, n, f, f2, f3, mps, nld, nad, nsk, tn, bm, hp;
	double s, s2, nl;
	cell* ct;

	if (floor(1 + t - P.TimeStep) != floor(1 + t))
	{
		bm = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));
		mps = 2 * ((int)P.PlaceTypeMeanSize[HOTEL_PLACE_TYPE]) - P.NumThreads - 1;
		k = (int)floor(t);
		d = k % MAX_TRAVEL_TIME;
		nad = nld = nsk = 0;
#pragma omp parallel for private(i,i2,j,k,l,d2,m,n,s,f,f2,f3,tn,hp) reduction(+:nad,nsk) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (i = tn; i < P.Nairports; i += P.NumThreads)
				if ((Airports[i].total_traffic > 0) && (Airports[i].num_mcell > 0))
				{
					s = Airports[i].total_traffic;
					if ((t > P.AirportCloseTimeStart) && (t < P.AirportCloseTimeStart + P.AirportCloseTimeStartBase))
						s *= P.AirportCloseEffectiveness;
					n = (s > 0) ? ((int)ignpoi_mt((double)s, tn)) : 0;
					f3 = 0;
					j = 0;
					while (j < n)
					{
						s = ranf_mt(tn);
						l = Airports[i].Inv_DestMcells[(int)floor(s * 1024)];
						while (Airports[i].DestMcells[l].prob < s) l++;
						l = Airports[i].DestMcells[l].id;
						k = (int)(ranf_mt(tn) * ((double)Mcells[l].n));
						i2 = Mcells[l].members[k];
						if ((abs(Hosts[i2].inf) < InfStat_InfectiousAsymptomaticNotCase) && (Hosts[i2].inf != InfStat_Case))
						{
							d2 = HOST_AGE_GROUP(i2);
							if ((P.RelativeTravelRate[d2] == 1) || (ranf_mt(tn) < P.RelativeTravelRate[d2]))
							{
								f2 = 1;
#pragma omp critical
								{
									if (Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] == -1)
									{
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -2;
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
										if (dist2_raw(Airports[i].loc_x, Airports[i].loc_y, Airports[k].loc_x, Airports[k].loc_y) > P.MoveRestrRadius2)
										{
											if (ranf_mt(tn) > P.MoveRestrEffect)
											{
												f2 = 1;
												nsk++;
												j++;
#pragma omp critical
												Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
											}
										}
									}
									if (!f2)
									{
										f = 1;
										do
										{
											s = ranf_mt(tn);
											m = Airports[k].Inv_DestPlaces[(int)floor(s * 1024)];
											while (Airports[k].DestPlaces[m].prob < s) m++;
											l = Airports[k].DestPlaces[m].id;
#pragma omp critical
											{
												if ((hp = Places[HOTEL_PLACE_TYPE][l].n) < mps)
												{
													f = 0;
													Places[HOTEL_PLACE_TYPE][l].n++;
												}
											}
											if (!f)
											{
												f3 = 0;
												Places[HOTEL_PLACE_TYPE][l].members[hp] = i2;
												d2 = (d + P.InvJourneyDurationDistrib[(int)(ranf_mt(tn) * 1024.0)]) % MAX_TRAVEL_TIME;
												Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = l;
												Hosts[i2].Travelling = 1 + d2;
												nad++;
												j++;
											}
											f2++;
										} while ((f) && (f2 < 300));
										if (f)
										{
#pragma omp critical
											Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
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
		nl = ((double)P.PlaceTypeMeanSize[HOTEL_PLACE_TYPE]) * P.HotelPropLocal / P.MeanLocalJourneyTime;
		nsk = 0;
#pragma omp parallel for private(c,i,i2,j,l,d2,m,n,s,s2,ct,f,f2,f3,tn,hp) reduction(+:nld,nsk) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (i = tn; i < P.Nplace[HOTEL_PLACE_TYPE]; i += P.NumThreads)
			{
				c = ((int)(Places[HOTEL_PLACE_TYPE][i].loc_x / P.cwidth)) * P.nch + ((int)(Places[HOTEL_PLACE_TYPE][i].loc_y / P.cheight));
				n = (int)ignpoi_mt(nl * Cells[c].tot_prob, tn);
				if (Places[HOTEL_PLACE_TYPE][i].n + n > mps)
				{
					nsk += (Places[HOTEL_PLACE_TYPE][i].n + n - mps);
					n = mps - Places[HOTEL_PLACE_TYPE][i].n;
				}
				for (j = 0; j < n; j++)
				{
					do
					{
						f = 0;
						s = ranf_mt(tn);
						l = Cells[c].InvCDF[(int)floor(s * 1024)];
						while (Cells[c].cum_trans[l] < s) l++;
						ct = CellLookup[l];
						m = (int)(ranf_mt(tn) * ((double)ct->S0));
						if (m < (ct->S + ct->L))
						{
							i2 = ct->susceptible[m];
							d2 = HOST_AGE_GROUP(i2);
							f3 = 0;
							if ((Hosts[i2].Travelling == 0) && ((P.RelativeTravelRate[d2] == 1) || (ranf_mt(tn) < P.RelativeTravelRate[d2])))
							{
#pragma omp critical
								{if (Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] == -1) { Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -2; f3 = 1; }}
							}
							if (f3)
							{
								s2 = dist2_raw(Households[Hosts[i2].hh].loc_x, Households[Hosts[i2].hh].loc_y, Places[HOTEL_PLACE_TYPE][i].loc_x, Places[HOTEL_PLACE_TYPE][i].loc_y);
								f2 = 1;
								if ((bm) && (s2 > P.MoveRestrRadius2))
								{
									if (ranf_mt(tn) >= P.MoveRestrEffect)
									{
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
										nsk++;
										f2 = 0;
									}
								}
								if (f2)
								{
									s = numKernel(s2) / Cells[c].max_trans[l];
									if (ranf_mt(tn) >= s)
									{
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = -1;
										f = 1;
									}
									else
									{
										d2 = (d + P.InvLocalJourneyDurationDistrib[(int)(ranf_mt(tn) * 1024.0)]) % MAX_TRAVEL_TIME;
										hp = Places[HOTEL_PLACE_TYPE][i].n;
										Places[HOTEL_PLACE_TYPE][i].n++;
										Places[HOTEL_PLACE_TYPE][i].members[hp] = i2;
										Hosts[i2].Travelling = 1 + d2;
										nld++;
#pragma omp critical
										Hosts[i2].PlaceLinks[HOTEL_PLACE_TYPE] = i;
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
	//// This function loops over infected people, and decides whom to infect. Structure is 1) #pragma loop over all cells then 1a) infectious people, which chooses who they will infect, adds them to a queue
	//// Next 2) #pragma loop infects those people from queue (using DoInfect function). This is to avoid race conditions. 
	//// Loop 1a) calculates the force of infection exerted by each infected person on (and therefore number of new infections to) i) their house; ii) their place(s); iii) other spatial cells. 
	//// Each force of infection includes infectiousness and susceptibility components. 
	//// Infectiousness is (broadly) a function of 1 person (their age, treatment status, places, no. people in their household etc.) 
	//// Susceptibility is (broadly) a function of 2 people (a person's susceptbility TO ANOTHER PERSON / potential infector) 
	//// After loop 1a) over infectious people, spatial infections are doled out. 

	int i, j, k, l, m;
	int n; //// number of people you could potentially infect in your place group, then number of potential spatial infections doled out by cell on other cells. 
	int i2, b, i3, f, f2, tn, cq /*cell queue*/, bm, ci /*person index*/;
	double seasonality, sbeta, hbeta;
	//// various quantities of force of infection, including "infectiousness" and "susceptibility" components
	double s; // household FOI on fellow household member, then place susceptibility, then random number for spatial infections allocation* / ;
	double s2; // spatial infectiousness, then distance in spatial infections allocation
	double s3; // household, then place infectiousness
	double s4; // place infectiousness (copy of s3 as some code commented out
	double s5; //// total spatial infectiousness summed over all infectious people in cell. 
	double fp; //// false positive
	double s3_scaled, s4_scaled; //scaled up versions of place infectiousness for digital contact tracing in places: ggilani 10/03/20
	int n_unscaled; //unscaled value of number of contacts for digital contact tracing in places: ggilani 10/03/20
	cell* c, *ct;
	microcell* mi, *mt, *mp;
	unsigned short int ts;
	person* si;

	if (!P.DoSeasonality)	seasonality = 1.0;
	else					seasonality = P.Seasonality[((int)t) % DAYS_PER_YEAR];

	ts = (unsigned short int) (P.TimeStepsPerDay * t);
	fp = P.TimeStep / (1 - P.FalsePositiveRate);
	sbeta = seasonality * fp * P.LocalBeta;
	hbeta = (P.DoHouseholds) ? (seasonality * fp * P.HouseholdTrans) : 0;
	bm = ((P.DoBlanketMoveRestr) && (t >= P.MoveRestrTimeStart) && (t < P.MoveRestrTimeStart + P.MoveRestrDuration));

#pragma omp parallel for private(j,k,l,m,n,i2,b,i3,f,f2,s,s2,s3,s4,s5,c,ct,mi,mt,mp,cq,ci,si) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)
		for (b = tn; b < P.NCP; b += P.NumThreads) //// loop over (in parallel) all populated cells. Loop 1)
		{
			c = CellLookup[b];
			s5 = 0; ///// spatial infectiousness summed over all infectious people in loop below
			for (j = 0; j < c->I; j++) //// Loop over all infectious people c->I in cell c. Loop 1a)
			{
				//// get person index ci
				ci = c->infected[j];
				//// get person si from Hosts (array of people), using pointer arithmetic.
				si = Hosts + ci;

				//// Household FOI component
				if (hbeta > 0)
				{
					if ((Households[si->hh].nh > 1) && (!si->Travelling))
					{
						l = Households[si->hh].FirstPerson;
						m = l + Households[si->hh].nh;
						s3 = hbeta * CalcHouseInf(ci, ts);
						f = 0;
						for (i3 = l; (i3 < m) && (!f); i3++) //// loop over people in household
							for (i2 = 0; (i2 < P.PlaceTypeNum) && (!f); i2++) //// loop over place types
								if (Hosts[i3].PlaceLinks[i2] >= 0) //// if person in household has any sort of link to place type
								{
									f = PLACE_CLOSED(i2, Hosts[i3].PlaceLinks[i2]);
								}
						if (f) { s3 *= P.PlaceCloseHouseholdRelContact; }/* NumPCD++;}*/ //// if people in your household are absent from places, person si/ci is more infectious to them, as they spend more time at home. 
						for (i3 = l; i3 < m; i3++) //// loop over all people in household (note goes from l to m - 1)
							if ((Hosts[i3].inf == InfStat_Susceptible) && (!Hosts[i3].Travelling)) //// if people in household uninfected/susceptible and not travelling
							{
								s = s3 * CalcHouseSusc(i3, ts, ci, tn);		//// FOI ( = infectiousness x susceptibility) from person ci/si on fellow household member i3
								if (ranf_mt(tn) < s) //// if household member i3 will be infected...
								{
									cq = Hosts[i3].pcell % P.NumThreads;
									if ((StateT[tn].n_queue[cq] < P.InfQueuePeakLength)) //(Hosts[i3].infector==-1)&&
									{
										if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
											Hosts[i3].infector = Hosts[i3].infect_type = -1;
										else
										{
											Hosts[i3].infector = ci; //// assign person ci as infector of peron i3
											Hosts[i3].infect_type = 1 + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
										}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3; //// .... add i3 to queue.
									}
								}
							}
					}
				}
				//// Place FOI component
				if (P.DoPlaces)
				{
					if (!HOST_ABSENT(ci))
					{
						mi = Mcells + si->mcell;
						for (k = 0; k < P.PlaceTypeNum; k++) //// loop over all place types
						{
							l = si->PlaceLinks[k];
							if ((l >= 0) && (!PLACE_CLOSED(k, l))) //// l>=0 means if place type k is relevant to person si. (And obviously if place isn't closed). 
							{
								s3 = fp * seasonality * CalcPlaceInf(ci, k, ts);
								//if doing digital contact tracing, scale up infectiousness s3 - but keep track of original value to make sure we scale down properly
								if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1) && (Hosts[ci].detected == 1))
								{
									s3_scaled = s3 * P.ScalingFactorPlaceDigitalContacts;
								}
								mp = Mcells + Places[k][l].mcell;
								if (bm)
								{
									if ((dist2_raw(Households[si->hh].loc_x, Households[si->hh].loc_y,
										Places[k][l].loc_x, Places[k][l].loc_y) > P.MoveRestrRadius2))
										s3 *= P.MoveRestrEffect;
									if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1) && (Hosts[ci].detected == 1))
									{
										s3_scaled *= P.MoveRestrEffect;
									}
								}
								else if ((mi->moverest != mp->moverest) && ((mi->moverest == 2) || (mp->moverest == 2)))
								{
									s3 *= P.MoveRestrEffect;
									if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1) && (Hosts[ci].detected == 1))
									{
										s3_scaled *= P.MoveRestrEffect;
									}
								}

								if ((k != HOTEL_PLACE_TYPE) && (!si->Travelling))
								{
									i2 = (si->PlaceGroupLinks[k]);
									s4 = s3;
									//s4=s3*(1-P.PlaceTypePropBetweenGroupLinks[k]*P.PlaceTypeGroupSizeParam1[k]/((double) Places[k][l].n));


									if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration))
									{
										if ((Hosts[ci].digitalContactTracingUser) && (Hosts[ci].detected))
										{
											s4_scaled = s3_scaled;
											if (s4 < 0)
											{
												ERR_CRITICAL_FMT("@@@ %lg\n", s4);
											}
											else if (s4_scaled >= 1)	//// if place infectiousness above threshold, consider everyone in group a potential infectee...
												n = Places[k][l].group_size[i2];
											else				//// ... otherwise randomly sample (from binomial distribution) number of potential infectees in this place. 
												n = (int)ignbin_mt((long)Places[k][l].group_size[i2], s4_scaled, tn);
											if (n > 0) SampleWithoutReplacement(tn, n, Places[k][l].group_size[i2]); //// changes thread-specific SamplingQueue.

											//also calculate unscaled versions of n, to ensure than we rescale correctly later
											if (s4 >= 1)
											{
												n_unscaled = Places[k][l].group_size[i2];
											}
											else
											{
												n_unscaled = (int)ignbin_mt((long)Places[k][l].group_size[i2], s4, tn);
											}
										}
									}
									else
									{

										if (s4 < 0)
										{
											ERR_CRITICAL_FMT("@@@ %lg\n", s4);
										}
										else if (s4 >= 1)	//// if place infectiousness above threshold, consider everyone in group a potential infectee...
											n = Places[k][l].group_size[i2];
										else				//// ... otherwise randomly sample (from binomial distribution) number of potential infectees in this place. 
											n = (int)ignbin_mt((long)Places[k][l].group_size[i2], s4, tn);
										if (n > 0) SampleWithoutReplacement(tn, n, Places[k][l].group_size[i2]); //// changes thread-specific SamplingQueue.
									}



									for (m = 0; m < n; m++)
									{
										i3 = Places[k][l].members[Places[k][l].group_start[i2] + SamplingQueue[tn][m]];
										if ((Hosts[i3].inf == InfStat_Susceptible) && (!HOST_ABSENT(i3))) //// if person i3 uninfected and not absent.
										{
											mt = Mcells + Hosts[i3].mcell;
											ct = Cells + Hosts[i3].pcell;
											s = CalcPlaceSusc(i3, k, ts, ci, tn);
											//These are all contacts of the infector, so we need to add them to list of contacts (if not already listed), and we need to scale down susceptibility by the appropriate amount so we don't over infect
											if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration))
											{
												if ((Hosts[ci].detected) && (Hosts[ci].digitalContactTracingUser))
												{
													//scale down susceptibility so we don't over accept
													s *= (((double)n_unscaled) / ((double)n)); //scale down relative to the maximum number of contacts - can probably do this better
													//if infectee is also a user, add them as a contact
													if (Hosts[i3].digitalContactTracingUser)
													{
														//to try and avoid adding duplicate contacts to the list
														bool dctflag = true;
														int ii = 0;
														while ((ii < Hosts[ci].ncontacts) && dctflag)
														{
															if (Hosts[ci].contacts[ii] == i3)
															{
																dctflag = false;
															}
															ii++;
														}
														if (dctflag)
														{
															if (Hosts[ci].ncontacts >= MAX_CONTACTS)
															{
																ERR_CRITICAL("Exceeded maximum number of contacts per host!\n");
															}
															else
															{
																Hosts[ci].contacts[Hosts[ci].ncontacts] = i3;
																Hosts[ci].ncontacts++;
															}
														}
													}
												}
											}

											if (bm)
											{
												if ((dist2_raw(Households[Hosts[i3].hh].loc_x, Households[Hosts[i3].hh].loc_y,
													Places[k][l].loc_x, Places[k][l].loc_y) > P.MoveRestrRadius2))
													s *= P.MoveRestrEffect;
											}
											else if ((mt->moverest != mp->moverest) && ((mt->moverest == 2) || (mp->moverest == 2)))
												s *= P.MoveRestrEffect;
											if ((s == 1) || (ranf_mt(tn) < s))
											{
												cq = Hosts[i3].pcell % P.NumThreads;
												if ((StateT[tn].n_queue[cq] < P.InfQueuePeakLength)) //(Hosts[i3].infector==-1)&&
												{
													if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
														Hosts[i3].infector = Hosts[i3].infect_type = -1;
													else
													{
														Hosts[i3].infector = ci;
														Hosts[i3].infect_type = 2 + k + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
													}
													StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
												}
											}
										}
									}
								}
								if ((k == HOTEL_PLACE_TYPE) || (!si->Travelling))
								{
									s3 *= P.PlaceTypePropBetweenGroupLinks[k] * P.PlaceTypeGroupSizeParam1[k] / ((double)Places[k][l].n);
									if (s3 < 0)
									{
										ERR_CRITICAL_FMT("@@@ %lg\n", s3);
									}
									else if (s3 >= 1)
										n = Places[k][l].n;
									else
										n = (int)ignbin_mt((long)Places[k][l].n, s3, tn);
									if (n > 0) SampleWithoutReplacement(tn, n, Places[k][l].n);
									for (m = 0; m < n; m++)
									{
										i3 = Places[k][l].members[SamplingQueue[tn][m]];
										if ((Hosts[i3].inf == InfStat_Susceptible) && (!HOST_ABSENT(i3)))
										{
											mt = Mcells + Hosts[i3].mcell;
											ct = Cells + Hosts[i3].pcell;
											s = CalcPlaceSusc(i3, k, ts, ci, tn);


											//These are all contacts of the infector, so we need to add them to list of contacts (if not already listed), and we need to scale down susceptibility by the appropriate amount so we don't over infect
											if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].detected == 1) && (Hosts[ci].digitalContactTracingUser == 1))
											{

												//scale down susceptibility so we don't over accept
												s *= (((double)n_unscaled) / ((double)n)); //scale down relative to the maximum number of contacts - can probably do this better
												//if infectee is also a user, add them as a contact
												if (Hosts[i3].digitalContactTracingUser == 1)
												{
													//to try and avoid adding duplicate contacts to the list
													bool dctflag = true;
													int ii = 0;
													while ((ii < Hosts[ci].ncontacts) && dctflag)
													{
														if (Hosts[ci].contacts[ii] == i3)
														{
															dctflag = false;
														}
														ii++;
													}
													if (dctflag)
													{
														if (Hosts[ci].ncontacts >= MAX_CONTACTS)
														{
															ERR_CRITICAL("Exceeded maximum number of contacts per host!\n");
														}
														else
														{
															Hosts[ci].contacts[Hosts[ci].ncontacts] = i3;
															Hosts[ci].ncontacts++;
														}
													}

												}
											}

											if (bm)
											{
												if ((dist2_raw(Households[Hosts[i3].hh].loc_x, Households[Hosts[i3].hh].loc_y,
													Places[k][l].loc_x, Places[k][l].loc_y) > P.MoveRestrRadius2))
													s *= P.MoveRestrEffect;
											}
											else if ((mt->moverest != mp->moverest) && ((mt->moverest == 2) || (mp->moverest == 2)))
												s *= P.MoveRestrEffect;
											if ((s == 1) || (ranf_mt(tn) < s))
											{
												cq = Hosts[i3].pcell % P.NumThreads;
												if ((StateT[tn].n_queue[cq] < P.InfQueuePeakLength))//(Hosts[i3].infector==-1)&&
												{
													if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
														Hosts[i3].infector = Hosts[i3].infect_type = -1;
													else
													{
														Hosts[i3].infector = ci;
														Hosts[i3].infect_type = 2 + k + NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
													}
													StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
												}
											}
										}
									}
								}
							}
						}
					}
				}
				//// Spatial FOI component
				if (sbeta > 0) //// sum spatial infectiousness over all infected people, the infections from which are allocated after loop over infected people. 
				{
					if (si->Travelling) //// if host currently away from their cell, they cannot add to their cell's spatial infectiousness.
					{
						s2 = 0; f = 0;
					}
					else
					{
						s2 = CalcSpatialInf(ci, ts);
						//if do digital contact tracing, scale up spatial infectiousness of infectives who are using the app and will be detected
						if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart)
							&& (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ci].digitalContactTracingUser == 1) && (Hosts[ci].detected == 1))
						{
							//scale up infectiousness so that we pick more spatial contacts
							s2 *= P.ScalingFactorSpatialDigitalContacts;
						}
						f = 0;
						if (P.DoPlaces)
							for (i3 = 0; (i3 < P.PlaceTypeNum) && (!f); i3++)
								if (si->PlaceLinks[i3] >= 0) //// if person has a link to place of type i3...
									f = PLACE_CLOSED(i3, si->PlaceLinks[i3]); //// find out if that place of type i3 is closed.
					}
					if (f) //// if place is closed then adjust the spatial infectiousness (similar logic to household infectiousness: place closure affects spatial infectiousness _
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
			}
			//// Allocate spatial infections
			if (s5 > 0) //// if spatial infectiousness positive
			{
				n = (int)ignpoi_mt(s5 * sbeta * ((double)c->tot_prob), tn); //// number people this cell's population might infect elsewhere. poisson random number based on spatial infectiousness s5, sbeta (seasonality) and this cell's "probability" (guessing this is a function of its population and geographical size). 
				i2 = c->I;
				if (n > 0) //// this block normalises cumulative infectiousness cell_inf by person. s5 is the total cumulative spatial infectiousness. Reason is so that infector can be chosen using ranf_mt, which returns random number between 0 and 1. 
				{
					//// normalise by cumulative spatial infectiousness. 
					for (j = 0; j < i2 - 1; j++) StateT[tn].cell_inf[j] = (float)(StateT[tn].cell_inf[j] / s5);
					//// does same as the above loop just a slightly faster calculation. i.e. StateT[tn].cell_inf[i2 - 1] / s5 would equal 1 or -1 anyway. 
					StateT[tn].cell_inf[i2 - 1] = (StateT[tn].cell_inf[i2 - 1] < 0) ? -1.0f : 1.0f;
				}
				for (k = 0; k < n; k++)  //// loop over infections to dole out. roughly speaking, this determines which infectious person in cell c infects which person elsewhere. 
				{
					//// decide on infector ci/si from cell c.
					if (i2 == 1)	j = 0;
					else				//// roughly speaking, this determines which infectious person in cell c infects which person elsewhere
					{
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
					ci = c->infected[j];
					si = Hosts + ci;

					//// decide on infectee outside cell c. 
					do
					{
						//// chooses which cell person will infect
						s = ranf_mt(tn);
						l = c->InvCDF[(int)floor(s * 1024)];
						while (c->cum_trans[l] < s) l++;
						ct = CellLookup[l];

						///// pick random person m within susceptibles of cell ct (S0 initial number susceptibles within cell). 
						m = (int)(ranf_mt(tn) * ((double)ct->S0));
						i3 = ct->susceptible[m];
						s2 = dist2(Hosts + i3, Hosts + ci); /// calculate distance squared between this susceptible person and person ci/si identified earlier
						s = numKernel(s2) / c->max_trans[l]; //// acceptance probability 
						f2 = 0;
						if ((ranf_mt(tn) >= s) || (abs(Hosts[i3].inf) == InfStat_Dead)) //// if rejected, or infectee i3/m already dead, ensure do-while evaluated again (i.e. choose a new infectee). 
						{
							f2 = 1;
						}
						//add cross border effect - relative contact assumed to be less if infector and infectee live in different countries: ggilani 09/12/14
						//else if(((AdUnits[Mcells[Hosts[ci].mcell].adunit].id/P.CountryDivisor)!=(AdUnits[Mcells[Hosts[i3].mcell].adunit].id/P.CountryDivisor))&&(ranf_mt(tn)>=P.PropCrossBorderInf)) //checking to see if they are in the same country
						//{
						//	f2=1;
						//}
						else if (m < ct->S) //// if person m susceptible. if not, then block of code not evaluated, we finish this do-while without having allocated infectee, or anybody to infection queue. 
						{
							if ((!Hosts[i3].Travelling) && ((c != ct) || (Hosts[i3].hh != si->hh))) //// if potential infectee not travelling, is not part of cell c and doesn't share a household with infector. 
							{

								mi = Mcells + si->mcell;
								mt = Mcells + Hosts[i3].mcell;
								s = CalcSpatialSusc(i3, ts, ci, tn);

								//so this person is a contact - but might not be infected. if we are doing digital contact tracing, we want to add the person to the contacts list, if both are users
								if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[si->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration))
								{
									if ((Hosts[ci].detected) && (Hosts[ci].digitalContactTracingUser))
									{
										//scale down susceptibility so we don't over accept
										s /= P.ScalingFactorSpatialDigitalContacts;
										//if infectee is also a user, add them as a contact
										if (Hosts[i3].digitalContactTracingUser)
										{
											//to try and avoid adding duplicate contacts to the list
											bool dctflag = true;
											int ii = 0;
											while ((ii < Hosts[ci].ncontacts) && dctflag)
											{
												if (Hosts[ci].contacts[ii] == i3)
												{
													dctflag = false;
												}
												ii++;
											}
											if (dctflag)
											{
												if (Hosts[ci].ncontacts >= MAX_CONTACTS)
												{
													ERR_CRITICAL("Exceeded maximum number of contacts per host!\n");
												}
												else
												{
													Hosts[ci].contacts[Hosts[ci].ncontacts] = i3;
													Hosts[ci].ncontacts++;
												}
											}
										}
									}

								}

								if (bm)
								{
									if ((dist2_raw(Households[si->hh].loc_x, Households[si->hh].loc_y,
										Households[Hosts[i3].hh].loc_x, Households[Hosts[i3].hh].loc_y) > P.MoveRestrRadius2))
										s *= P.MoveRestrEffect;
								}
								else if ((mt->moverest != mi->moverest) && ((mt->moverest == 2) || (mi->moverest == 2)))
									s *= P.MoveRestrEffect;
								if (!f) //// if infector did not have place closed, loop over place types of infectee i3 to see if their places had closed. If they had, amend their susceptibility. 
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
									if ((Hosts[i3].inf == InfStat_Susceptible) && (StateT[tn].n_queue[cq] < P.InfQueuePeakLength)) //Hosts[i3].infector==-1
									{
										if ((P.FalsePositiveRate > 0) && (ranf_mt(tn) < P.FalsePositiveRate))
											Hosts[i3].infector = Hosts[i3].infect_type = -1;
										else
										{
											Hosts[i3].infector = ci;
											Hosts[i3].infect_type = 2 + 2 * NUM_PLACE_TYPES + INFECT_TYPE_MASK * (1 + si->infect_type / INFECT_TYPE_MASK);
										}
										StateT[tn].inf_queue[cq][StateT[tn].n_queue[cq]++] = i3;
									}
								}
							}
						}
					} while (f2);
				}
			}
		}


//#pragma omp parallel for private(i,k) schedule(static,1)
	for (j = 0; j < P.NumThreads; j++)
	{
		for (k = 0; k < P.NumThreads; k++)
		{
			for (i = 0; i < StateT[k].n_queue[j]; i++)
			{
				if (Hosts[StateT[k].inf_queue[j][i]].infect_type == -1) //// i.e. if host doesn't have an infector
					DoFalseCase(StateT[k].inf_queue[j][i], t, ts, j);
				else
					DoInfect(StateT[k].inf_queue[j][i], t, j, run);
			}
			StateT[k].n_queue[j] = 0;
		}
	}
}

void IncubRecoverySweep(double t, int run)
{
	int i, j, k, l, b, tn, ci;
	cell* c;
	person* si;
	unsigned short int ts; //// this timestep
	unsigned short int tc; //// time at which person becomes case (i.e. moves from infectious and asymptomatic to infectious and symptomatic). 

	ts = (unsigned short int) (P.TimeStepsPerDay * t);

	for (i = 0; i < P.NumHolidays; i++)
	{
		if ((t + P.TimeStep >= P.HolidayStartTime[i]+P.PreControlClusterIdHolOffset) && (t < P.HolidayStartTime[i]+ P.PreControlClusterIdHolOffset))
		{
#pragma omp parallel for private(j,k,l,b,tn) schedule(static,1)
			for (tn = 0; tn < P.NumThreads; tn++)
			{
				for (b = tn; b < P.N; b += P.NumThreads)
				{
					for (j = 0; j < P.PlaceTypeNum; j++)
					{
						l = Hosts[b].PlaceLinks[j];
						if (l >= 0)
						{
							if ((P.HolidayEffect[j] < 1) && ((P.HolidayEffect[j] == 0) || (ranf_mt(tn) >= P.HolidayEffect[j])))
							{
								k = (int)(P.HolidayStartTime[i] * P.TimeStepsPerDay);
								if (Places[j][l].close_start_time > k)	Places[j][l].close_start_time = (unsigned short)k;
								if (Hosts[b].absent_start_time > k)		Hosts[b].absent_start_time = (unsigned short)k;

								k = (int)((P.HolidayStartTime[i] + P.HolidayDuration[i]) * P.TimeStepsPerDay);
								if (Places[j][l].close_end_time < k)	Places[j][l].close_end_time = (unsigned short)k;
								if (Hosts[b].absent_stop_time < k)		Hosts[b].absent_stop_time = (unsigned short)k;
							}
						}
					}
				}
			}
		}
	}

#pragma omp parallel for private(j,b,c,tn,tc,ci,si) schedule(static,1)
	for (tn = 0; tn < P.NumThreads; tn++)	//// loop over threads
		for (b = tn; b < P.NCP; b += P.NumThreads)	//// loop/step over populated cells
		{
			c = CellLookup[b]; //// find (pointer-to) cell. 
			for (j = ((int)c->L - 1); j >= 0; j--) //// loop backwards over latently infected people, hence it starts from L - 1 and goes to zero. Runs backwards because of pointer swapping?
				if (ts >= Hosts[c->latent[j]].latent_time) //// if now after time at which person became infectious (latent_time a slight misnomer). 
					DoIncub(c->latent[j], ts, tn, run); //// move infected person from latently infected (L) to infectious (I), but not symptomatic
			//StateT[tn].n_queue[0] = StateT[tn].n_queue[1] = 0;
			for (j = c->I - 1; j >= 0; j--) ///// loop backwards over Infectious people. Runs backwards because of pointer swapping?
			{
				ci = c->infected[j];	//// person index
				si = Hosts + ci;		//// person
				tc = si->latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep)); //// time that person si/ci becomes case (symptomatic)...
				if ((P.DoSymptoms) && (ts == tc)) //// ... if now is that time...
					DoCase(ci, t, ts, tn);  //// ... change infectious (but asymptomatic) person to infectous and symptomatic. If doing severity, this contains DoMild and DoILI. 

				if (P.DoSeverity)
				{
					if (ts >= si->SARI_time)					DoSARI(ci, tn);	//// see if you can dispense with inequalities by initializing SARI_time, Critical_time etc. to USHRT_MAX
					if (ts >= si->Critical_time)				DoCritical(ci, tn);
					if (ts >= si->RecoveringFromCritical_time)	DoRecoveringFromCritical(ci, tn);
					if (ts >= si->recovery_time)
					{
						if (si->to_die)
							DoDeath_FromCriticalorSARI(ci, tn);
						else
							DoRecover_FromSeverity(ci, tn);
					}
				}

				//Adding code to assign recovery or death when leaving the infectious class: ggilani - 22/10/14
				if (ts >= si->recovery_time)
				{
					if (!si->to_die) //// if person si recovers and this timestep is after they've recovered
						DoRecover(ci, run, tn);
						//StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++] = ci; //// add them to end of 0th thread of inf queue. Don't get why 0 here.
					else if (si->to_die) /// if they die and this timestep is after they've died.
					{
						if (HOST_TREATED(ci) && (ranf_mt(tn) < P.TreatDeathDrop))
							DoRecover(ci, run, tn);
							//StateT[tn].inf_queue[0][StateT[tn].n_queue[0]++] = ci; //// add them to end of 0th thread of inf queue. Don't get why 0 here.
						else
							DoDeath(ci, run, tn);
							//StateT[tn].inf_queue[0][StateT[tn].n_queue[1]++] = ci;
					}
				}
			}

			//for (j = 0; j < StateT[tn].n_queue[0]; j++) DoRecover(StateT[tn].inf_queue[0][j], run, tn);
			//for (j = 0; j < StateT[tn].n_queue[1]; j++) DoDeath(StateT[tn].inf_queue[1][j], run, tn);
			//StateT[tn].n_queue[0] = StateT[tn].n_queue[1] = 0;
		}
}

void DigitalContactTracingSweep(double t)
{
	/*
	 * Function: Update DigitalContactTracing
	 *
	 * Purpose: to count the number of people in each admin unit who are being digitally contact traced each day
	 * Parameters: none
	 * Returns: none
	 *
	 * Author: ggilani, 10/03/20
	 */
	int i, j, k;
	unsigned short int ts;
	//int numCT;

	//find current time step
	ts = (unsigned short int) (P.TimeStepsPerDay * t);


	for (i = 0; i < P.NumAdunits; i++)
	{
		if ((t >= AdUnits[i].DigitalContactTracingTimeStart) && (t < AdUnits[i].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration))
		{
			//take care of people no longer requiring contact tracing first
			for (j = 0; j < AdUnits[i].ndct;)
			{
				if (ts >= (Hosts[AdUnits[i].dct[j]].dct_end_time))
				{
					//stop contact tracing this host
					Hosts[AdUnits[i].dct[j]].digitalContactTraced = 0;
					k = AdUnits[i].dct[j];
					AdUnits[i].dct[j] = AdUnits[i].dct[AdUnits[i].ndct - 1];
					AdUnits[i].dct[AdUnits[i].ndct - 1] = k;
					AdUnits[i].ndct--;
				}
				else
				{
					j++;
				}
			}

			//next transfer people who are waiting in the queue to be traced to be contact traced, providing that they have passed the start time
			for (j = 0; j < AdUnits[i].ndct_queue;)
			{

				if (ts >= (Hosts[AdUnits[i].dct_queue[j]].dct_start_time))
				{
					//if this host is already under isolation, we don't want to re-add them to the list, only update their end time
					if (Hosts[AdUnits[i].dct[AdUnits[i].ndct]].digitalContactTraced == 1)
					{
						//update end time of isolation, based on list position of host
						Hosts[AdUnits[i].dct[AdUnits[i].ndct]].dct_end_time = ts + (unsigned short int) (P.LengthDigitalContactIsolation * P.TimeStepsPerDay);
					}
					else
					{
						//move this contact from the queue to the main contact tracing list, and mark them as digitally contact traced, and set end time
						AdUnits[i].dct[AdUnits[i].ndct] = AdUnits[i].dct_queue[j];
						Hosts[AdUnits[i].dct[AdUnits[i].ndct]].digitalContactTraced = 1;
						Hosts[AdUnits[i].dct[AdUnits[i].ndct]].dct_end_time = Hosts[AdUnits[i].dct[AdUnits[i].ndct]].dct_start_time + (unsigned short int) (P.LengthDigitalContactIsolation * P.TimeStepsPerDay);
						//update number being traced in AdUnit
						AdUnits[i].ndct++;
						//increment state variables
						StateT[0].cumDCT_adunit[i]++; //NEED TO CHECK THESE - but storing everyone in queue 0 at the moment as this part doesn't rely on threads
						StateT[0].cumDCT++;
					}

					//now remove this case from the queue
					k = AdUnits[i].dct_queue[j];
					AdUnits[i].dct_queue[j] = AdUnits[i].dct_queue[AdUnits[i].ndct_queue - 1];
					AdUnits[i].dct_queue[AdUnits[i].ndct_queue - 1] = k;
					AdUnits[i].ndct_queue--;
				}
				else
				{
					j++;
				}
			}

			//finally add people from the State queues to the AdUnit queues and set the time at which they start being digitally contact traced.
			for (j = 0; j < P.NumThreads; j++)
			{
				for (k = 0; k < StateT[j].ndct_queue[i]; k++)
				{
					AdUnits[i].dct_queue[k + AdUnits[i].ndct_queue] = StateT[j].dct_queue[i][k];
					//set contact tracing start time
					Hosts[AdUnits[i].dct_queue[k + AdUnits[i].ndct_queue]].dct_start_time = ts + (unsigned short int) (P.DigitalContactTracingDelay * P.TimeStepsPerDay);
					AdUnits[i].ndct++;

				}
				AdUnits[i].ndct_queue += StateT[j].ndct_queue[i];
			}

		}
	}

	//reset State variable queues
	for (i = 0; i < P.NumAdunits; i++)
	{
		//AdUnits[i].nct_queue = 0;
		for (j = 0; j < P.NumThreads; j++)
		{
			StateT[j].ndct_queue[i] = 0;
		}
	}


}

int TreatSweep(double t)
{
	///// function loops over microcells to decide which cells are treated (either with treatment, vaccine, social distancing, movement restrictions etc.)

	int i, j, i2, j2, k, l, m, b;
	int rd, prd, mrd; /// rd = radius; prd = place close radium; mrd = movement restriction radius in terms of microcells. 
	int f, f1, f2, f3, f4; //// various fail conditions. Used for other things 
	int tn, bs, adi, ad, ad2;
	int minx, maxx, miny, trig_thresh, nckwp;
	int interventionFlag; //added this as a way to help filter out when interventions start

	//// time steps
	unsigned short int ts;		////  time-step now. 
	unsigned short int tstf;	////  time-step treatment finish
	unsigned short int tstb;	////  time-step treatment begin
	unsigned short int tsvb;	////  time-step vaccination begin
	unsigned short int tspf;	////  time-step place closure finish
	unsigned short int tsmb;	////  time-step movement restriction begin
	unsigned short int tsmf;	////  time-step movement restriction finish
	unsigned short int tssdf;	////  time-step social distancing finish (e.g. household quarantine)
	unsigned short int tskwpf;	////  time-step key worker place closure finish
	int global_trig;
	double r;

	ts = (unsigned short int) (P.TimeStepsPerDay * t);
	rd = (int)ceil(P.TreatRadius / P.mcwidth);
	prd = (int)ceil(P.PlaceCloseRadius / P.mcwidth);
	mrd = (int)ceil(P.MoveRestrRadius / P.mcwidth);
	f = f1 = 0;
	if (P.DoGlobalTriggers)
	{
		if (P.DoPerCapitaTriggers)
			global_trig = (int)floor(((double)State.trigDC) * P.GlobalIncThreshPop / ((double)P.N));
		else
			global_trig = State.trigDC;
	}
	else
		global_trig = 0;

	///// block loops over places and determines whom to prophylactically treat
	if ((P.DoPlaces) && (t >= P.TreatTimeStart) && (t < P.TreatTimeStart + P.TreatPlaceGeogDuration) && (State.cumT < P.TreatMaxCourses))
	{
		tstf = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatDelayMean + P.TreatProphCourseLength));
#pragma omp parallel for private(i,j,k,l,m,f) reduction(+:f1) schedule(static,1)
		for (i = 0; i < P.NumThreads; i++)
			for (j = 0; j < P.PlaceTypeNum; j++)
			{
				for (k = 0; k < StateT[i].np_queue[j]; k++)
				{
					l = StateT[i].p_queue[j][k];
					if (P.DoPlaceGroupTreat)
					{
						f = StateT[i].pg_queue[j][k];
						for (m = ((int)Places[j][l].group_start[f]); m < ((int)(Places[j][l].group_start[f] + Places[j][l].group_size[f])); m++)
						{
							/*							if((Places[j][l].members[m]<0)||(Places[j][l].members[m]>P.N-1))
															fprintf(stderr,"\n*** npq=%i gn=%i h=%i m=%i j=%i l=%i f=%i s=%i n=%i ***\n",
																StateT[i].np_queue[j],
																Places[j][l].n,
																Places[j][l].members[m],
																m,j,l,f,
																(int) Places[j][l].group_start[f],
																(int) Places[j][l].group_size[f]);
														else
							*/
							if ((!HOST_TO_BE_TREATED(Places[j][l].members[m])) && ((P.TreatPlaceTotalProp[j] == 1) || (ranf_mt(i) < P.TreatPlaceTotalProp[j])))
								DoProph(Places[j][l].members[m], ts, i);
						}
					}
					else
					{
						if ((Places[j][l].treat) && (!PLACE_TREATED(j, l)))
						{
							f1 = 1;
							Places[j][l].treat_end_time = tstf;
							for (m = 0; m < Places[j][l].n; m++)
								if (!HOST_TO_BE_TREATED(Places[j][l].members[m]))
								{
									if ((P.TreatPlaceTotalProp[j] == 1) || (ranf_mt(i) < P.TreatPlaceTotalProp[j]))
										DoProph(Places[j][l].members[m], ts, i);
								}
						}
						Places[j][l].treat = 0;
					}
				}
				StateT[i].np_queue[j] = 0;
			}
	}

	///// block vaccinates everyone in mass vaccination queue. Don't know why loop is done twice (although State.mvacc_cum is reset at end so will relate to that)
	if ((P.DoMassVacc) && (t >= P.VaccTimeStart))
		for (j = 0; j < 2; j++)
		{
			m = (int)P.VaccMaxCourses;
			if (m > State.n_mvacc) m = State.n_mvacc;
#pragma omp parallel for private(i) schedule(static,1000)
			for (i = State.mvacc_cum; i < m; i++)
				DoVacc(State.mvacc_queue[i], ts);
			State.mvacc_cum = m;
		}

	if ((t >= P.TreatTimeStart) || (t >= P.VaccTimeStartGeo) || (t >= P.PlaceCloseTimeStart) || (t >= P.MoveRestrTimeStart) || (t >= P.SocDistTimeStart) || (t >= P.KeyWorkerProphTimeStart)) //changed this to start time geo
	{
		tstf = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatProphCourseLength) - 1);
		tstb = (unsigned short int) (P.TimeStepsPerDay * (t + P.TreatDelayMean));
		tsvb = (unsigned short int) (P.TimeStepsPerDay * (t + P.VaccDelayMean));
		tspf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.PlaceCloseDelayMean + P.PlaceCloseDuration));
		tsmf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.MoveRestrDuration));
		tsmb = (unsigned short int) floor(P.TimeStepsPerDay * (t + P.MoveDelayMean));
		tssdf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.SocDistDuration));
		tskwpf = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.KeyWorkerProphRenewalDuration));
		nckwp = (int)ceil(P.KeyWorkerProphDuration / P.TreatProphCourseLength);

#pragma omp parallel for private(tn,i,i2,j,j2,k,l,m,b,bs,minx,maxx,miny,f2,f3,f4,trig_thresh,r,ad,ad2,adi,interventionFlag) reduction(+:f) schedule(static,1)
		for (tn = 0; tn < P.NumThreads; tn++)
			for (bs = tn; bs < P.NMCP; bs += P.NumThreads) //// loop over populated microcells
			{
				b = (int)(McellLookup[bs] - Mcells); //// microcell number
				adi = (P.DoAdUnits) ? Mcells[b].adunit : -1;
				ad = (P.DoAdUnits) ? AdUnits[adi].id : 0;

					//// Code block goes through various types of treatments/interventions (vaccination/movement restrictions etc.), 
					//// assesses whether various triggers (counts) are over a certain threshold, (specified in ReadParams)
					//// and then implements those treatments by setting various flags (i.e. .treat/ .vacc etc.) by microcell.
					//// Further, this block assigns all microcells that are within this admin unit (and around this microcell) to be treated, using the flags set to avoid duplication. 

					if ((Mcells[b].treat == 2) && (ts >= Mcells[b].treat_end_time))
					{
						f = 1;
						Mcells[b].treat = 0;
					}
					if ((Mcells[b].treat == 1) && (ts >= Mcells[b].treat_start_time))
					{
						f = 1;
						Mcells[b].treat = 2;
						Mcells[b].treat_trig = 0;
						Mcells[b].treat_end_time = tstf;
						for (i = 0; i < Mcells[b].n; i++)
						{
							l = Mcells[b].members[i];
							if ((!HOST_TO_BE_TREATED(l)) && ((P.TreatPropRadial == 1) || (ranf_mt(tn) < P.TreatPropRadial)))
								DoProphNoDelay(l, ts, tn, 1);
						}
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.TreatCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : (int)P.TreatCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.TreatCellIncThresh)) / P.IncThreshPop)) : (int)P.TreatCellIncThresh;
						f2 = (Mcells[b].treat_trig >= trig_thresh);
					}
					if ((t >= P.TreatTimeStart) && (Mcells[b].treat == 0) && (f2) && (P.TreatRadius2 > 0))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						maxx = 0;
						i = j = m = f2 = 0;
						l = f3 = 1;
						if ((!P.TreatByAdminUnit) || (ad > 0))
						{
							ad2 = ad / P.TreatAdminUnitDivisor;
							do
							{
								if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								{
									if (P.TreatByAdminUnit)
										f4 = (AdUnits[Mcells[k].adunit].id / P.TreatAdminUnitDivisor == ad2);
									else
										f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.TreatRadius2);
									if (f4)
									{
										f = f2 = 1;
										if ((Mcells[k].n > 0) && (Mcells[k].treat == 0))
										{
											Mcells[k].treat_start_time = tstb;
											Mcells[k].treat = 1;
											maxx += Mcells[k].n;
										}
									}
								}
								if (j == 0)
									minx = minx + 1;
								else if (j == 1)
									miny = miny - 1;
								else if (j == 2)
									minx = minx - 1;
								else if (j == 3)
									miny = miny + 1;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = (j + 1) % 4;
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == 1) { f3 = f2; f2 = 0; }
								}
								k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
							} while ((f3) && (maxx < P.TreatMaxCoursesPerCase));
						}
					}

					//// vaccinates proportion VaccProp of people in microcell (or at least adds them to geovacc_queue). 
					if ((Mcells[b].vacc == 1) && (ts >= Mcells[b].vacc_start_time))
					{
						f = 1;
						Mcells[b].vacc_trig = 0;
						//if(State.cumVG+P.NumThreads*Mcells[b].n<P.VaccMaxCourses) //changed to VG - commented this out for now, we'll add everyone to queues and deal with the number of doses available in the vaccination function
						{
							for (i = 0; i < Mcells[b].n; i++)
							{
								l = Mcells[b].members[i];
								//#pragma omp critical (state_cumV_daily) //added this
								if (((P.VaccProp == 1) || (ranf_mt(tn) < P.VaccProp)))
								{
									//add to the queue
									DoVaccNoDelay(l,ts);
								}
							}
							Mcells[b].vacc = 2;
						}
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.VaccCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : (int)P.VaccCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.VaccCellIncThresh)) / P.IncThreshPop)) : (int)P.VaccCellIncThresh;
						f2 = (Mcells[b].treat_trig >= trig_thresh);
					}
					if ((!P.DoMassVacc) && (P.VaccRadius2 > 0) && (t >= P.VaccTimeStartGeo) && (Mcells[b].vacc == 0) && (f2)) //changed from VaccTimeStart to VaccTimeStarGeo
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						i = j = m = f2 = 0;
						l = f3 = 1;
						if ((!P.VaccByAdminUnit) || (ad > 0))
						{
							ad2 = ad / P.VaccAdminUnitDivisor;
							do
							{
								if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								{
									if (P.VaccByAdminUnit)
									{
										f4 = (AdUnits[Mcells[k].adunit].id / P.VaccAdminUnitDivisor == ad2);
										r = 1e20;
									}
									else
										f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.VaccRadius2);
									if (f4)
									{
										f = f2 = 1;
										if (r < P.VaccMinRadius2)
											Mcells[k].vacc = 3;
										else if ((Mcells[k].n > 0) && (Mcells[k].vacc == 0))
										{
											Mcells[k].vacc_start_time = tsvb;
											Mcells[k].vacc = 1;
										}
									}
								}
								if (j == 0)
									minx = minx + 1;
								else if (j == 1)
									miny = miny - 1;
								else if (j == 2)
									minx = minx - 1;
								else if (j == 3)
									miny = miny + 1;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = (j + 1) % 4;
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == 1) { f3 = f2; f2 = 0; }
								}
								k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
							} while (f3);
						}
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig < P.PlaceCloseCellIncStopThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.PlaceCloseCellIncStopThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncStopThresh;
						f2 = (State.trigDC_adunit[adi] < trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.PlaceCloseCellIncStopThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncStopThresh;
						f2 = (Mcells[b].treat_trig < trig_thresh);
					}
					if ((t >= P.PlaceCloseTimeStart) && (Mcells[b].placeclose == 2) && ((f2) || (ts >= Mcells[b].place_end_time)))
					{
						f = 1;
						Mcells[b].placeclose = P.DoPlaceCloseOnceOnly;
						Mcells[b].place_end_time = ts;
						Mcells[b].place_trig = 0;
						if (f2)
						{
							for (j2 = 0; j2 < P.PlaceTypeNum; j2++)
								if (j2 != HOTEL_PLACE_TYPE)
									for (i2 = 0; i2 < Mcells[b].np[j2]; i2++)
										DoPlaceOpen(j2, Mcells[b].places[j2][i2], ts, tn);
						}
					}
					if ((P.DoPlaces) && (t >= P.PlaceCloseTimeStart) && (Mcells[b].placeclose == 0))
					{
						if (P.DoGlobalTriggers)
							f2 = (global_trig >= P.PlaceCloseCellIncThresh);
						else if (P.DoAdminTriggers)
						{
							trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
							f2 = (State.trigDC_adunit[adi] > trig_thresh);
						}
						else
						{
							trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.PlaceCloseCellIncThresh)) / P.IncThreshPop)) : P.PlaceCloseCellIncThresh;
							f2 = (Mcells[b].treat_trig >= trig_thresh);
						}
						if (((P.PlaceCloseByAdminUnit) && (AdUnits[Mcells[b].adunit].place_close_trig < USHRT_MAX - 1)
							&& (((double)AdUnits[Mcells[b].adunit].place_close_trig) / ((double)AdUnits[Mcells[b].adunit].NP) > P.PlaceCloseAdunitPropThresh))
							|| ((!P.PlaceCloseByAdminUnit) && (f2)))
						{
							//							if(P.PlaceCloseByAdminUnit) AdUnits[Mcells[b].adunit].place_close_trig=USHRT_MAX-1; // This means schools only close once
							interventionFlag = 1;
							if (P.DoInterventionDelaysByAdUnit)
								if ((t <= AdUnits[Mcells[b].adunit].PlaceCloseTimeStart) || (t >= (AdUnits[Mcells[b].adunit].PlaceCloseTimeStart + AdUnits[Mcells[b].adunit].PlaceCloseDuration)))
									interventionFlag = 0;

							if (interventionFlag == 1)
								if ((!P.PlaceCloseByAdminUnit) || (ad > 0))
								{
									ad2 = ad / P.PlaceCloseAdminUnitDivisor;
									if ((Mcells[b].n > 0) && (Mcells[b].placeclose == 0))
									{
										//if doing intervention delays and durations by admin unit based on global triggers
										if (P.DoInterventionDelaysByAdUnit)
											Mcells[b].place_end_time = (unsigned short int) ceil(P.TimeStepsPerDay * (t + P.PlaceCloseDelayMean + AdUnits[Mcells[b].adunit].PlaceCloseDuration));
										else
											Mcells[b].place_end_time = tspf;
										Mcells[b].place_trig = 0;
										Mcells[b].placeclose = 2;
										for (j2 = 0; j2 < P.PlaceTypeNum; j2++)
											if (j2 != HOTEL_PLACE_TYPE)
												for (i2 = 0; i2 < Mcells[b].np[j2]; i2++)
													DoPlaceClose(j2, Mcells[b].places[j2][i2], ts, tn, 1);
									}
								}
						}
					}
					if ((Mcells[b].moverest == 2) && (ts >= Mcells[b].move_end_time))
					{
						f = 1;
						Mcells[b].moverest = P.DoMoveRestrOnceOnly;
					}
					if ((Mcells[b].moverest == 1) && (ts >= Mcells[b].move_start_time))
					{
						f = 1;
						Mcells[b].moverest = 2;
						Mcells[b].move_trig = 0;
						Mcells[b].move_end_time = tsmf;
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.MoveRestrCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.MoveRestrCellIncThresh)) / P.IncThreshPop)) : P.MoveRestrCellIncThresh;
						f2 = (Mcells[b].treat_trig >= trig_thresh);
					}

					if ((t >= P.MoveRestrTimeStart) && (Mcells[b].moverest == 0) && (f2))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						i = j = m = f2 = 0;
						l = f3 = 1;
						if ((!P.MoveRestrByAdminUnit) || (ad > 0))
						{
							ad2 = ad / P.MoveRestrAdminUnitDivisor;
							do
							{
								if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								{
									if (P.MoveRestrByAdminUnit)
										f4 = (AdUnits[Mcells[k].adunit].id / P.MoveRestrAdminUnitDivisor == ad2);
									else
										f4 = ((r = dist2_mm(Mcells + b, Mcells + k)) < P.MoveRestrRadius2);
									if (f4)
									{
										f = f2 = 1;
										if ((Mcells[k].n > 0) && (Mcells[k].moverest == 0))
										{
											Mcells[k].move_start_time = tsmb;
											Mcells[k].moverest = 1;
										}
									}
								}
								if (j == 0)
									minx = minx + 1;
								else if (j == 1)
									miny = miny - 1;
								else if (j == 2)
									minx = minx - 1;
								else if (j == 3)
									miny = miny + 1;
								m = (m + 1) % l;
								if (m == 0)
								{
									j = (j + 1) % 4;
									i = (i + 1) % 2;
									if (i == 0) l++;
									if (j == 1) { f3 = f2; f2 = 0; }
								}
								k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
							} while (f3);
						}
					}

					if (P.DoGlobalTriggers)
						f2 = (global_trig < P.SocDistCellIncStopThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.SocDistCellIncStopThresh)) / P.IncThreshPop)) : P.SocDistCellIncStopThresh;
						f2 = (State.trigDC_adunit[adi] < trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.SocDistCellIncStopThresh)) / P.IncThreshPop)) : P.SocDistCellIncStopThresh;
						f2 = (Mcells[b].treat_trig < trig_thresh);
					}
					if ((t >= P.SocDistTimeStart) && (Mcells[b].socdist == 2) && ((f2) || (ts >= Mcells[b].socdist_end_time)))
					{
						f = 1;
						Mcells[b].socdist = P.DoSocDistOnceOnly;
						Mcells[b].socdist_trig = 0;
						Mcells[b].socdist_end_time = ts;
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.SocDistCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
						f2 = (State.trigDC_adunit[adi] >= trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.SocDistCellIncThresh)) / P.IncThreshPop)) : P.SocDistCellIncThresh;
						f2 = (Mcells[b].treat_trig >= trig_thresh);
					}
					if ((t >= P.SocDistTimeStart) && (Mcells[b].socdist == 0) && (f2))
					{
						//some code to try and deal with intervention delays and durations by admin unit based on global triggers
						interventionFlag = 1;

						if (P.DoInterventionDelaysByAdUnit)
							if ((t <= AdUnits[Mcells[b].adunit].SocialDistanceTimeStart) || (t >= (AdUnits[Mcells[b].adunit].SocialDistanceTimeStart + AdUnits[Mcells[b].adunit].SocialDistanceDuration)))
								interventionFlag = 0;

						if (interventionFlag == 1)
							if ((Mcells[b].n > 0) && (Mcells[b].socdist == 0))
							{
								Mcells[b].socdist = 2;
								Mcells[b].socdist_trig = 0;
								if (P.DoInterventionDelaysByAdUnit)
									Mcells[b].socdist_end_time = (unsigned short int) ceil(P.TimeStepsPerDay * (t + AdUnits[Mcells[b].adunit].SocialDistanceDuration));
								else
									Mcells[b].socdist_end_time = tssdf;
							}
					}
					if ((Mcells[b].keyworkerproph == 2) && (ts >= Mcells[b].keyworkerproph_end_time))
					{
						f = 1;
						Mcells[b].keyworkerproph = P.DoKeyWorkerProphOnceOnly;
					}
					if (P.DoGlobalTriggers)
						f2 = (global_trig >= P.KeyWorkerProphCellIncThresh);
					else if (P.DoAdminTriggers)
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(AdUnits[adi].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
						f2 = (State.trigDC_adunit[adi] > trig_thresh);
					}
					else
					{
						trig_thresh = (P.DoPerCapitaTriggers) ? ((int)ceil(((double)(Mcells[b].n * P.KeyWorkerProphCellIncThresh)) / P.IncThreshPop)) : P.KeyWorkerProphCellIncThresh;
						f2 = (Mcells[b].treat_trig >= trig_thresh);
					}
					if ((P.DoPlaces) && (t >= P.KeyWorkerProphTimeStart) && (Mcells[b].keyworkerproph == 0) && (f2))
					{
						minx = (b / P.nmch); miny = (b % P.nmch);
						k = b;
						i = j = m = f2 = 0;
						l = f3 = 1;
						do
						{
							if ((minx >= 0) && (minx < P.nmcw) && (miny >= 0) && (miny < P.nmch))
								if (dist2_mm(Mcells + b, Mcells + k) < P.KeyWorkerProphRadius2)
								{
									f = f2 = 1;
									if ((Mcells[k].n > 0) && (Mcells[k].keyworkerproph == 0))
									{
										Mcells[k].keyworkerproph = 2;
										Mcells[k].keyworkerproph_trig = 0;
										Mcells[k].keyworkerproph_end_time = tskwpf;
										for (i2 = 0; i2 < Mcells[k].n; i2++)
										{
											j2 = Mcells[k].members[i2];
											if ((Hosts[j2].keyworker) && (!HOST_TO_BE_TREATED(j2)))
												DoProphNoDelay(j2, ts, tn, nckwp);
										}
									}
								}
							if (j == 0)
								minx = minx + 1;
							else if (j == 1)
								miny = miny - 1;
							else if (j == 2)
								minx = minx - 1;
							else if (j == 3)
								miny = miny + 1;
							m = (m + 1) % l;
							if (m == 0)
							{
								j = (j + 1) % 4;
								i = (i + 1) % 2;
								if (i == 0) l++;
								if (j == 1) { f3 = f2; f2 = 0; }
							}
							k = ((minx + P.nmcw) % P.nmcw) * P.nmch + (miny + P.nmch) % P.nmch;
						} while (f3);
					}

			}
		for (i = 0; i < P.NumThreads; i++)
		{
			State.cumT += StateT[i].cumT;
			State.cumTP += StateT[i].cumTP;
			State.cumUT += StateT[i].cumUT;
			//State.cumV+=StateT[i].cumV;
			StateT[i].cumT = StateT[i].cumUT = StateT[i].cumTP = 0; //StateT[i].cumV=0;
		}
	}
	f += f1;


	return (f > 0);
}
