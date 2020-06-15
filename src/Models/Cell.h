#pragma once

/**
 * @brief Holds microcells.
 *
 * Keeps track of susceptible, latent and infected people (in addition to details like who
 * is vaccinated, treated etc.) Also contains data for the spatial gravity model for social
 * interactions (probability distributions).
*/
struct Cell
{
	int n, S, L, I, R, D, cumTC, S0, tot_treat, tot_vacc;
	int* members, *susceptible, *latent, *infected; //// pointers to people in cell. e.g. *susceptible identifies where the final susceptible member of cel is.
	int* InvCDF;
	float tot_prob, *cum_trans, *max_trans;
	short int CurInterv[MAX_INTERVENTION_TYPES];
};
