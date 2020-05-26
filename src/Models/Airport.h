#pragma once

#include "IndexList.h"
#include "Coordinates/Vector2.h"
#include "Place.h"

/**
 * @brief Airport state.
 *
 * Not used for COVID-19 right now. Might be more relevant for USA and
 * other countries that have lots of internal flights. Slows the simulation.
 */
struct Airport
{
	int num_mcell, num_place, Inv_prop_traffic[129], Inv_DestMcells[1025], Inv_DestPlaces[1025];
	unsigned short int num_connected, *conn_airports;
	float total_traffic;
	Vector2<float> loc;
	float* prop_traffic;
	IndexList* DestMcells, *DestPlaces;

	double distance_to(const Airport &other) const;
	double distance_to_squared(const Airport &other) const;
	double distance_to(const Place &other) const;
	double distance_to_squared(const Place &other) const;
};

