#pragma once

#include <Coordinates/Location.h>
#include "IndexList.h"

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
	Location loc;
	float* prop_traffic;
	IndexList* DestMcells, *DestPlaces;

	double distance_to(Airport* other) const;
	double distance_squared_to(Airport* other) const;
};

