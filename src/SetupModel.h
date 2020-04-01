#pragma once

#ifndef SPATIALSIM_SETUPMODEL_H_INCLUDED_
#define SPATIALSIM_SETUPMODEL_H_INCLUDED_

void SetupModel(char*, char*, char*, char*);
void SetupPopulation(char*, char*, char*);
void SetupAirports(void);
void SetupRoads(void); //added new function to take care of roads: ggilani - 12/02/15
void AssignHouseholdAges(int, int, int);
void AssignPeopleToPlaces(void);
void StratifyPlaces(void);
void LoadPeopleToPlaces(char*);
void SavePeopleToPlaces(char*);
void SaveAgeDistrib(void);

// network file format version; update this number when you make changes to the format of the
// network file, to ensure old/incompatible files are not loaded.
#define NETWORK_FILE_VERSION (1)

typedef struct BIN_FILE {
	double x, y, pop;
	int cnt, ad;
} bin_file;

extern char OutFile[1024],OutDensFile[1024], RoadNetworkFile[1024];

#endif // SPATIALSIM_SETUPMODEL_H_INCLUDED_
