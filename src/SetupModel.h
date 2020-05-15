#ifndef COVIDSIM_SETUPMODEL_H_INCLUDED_
#define COVIDSIM_SETUPMODEL_H_INCLUDED_

void SetupModel(char*, char*, char*, char*);
void SetupPopulation(char*, char*, char*);
void SetupAirports(void);
void AssignHouseholdAges(int, int, int);
void AssignPeopleToPlaces(void);
void StratifyPlaces(void);
void LoadPeopleToPlaces(char*);
void SavePeopleToPlaces(char*);
void SaveAgeDistrib(void);
void UpdateProbs(int);

// network file format version; update this number when you make changes to the format of the
// network file, to ensure old/incompatible files are not loaded.
const int NETWORK_FILE_VERSION = 1;

struct BinFile
{
	double x, y, pop;
	int cnt, ad;
};

extern char OutFile[1024],OutDensFile[1024];

#endif // COVIDSIM_SETUPMODEL_H_INCLUDED_
