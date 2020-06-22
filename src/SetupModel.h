#ifndef COVIDSIM_SETUPMODEL_H_INCLUDED_
#define COVIDSIM_SETUPMODEL_H_INCLUDED_

#include <string>

void SetupModel(std::string const&, std::string const&, std::string const&, std::string const&, std::string const&, std::string const&, std::string const&);
void SetupPopulation(std::string const&, std::string const&, std::string const&, std::string const&);
void SetupAirports(void);
void AssignHouseholdAges(int, int, int, std::string const&);
void AssignPeopleToPlaces(void);
void StratifyPlaces(void);
void LoadPeopleToPlaces(std::string const&);
void SavePeopleToPlaces(std::string const&);
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

extern std::string OutFile;

#endif // COVIDSIM_SETUPMODEL_H_INCLUDED_
