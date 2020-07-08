#ifndef COVIDSIM_SETUPMODEL_H_INCLUDED_
#define COVIDSIM_SETUPMODEL_H_INCLUDED_

#include <string>

int ReadFitIter(std::string const&);
void ResetTimeSeries(void);
void InitTransmissionCoeffs(void);

/**
 * Initialize and set up model functions before running the simulation.
 *
 * @param density_file		Input population density file path. An empty string will cause
 * 							an even population distribution to be used. A `.txt` file will
 * 							be converted to binary and can be saved with `out_density_file`.
 * 							A binary file saved via `out_density_file` from a previous run
 * 							can be also used here.
 * @param out_density_file	Output population density file path. An empty string won't save
 * 							the binary form of the `density_file` contents.
 * @param load_network_file Population model file path to load from a previous run. An empty
 * 							string will cause a new population model to be generated.
 * @param save_network_file Population model file path to save for subsequent runs. An empty
 * 							string will prevent the model generated from being saved.
 * @param school_file		School file path. An empty string will prevent schools from being
 * 							assigned during population model generation.
 * @param reg_demog_file	Regional demography file path. An empty string will cause static
 * 							administrative units to be used for population age distribution
 * @param out_file_base		Output file path prefix
 */
void SetupModel(std::string const& density_file, std::string const& out_density_file, std::string const& load_network_file,
				std::string const& save_network_file, std::string const& school_file, std::string const& reg_demog_file,
				std::string const& out_file_base);

/**
 * Initialize the microcell and population arrays.
 *
 * @param density_file		Input population density file path. An empty string will cause
 * 							an even population distribution to be used. A `.txt` file will
 * 							be converted to binary and can be saved with `out_density_file`.
 * 							A binary file saved via `out_density_file` from a previous run
 * 							can be also used here.
 * @param out_density_file	Output population density file path. An empty string won't save
 * 							the binary form of the `density_file` contents.
 * @param school_file		School file path. An empty string will prevent schools from being
 * 							assigned during population model generation.
 * @param reg_demog_file	Regional demography file path. An empty string will cause static
 * 							administrative units to be used for population age distribution
 */
void SetupPopulation(std::string const& density_file, std::string const& out_density_file, std::string const& school_file, std::string const& reg_demog_file);

/**
 * Initialize airports by assigning them to microcells and mapping to hotels.
 */

void SetupAirports(void);

/**
 * Populate ages for each host household randomly using distribution.
 *
 * @param n					Number of people in household
 * @param pers				Person to assign
 * @param tn				Thread number
 * @param do_adunit_demog	Flag if regional demography file exists
 */
void AssignHouseholdAges(int n, int pers, int tn, bool do_adunit_demog);
void AssignPeopleToPlaces(void);
void StratifyPlaces(void);

/**
 * Load a population model from a previous simulation instead of generating one.
 *
 * @param load_network_file Population model file path to load from a previous run. An empty
 * 							string will cause a new population model to be generated.
 */
void LoadPeopleToPlaces(std::string const& load_network_file);

/**
 * Save a population model so that subsequent simulations do not have to generate it.
 *
 * @param save_network_file Population model file path to save for subsequent runs. An empty
 * 							string will prevent the model generated from being saved.
 */
void SavePeopleToPlaces(std::string const& save_network_file);

/**
 * Save the age distribution to a file.
 * 
 * @param output_file_base	Output file path prefix
 */
void SaveAgeDistrib(std::string const& output_file_base);
void UpdateProbs(int);

// network file format version; update this number when you make changes to the format of the
// network file, to ensure old/incompatible files are not loaded.
const int NETWORK_FILE_VERSION = 1;

struct BinFile
{
	double x, y, pop;
	int cnt, ad;
};

#endif // COVIDSIM_SETUPMODEL_H_INCLUDED_
