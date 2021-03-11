/** \file  Params.cpp
 *  \brief Support for parsing parameter files
 */

#include "ReadParams.h"

const double ICDF_START = 100.0;


ParamMap Params::read_params_map(const char* file)
{
	std::ifstream fin(file);  // The file to be read
	std::string buf;       // Buffer for a line of text

	ParamMap param_map;   // Resulting map of (key, value)
	ParamIter iter;       // Iterator, for checking duplicate keys
	std::string key;      // Current key being processed
	std::string value;    // Current value being built
	bool skipping = true;
	while (std::getline(fin, buf))    // Each line. (buffer includes new-line)
	{
		buf.erase(0, buf.find_first_not_of(" \t\n\r"));     // Remove lead space
		buf.erase(buf.find_last_not_of(" \t\n\r") + 1);     // ... and trailing.

		if (skipping) 
		{   // Currently, we're not doing anything interesting.
			if ((buf.length() > 1) && (buf.compare(0, 1, "[") == 0)) // We found a key
			{
				skipping = false;
				key = buf.substr(1, buf.length() - 2);
			}
		}
		else
		{    // Not skipping...
			if ((buf.length() == 0) || (buf.compare(0, 1, "*") == 0) ||
				(buf.compare(0, 1, "^") == 0) || (buf.compare(0, 1, "=") == 0) || (buf.compare(0, 1, "[") == 0)) 
			{
				value.erase(0, value.find_first_not_of(" \t\n\r"));
				value.erase(value.find_last_not_of(" \t\n\r") + 1);

				// ERR if the key already exists.
				iter = param_map.find(key);
				if (iter != param_map.end()) {
					ERR_CRITICAL_FMT("Duplicate parameter values for %s\n(1):%s\n(2):%s\n",
						key.c_str(), iter->second.c_str(), value.c_str());
				}
				param_map.insert(ParamPair(key, value));
				value.clear();
				key.clear();
				if (buf.compare(0, 1, "[") == 0) 
				{
					key = buf.substr(1, buf.length() - 2);
				}
				else
				{
					skipping = true;
				}
			}
			else 
			{
				value.append(buf);
				value.append("\n");
			}
		}
	}

	if (!key.empty()) 
	{
		value.erase(0, value.find_first_not_of(" \t\n\r"));
		value.erase(value.find_last_not_of(" \t\n\r") + 1);

		iter = param_map.find(key);
		if (iter != param_map.end()) {
			ERR_CRITICAL_FMT("Duplicate parameter values for %s\n(1):%s\n(2):%s\n",
				key.c_str(), iter->second.c_str(), value.c_str());
		}

		param_map.insert(ParamPair(key, value));
	}

	fin.close();
	return param_map;
}

/*****************************************************************************/

std::string Params::clp_overwrite(std::string value, Param* P) {
	if (value.at(0) != '#')
	{
		return value;
	} 
	else 
	{
	  int clp_no = std::stoi(value.substr(1, std::string::npos));
		if ((clp_no < 0) || (clp_no > 99)) 
		{
			ERR_CRITICAL_FMT("CLP %d is out of bounds reading parameters\n", clp_no);
		}
		return std::to_string(P->clP[clp_no]);
	}
}

/*****************************************************************************/

std::string Params::lookup_param(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, Param* P)
{

	ParamIter iter = params.find(param_name);
	if (iter != params.end())
	{
		return Params::clp_overwrite(iter->second, P);
	}
	else {
		iter = fallback.find(param_name);
		if (iter != fallback.end()) {
			return Params::clp_overwrite(iter->second, P);
		}
		else if (base == fallback)
		{
			return "NULL";
		}
		else 
		{
			iter = base.find(param_name);
			if (iter != base.end())
			{
				return Params::clp_overwrite(iter->second, P);
			}
			else
			{
				return "NULL";
			}
		}
	}
}

std::string Params::lookup_param(ParamMap &fallback, ParamMap &params, std::string param_name, Param* P)
{
  return Params::lookup_param(fallback, fallback, params, param_name, P);
}

/*****************************************************************************/

bool Params::param_found(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name)
{

	ParamIter iter = params.find(param_name);
	if (iter != params.end())
	{
		return true;
	}
	else
	{
		iter = fallback.find(param_name);
		if (iter != fallback.end())
		{
			return true;
		}
		else if (base == fallback)
		{
			return false;
		}
		else
		{
			iter = base.find(param_name);
			if (iter != base.end())
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
}

bool Params::param_found(ParamMap &fallback, ParamMap &params, std::string param_name)
{
	return Params::param_found(fallback, fallback, params, param_name);
}

/*****************************************************************************************/

double Params::get_double(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, double default_value, bool err_on_missing, Param* P)
{
	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0)
	{
		std::string::size_type idx;
		return std::stod(str_value, &idx);
	}
	if (err_on_missing)
	{
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
	} else {
		return default_value;
	}
}

double Params::get_double(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, double default_value, Param* P)
{
	return Params::get_double(base, fallback, params, param_name, default_value, false, P);
}

double Params::get_double(ParamMap &fallback, ParamMap &params, std::string param_name, double default_value, Param* P)
{
	return Params::get_double(fallback, fallback, params, param_name, default_value, false, P);
}

double Params::req_double(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, Param* P)
{
	return Params::get_double(base, fallback, params, param_name, 0, true, P);
}

double Params::req_double(ParamMap &fallback, ParamMap &params, std::string param_name, Param* P)
{
	return Params::req_double(fallback, fallback, params, param_name, P);
}


/*****************************************************************************************/

int Params::get_int(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, int default_value, bool err_on_missing, Param* P)
{
	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);
	
	if (str_value.compare("NULL") != 0)
	{
		std::string::size_type idx;
		uint64_t lval = std::stoll(str_value, &idx);
		int result;
		if ((lval > (double) INT32_MAX) || (lval < (double) INT32_MIN))
		{
			Files::xfprintf_stderr("Warning: int overflow param %s (%s)\n", param_name.c_str(), str_value.c_str());
			result = (int) lval;
		}
		else
		{
			result = std::stoi(str_value, &idx);
		}
		return result;
	}
	if (err_on_missing)
	{
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
	}
	else
	{
		return default_value;
	}
}

int Params::get_int(ParamMap &fallback, ParamMap &params, std::string param_name, int default_value, Param* P)
{
	return Params::get_int(fallback, fallback, params, param_name, default_value, false, P);
}

int Params::get_int(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, int default_value, Param* P)
{
	return Params::get_int(base, fallback, params, param_name, default_value, false, P);
}

int Params::req_int(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, Param* P)
{
	return Params::get_int(base, fallback, params, param_name, 0, true, P);
}

int Params::req_int(ParamMap &fallback, ParamMap &params, std::string param_name, Param* P)
{
	return Params::req_int(fallback, fallback, params, param_name, P);
}

/***********************************************************************************/

void Params::get_double_vec(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, double* array, int expected, double default_value, int default_size, bool err_on_missing, Param* P, bool force_fail)
{

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if ((str_value.compare("NULL") != 0) && (!force_fail))
	{
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int index = 0;
		while (str_stream >> buffer)
		{
			if (!buffer.empty())
			{
			  if (index < expected) array[index] = std::stod(buffer, &idx);
			  index++;
			}
		}
		if (index != expected)
		{
			Files::xfprintf_stderr("Warning - Extra elements for %s (%d - only needed %d)\n", param_name.c_str(), index, expected);
		}
	}
	else
	{
		for (int i = 0; i < default_size; i++) array[i] = default_value;
		if (err_on_missing)
		{
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
		}
	}
}

void Params::get_double_vec(ParamMap &fallback, ParamMap &params, std::string param_name, double* array, int expected, double default_value, int default_size, Param* P)
{
	Params::get_double_vec(fallback, fallback, params, param_name, array, expected, default_value, default_size, false, P, false);
}

void Params::req_double_vec(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, double* array, int expected, Param* P)
{
	Params::get_double_vec(base, fallback, params, param_name, array, expected, 0, 0, true, P, false);
}

void Params::req_double_vec(ParamMap &fallback, ParamMap &params, std::string param_name, double* array, int expected, Param* P)
{
	Params::get_double_vec(fallback, fallback, params, param_name, array, expected, 0, 0, true, P, false);
}

void Params::get_double_vec_ff(bool force_fail, ParamMap &fallback, ParamMap &params, std::string param_name, double* array, int expected, double default_value, Param* P)
{
	Params::get_double_vec(fallback, fallback, params, param_name, array, expected, default_value, expected, false, P, force_fail);
}

/***********************************************************************************/

void Params::get_int_vec(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, int* array, int expected, int default_value, int default_size, bool err_on_missing, Param* P, bool force_fail)
{

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if ((str_value.compare("NULL") != 0) && (!force_fail))
	{
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int index = 0;
		while (str_stream >> buffer)
		{
			if (!buffer.empty())
			{
				uint64_t lval = std::stoll(buffer, &idx);
				int result;
				if ((lval > (double)INT32_MAX) || (lval < (double)INT32_MIN))
				{
					Files::xfprintf_stderr("Warning: int overflow param %s (%s)\n", param_name.c_str(), buffer.c_str());
					result = (int) lval;
				}
				else
				{
					result = std::stoi(buffer, &idx);
				}

				array[index] = result;
				index++;
			}
		}
		if (index != expected)
		{
			Files::xfprintf_stderr("Warning - Extra elements for %s (%d - only needed %d)\n", param_name.c_str(), index, expected);
		}
	}
	else
	{
		for (int i = 0; i < default_size; i++) array[i] = default_value;
		if (err_on_missing)
		{
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
		}
	}
}

void Params::get_int_vec(ParamMap &fallback, ParamMap &params, std::string param_name, int* array, int expected, int default_value, int default_size, bool err_on_missing, Param* P, bool force_fail)
{
	Params::get_int_vec(fallback, fallback, params, param_name, array, expected, default_value, default_size, err_on_missing, P, force_fail);
}

void Params::get_int_vec(ParamMap &fallback, ParamMap &params, std::string param_name, int* array, int expected, int default_value, int default_size, Param* P)
{
	Params::get_int_vec(fallback, fallback, params, param_name, array, expected, default_value, default_size, false, P, false);
}

void Params::req_int_vec(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, int* array, int expected, Param* P)
{
	Params::get_int_vec(base, fallback, params, param_name, array, expected, 0, 0, false, P, false);
}

void Params::req_int_vec(ParamMap &fallback, ParamMap &params, std::string param_name, int* array, int expected, Param* P)
{
	Params::get_int_vec(fallback, fallback, params, param_name, array, expected, 0, 0, false, P, false);
}

void Params::get_int_vec_ff(bool force_fail, ParamMap &fallback, ParamMap &params, std::string param_name, int* array, int expected, int default_value, Param* P)
{
	Params::get_int_vec(fallback, fallback, params, param_name, array, expected, default_value, expected, false, P, force_fail);
}

/***********************************************************************************/

int Params::req_string_vec(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, char** array, int expected, Param* P)
{

	std::string str_value = Params::lookup_param(fallback, params, param_name, P);
	if (str_value.compare("NULL") != 0)
	{
		std::stringstream str_stream(str_value);
		std::string buffer;
		int index = 0;
		while (str_stream >> buffer)
		{
			if (!buffer.empty())
			{
				array[index] = new char[buffer.length() + 1];
				strcpy(array[index], buffer.c_str());
				index++;
			}
		}
		return index;
	}
	else
	{
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
	}
}

int Params::req_string_vec(ParamMap &fallback, ParamMap &params, std::string param_name, char** array, int expected, Param* P)
{
	return Params::req_string_vec(fallback, fallback, params, param_name, array, expected, P);
}

/***********************************************************************************/

void Params::get_double_matrix(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name, double** array, int sizex, int sizey, double default_value, bool err_on_missing, Param* P)
{
	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0)
	{
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int x = 0;
		int y = 0;
		while (str_stream >> buffer)
		{
			if (!buffer.empty())
			{
				if ((y < sizey) && (x < sizex)) array[x][y] = std::stod(buffer, &idx);
				x++;
				if (x == sizex)
				{
					y++;
					x = 0;
				}
			}
		}
		if ((x >= sizex) || (y >= sizey))
		{
			Files::xfprintf_stderr("Warning: more data available for %s than was read.\n", param_name.c_str());
		}
	}
	else 
	{
		if (err_on_missing)
		{
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
		}
		else
		{
	    for (int x = 0; x < sizex; x++)
			{
		    for (int y = 0; y < sizey; y++)
				{
			    array[x][y] = default_value;
				}
			}
		}
	}
}

void Params::get_double_matrix(ParamMap &fallback, ParamMap &params, std::string param_name, double** array, int sizex, int sizey, double default_value, Param* P)
{
	Params::get_double_matrix(fallback, fallback, params, param_name, array, sizex, sizey, default_value, false, P);
}

void Params::get_inverse_cdf(ParamMap fallback, ParamMap params, const char* icdf_name, InverseCdf* inverseCdf, Param* P, double start_value)
{
	Params::get_double_vec(fallback, params, icdf_name, inverseCdf->get_values(), CDF_RES + 1, 0, CDF_RES + 1, P);
	if (!Params::param_found(fallback, params, icdf_name))
	{
		inverseCdf->set_neg_log(start_value);
	}
	inverseCdf->assign_exponent();
}


double** create_2d_double(int sizex, int sizey)
{
	double** arr = new double* [sizex];
	for (int i = 0; i < sizex; i++)
	{
		arr[i] = new double[sizey];
	}
	return arr;
}

void Params::alloc_params(Param* P)
{
	P->LocationInitialInfection = create_2d_double(MAX_NUM_SEED_LOCATIONS, 2);
	P->WAIFW_Matrix = create_2d_double(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
	P->WAIFW_Matrix_SpatialOnly = create_2d_double(NUM_AGE_GROUPS, NUM_AGE_GROUPS);
	P->SD_PlaceEffects_OverTime = create_2d_double(MAX_NUM_INTERVENTION_CHANGE_TIMES, NUM_PLACE_TYPES);
	P->Enhanced_SD_PlaceEffects_OverTime = create_2d_double(MAX_NUM_INTERVENTION_CHANGE_TIMES, NUM_PLACE_TYPES);
	P->HQ_PlaceEffects_OverTime = create_2d_double(MAX_NUM_INTERVENTION_CHANGE_TIMES, NUM_PLACE_TYPES);
	P->PC_PlaceEffects_OverTime = create_2d_double(MAX_NUM_INTERVENTION_CHANGE_TIMES, NUM_PLACE_TYPES);
	P->PC_PropAttending_OverTime = create_2d_double(MAX_NUM_INTERVENTION_CHANGE_TIMES, NUM_PLACE_TYPES);
	P->HouseholdSizeDistrib = create_2d_double(MAX_ADUNITS, MAX_HOUSEHOLD_SIZE);
	P->PropAgeGroup = create_2d_double(MAX_ADUNITS, NUM_AGE_GROUPS);
	P->PopByAdunit = create_2d_double(MAX_ADUNITS, 2);
	P->InvLifeExpecDist = create_2d_double(MAX_ADUNITS, 1001);
}
/**************************************************************************************************************/

void Params::airport_params(ParamMap adm_params, ParamMap pre_params, ParamMap params, Param* P)
{
	P->DoAirports = Params::get_int(params, pre_params, "Include air travel", 0, P);
	if (P->DoAirports == 0)  // Airports disabled => all places are not to do with airports, and we have no hotels
	{
		P->PlaceTypeNoAirNum = P->PlaceTypeNum;
		P->HotelPlaceType = P->PlaceTypeNum;
		return;
	}
	
	// When airports are activated we must have at least one airport place
	// // and a hotel type.
	P->PlaceTypeNoAirNum = Params::req_int(pre_params, adm_params, "Number of non-airport places", P);
	P->HotelPlaceType = Params::req_int(pre_params, adm_params, "Hotel place type", P);
	if (P->PlaceTypeNoAirNum >= P->PlaceTypeNum)
	{
		ERR_CRITICAL_FMT("[Number of non-airport places] parameter (%d) is greater than number of places (%d).\n", P->PlaceTypeNoAirNum, P->PlaceTypeNum);
	}
	if (P->HotelPlaceType < P->PlaceTypeNoAirNum || P->HotelPlaceType >= P->PlaceTypeNum)
	{
		ERR_CRITICAL_FMT("[Hotel place type] parameter (%d) not in the range [%d, %d)\n", P->HotelPlaceType, P->PlaceTypeNoAirNum, P->PlaceTypeNum);
	}

	P->AirportTrafficScale = Params::get_double(params, pre_params, "Scaling factor for input file to convert to daily traffic", 1.0, P);
	P->HotelPropLocal = Params::get_double(params, pre_params, "Proportion of hotel attendees who are local", 0, P);

	// If params are not specified, get_double_vec here fills with zeroes.

	Params::get_double_vec(params, pre_params, "Distribution of duration of air journeys", P->JourneyDurationDistrib, MAX_TRAVEL_TIME, 0, MAX_TRAVEL_TIME, P);
	if (!Params::param_found(params, pre_params, "Distribution of duration of air journeys"))
	{
		P->JourneyDurationDistrib[0] = 1;
	}
	Params::get_double_vec(params, pre_params, "Distribution of duration of local journeys", P->LocalJourneyDurationDistrib, MAX_TRAVEL_TIME, 0, MAX_TRAVEL_TIME, P);
	if (!Params::param_found(params, pre_params, "Distribution of duration of local journeys"))
	{
		P->LocalJourneyDurationDistrib[0] = 1;
	}
	P->MeanJourneyTime = 0;
	P->MeanLocalJourneyTime = 0;
	for (int i = 0; i < MAX_TRAVEL_TIME; i++)
	{
		P->MeanJourneyTime += ((double)(i)) * P->JourneyDurationDistrib[i];
		P->MeanLocalJourneyTime += ((double)(i)) * P->LocalJourneyDurationDistrib[i];
	}
	Files::xfprintf_stderr("Mean duration of local journeys = %lf days\n", P->MeanLocalJourneyTime);
	for (int i = 1; i < MAX_TRAVEL_TIME; i++)
	{
		P->JourneyDurationDistrib[i] += P->JourneyDurationDistrib[i - 1];
		P->LocalJourneyDurationDistrib[i] += P->LocalJourneyDurationDistrib[i - 1];
	}
	int j1 = 0;
	int j2 = 0;
	for (int i = 0; i <= 1024; i++)
	{
		double s = ((double) i) / 1024;
		while (P->JourneyDurationDistrib[j1] < s) j1++;
		P->InvJourneyDurationDistrib[i] = j1;
		while (P->LocalJourneyDurationDistrib[j2] < s) j2++;
		P->InvLocalJourneyDurationDistrib[i] = j2;
	}
}

/**************************************************************************************************************/

void Params::ReadParams(std::string const& ParamFile, std::string const& PreParamFile, std::string const& AdUnitFile, Param* P, AdminUnit* AdUnits)
{
	double s, t, AgeSuscScale;
	int i, j, k, f, nc, na;

	char* CountryNameBuf = new char[128 * MAX_COUNTRIES];
	char* AdunitListNamesBuf = new char[128 * MAX_ADUNITS];
	char** CountryNames = new char* [MAX_COUNTRIES];
	char** AdunitListNames = new char* [MAX_ADUNITS];

	AgeSuscScale = 1.0;
	ParamMap params = Params::read_params_map(ParamFile.c_str());
	ParamMap pre_params = Params::read_params_map(PreParamFile.c_str());
	ParamMap adm_params = Params::read_params_map(AdUnitFile.c_str());

	if (P->FitIter == 0)
	{
		for (i = 0; i < MAX_COUNTRIES; i++) {
			CountryNames[i] = CountryNameBuf + INT64_C(128) * i;
			CountryNames[i][0] = 0;
		}
		for (i = 0; i < MAX_ADUNITS; i++) {
			AdunitListNames[i] = AdunitListNamesBuf + INT64_C(128) * i;
			AdunitListNames[i][0] = 0;
		}
		for (i = 0; i < 100; i++) P->clP_copies[i] = 0;
		P->LongitudeCutLine = Params::get_double(params, adm_params, "Longitude cut line", -360.0, P);
		P->ModelTimeStep = Params::req_double(params, pre_params, "Update timestep", P);
		P->OutputTimeStep = Params::req_double(params, pre_params, "Sampling timestep", P);
		if (P->ModelTimeStep > P->OutputTimeStep) ERR_CRITICAL("Update step must be smaller than sampling step\n");
		t = ceil(P->OutputTimeStep / P->ModelTimeStep - 1e-6);
		P->NumModelTimeStepsPerOutputTimeStep = (int)t;
		P->ModelTimeStep = P->OutputTimeStep / t;
		P->TimeStepsPerDay = ceil(1.0 / P->ModelTimeStep - 1e-6);
		Files::xfprintf_stderr("Update step = %lf\nSampling step = %lf\nUpdates per sample=%i\nTimeStepsPerDay=%lf\n", P->ModelTimeStep, P->OutputTimeStep, P->NumModelTimeStepsPerOutputTimeStep, P->TimeStepsPerDay);
		P->SimulationDuration = Params::req_double(params, pre_params, "Sampling time", P);
		P->NumOutputTimeSteps = 1 + (int)ceil(P->SimulationDuration / P->OutputTimeStep);
		P->PopSize = Params::req_int(pre_params, adm_params, "Population size", P);
		if (P->NumRealisations == 0)
		{
			P->NumRealisations = Params::req_int(params, pre_params, "Number of realisations", P);
			P->NumNonExtinctRealisations = Params::get_int(params, pre_params, "Number of non-extinct realisations", P->NumRealisations, P);
		}
		else
		{
			P->NumNonExtinctRealisations = P->NumRealisations;
		}
		P->SmallEpidemicCases = Params::get_int(params, pre_params, "Maximum number of cases defining small outbreak", -1, P);

		P->NumCells = -1;
		P->NMCL = Params::req_int(params, pre_params, "Number of micro-cells per spatial cell width", P);
		//added parameter to reset seeds after every run
		P->ResetSeeds = Params::get_int(params, pre_params, "Reset seeds for every run", 0, P);
		if (P->ResetSeeds != 0)
		{
			P->KeepSameSeeds = Params::get_int(params, pre_params, "Keep same seeds for every run", 0, P); //added this to control which seeds are used: ggilani 27/11/19
		}
		P->ResetSeedsPostIntervention = Params::get_int(params, pre_params, "Reset seeds after intervention", 0, P);
		if (P->ResetSeedsPostIntervention != 0)
		{
			P->TimeToResetSeeds = Params::get_int(params, pre_params, "Time to reset seeds after intervention", 1000000, P);
		}
		P->DoHouseholds = Params::get_int(pre_params, adm_params, "Include households", 1, P);

		P->KernelLookup.size_ = Params::get_int(pre_params, adm_params, "Kernel resolution", 4000000, P);
		if (P->KernelLookup.size_ < 2000000)
		{
			ERR_CRITICAL_FMT("[Kernel resolution] needs to be at least 2000000 - not %d", P->KernelLookup.size_);
		}
		P->KernelLookup.expansion_factor_ = Params::get_int(pre_params, adm_params, "Kernel higher resolution factor", P->KernelLookup.size_ / 1600, P);
		if (P->KernelLookup.expansion_factor_ < 1 || P->KernelLookup.expansion_factor_ >= P->KernelLookup.size_)
		{
			ERR_CRITICAL_FMT("[Kernel higher resolution factor] needs to be in range [1, P->NKR = %d) - not %d", P->KernelLookup.size_, P->KernelLookup.expansion_factor_);
		}
	}
	P->OutputAge = Params::get_int(params, pre_params, "OutputAge", 1, P);
	P->OutputSeverity = Params::get_int(params, pre_params, "OutputSeverity", 1, P);
	P->OutputSeverityAdminUnit = Params::get_int(params, pre_params, "OutputSeverityAdminUnit", 1, P);
	P->OutputSeverityAge = Params::get_int(params, pre_params, "OutputSeverityAge", 1, P);
	P->OutputAdUnitAge = Params::get_int(params, pre_params, "OutputAdUnitAge", 0, P);
	P->OutputR0 = Params::get_int(params, pre_params, "OutputR0", 0, P);
	P->OutputControls = Params::get_int(params, pre_params, "OutputControls", 0, P);
	P->OutputCountry = Params::get_int(params, pre_params, "OutputCountry", 0, P);
	P->OutputAdUnitVar = Params::get_int(params, pre_params, "OutputAdUnitVar", 0, P);
	P->OutputHousehold = Params::get_int(params, pre_params, "OutputHousehold", 0, P);
	P->OutputInfType = Params::get_int(params, pre_params, "OutputInfType", 0, P);
	P->OutputNonSeverity = Params::get_int(params, pre_params, "OutputNonSeverity", 0, P);
	P->OutputNonSummaryResults = Params::get_int(params, pre_params, "OutputNonSummaryResults", 0, P);

	if (P->DoHouseholds != 0)
	{
		Params::req_double_vec(pre_params, adm_params, "Household size distribution", P->HouseholdSizeDistrib[0], MAX_HOUSEHOLD_SIZE, P);
		P->HouseholdTrans = Params::req_double(params, pre_params, "Household attack rate", P);
		P->HouseholdTransPow = Params::req_double(params, pre_params, "Household transmission denominator power", P);
		P->DoCorrectAgeDist = Params::get_int(pre_params, adm_params, "Correct age distribution after household allocation to exactly match specified demography", 0, P);
	}
	else
	{
		P->HouseholdTrans = 0.0;
		P->HouseholdTransPow = 1.0;
		P->HouseholdSizeDistrib[0][0] = 1.0;
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P->HouseholdSizeDistrib[0][i] = 0;
	}
	if (P->FitIter == 0)
	{
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P->HouseholdSizeDistrib[0][i] = P->HouseholdSizeDistrib[0][i] + P->HouseholdSizeDistrib[0][i - 1];
		P->HouseholdDenomLookup[0] = 1.0;
		for (i = 1; i < MAX_HOUSEHOLD_SIZE; i++)
			P->HouseholdDenomLookup[i] = 1 / pow(((double)(INT64_C(1) + i)), P->HouseholdTransPow);
		P->DoAdUnits = Params::get_int(pre_params, adm_params, "Include administrative units within countries", 1, P);
		P->CountryDivisor = Params::get_int(pre_params, adm_params, "Divisor for countries", 1, P);
		if (P->DoAdUnits != 0)
		{

			char** AdunitNames = (char**)Memory::xmalloc(3 * ADUNIT_LOOKUP_SIZE * sizeof(char*));
			char* AdunitNamesBuf = (char*)Memory::xmalloc(3 * ADUNIT_LOOKUP_SIZE * 360 * sizeof(char));

			for (i = 0; i < ADUNIT_LOOKUP_SIZE; i++)
			{
				P->AdunitLevel1Lookup[i] = -1;
				AdunitNames[3 * i] = AdunitNamesBuf + INT64_C(3) * i * 360;
				AdunitNames[3 * i + 1] = AdunitNamesBuf + INT64_C(3) * i * 360 + 60;
				AdunitNames[3 * i + 2] = AdunitNamesBuf + INT64_C(3) * i * 360 + 160;
			}
			P->AdunitLevel1Divisor = Params::get_int(pre_params, adm_params, "Divisor for level 1 administrative units", 1, P);
			P->AdunitLevel1Mask = Params::get_int(pre_params, adm_params, "Mask for level 1 administrative units", 1000000000, P);
			na = Params::req_string_vec(pre_params, adm_params, "Codes and country/province names for admin units", AdunitNames, 3 * ADUNIT_LOOKUP_SIZE, P) / 3;
			nc = Params::get_int(pre_params, adm_params, "Number of countries to include", 0, P);
			if ((na > 0) && (nc > 0))
			{
				P->DoAdunitBoundaries = (nc > 0);
				nc = abs(nc);
				Params::req_string_vec(pre_params, adm_params, "List of names of countries to include", CountryNames, nc, P);
				P->NumAdunits = 0;
				for (i = 0; i < na; i++)
					for (j = 0; j < nc; j++)
						if ((AdunitNames[3 * i + 1][0]) && (!strcmp(AdunitNames[3 * i + 1], CountryNames[j])) && (atoi(AdunitNames[3 * i]) != 0))
						{
							AdUnits[P->NumAdunits].id = atoi(AdunitNames[3 * i]);
							P->AdunitLevel1Lookup[(AdUnits[P->NumAdunits].id % P->AdunitLevel1Mask) / P->AdunitLevel1Divisor] = P->NumAdunits;
							if (strlen(AdunitNames[3 * i + 1]) < 100) strcpy(AdUnits[P->NumAdunits].cnt_name, AdunitNames[3 * i + 1]);
							if (strlen(AdunitNames[3 * i + 2]) < 200) strcpy(AdUnits[P->NumAdunits].ad_name, AdunitNames[3 * i + 2]);
							//	Files::xfprintf_stderr("%i %s %s ## ",AdUnits[P->NumAdunits].id,AdUnits[P->NumAdunits].cnt_name,AdUnits[P->NumAdunits].ad_name);
							P->NumAdunits++;
						}
			}
			else
			{
				P->NumAdunits = Params::get_int(pre_params, adm_params, "Number of level 1 administrative units to include", 0, P);
				if (P->NumAdunits > 0)
				{
					P->DoAdunitBoundaries = 1;
					if (P->NumAdunits > MAX_ADUNITS) ERR_CRITICAL("MAX_ADUNITS too small.\n");
					Params::req_string_vec(pre_params, adm_params, "List of level 1 administrative units to include", AdunitListNames, P->NumAdunits, P);
					na = P->NumAdunits;
					for (i = 0; i < P->NumAdunits; i++)
					{
						f = 0;
						if (na > 0)
						{
							for (j = 0; (j < na) && (!f); j++) f = (!strcmp(AdunitNames[3 * j + 2], AdunitListNames[i]));
							if (f) k = atoi(AdunitNames[3 * (j - 1)]);
						}
						if ((na == 0) || (!f)) k = atoi(AdunitListNames[i]);
						AdUnits[i].id = k;
						P->AdunitLevel1Lookup[(k % P->AdunitLevel1Mask) / P->AdunitLevel1Divisor] = i;
						for (j = 0; j < na; j++)
							if (atoi(AdunitNames[3 * j]) == k)
							{
								if (strlen(AdunitNames[3 * j + 1]) < 100) strcpy(AdUnits[i].cnt_name, AdunitNames[3 * j + 1]);
								if (strlen(AdunitNames[3 * j + 2]) < 200) strcpy(AdUnits[i].ad_name, AdunitNames[3 * j + 2]);
								j = na;
							}
					}
				}
				else
					P->DoAdunitBoundaries = 0;
			}
			Memory::xfree(AdunitNames);
			Memory::xfree(AdunitNamesBuf);

			P->DoAdunitOutput = Params::get_int(params, pre_params, "Output incidence by administrative unit", 0, P);
			P->DoAdunitBoundaryOutput = Params::get_int(pre_params, adm_params, "Draw administrative unit boundaries on maps", 0, P);
			P->DoCorrectAdunitPop = Params::get_int(pre_params, adm_params, "Correct administrative unit populations", 0, P);
			P->DoSpecifyPop = Params::get_int(pre_params, adm_params, "Fix population size at specified value", 0, P);
			Files::xfprintf_stderr("Using %i administrative units\n", P->NumAdunits);
			P->AdunitBitmapDivisor = Params::get_int(pre_params, adm_params, "Divisor for administrative unit codes for boundary plotting on bitmaps", 1, P);
			P->DoOutputPlaceDistForOneAdunit = Params::get_int(params, pre_params, "Only output household to place distance distribution for one administrative unit", 0, P);
			if (P->DoOutputPlaceDistForOneAdunit != 0)
			{
				P->DoOutputPlaceDistForOneAdunit = Params::get_int(params, pre_params, "Administrative unit for which household to place distance distribution to be output", 0, P);
			}
		}
		else
		{
			P->DoAdunitBoundaries = P->DoAdunitBoundaryOutput = P->DoAdunitOutput = P->DoCorrectAdunitPop = P->DoSpecifyPop = 0;
			P->AdunitLevel1Divisor = 1; P->AdunitLevel1Mask = 1000000000;
			P->AdunitBitmapDivisor = P->AdunitLevel1Divisor;
		}
	}

	P->DoAge = Params::get_int(pre_params, adm_params, "Include age", 1, P);
	if (P->DoAge == 0)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++) {
			P->PropAgeGroup[0][i] = 1.0 / NUM_AGE_GROUPS;
			P->InitialImmunity[i] = 0;
			P->AgeInfectiousness[i] = P->AgeSusceptibility[i] = 1;
			P->RelativeSpatialContact[i] = P->RelativeTravelRate[i] = 1.0;
		}
	}
	else
	{

		P->DoPartialImmunity = Params::get_int(params, pre_params, "Initial immunity acts as partial immunity", 1, P);
		if ((P->DoHouseholds != 0) && (P->DoPartialImmunity == 0))
		{
			P->DoWholeHouseholdImmunity = Params::get_int(params, pre_params, "Initial immunity applied to all household members", 0, P);
		}
		else
			P->DoWholeHouseholdImmunity = 0;
		Params::get_double_vec(params, pre_params, "Initial immunity profile by age", P->InitialImmunity, NUM_AGE_GROUPS, 0, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "Relative spatial contact rates by age", P->RelativeSpatialContact, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, P);
		if (Params::get_int(params, pre_params, "Apply spatial contact rates by age to susceptibles as well as infecteds", 0, P) != 0)
			for (i = 0; i < NUM_AGE_GROUPS; i++)	P->RelativeSpatialContactSusc[i] = P->RelativeSpatialContact[i];
		else
			for (i = 0; i < NUM_AGE_GROUPS; i++)	P->RelativeSpatialContactSusc[i] = 1.0;
		Params::get_double_vec(params, pre_params, "Age-dependent infectiousness", P->AgeInfectiousness, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "Age-dependent susceptibility", P->AgeSusceptibility, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, P);
		Params::req_double_vec(pre_params, adm_params, "Age distribution of population", P->PropAgeGroup[0], NUM_AGE_GROUPS, P);
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			t += P->PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P->PropAgeGroup[0][i] /= t;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			if (P->AgeSusceptibility[i] > t) t = P->AgeSusceptibility[i]; //peak susc has to be 1
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P->AgeSusceptibility[i] /= t;
		AgeSuscScale = t;
		if (P->DoHouseholds) P->HouseholdTrans *= AgeSuscScale;
		Params::get_double_vec(pre_params, adm_params, "Relative travel rates by age", P->RelativeTravelRate, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, P);

		//if (!GetInputParameter2(params, pre_params, "WAIFW matrix", "%lf", (void*)P->WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0))
		if (!Params::param_found(params, pre_params, "WAIFW matrix"))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					P->WAIFW_Matrix[i][j] = 1.0;
		}
		else
		{
			Params::get_double_matrix(params, pre_params, "WAIFW matrix", P->WAIFW_Matrix, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 1.0, P);

			/* WAIFW matrix needs to be scaled to have max value of 1.
			1st index of matrix specifies host being infected, second the infector.
			Overall age variation in infectiousness/contact rates/susceptibility should be factored
			out of WAIFW_matrix and put in Age dep infectiousness/susceptibility for efficiency. */
			t = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					if (P->WAIFW_Matrix[i][j] > t) t = P->WAIFW_Matrix[i][j];
			if (t > 0)
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P->WAIFW_Matrix[i][j] /= t;
			}
			else
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P->WAIFW_Matrix[i][j] = 1.0;
			}
		}
		if (!Params::param_found(params, pre_params, "WAIFW matrix spatial infections only"))
		{
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					P->WAIFW_Matrix_SpatialOnly[i][j] = 1.0;

			P->Got_WAIFW_Matrix_Spatial = 0;
		}
		else
		{
			Params::get_double_matrix(params, pre_params, "WAIFW matrix spatial infections only", P->WAIFW_Matrix_SpatialOnly, NUM_AGE_GROUPS, NUM_AGE_GROUPS, 0, P);
			P->Got_WAIFW_Matrix_Spatial = 1;
			/* WAIFW matrix needs to be scaled to have max value of 1.
			1st index of matrix specifies host being infected, second the infector.
			Overall age variation in infectiousness/contact rates/susceptibility should be factored
			out of WAIFW_matrix and put in Age dep infectiousness/susceptibility for efficiency. */

			double Maximum = 0;
			for (i = 0; i < NUM_AGE_GROUPS; i++)
				for (j = 0; j < NUM_AGE_GROUPS; j++)
					if (P->WAIFW_Matrix_SpatialOnly[i][j] > Maximum) Maximum = P->WAIFW_Matrix_SpatialOnly[i][j];
			if (Maximum > 0)
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P->WAIFW_Matrix_SpatialOnly[i][j] /= Maximum;
			}
			else
			{
				for (i = 0; i < NUM_AGE_GROUPS; i++)
					for (j = 0; j < NUM_AGE_GROUPS; j++)
						P->WAIFW_Matrix_SpatialOnly[i][j] = 1.0;
			}

		}

		P->DoDeath = 0;
		t = 0;
		for (i = 0; i < NUM_AGE_GROUPS; i++)	t += P->AgeInfectiousness[i] * P->PropAgeGroup[0][i];
		for (i = 0; i < NUM_AGE_GROUPS; i++)	P->AgeInfectiousness[i] /= t;
	}

	if (P->FitIter == 0)
	{
		P->DoSpatial = Params::get_int(pre_params, adm_params, "Include spatial transmission", 1, P);
		P->MoveKernel.type_ = Params::req_int(pre_params, adm_params, "Kernel type", P);
		P->MoveKernel.scale_ = Params::req_double(pre_params, adm_params, "Kernel scale", P);
		if (P->KernelOffsetScale != 1)
		{
			P->MoveKernel.scale_ *= P->KernelOffsetScale;
		}
		P->MoveKernel.p3_ = Params::get_double(pre_params, adm_params, "Kernel 3rd param", 0, P);
		P->MoveKernel.p4_ = Params::get_double(pre_params, adm_params, "Kernel 4th param", 0, P);
		P->MoveKernel.shape_ = Params::get_double(pre_params, adm_params, "Kernel Shape", 1.0, P);
		if (P->KernelPowerScale != 1)
		{
			P->MoveKernel.shape_ *= P->KernelPowerScale;
		}
		P->AirportKernel.type_ = Params::get_int(pre_params, adm_params, "Airport Kernel Type", P->MoveKernel.type_, P);
		P->AirportKernel.scale_ = Params::get_double(pre_params, adm_params, "Airport Kernel Scale", P->MoveKernel.scale_, P);
		P->AirportKernel.shape_ = Params::get_double(pre_params, adm_params, "Airport Kernel Shape", P->MoveKernel.shape_, P);
		P->AirportKernel.p3_ = Params::get_double(pre_params, adm_params, "Airport Kernel 3rd param", P->MoveKernel.p3_, P);
		P->AirportKernel.p4_ = Params::get_double(pre_params, adm_params, "Airport Kernel 4th param", P->MoveKernel.p4_, P);

		P->DoPlaces = Params::get_int(pre_params, adm_params, "Include places", 1, P);
		if (P->DoPlaces != 0)
		{
			P->PlaceTypeNum = Params::get_int(pre_params, adm_params, "Number of types of places", 0, P);
			if (P->PlaceTypeNum == 0) P->DoPlaces = P->DoAirports = 0;
		}
		else
			P->PlaceTypeNum = P->DoAirports = 0;
	}
	if (P->DoPlaces != 0)
	{
		P->CareHomeResidentHouseholdScaling = Params::get_double(pre_params, adm_params, "Scaling of household contacts for care home residents", 1.0, P);
		P->CareHomeResidentSpatialScaling = Params::get_double(pre_params, adm_params, "Scaling of spatial contacts for care home residents", 1.0, P);
		P->CareHomeResidentPlaceScaling = Params::get_double(pre_params, adm_params, "Scaling of between group (home) contacts for care home residents", 1.0, P);
		P->CareHomeWorkerGroupScaling = Params::get_double(pre_params, adm_params, "Scaling of within group (home) contacts for care home workers", 1.0, P);
		P->CareHomeRelProbHosp = Params::get_double(pre_params, adm_params, "Relative probability that care home residents are hospitalised", 1.0, P);

		if (P->FitIter == 0)
		{
			if (P->PlaceTypeNum > NUM_PLACE_TYPES) ERR_CRITICAL("Too many place types\n");
			P->CareHomePlaceType = Params::get_int(pre_params, adm_params, "Place type number for care homes", -1, P);
			P->CareHomeAllowInitialInfections = Params::get_int(pre_params, adm_params, "Allow initial infections to be in care homes", 0, P);
			P->CareHomeResidentMinimumAge = Params::get_int(pre_params, adm_params, "Minimum age of care home residents", 1000, P);

			Params::req_int_vec(pre_params, adm_params, "Minimum age for age group 1 in place types", P->PlaceTypeAgeMin, P->PlaceTypeNum, P);
			Params::req_int_vec(pre_params, adm_params, "Maximum age for age group 1 in place types", P->PlaceTypeAgeMax, P->PlaceTypeNum, P);
			Params::req_double_vec(pre_params, adm_params, "Proportion of age group 1 in place types", P->PlaceTypePropAgeGroup, P->PlaceTypeNum, P);

			if (!Params::param_found(pre_params, adm_params, "Proportion of age group 2 in place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P->PlaceTypePropAgeGroup2[i] = 0;
					P->PlaceTypeAgeMin2[i] = 0;
					P->PlaceTypeAgeMax2[i] = 1000;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Proportion of age group 2 in place types", P->PlaceTypePropAgeGroup2, P->PlaceTypeNum, P);
				Params::req_int_vec(pre_params, adm_params, "Minimum age for age group 2 in place types", P->PlaceTypeAgeMin2, P->PlaceTypeNum, P);
				Params::req_int_vec(pre_params, adm_params, "Maximum age for age group 2 in place types", P->PlaceTypeAgeMax2, P->PlaceTypeNum, P);
			}
			if (!Params::param_found(pre_params, adm_params, "Proportion of age group 3 in place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P->PlaceTypePropAgeGroup3[i] = 0;
					P->PlaceTypeAgeMin3[i] = 0;
					P->PlaceTypeAgeMax3[i] = 1000;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Proportion of age group 3 in place types", P->PlaceTypePropAgeGroup3, P->PlaceTypeNum, P);
				Params::req_int_vec(pre_params, adm_params, "Minimum age for age group 3 in place types", P->PlaceTypeAgeMin3, P->PlaceTypeNum, P);
				Params::req_int_vec(pre_params, adm_params, "Maximum age for age group 3 in place types", P->PlaceTypeAgeMax3, P->PlaceTypeNum, P);
			}
			if (!Params::param_found(pre_params, adm_params, "Kernel shape params for place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P->PlaceTypeKernelShape[i] = P->MoveKernel.shape_;
					P->PlaceTypeKernelScale[i] = P->MoveKernel.scale_;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Kernel shape params for place types", P->PlaceTypeKernelShape, P->PlaceTypeNum, P);
				Params::req_double_vec(pre_params, adm_params, "Kernel scale params for place types", P->PlaceTypeKernelScale, P->PlaceTypeNum, P);
			}
			if (!Params::param_found(pre_params, adm_params, "Kernel 3rd param for place types"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)
				{
					P->PlaceTypeKernelP3[i] = P->MoveKernel.p3_;
					P->PlaceTypeKernelP4[i] = P->MoveKernel.p4_;
				}
			}
			else
			{
				Params::req_double_vec(pre_params, adm_params, "Kernel 3rd param for place types", P->PlaceTypeKernelP3, P->PlaceTypeNum, P);
				Params::req_double_vec(pre_params, adm_params, "Kernel 4th param for place types", P->PlaceTypeKernelP4, P->PlaceTypeNum, P);
			}
			Params::get_int_vec(pre_params, adm_params, "Number of closest places people pick from (0=all) for place types", P->PlaceTypeNearestNeighb, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
			if (P->DoAdUnits != 0)
			{
				Params::get_double_vec(params, pre_params, "Degree to which crossing administrative unit boundaries to go to places is inhibited", P->InhibitInterAdunitPlaceAssignment, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
			}

			Params::airport_params(adm_params, pre_params, params, P);

			Params::req_double_vec(pre_params, adm_params, "Mean size of place types", P->PlaceTypeMeanSize, P->PlaceTypeNum, P);
			Params::req_double_vec(pre_params, adm_params, "Param 1 of place group size distribution", P->PlaceTypeGroupSizeParam1, P->PlaceTypeNum, P);
			Params::get_double_vec(pre_params, adm_params, "Power of place size distribution", P->PlaceTypeSizePower, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);

			//added to enable lognormal distribution - ggilani 09/02/17
			Params::get_double_vec(pre_params, adm_params, "Standard deviation of place size distribution", P->PlaceTypeSizeSD, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
			Params::get_double_vec(pre_params, adm_params, "Offset of place size distribution", P->PlaceTypeSizeOffset, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
			Params::get_double_vec(pre_params, adm_params, "Maximum of place size distribution", P->PlaceTypeSizeMax, P->PlaceTypeNum, 1e20, NUM_PLACE_TYPES, P);
			Params::get_double_vec(pre_params, adm_params, "Minimum of place size distribution", P->PlaceTypeSizeMin, P->PlaceTypeNum, 1.0, NUM_PLACE_TYPES, P);
			Params::get_int_vec(pre_params, adm_params, "Kernel type for place types", P->PlaceTypeKernelType, P->PlaceTypeNum, P->MoveKernel.type_, NUM_PLACE_TYPES, P);
			Params::get_double_vec(pre_params, adm_params, "Place overlap matrix", P->PlaceExclusivityMatrix, P->PlaceTypeNum * P->PlaceTypeNum, 0, P->PlaceTypeNum * P->PlaceTypeNum, P);
			if (!Params::param_found(pre_params, adm_params, "Place overlap matrix"))
			{
				for (i = 0; i < NUM_PLACE_TYPES; i++)  // get_double_vec will set the zeroes if missing;
					P->PlaceExclusivityMatrix[i * (NUM_PLACE_TYPES + 1)] = 1; // this line sets the diagonal to 1 (identity matrix)
			}
		}
		/* Note P->PlaceExclusivityMatrix not used at present - places assumed exclusive (each person belongs to 0 or 1 place) */

		Params::req_double_vec(params, pre_params, "Proportion of between group place links", P->PlaceTypePropBetweenGroupLinks, P->PlaceTypeNum, P);
		Params::req_double_vec(params, pre_params, "Relative transmission rates for place types", P->PlaceTypeTrans, P->PlaceTypeNum, P);
		for (i = 0; i < P->PlaceTypeNum; i++) P->PlaceTypeTrans[i] *= AgeSuscScale;
	}

	Params::get_double_vec(params, pre_params, "Daily seasonality coefficients", P->Seasonality, DAYS_PER_YEAR, 1, DAYS_PER_YEAR, P);
	if (!Params::param_found(params, pre_params, "Daily seasonality coefficients"))
	{
		P->DoSeasonality = 0;
	}
	else
	{
		P->DoSeasonality = 1;
		s = 0;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			s += P->Seasonality[i];
		s += 1e-20;
		s /= DAYS_PER_YEAR;
		for (i = 0; i < DAYS_PER_YEAR; i++)
			P->Seasonality[i] /= s;
	}
	P->NumSeedLocations = Params::get_int(pre_params, adm_params, "Number of seed locations", 1, P);
	if (P->NumSeedLocations > MAX_NUM_SEED_LOCATIONS)
	{
		Files::xfprintf_stderr("Too many seed locations\n");
		P->NumSeedLocations = MAX_NUM_SEED_LOCATIONS;
	}
	Params::req_int_vec(pre_params, adm_params, "Initial number of infecteds", P->NumInitialInfections, P->NumSeedLocations, P);
	Params::get_double_matrix(pre_params, adm_params, "Location of initial infecteds", P->LocationInitialInfection, P->NumSeedLocations, 2, 0, P);
	P->MinPopDensForInitialInfection = Params::get_int(pre_params, adm_params, "Minimum population in microcell of initial infection", 0, P);
	P->MaxPopDensForInitialInfection = Params::get_int(pre_params, adm_params, "Maximum population in microcell of initial infection", 10000000, P);
	P->MaxAgeForInitialInfection = Params::get_int(pre_params, adm_params, "Maximum age of initial infections", 1000, P);
	P->DoRandomInitialInfectionLoc = Params::get_int(pre_params, adm_params, "Randomise initial infection location", 1, P);
	P->DoAllInitialInfectioninSameLoc = Params::get_int(pre_params, adm_params, "All initial infections located in same microcell", 0, P);

	if (Params::param_found(pre_params, adm_params, "Day of year of start of seeding"))
	{
		P->InitialInfectionCalTime = Params::req_double(pre_params, adm_params, "Day of year of start of seeding", P);
		if (Params::param_found(params, pre_params, "Scaling of infection seeding"))
		{
			P->SeedingScaling = Params::req_double(params, pre_params, "Scaling of infection seeding", P);
			P->DoNoCalibration = 1;
		}
		else
		{
			P->SeedingScaling = 1.0;
			P->DoNoCalibration = 0;
		}
	}
	else
	{
		P->SeedingScaling = 1.0;
		P->DoNoCalibration = 0;
		P->InitialInfectionCalTime = -1;
	}
	if (P->FitIter == 0)
	{
		if (P->DoAdUnits != 0)
		{
			if (!Params::param_found(pre_params, adm_params, "Administrative unit to seed initial infection into"))
			{
				Params::req_string_vec(pre_params, adm_params, "Administrative unit to seed initial infection into", AdunitListNames, P->NumSeedLocations, P);
				for (i = 0; i < P->NumSeedLocations; i++) P->InitialInfectionsAdminUnit[i] = 0;
			}
			else
				for (i = 0; i < P->NumSeedLocations; i++)
				{
					f = 0;
					if (P->NumAdunits > 0)
					{
						for (j = 0; (j < P->NumAdunits) && (!f); j++) f = (strcmp(AdUnits[j].ad_name, AdunitListNames[i]) == 0);
						if (f) k = AdUnits[j - 1].id;
					}
					if (!f) k = atoi(AdunitListNames[i]);
					P->InitialInfectionsAdminUnit[i] = k;
					P->InitialInfectionsAdminUnitId[i] = P->AdunitLevel1Lookup[(k % P->AdunitLevel1Mask) / P->AdunitLevel1Divisor];
				}
			Params::get_double_vec(pre_params, adm_params, "Administrative unit seeding weights", P->InitialInfectionsAdminUnitWeight, P->NumSeedLocations, 1.0, P->NumSeedLocations, P);
			s = 0;
			for (i = 0; i < P->NumSeedLocations; i++) s += P->InitialInfectionsAdminUnitWeight[i];
			for (i = 0; i < P->NumSeedLocations; i++) P->InitialInfectionsAdminUnitWeight[i] /= s;
		}
		else
		{
			for (i = 0; i < P->NumSeedLocations; i++) P->InitialInfectionsAdminUnit[i] = 0;
		}
	}
	P->InfectionImportRate1 = Params::get_double(params, pre_params, "Initial rate of importation of infections", 0, P);
	P->InfectionImportRate2 = Params::get_double(params, pre_params, "Changed rate of importation of infections", 0, P);
	P->InfectionImportChangeTime = Params::get_double(params, pre_params, "Time when infection rate changes", 1e10, P);
	P->DoImportsViaAirports = Params::get_int(params, pre_params, "Imports via air travel", 0, P);
	P->DurImportTimeProfile = Params::get_int(params, pre_params, "Length of importation time profile provided", 0, P);
	if (P->DurImportTimeProfile > 0)
	{
		if (P->DurImportTimeProfile >= MAX_DUR_IMPORT_PROFILE) ERR_CRITICAL("MAX_DUR_IMPORT_PROFILE too small\n");
		Params::req_double_vec(params, pre_params, "Daily importation time profile", P->ImportInfectionTimeProfile, P->DurImportTimeProfile, P);
	}

	P->R0 = Params::req_double(params, pre_params, "Reproduction number", P);
	if (Params::param_found(params, pre_params, "Beta for spatial transmission"))
	{
		P->LocalBeta = Params::req_double(params, pre_params, "Beta for spatial transmission", P);
		P->FixLocalBeta = 1;
	}
	else
	{
		P->LocalBeta = -1.0;
		P->FixLocalBeta = 0;
	}
	P->InfectiousPeriod = Params::req_double(params, pre_params, "Infectious period", P);
	P->SusceptibilitySD = Params::get_double(params, pre_params, "SD of individual variation in susceptibility", 0, P);
	P->InfectiousnessSD = Params::get_double(params, pre_params, "SD of individual variation in infectiousness", 0, P);
	if (Params::param_found(params, pre_params, "k of individual variation in infectiousness"))
		P->InfectiousnessSD = 1.0 / sqrt(Params::req_double(params, pre_params, "k of individual variation in infectiousness", P));
	P->NoInfectiousnessSDinHH = Params::get_int(params, pre_params, "k does not apply in households", 0, P);
	P->DoInfectiousnessProfile = Params::get_int(params, pre_params, "Model time varying infectiousness", 0, P);
	P->R0DensityScalePower = Params::get_double(params, pre_params, "Power of scaling of spatial R0 with density", 0, P);
	if (P->DoInfectiousnessProfile != 0)
	{
		Params::get_double_vec(params, pre_params, "Infectiousness profile", P->infectious_prof, INFPROF_RES, 1, INFPROF_RES, P);
		k = (int)ceil(P->InfectiousPeriod / P->ModelTimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		s = 0;
		P->infectious_prof[INFPROF_RES] = 0;
		for (i = 0; i < MAX_INFECTIOUS_STEPS; i++)	P->infectiousness[i] = 0;
		for (i = 0; i < k; i++)
		{
			t = (((double)i) * P->ModelTimeStep / P->InfectiousPeriod * INFPROF_RES);
			j = (int)t;
			t -= (double)j;
			if (j < INFPROF_RES)
				s += (P->infectiousness[i] = P->infectious_prof[j] * (1 - t) + P->infectious_prof[j + 1] * t);
			else
				s += (P->infectiousness[i] = P->infectious_prof[INFPROF_RES]);
		}
		s /= ((double)k);
		for (i = 0; i <= k; i++) P->infectiousness[i] /= s;
		P->infectious_icdf.assign_exponent(-1.0);
	}
	else
	{
		if (Params::param_found(params, pre_params, "Infectious period inverse CDF"))
		{
			Params::req_double_vec(params, pre_params, "Infectious period inverse CDF", P->infectious_icdf.get_values(), CDF_RES + 1, P);
			P->infectious_icdf.set_neg_log(ICDF_START);
		}
		k = (int)ceil(P->InfectiousPeriod * P->infectious_icdf[CDF_RES] / P->ModelTimeStep);
		if (k >= MAX_INFECTIOUS_STEPS) ERR_CRITICAL("MAX_INFECTIOUS_STEPS not big enough\n");
		for (i = 0; i < k; i++) P->infectiousness[i] = 1.0;
		P->infectiousness[k] = 0;
		P->infectious_icdf.assign_exponent();
	}
	P->DoLatent = Params::get_int(params, pre_params, "Include latent period", 0, P);
	if (P->DoLatent != 0)
	{
		P->LatentPeriod = Params::req_double(params, pre_params, "Latent period", P);
		Params::get_inverse_cdf(params, pre_params, "Latent period inverse CDF", &P->latent_icdf, P, 1e10);
	}

	P->DoSymptoms = Params::get_int(params, pre_params, "Include symptoms", 0, P);
	if (P->DoSymptoms == 0)
	{
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P->ProportionSymptomatic[i] = 0;
		P->FalsePositiveRate = 0;
		P->SymptInfectiousness = P->AsymptInfectiousness = 1.0;
		P->LatentToSymptDelay = 0;
	}
	else
	{
		if (P->DoAge != 0)
			Params::req_double_vec(params, pre_params, "Proportion symptomatic by age group", P->ProportionSymptomatic, NUM_AGE_GROUPS, P);
		else
		{
			P->ProportionSymptomatic[0] = Params::req_double(params, pre_params, "Proportion symptomatic", P);
			for (i = 1; i < NUM_AGE_GROUPS; i++)
				P->ProportionSymptomatic[i] = P->ProportionSymptomatic[0];
		}
		P->LatentToSymptDelay = Params::req_double(params, pre_params, "Delay from end of latent period to start of symptoms", P);
		P->SymptSpatialContactRate = Params::req_double(params, pre_params, "Relative rate of random contacts if symptomatic", P);
		P->SymptInfectiousness = Params::get_double(params, pre_params, "Symptomatic infectiousness relative to asymptomatic", 1.0, P);
		P->AsymptInfectiousness = Params::get_double(params, pre_params, "Asymptomatic infectiousness relative to symptomatic", 1.0, P);
		P->DoRealSymptWithdrawal = Params::get_int(params, pre_params, "Model symptomatic withdrawal to home as true absenteeism", 0, P);
		if (P->DoPlaces != 0)
		{
			Params::req_double_vec(params, pre_params, "Relative level of place attendance if symptomatic", P->SymptPlaceTypeContactRate, P->PlaceTypeNum, P);
			if (P->DoRealSymptWithdrawal != 0)
			{
				for (j = 0; j < NUM_PLACE_TYPES; j++)
				{
					P->SymptPlaceTypeWithdrawalProp[j] = 1.0 - P->SymptPlaceTypeContactRate[j];
					P->SymptPlaceTypeContactRate[j] = 1.0;
				}
			}
			else
				for (j = 0; j < NUM_PLACE_TYPES; j++) P->SymptPlaceTypeWithdrawalProp[j] = 0.0;
		}
		P->CaseAbsentChildAgeCutoff = Params::get_int(params, pre_params, "Maximum age of child at home for whom one adult also stays at home", 0, P);
		P->CaseAbsentChildPropAdultCarers = Params::get_double(params, pre_params, "Proportion of children at home for whom one adult also stays at home", 0, P);
		P->PlaceCloseRoundHousehold = Params::get_int(params, pre_params, "Place close round household", 1, P);
		P->AbsenteeismPlaceClosure = Params::get_int(params, pre_params, "Absenteeism place closure", 0, P);
		if (P->AbsenteeismPlaceClosure != 0)
		{
			P->CaseAbsenteeismDelay = 0; // Set to zero for tracking absenteeism
			P->MaxAbsentTime = Params::get_int(params, pre_params, "Max absent time", MAX_ABSENT_TIME, P);
			if (P->MaxAbsentTime > MAX_ABSENT_TIME || P->MaxAbsentTime < 0)
			{
				ERR_CRITICAL_FMT("[Max absent time] out of range (%d), should be in range [0, %d]", P->MaxAbsentTime, MAX_ABSENT_TIME);
			}
		}
		else
		{
			P->CaseAbsenteeismDelay = Params::get_double(params, pre_params, "Delay in starting place absenteeism for cases who withdraw", 0, P);
			P->MaxAbsentTime = 0; // Not used when !P->AbsenteeismPlaceClosure
		}
		P->CaseAbsenteeismDuration = Params::get_double(params, pre_params, "Duration of place absenteeism for cases who withdraw", 7, P);

		P->FalsePositiveRate = Params::get_double(params, pre_params, "False positive rate", 0.0, P);
		P->FalsePositivePerCapitaIncidence = Params::get_double(params, pre_params, "False positive per capita incidence", 0.0, P);
		Params::get_double_vec(params, pre_params, "False positive relative incidence by age", P->FalsePositiveAgeRate, NUM_AGE_GROUPS, 1.0, NUM_AGE_GROUPS, P);
	}

	P->SeroConvMaxSens = Params::get_double(params, pre_params, "Maximum sensitivity of serology assay", 1.0, P);
	P->SeroConvP1 = Params::get_double(params, pre_params, "Seroconversion model parameter 1", 14.0, P);
	P->SeroConvP2 = Params::get_double(params, pre_params, "Seroconversion model parameter 2", 3.0, P);
	P->SeroConvSpec = Params::get_double(params, pre_params, "Specificity of serology assay", 1.0, P);
	P->InfPrevSurveyScale = Params::get_double(params, pre_params, "Scaling of modelled infection prevalence to match surveys", 1.0, P);

	P->DoSeverity = Params::get_int(params, pre_params, "Do Severity Analysis", 0, P);
	if (P->DoSeverity != 0)
	{
		P->ScaleSymptProportions = Params::get_double(params, pre_params, "Factor to scale IFR", 1.0, P);
		//// Means for icdf's.
		P->Mean_TimeToTest = Params::get_double(params, pre_params, "MeanTimeToTest", 0.0, P);
		P->Mean_TimeToTestOffset = Params::get_double(params, pre_params, "MeanTimeToTestOffset", 1.0, P);
		P->Mean_TimeToTestCriticalOffset = Params::get_double(params, pre_params, "MeanTimeToTestCriticalOffset", 1.0, P);
		P->Mean_TimeToTestCritRecovOffset = Params::get_double(params, pre_params, "MeanTimeToTestCritRecovOffset", 1.0, P);
		if (Params::get_int(params, pre_params, "Age dependent severity delays", 0, P) == 0)
		{
			P->Mean_MildToRecovery[0] = Params::req_double(params, pre_params, "Mean_MildToRecovery", P);
			P->Mean_ILIToRecovery[0] = Params::req_double(params, pre_params, "Mean_ILIToRecovery", P);
			P->Mean_SARIToRecovery[0] = Params::req_double(params, pre_params, "Mean_SARIToRecovery", P);
			P->Mean_CriticalToCritRecov[0] = Params::req_double(params, pre_params, "Mean_CriticalToCritRecov", P);
			P->Mean_CritRecovToRecov[0] = Params::req_double(params, pre_params, "Mean_CritRecovToRecov", P);
			P->Mean_ILIToSARI[0] = Params::req_double(params, pre_params, "Mean_ILIToSARI", P);
			P->Mean_ILIToDeath[0] = Params::get_double(params, pre_params, "Mean_ILIToDeath", 7.0, P);
			P->Mean_SARIToCritical[0] = Params::req_double(params, pre_params, "Mean_SARIToCritical", P);
			P->Mean_SARIToDeath[0] = Params::req_double(params, pre_params, "Mean_SARIToDeath", P);
			P->Mean_CriticalToDeath[0] = Params::req_double(params, pre_params, "Mean_CriticalToDeath", P);
			for (int AgeGroup = 1; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
			{
				P->Mean_MildToRecovery[AgeGroup] = P->Mean_MildToRecovery[0];
				P->Mean_ILIToRecovery[AgeGroup] = P->Mean_ILIToRecovery[0];
				P->Mean_SARIToRecovery[AgeGroup] = P->Mean_SARIToRecovery[0];
				P->Mean_CriticalToCritRecov[AgeGroup] = P->Mean_CriticalToCritRecov[0];
				P->Mean_CritRecovToRecov[AgeGroup] = P->Mean_CritRecovToRecov[0];
				P->Mean_ILIToSARI[AgeGroup] = P->Mean_ILIToSARI[0];
				P->Mean_ILIToDeath[AgeGroup] = P->Mean_ILIToDeath[0];
				P->Mean_SARIToCritical[AgeGroup] = P->Mean_SARIToCritical[0];
				P->Mean_SARIToDeath[AgeGroup] = P->Mean_SARIToDeath[0];
				P->Mean_CriticalToDeath[AgeGroup] = P->Mean_CriticalToDeath[0];
			}
		}
		else
		{
			Params::req_double_vec(params, pre_params, "Mean_MildToRecovery", P->Mean_MildToRecovery, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_ILIToRecovery", P->Mean_ILIToRecovery, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_SARIToRecovery", P->Mean_SARIToRecovery, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_CriticalToCritRecov", P->Mean_CriticalToCritRecov, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_CritRecovToRecov", P->Mean_CritRecovToRecov, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_ILIToSARI", P->Mean_ILIToSARI, NUM_AGE_GROUPS, P);
			Params::get_double_vec(params, pre_params, "Mean_ILIToDeath", P->Mean_ILIToDeath, NUM_AGE_GROUPS, 7.0, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_SARIToCritical", P->Mean_SARIToCritical, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_SARIToDeath", P->Mean_SARIToDeath, NUM_AGE_GROUPS, P);
			Params::req_double_vec(params, pre_params, "Mean_CriticalToDeath", P->Mean_CriticalToDeath, NUM_AGE_GROUPS, P);
		}

		//// Get InverseCDFs
		Params::get_inverse_cdf(params, pre_params, "MildToRecovery_icdf", &P->MildToRecovery_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "ILIToRecovery_icdf", &P->ILIToRecovery_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "ILIToDeath_icdf", &P->ILIToDeath_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "SARIToRecovery_icdf", &P->SARIToRecovery_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "CriticalToCritRecov_icdf", &P->CriticalToCritRecov_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "CritRecovToRecov_icdf", &P->CritRecovToRecov_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "ILIToSARI_icdf", &P->ILIToSARI_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "SARIToCritical_icdf", &P->SARIToCritical_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "SARIToDeath_icdf", &P->SARIToDeath_icdf, P, ICDF_START);
		Params::get_inverse_cdf(params, pre_params, "CriticalToDeath_icdf", &P->CriticalToDeath_icdf, P, ICDF_START);

		// If you decide to decompose Critical -> Death transition into Critical -> Stepdown and Stepdown -> Death, use the block below.
		P->IncludeStepDownToDeath = Params::get_int(params, pre_params, "IncludeStepDownToDeath", 0, P);
		if (P->IncludeStepDownToDeath == 0) /// for backwards compatibility. If Stepdown to death not included (or if unspecified), set stepdown->death = stepdown->recovery.
		{
			for (int quantile = 0; quantile <= CDF_RES; quantile++)
				P->StepdownToDeath_icdf[quantile] = P->CritRecovToRecov_icdf[quantile];
			for (int AgeGroup = 0; AgeGroup < NUM_AGE_GROUPS; AgeGroup++)
				P->Mean_StepdownToDeath[AgeGroup] = P->Mean_CritRecovToRecov[AgeGroup];
		}
		else
		{
			Params::req_double_vec(params, pre_params, "Mean_StepdownToDeath", P->Mean_StepdownToDeath, NUM_AGE_GROUPS, P);
			Params::get_inverse_cdf(params, pre_params, "StepdownToDeath_icdf", &P->StepdownToDeath_icdf, P, ICDF_START);
		}

		Params::get_double_vec(params, pre_params, "Prop_Mild_ByAge", P->Prop_Mild_ByAge, NUM_AGE_GROUPS, 0.5, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "Prop_ILI_ByAge", P->Prop_ILI_ByAge, NUM_AGE_GROUPS, 0.3, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "Prop_SARI_ByAge", P->Prop_SARI_ByAge, NUM_AGE_GROUPS, 0.15, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "Prop_Critical_ByAge", P->Prop_Critical_ByAge, NUM_AGE_GROUPS, 0.05, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "CFR_SARI_ByAge", P->CFR_SARI_ByAge, NUM_AGE_GROUPS, 0.5, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "CFR_Critical_ByAge", P->CFR_Critical_ByAge, NUM_AGE_GROUPS, 0.5, NUM_AGE_GROUPS, P);
		Params::get_double_vec(params, pre_params, "CFR_ILI_ByAge", P->CFR_ILI_ByAge, NUM_AGE_GROUPS, 0, NUM_AGE_GROUPS, P);

		//Add param to allow severity to be uniformly scaled up or down.
		for (i = 0; i < NUM_AGE_GROUPS; i++)
		{
			P->Prop_SARI_ByAge[i] *= P->ScaleSymptProportions;
			P->Prop_Critical_ByAge[i] *= P->ScaleSymptProportions;
			P->Prop_ILI_ByAge[i] = 1.0 - P->Prop_Mild_ByAge[i] - P->Prop_SARI_ByAge[i] - P->Prop_Critical_ByAge[i];
		}
	}
	if (P->FitIter == 0)
	{
		if (Params::param_found(params, pre_params, "Bounding box for bitmap"))
		{
			double* buf = new double[4];
			Params::req_double_vec(params, pre_params, "Bounding box for bitmap", buf, 4, P);
			P->BoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(buf[0], buf[1]);
			P->BoundingBox.top_right() = CovidSim::Geometry::Vector2d(buf[2], buf[3]);
			delete[] buf;
		}
		else
		{
			P->BoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(0.0, 0.0);
			P->BoundingBox.top_right() = CovidSim::Geometry::Vector2d(1.0, 1.0);
		}
		if (Params::param_found(params, pre_params, "Spatial domain for simulation"))
		{
			double* buf = new double[4];
			Params::req_double_vec(params, pre_params, "Spatial domain for simulation", buf, 4, P);
			P->SpatialBoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(buf[0], buf[1]);
			P->SpatialBoundingBox.top_right() = CovidSim::Geometry::Vector2d(buf[2], buf[3]);
			delete[] buf;
		}
		else
		{
			P->SpatialBoundingBox.bottom_left() = CovidSim::Geometry::Vector2d(0.0, 0.0);
			P->SpatialBoundingBox.top_right() = CovidSim::Geometry::Vector2d(1.0, 1.0);
		}
		P->in_cells_.width = Params::get_double(params, pre_params, "Grid size", 1.0 / 120.0, P);
		P->DoUTM_coords = Params::get_int(params, pre_params, "Use long/lat coord system", 1, P);
		P->BitmapScale = Params::get_double(params, pre_params, "Bitmap scale", 1, P);
		P->BitmapAspectScale = Params::get_double(params, pre_params, "Bitmap y:x aspect scaling", 1, P);
		P->BitmapMovieFrame = Params::get_int(params, pre_params, "Bitmap movie frame interval", 250, P);
		P->OutputBitmap = Params::get_int(params, pre_params, "Output bitmap", 0, P);
		P->OutputBitmapDetected = Params::get_int(params, pre_params, "Output bitmap detected", 0, P);
		P->DoImmuneBitmap = Params::get_int(params, pre_params, "Output immunity on bitmap", 0, P);
		P->DoInfectionTree = Params::get_int(params, pre_params, "Output infection tree", 0, P);
		P->DoOneGen = Params::get_int(params, pre_params, "Do one generation", 0, P);
		P->OutputEveryRealisation = Params::get_int(params, pre_params, "Output every realisation", 0, P);
		P->MaxCorrSample = Params::get_int(params, pre_params, "Maximum number to sample for correlations", 1000000000, P);
		P->DoSI = Params::get_int(params, pre_params, "Assume SI model", 0, P);
		P->DoPeriodicBoundaries = Params::get_int(params, pre_params, "Assume periodic boundary conditions", 0, P);
		P->OutputOnlyNonExtinct = Params::get_int(params, pre_params, "Only output non-extinct realisations", 0, P);

		P->DoPerCapitaTriggers = Params::get_int(params, pre_params, "Use cases per thousand threshold for area controls", 0, P);
		P->DoGlobalTriggers = Params::get_int(params, pre_params, "Use global triggers for interventions", 0, P);
		P->DoAdminTriggers = Params::get_int(params, pre_params, "Use admin unit triggers for interventions", 0, P);
		P->DoICUTriggers = Params::get_int(params, pre_params, "Use ICU case triggers for interventions", 0, P);
		if (P->DoGlobalTriggers != 0) P->DoAdminTriggers = 0;
		P->IncThreshPop = Params::get_int(params, pre_params, "Divisor for per-capita area threshold (default 1000)", 1000, P);
		P->GlobalIncThreshPop = Params::get_int(params, pre_params, "Divisor for per-capita global threshold (default 1000)", 1000, P);

		P->TriggersSamplingInterval = Params::get_int(params, pre_params, "Number of sampling intervals over which cumulative incidence measured for global trigger", 10000000, P);
		P->PostAlertControlPropCasesId = Params::get_double(params, pre_params, "Proportion of cases detected for treatment", 1, P);
		P->PreAlertControlPropCasesId = Params::get_double(params, pre_params, "Proportion of cases detected before outbreak alert", 1, P);
		P->TriggerAlertOnDeaths = Params::get_int(params, pre_params, "Trigger alert on deaths", 0, P);
	}
	if (P->TriggerAlertOnDeaths != 0)
	{
		P->CaseOrDeathThresholdBeforeAlert = Params::get_int(params, pre_params, "Number of deaths accummulated before alert", 0, P);
	}
	else
	{
		P->CaseOrDeathThresholdBeforeAlert = Params::get_int(params, pre_params, "Number of detected cases needed before outbreak alert triggered", 0, P);
	}

	if (P->CaseOrDeathThresholdBeforeAlert_CommandLine > 0) P->CaseOrDeathThresholdBeforeAlert = P->CaseOrDeathThresholdBeforeAlert_CommandLine;
	P->DoAlertTriggerAfterInterv = Params::get_int(params, pre_params, "Alert trigger starts after interventions", 0, P);
	P->DateTriggerReached_CalTime = Params::get_double(params, pre_params, "Day of year trigger is reached", -1, P);
	if (P->DoAlertTriggerAfterInterv != 0)
	{
		P->Interventions_StartDate_CalTime = Params::req_double(params, pre_params, "Day of year interventions start", P);
		if (P->DateTriggerReached_CalTime <= P->Interventions_StartDate_CalTime)
			P->DoAlertTriggerAfterInterv = 0;
		else
		{
			P->AlertTriggerAfterIntervThreshold = P->CaseOrDeathThresholdBeforeAlert;
			P->CaseOrDeathThresholdBeforeAlert = 1000;
			Files::xfprintf_stderr("Threshold of %i deaths by day %lg\n", P->AlertTriggerAfterIntervThreshold, P->DateTriggerReached_CalTime);
		}
	}
	else
	{
		P->Interventions_StartDate_CalTime = P->DateTriggerReached_CalTime;
	}
	if (P->FitIter == 0)
	{
		P->CaseOrDeathThresholdBeforeAlert_Fixed = P->CaseOrDeathThresholdBeforeAlert;
	}

	P->WindowToEvaluateTriggerAlert = Params::get_int(params, pre_params, "Number of days to accummulate cases/deaths before alert", 1000, P);

	P->DoPlaceGroupTreat = Params::get_int(params, pre_params, "Only treat mixing groups within places", 0, P);

	P->TreatCellIncThresh = Params::get_double(params, pre_params, "Treatment trigger incidence per cell", INT32_MAX, P);
	P->CaseIsolation_CellIncThresh = Params::get_double(params, pre_params, "Case isolation trigger incidence per cell", P->TreatCellIncThresh, P);
	P->HHQuar_CellIncThresh = Params::get_double(params, pre_params, "Household quarantine trigger incidence per cell", P->TreatCellIncThresh, P);

	P->TreatSuscDrop = Params::get_double(params, pre_params, "Relative susceptibility of treated individual", 1, P);
	P->TreatInfDrop = Params::get_double(params, pre_params, "Relative infectiousness of treated individual", 1, P);
	P->TreatDeathDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases resulting in death prevented by treatment", 0, P);
	P->TreatSympDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases prevented by treatment", 0, P);
	P->TreatDelayMean = Params::get_double(params, pre_params, "Delay to treat cell", 0, P);
	P->TreatCaseCourseLength = Params::get_double(params, pre_params, "Duration of course of treatment", 5, P);
	P->TreatProphCourseLength = Params::get_double(params, pre_params, "Duration of course of prophylaxis", 10, P);
	P->TreatPropCases = Params::get_double(params, pre_params, "Proportion of detected cases treated", 1, P);
	if (P->DoHouseholds != 0)
	{
		P->TreatPropCaseHouseholds = Params::get_double(params, pre_params, "Proportion of households of cases treated", 0, P);
		P->TreatHouseholdsDuration = Params::get_double(params, pre_params, "Duration of household prophylaxis policy", USHRT_MAX / P->TimeStepsPerDay, P);
	}
	// Check below - "Proportional treated" will always be ignored.
	//if (!GetInputParameter2(params, pre_params, "Proportion treated", "%lf", (void*) & (P->TreatPropRadial), 1, 1, 0)) P->TreatPropRadial = 1.0;
	//if (!GetInputParameter2(params, pre_params, "Proportion treated in radial prophylaxis", "%lf", (void*) & (P->TreatPropRadial), 1, 1, 0)) P->TreatPropRadial = 1.0;

	P->TreatPropRadial = Params::get_double(params, pre_params, "Proportion treated in radial prophylaxis", 1.0, P);
	P->TreatRadius = Params::get_double(params, pre_params, "Treatment radius", 0, P);
	P->TreatPlaceGeogDuration = Params::get_double(params, pre_params, "Duration of place/geographic prophylaxis policy", USHRT_MAX / P->TimeStepsPerDay, P);
	P->TreatTimeStartBase = Params::get_double(params, pre_params, "Treatment start time", USHRT_MAX / P->TimeStepsPerDay, P);
	if (P->DoPlaces != 0)
	{
		Params::get_double_vec(params, pre_params, "Proportion of places treated after case detected", P->TreatPlaceProbCaseId, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
		Params::get_double_vec(params, pre_params, "Proportion of people treated in targeted places", P->TreatPlaceTotalProp, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
	}
	P->TreatMaxCoursesBase = Params::get_double(params, pre_params, "Maximum number of doses available", 1e20, P);
	P->TreatNewCoursesStartTime = Params::get_double(params, pre_params, "Start time of additional treatment production", USHRT_MAX / P->TimeStepsPerDay, P);
	P->TreatNewCoursesRate = Params::get_double(params, pre_params, "Rate of additional treatment production (courses per day)", 0, P);
	P->TreatMaxCoursesPerCase = Params::get_int(params, pre_params, "Maximum number of people targeted with radial prophylaxis per case", INT32_MAX, P);

	if (P->DoAdUnits != 0)
	{
		P->TreatByAdminUnit = Params::get_int(params, pre_params, "Treat administrative units rather than rings", 0, P);
		P->TreatAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for treatment", 1, P);
		if ((P->TreatAdminUnitDivisor == 0) || (P->TreatByAdminUnit == 0)) { P->TreatByAdminUnit = 0; P->TreatAdminUnitDivisor = 1; }
	}
	else
	{
		P->TreatAdminUnitDivisor = 1; P->TreatByAdminUnit = 0;
	}

	P->VaccCellIncThresh = Params::get_double(params, pre_params, "Vaccination trigger incidence per cell", 1000000000, P);
	P->VaccSuscDrop = Params::get_double(params, pre_params, "Relative susceptibility of vaccinated individual", 1, P);
	P->VaccSuscDrop2 = Params::get_double(params, pre_params, "Relative susceptibility of individual vaccinated after switch time", 1, P);
	P->VaccTimeEfficacySwitch = Params::get_double(params, pre_params, "Switch time at which vaccine efficacy increases", USHRT_MAX / P->TimeStepsPerDay, P);
	P->VaccEfficacyDecay = Params::get_double(params, pre_params, "Decay rate of vaccine efficacy (per year)", 0, P);
	P->VaccEfficacyDecay /= DAYS_PER_YEAR;
	P->VaccInfDrop = Params::get_double(params, pre_params, "Relative infectiousness of vaccinated individual", 1, P);
	P->VaccMortDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases resulting in death prevented by vaccination", 0, P);
	P->VaccSympDrop = Params::get_double(params, pre_params, "Proportion of symptomatic cases prevented by vaccination", 0, P);
	P->VaccDelayMean = Params::get_double(params, pre_params, "Delay to vaccinate", 0, P);

	P->VaccTimeToEfficacy = Params::get_double(params, pre_params, "Delay from vaccination to full protection", 0, P);

	P->VaccCampaignInterval = Params::get_double(params, pre_params, "Years between rounds of vaccination", 1e10, P);
	P->VaccDosePerDay = Params::get_int(params, pre_params, "Max vaccine doses per day", -1, P);
	P->VaccCampaignInterval *= DAYS_PER_YEAR;
	P->VaccMaxRounds = Params::get_int(params, pre_params, "Maximum number of rounds of vaccination", 1, P);
	if (P->DoHouseholds != 0)
	{
		P->VaccPropCaseHouseholds = Params::get_double(params, pre_params, "Proportion of households of cases vaccinated", 0, P);
		P->VaccHouseholdsDuration = Params::get_double(params, pre_params, "Duration of household vaccination policy", USHRT_MAX / P->TimeStepsPerDay, P);
	}

	P->VaccTimeStartBase = Params::get_double(params, pre_params, "Vaccination start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->VaccProp = Params::get_double(params, pre_params, "Proportion of population vaccinated", 0, P);
	P->VaccCoverageIncreasePeriod = Params::get_double(params, pre_params, "Time taken to reach max vaccination coverage (in years)", 0, P);
	P->VaccCoverageIncreasePeriod *= DAYS_PER_YEAR;
	P->VaccTimeStartGeo = Params::get_double(params, pre_params, "Time to start geographic vaccination", 1e10, P);
	P->VaccRadius = Params::get_double(params, pre_params, "Vaccination radius", 0, P);
	P->VaccMinRadius = Params::get_double(params, pre_params, "Minimum radius from case to vaccinate", 0, P);
	P->VaccMaxCoursesBase = Params::get_double(params, pre_params, "Maximum number of vaccine courses available", 1e20, P);
	P->VaccNewCoursesStartTime = Params::get_double(params, pre_params, "Start time of additional vaccine production", USHRT_MAX / P->TimeStepsPerDay, P);
	P->VaccNewCoursesEndTime = Params::get_double(params, pre_params, "End time of additional vaccine production", USHRT_MAX / P->TimeStepsPerDay, P);
	P->VaccNewCoursesRate = Params::get_double(params, pre_params, "Rate of additional vaccine production (courses per day)", 0, P);
	P->DoMassVacc = Params::get_int(params, pre_params, "Apply mass rather than reactive vaccination", 0, P);
	if (Params::param_found(params, pre_params, "Priority age range for mass vaccination")) {
		Params::req_int_vec(params, pre_params, "Priority age range for mass vaccination", P->VaccPriorityGroupAge, 2, P);
	}
	else {
		P->VaccPriorityGroupAge[0] = 1; P->VaccPriorityGroupAge[1] = 0;
	}

	if (P->DoAdUnits != 0)
	{
		P->VaccByAdminUnit = Params::get_int(params, pre_params, "Vaccinate administrative units rather than rings", 0, P);
		P->VaccAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for vaccination", 1, P);
		if ((P->VaccAdminUnitDivisor == 0) || (P->VaccByAdminUnit == 0)) P->VaccAdminUnitDivisor = 1;
	}
	else
	{
		P->VaccAdminUnitDivisor = 1; P->VaccByAdminUnit = 0;
	}

	P->MoveRestrCellIncThresh = Params::get_int(params, pre_params, "Movement restrictions trigger incidence per cell", INT32_MAX, P);
	P->MoveDelayMean = Params::get_double(params, pre_params, "Delay to start movement restrictions", 0, P);
	P->MoveRestrDuration = Params::get_double(params, pre_params, "Duration of movement restrictions", 7, P);
	P->MoveRestrEffect = Params::get_double(params, pre_params, "Residual movements after restrictions", 0, P);
	P->MoveRestrRadius = Params::get_double(params, pre_params, "Minimum radius of movement restrictions", 0, P);
	P->MoveRestrTimeStartBase = Params::get_double(params, pre_params, "Movement restrictions start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->DoBlanketMoveRestr = Params::get_int(params, pre_params, "Impose blanket movement restrictions", 0, P);
	P->DoMoveRestrOnceOnly = Params::get_int(params, pre_params, "Movement restrictions only once", 0, P);
	//if (P->DoMoveRestrOnceOnly) P->DoMoveRestrOnceOnly = 4; //// don't need this anymore with TreatStat option. Keep it as a boolean.
	if (P->DoAdUnits != 0)
	{
		P->MoveRestrByAdminUnit = Params::get_int(params, pre_params, "Movement restrictions in administrative units rather than rings", 0, P);
		P->MoveRestrAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for movement restrictions", 1, P);
		if ((P->MoveRestrAdminUnitDivisor == 0) || (P->MoveRestrByAdminUnit == 0)) P->MoveRestrAdminUnitDivisor = 1;
	}
	else
	{
		P->MoveRestrAdminUnitDivisor = 1; P->MoveRestrByAdminUnit = 0;
	}

	//Intervention delays and durations by admin unit: ggilani 16/03/20
	P->DoInterventionDelaysByAdUnit = Params::get_int(params, pre_params, "Include intervention delays by admin unit", 0, P);
	if (P->DoInterventionDelaysByAdUnit)
	{
		//Set up arrays to temporarily store parameters per admin unit
		double AdunitDelayToSocialDistance[MAX_ADUNITS];
		double AdunitDelayToHQuarantine[MAX_ADUNITS];
		double AdunitDelayToCaseIsolation[MAX_ADUNITS];
		double AdunitDelayToPlaceClose[MAX_ADUNITS];
		double AdunitDurationSocialDistance[MAX_ADUNITS];
		double AdunitDurationHQuarantine[MAX_ADUNITS];
		double AdunitDurationCaseIsolation[MAX_ADUNITS];
		double AdunitDurationPlaceClose[MAX_ADUNITS];

		Params::get_double_vec(params, pre_params, "Delay to social distancing by admin unit", AdunitDelayToSocialDistance, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Delay to household quarantine by admin unit", AdunitDelayToHQuarantine, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Delay to case isolation by admin unit", AdunitDelayToCaseIsolation, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Delay to place closure by admin unit", AdunitDelayToPlaceClose, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Duration of social distancing by admin unit", AdunitDurationSocialDistance, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Duration of household quarantine by admin unit", AdunitDurationHQuarantine, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Duration of case isolation by admin unit", AdunitDurationCaseIsolation, P->NumAdunits, 0, P->NumAdunits, P);
		Params::get_double_vec(params, pre_params, "Duration of place closure by admin unit", AdunitDurationPlaceClose, P->NumAdunits, 0, P->NumAdunits, P);

		for (i = 0; i < P->NumAdunits; i++)
		{
			AdUnits[i].SocialDistanceDelay = AdunitDelayToSocialDistance[i];
			AdUnits[i].SocialDistanceDuration = AdunitDurationSocialDistance[i];
			AdUnits[i].HQuarantineDelay = AdunitDelayToHQuarantine[i];
			AdUnits[i].HQuarantineDuration = AdunitDurationHQuarantine[i];
			AdUnits[i].CaseIsolationDelay = AdunitDelayToCaseIsolation[i];
			AdUnits[i].CaseIsolationPolicyDuration = AdunitDurationCaseIsolation[i];
			AdUnits[i].PlaceCloseDelay = AdunitDelayToPlaceClose[i];
			AdUnits[i].PlaceCloseDuration = AdunitDurationPlaceClose[i];
		}
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** DIGITAL CONTACT TRACING
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	//New code for digital contact tracing - ggilani: 09/03/20
	P->DoDigitalContactTracing = Params::get_int(params, pre_params, "Include digital contact tracing", 0, P);
	if (P->DoDigitalContactTracing != 0)
	{
		P->DigitalContactTracing_CellIncThresh = Params::get_double(params, pre_params, "Digital contact tracing trigger incidence per cell", 1000000000, P);

		P->PropPopUsingDigitalContactTracing = Params::get_double(params, pre_params, "Proportion of population or households covered by digital contact tracing", 1, P);
		Params::get_double_vec(params, pre_params, "Proportion of smartphone users by age", P->ProportionSmartphoneUsersByAge, NUM_AGE_GROUPS, 1, NUM_AGE_GROUPS, P);
		if (P->DoPlaces != 0)
		{
			P->ClusterDigitalContactUsers = Params::get_int(params, pre_params, "Cluster digital app clusters by household", 0, P); // by default, don't cluster by location
		}
		else
		{
			P->ClusterDigitalContactUsers = 0;
		}
		P->ProportionDigitalContactsIsolate = Params::get_double(params, pre_params, "Proportion of digital contacts who self-isolate", 0, P);
		P->MaxDigitalContactsToTrace = Params::get_int(params, pre_params, "Maximum number of contacts to trace per index case", MAX_CONTACTS, P);
		P->DigitalContactTracingDelay = Params::get_double(params, pre_params, "Delay between isolation of index case and contacts", P->ModelTimeStep, P);
		//we really need one timestep between to make sure contact is not processed before index
		if (P->DigitalContactTracingDelay == 0) P->DigitalContactTracingDelay = P->ModelTimeStep;
		P->LengthDigitalContactIsolation = Params::get_double(params, pre_params, "Length of self-isolation for digital contacts", 0, P);
		P->ScalingFactorSpatialDigitalContacts = Params::get_double(params, pre_params, "Spatial scaling factor - digital contact tracing", 1, P);
		P->ScalingFactorPlaceDigitalContacts = Params::get_double(params, pre_params, "Place scaling factor - digital contact tracing", 1, P);
		P->DigitalContactTracingTimeStartBase = Params::get_double(params, pre_params, "Digital contact tracing start time", USHRT_MAX / P->TimeStepsPerDay, P);
		P->DigitalContactTracingPolicyDuration = Params::get_double(params, pre_params, "Duration of digital contact tracing policy", 7, P);
		P->OutputDigitalContactTracing = Params::get_int(params, pre_params, "Output digital contact tracing", 0, P);
		P->OutputDigitalContactDist = Params::get_int(params, pre_params, "Output digital contact distribution", 0, P);

		if (P->DoInterventionDelaysByAdUnit)
		{
			double AdunitDelayToDCT[MAX_ADUNITS];
			double AdunitDurationDCT[MAX_ADUNITS];

			Params::get_double_vec(params, pre_params, "Delay to digital contact tracing by admin unit", AdunitDelayToDCT, P->NumAdunits, 0, P->NumAdunits, P);
			Params::get_double_vec(params, pre_params, "Duration of digital contact tracing by admin unit", AdunitDurationDCT, P->NumAdunits, 0, P->NumAdunits, P);
			for (i = 0; i < P->NumAdunits; i++)
			{
				AdUnits[i].DCTDelay = AdunitDelayToDCT[i];
				AdUnits[i].DCTDuration = AdunitDurationDCT[i];
			}
		}
		P->DCTIsolateIndexCases = Params::get_int(params, pre_params, "Isolate index cases in digital contact tracing", 1, P);
		P->DCTCaseIsolationEffectiveness = Params::get_double(params, pre_params, "Residual contacts after digital contact tracing isolation", P->CaseIsolationEffectiveness, P);
		P->DCTCaseIsolationHouseEffectiveness = Params::get_double(params, pre_params, "Residual household contacts after digital contact tracing isolation", P->CaseIsolationHouseEffectiveness, P);
		//initialise total number of users to 0
		P->NDigitalContactUsers = 0;
		P->NDigitalHouseholdUsers = 0;

		P->DelayFromIndexCaseDetectionToDCTIsolation = Params::get_double(params, pre_params, "Delay between symptom onset and isolation for index case", 0, P);
		P->DoDCTTest = Params::get_int(params, pre_params, "Test index cases and contacts", 0, P);
		P->DelayToTestIndexCase = Params::get_double(params, pre_params, "Delay to test index case", 1, P);
		P->DelayToTestDCTContacts = Params::get_double(params, pre_params, "Delay to test DCT contacts", 7, P);
		P->SpecificityDCT = Params::get_double(params, pre_params, "Testing specificity - DCT", 1, P);
		P->SensitivityDCT = Params::get_double(params, pre_params, "Testing sensitivity - DCT", 1, P);
		P->FindContactsOfDCTContacts = Params::get_int(params, pre_params, "Find contacts of digital contacts", 0, P);
		P->RemoveContactsOfNegativeIndexCase = Params::get_int(params, pre_params, "Remove contacts of a negative index case", 0, P);
	}
	else
	{
		//Set these to 1 so it doesn't interfere with code if we aren't using digital contact tracing.

		P->ScalingFactorSpatialDigitalContacts = 1;
		P->ScalingFactorPlaceDigitalContacts = 1;
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** PLACE CLOSURE
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****


	P->PlaceCloseCellIncThresh1 = Params::get_int(params, pre_params, "Trigger incidence per cell for place closure", 1000000000, P);
	P->PlaceCloseCellIncThresh2 = Params::get_int(params, pre_params, "Trigger incidence per cell for second place closure", 1000000000, P);
	if (P->PlaceCloseCellIncThresh1 < 0) P->PlaceCloseCellIncThresh1 = 1000000000;
	if (P->PlaceCloseCellIncThresh2 < 0) P->PlaceCloseCellIncThresh2 = 1000000000;
	P->PlaceCloseCellIncStopThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for end of place closure", 0, P);
	P->PlaceCloseDelayMean = Params::get_double(params, pre_params, "Delay to start place closure", 0, P);
	P->PlaceCloseDurationBase = Params::get_double(params, pre_params, "Duration of place closure", 7, P);
	P->PlaceCloseDuration2 = Params::get_double(params, pre_params, "Duration of second place closure", 7, P);
	if (P->DoPlaces != 0)
	{
		Params::get_double_vec(params, pre_params, "Proportion of places remaining open after closure by place type", P->PlaceCloseEffect, P->PlaceTypeNum, 1, NUM_PLACE_TYPES, P);
		Params::get_double_vec(params, pre_params, "Proportional attendance after closure by place type", P->PlaceClosePropAttending, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
	}
	if (P->DoHouseholds != 0)
		P->PlaceCloseHouseholdRelContact = Params::get_double(params, pre_params, "Relative household contact rate after closure", 1, P);
	P->PlaceCloseSpatialRelContact = Params::get_double(params, pre_params, "Relative spatial contact rate after closure", 1, P);

	P->DoHolidays = Params::get_int(pre_params, adm_params, "Include holidays", 0, P);
	if (P->DoHolidays != 0)
	{
		Params::get_double_vec(pre_params, adm_params, "Proportion of places remaining open during holidays by place type", P->HolidayEffect, P->PlaceTypeNum, 1, NUM_PLACE_TYPES, P);
		P->NumHolidays = Params::get_int(pre_params, adm_params, "Number of holidays", 0, P);
		if (P->NumHolidays > DAYS_PER_YEAR) P->NumHolidays = DAYS_PER_YEAR;
		if (P->NumHolidays > 0)
		{
			Params::req_double_vec(pre_params, adm_params, "Holiday start times", P->HolidayStartTime, P->NumHolidays, P);
			Params::req_double_vec(pre_params, adm_params, "Holiday durations", P->HolidayDuration, P->NumHolidays, P);
		}
	}
	else
	{
		P->NumHolidays = 0;
	}
	P->PlaceCloseRadius = Params::get_double(params, pre_params, "Minimum radius for place closure", 0, P);
	P->PlaceCloseTimeStartBase = Params::get_double(params, pre_params, "Place closure start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->PlaceCloseTimeStartBase2 = Params::get_double(params, pre_params, "Place closure second start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->DoPlaceCloseOnceOnly = Params::get_int(params, pre_params, "Places close only once", 0, P);
	//if (P->DoPlaceCloseOnceOnly) P->DoPlaceCloseOnceOnly = 4; //// don't need this anymore with TreatStat option. Keep it as a boolean.
	P->PlaceCloseIncTrig1 = Params::get_int(params, pre_params, "Place closure incidence threshold", 1, P);
	P->PlaceCloseIncTrig2 = Params::get_int(params, pre_params, "Place closure second incidence threshold", P->PlaceCloseIncTrig1, P);
	P->PlaceCloseFracIncTrig = Params::get_double(params, pre_params, "Place closure fractional incidence threshold", 0, P);
	if ((P->DoAdUnits != 0) && (P->DoPlaces != 0))
	{
		P->PlaceCloseByAdminUnit = Params::get_int(params, pre_params, "Place closure in administrative units rather than rings", 0, P);
		P->PlaceCloseAdminUnitDivisor = Params::get_int(params, pre_params, "Administrative unit divisor for place closure", 1, P);
		Params::get_int_vec(params, pre_params, "Place types to close for admin unit closure (0/1 array)", P->PlaceCloseAdunitPlaceTypes, P->PlaceTypeNum, 0, P->PlaceTypeNum, P);
		P->PlaceCloseCasePropThresh = Params::get_double(params, pre_params, "Cumulative proportion of place members needing to become sick for admin unit closure", 2, P);
		P->PlaceCloseAdunitPropThresh = Params::get_double(params, pre_params, "Proportion of places in admin unit needing to pass threshold for place closure", 2, P);
		if ((P->PlaceCloseAdminUnitDivisor < 1) || (P->PlaceCloseByAdminUnit == 0)) P->PlaceCloseAdminUnitDivisor = 1;
	}
	else
	{
		P->PlaceCloseAdminUnitDivisor = 1; P->PlaceCloseByAdminUnit = 0;
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** SOCIAL DISTANCING
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	P->SocDistCellIncThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for social distancing", 1000000000, P);
	P->SocDistCellIncStopThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for end of social distancing", 0, P);
	P->SocDistDuration = Params::get_double(params, pre_params, "Duration of social distancing", 7, P);
	P->SocDistDuration2 = Params::get_double(params, pre_params, "Duration of social distancing after change", 7, P);
	if (P->DoPlaces != 0)
	{
		Params::get_double_vec(params, pre_params, "Relative place contact rate given social distancing by place type", P->SocDistPlaceEffect, P->PlaceTypeNum, 1, NUM_PLACE_TYPES, P);
		Params::get_double_vec(params, pre_params, "Relative place contact rate given enhanced social distancing by place type", P->EnhancedSocDistPlaceEffect, P->PlaceTypeNum, 1, NUM_PLACE_TYPES, P);
		if (Params::param_found(params, pre_params, "Relative place contact rate given social distancing by place type after change"))
		{
			Params::req_double_vec(params, pre_params, "Relative place contact rate given social distancing by place type after change", P->SocDistPlaceEffect2, P->PlaceTypeNum, P);
		}
		else {
			for (i = 0; i < NUM_PLACE_TYPES; i++) P->SocDistPlaceEffect2[i] = P->SocDistPlaceEffect[i];
		}

		if (Params::param_found(params, pre_params, "Relative place contact rate given enhanced social distancing by place type after change"))
		{
			Params::req_double_vec(params, pre_params, "Relative place contact rate given enhanced social distancing by place type after change", P->EnhancedSocDistPlaceEffect2, P->PlaceTypeNum, P);
		}
		else
		{
			for (i = 0; i < NUM_PLACE_TYPES; i++) P->EnhancedSocDistPlaceEffect2[i] = P->EnhancedSocDistPlaceEffect[i];
		}
	}
	if (P->DoHouseholds != 0)
	{
		P->SocDistHouseholdEffect = Params::get_double(params, pre_params, "Relative household contact rate given social distancing", 1, P);
		P->EnhancedSocDistHouseholdEffect = Params::get_double(params, pre_params, "Relative household contact rate given enhanced social distancing", 1, P);
		P->SocDistHouseholdEffect2 = Params::get_double(params, pre_params, "Relative household contact rate given social distancing after change", P->SocDistHouseholdEffect, P);
		P->EnhancedSocDistHouseholdEffect2 = Params::get_double(params, pre_params, "Relative household contact rate given enhanced social distancing after change", P->EnhancedSocDistHouseholdEffect, P);
		P->EnhancedSocDistClusterByHousehold = Params::get_int(params, pre_params, "Cluster compliance with enhanced social distancing by household", 0, P);
	}
	else
	{
		P->EnhancedSocDistClusterByHousehold = 0;
	}
	P->SocDistSpatialEffect = Params::get_double(params, pre_params, "Relative spatial contact rate given social distancing", 1, P);
	P->SocDistSpatialEffect2 = Params::get_double(params, pre_params, "Relative spatial contact rate given social distancing after change", P->SocDistSpatialEffect, P);
	P->SocDistRadius = Params::get_double(params, pre_params, "Minimum radius for social distancing", 0, P);
	P->SocDistTimeStartBase = Params::get_double(params, pre_params, "Social distancing start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->SocDistChangeDelay = Params::get_double(params, pre_params, "Delay for change in effectiveness of social distancing", USHRT_MAX / P->TimeStepsPerDay, P);
	if (Params::param_found(params, pre_params, "Proportion compliant with enhanced social distancing by age group"))
	{
		Params::req_double_vec(params, pre_params, "Proportion compliant with enhanced social distancing by age group", P->EnhancedSocDistProportionCompliant, NUM_AGE_GROUPS, P);

	}
	else
	{
		t = Params::get_double(params, pre_params, "Proportion compliant with enhanced social distancing", 0, P);
		for (i = 0; i < NUM_AGE_GROUPS; i++)
			P->EnhancedSocDistProportionCompliant[i] = t;
	}

	P->EnhancedSocDistSpatialEffect = Params::get_double(params, pre_params, "Relative spatial contact rate given enhanced social distancing", 1, P);
	P->EnhancedSocDistSpatialEffect2 = Params::get_double(params, pre_params, "Relative spatial contact rate given enhanced social distancing after change", P->EnhancedSocDistSpatialEffect, P);

	P->DoSocDistOnceOnly = Params::get_int(params, pre_params, "Social distancing only once", 0, P);
	//if (P->DoSocDistOnceOnly) P->DoSocDistOnceOnly = 4; //// don't need this anymore with TreatStat option. Keep it as a boolean.

	P->AirportCloseEffectiveness = Params::get_double(params, pre_params, "Airport closure effectiveness", 0, P);
	P->AirportCloseEffectiveness = 1.0 - P->AirportCloseEffectiveness;
	P->AirportCloseTimeStartBase = Params::get_double(params, pre_params, "Airport closure start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->AirportCloseDuration = Params::get_double(params, pre_params, "Airport closure duration", USHRT_MAX / P->TimeStepsPerDay, P);

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** HOUSEHOLD QUARANTINE
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	if (P->DoHouseholds != 0)
	{
		P->DoHQretrigger = Params::get_int(params, pre_params, "Retrigger household quarantine with each new case in quarantine window", 0, P);
		P->HQuarantineTimeStartBase = Params::get_double(params, pre_params, "Household quarantine start time", USHRT_MAX / P->TimeStepsPerDay, P);
		P->HQuarantineDelay = Params::get_double(params, pre_params, "Delay to start household quarantine", 0, P);
		P->HQuarantineHouseDuration = Params::get_double(params, pre_params, "Length of time households are quarantined", 0, P);
		P->HQuarantinePolicyDuration = Params::get_double(params, pre_params, "Duration of household quarantine policy", USHRT_MAX / P->TimeStepsPerDay, P);
		P->HQuarantineHouseEffect = Params::get_double(params, pre_params, "Relative household contact rate after quarantine", 1, P);
		if (P->DoPlaces != 0)
		{
			Params::get_double_vec(params, pre_params, "Residual place contacts after household quarantine by place type", P->HQuarantinePlaceEffect, P->PlaceTypeNum, 1, NUM_PLACE_TYPES, P);
		}
		P->HQuarantineSpatialEffect = Params::get_double(params, pre_params, "Residual spatial contacts after household quarantine", 1, P);
		P->HQuarantinePropHouseCompliant = Params::get_double(params, pre_params, "Household level compliance with quarantine", 1, P);
		P->HQuarantinePropIndivCompliant = Params::get_double(params, pre_params, "Individual level compliance with quarantine", 1, P);
	}
	else
	{
		P->HQuarantineTimeStartBase = 1e10;
	}
	P->CaseIsolationTimeStartBase = Params::get_double(params, pre_params, "Case isolation start time", USHRT_MAX / P->TimeStepsPerDay, P);
	P->CaseIsolationProp = Params::get_double(params, pre_params, "Proportion of detected cases isolated", 0, P);
	P->CaseIsolationDelay = Params::get_double(params, pre_params, "Delay to start case isolation", 0, P);
	P->CaseIsolationDuration = Params::get_double(params, pre_params, "Duration of case isolation", 0, P);
	P->CaseIsolationPolicyDuration = Params::get_double(params, pre_params, "Duration of case isolation policy", 1e10, P);
	P->CaseIsolationEffectiveness = Params::get_double(params, pre_params, "Residual contacts after case isolation", 1, P);
	if (P->DoHouseholds != 0)
	{
		P->CaseIsolationHouseEffectiveness = Params::get_double(params, pre_params, "Residual household contacts after case isolation", P->CaseIsolationEffectiveness, P);
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** VARIABLE EFFICACIES OVER TIME
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****

	P->VaryEfficaciesOverTime = Params::get_int(params, pre_params, "Vary efficacies over time", 0, P);
	//// **** number of change times
	if (P->VaryEfficaciesOverTime == 0)
	{
		P->Num_SD_ChangeTimes = 1;
		P->Num_CI_ChangeTimes = 1;
		P->Num_HQ_ChangeTimes = 1;
		P->Num_PC_ChangeTimes = 1;
		P->Num_DCT_ChangeTimes = 1;
	}
	else
	{
		P->Num_SD_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of social distancing", 1, P);
		P->Num_CI_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of case isolation", 1, P);
		P->Num_HQ_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of household quarantine", 1, P);
		P->Num_PC_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of place closure", 1, P);
		P->Num_DCT_ChangeTimes = Params::get_int(params, pre_params, "Number of change times for levels of digital contact tracing", 1, P);
	}

	//// **** change times:
	//// By default, initialize first change time to zero and all subsequent change times to occur after simulation time, i.e. single value e.g. of efficacy for social distancing.
	P->SD_ChangeTimes[0] = 0;
	P->CI_ChangeTimes[0] = 0;
	P->HQ_ChangeTimes[0] = 0;
	P->PC_ChangeTimes[0] = 0;
	P->DCT_ChangeTimes[0] = 0;
	for (int ChangeTime = 1; ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES; ChangeTime++)
	{
		P->SD_ChangeTimes[ChangeTime] = 1e10;
		P->CI_ChangeTimes[ChangeTime] = 1e10;
		P->HQ_ChangeTimes[ChangeTime] = 1e10;
		P->PC_ChangeTimes[ChangeTime] = 1e10;
		P->DCT_ChangeTimes[ChangeTime] = 1e10;
		P->CFR_ChangeTimes_CalTime[ChangeTime] = INT32_MAX; // Out of bounds for int
	}
	//// Get real values from (pre)param file

	Params::get_double_vec(params, pre_params, "Change times for levels of social distancing", P->SD_ChangeTimes, P->Num_SD_ChangeTimes, 0, P->Num_SD_ChangeTimes, P);
	Params::get_double_vec(params, pre_params, "Change times for levels of case isolation", P->CI_ChangeTimes, P->Num_CI_ChangeTimes, 0, P->Num_CI_ChangeTimes, P);
	Params::get_double_vec(params, pre_params, "Change times for levels of household quarantine", P->HQ_ChangeTimes, P->Num_HQ_ChangeTimes, 0, P->Num_HQ_ChangeTimes, P);
	Params::get_double_vec(params, pre_params, "Change times for levels of place closure", P->PC_ChangeTimes, P->Num_PC_ChangeTimes, 0, P->Num_PC_ChangeTimes, P);
	Params::get_double_vec(params, pre_params, "Change times for levels of digital contact tracing", P->DCT_ChangeTimes, P->Num_DCT_ChangeTimes, 0, P->Num_DCT_ChangeTimes, P);

	// initialize to zero (regardless of whether doing places or households).
	for (int ChangeTime = 0; ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES; ChangeTime++)
	{
		//// **** "efficacies"
		//// spatial
		P->SD_SpatialEffects_OverTime[ChangeTime] = 0;
		P->Enhanced_SD_SpatialEffects_OverTime[ChangeTime] = 0;
		P->CI_SpatialAndPlaceEffects_OverTime[ChangeTime] = 0;
		P->HQ_SpatialEffects_OverTime[ChangeTime] = 0;
		P->PC_SpatialEffects_OverTime[ChangeTime] = 0;
		P->DCT_SpatialAndPlaceEffects_OverTime[ChangeTime] = 0;

		//// Household
		P->SD_HouseholdEffects_OverTime[ChangeTime] = 0;
		P->Enhanced_SD_HouseholdEffects_OverTime[ChangeTime] = 0;
		P->CI_HouseholdEffects_OverTime[ChangeTime] = 0;
		P->HQ_HouseholdEffects_OverTime[ChangeTime] = 0;
		P->PC_HouseholdEffects_OverTime[ChangeTime] = 0;
		P->DCT_HouseholdEffects_OverTime[ChangeTime] = 0;

		//// place
		for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
		{
			P->SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = 0;
			P->Enhanced_SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = 0;
			P->HQ_PlaceEffects_OverTime[ChangeTime][PlaceType] = 0;
			P->PC_PlaceEffects_OverTime[ChangeTime][PlaceType] = 0;
		}
		P->PC_Durs_OverTime[ChangeTime] = 0;

		//// **** compliance
		P->CI_Prop_OverTime[ChangeTime] = 0;
		P->HQ_Individual_PropComply_OverTime[ChangeTime] = 0;
		P->HQ_Household_PropComply_OverTime[ChangeTime] = 0;
		P->DCT_Prop_OverTime[ChangeTime] = 0;
	}


	//// **** "efficacies": by default, initialize to values read in previously.
	///// spatial contact rates rates over time (and place too for CI and DCT)
	bool v = (P->VaryEfficaciesOverTime == 0);
	Params::get_double_vec_ff(v, params, pre_params, "Relative spatial contact rates over time given social distancing", P->SD_SpatialEffects_OverTime, P->Num_SD_ChangeTimes, P->SocDistSpatialEffect, P);
	Params::get_double_vec_ff(v, params, pre_params, "Relative spatial contact rates over time given enhanced social distancing", P->Enhanced_SD_SpatialEffects_OverTime, P->Num_SD_ChangeTimes, P->EnhancedSocDistSpatialEffect, P);
	Params::get_double_vec_ff(v, params, pre_params, "Residual contacts after case isolation over time", P->CI_SpatialAndPlaceEffects_OverTime, P->Num_CI_ChangeTimes, P->CaseIsolationEffectiveness, P);
	Params::get_double_vec_ff(v, params, pre_params, "Residual spatial contacts over time after household quarantine", P->HQ_SpatialEffects_OverTime, P->Num_HQ_ChangeTimes, P->HQuarantineSpatialEffect, P);
	Params::get_double_vec_ff(v, params, pre_params, "Relative spatial contact rates over time after place closure", P->PC_SpatialEffects_OverTime, P->Num_PC_ChangeTimes, P->PlaceCloseSpatialRelContact, P);
	Params::get_double_vec_ff(v, params, pre_params, "Residual contacts after digital contact tracing isolation over time", P->DCT_SpatialAndPlaceEffects_OverTime, P->Num_DCT_ChangeTimes, P->DCTCaseIsolationEffectiveness, P);

	///// Household contact rates over time
	if (P->DoHouseholds != 0)
	{
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time given social distancing", P->SD_HouseholdEffects_OverTime, P->Num_SD_ChangeTimes, P->SocDistHouseholdEffect, P);
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time given enhanced social distancing", P->Enhanced_SD_HouseholdEffects_OverTime, P->Num_SD_ChangeTimes, P->EnhancedSocDistHouseholdEffect, P);
		Params::get_double_vec_ff(v, params, pre_params, "Residual household contacts after case isolation over time", P->CI_HouseholdEffects_OverTime, P->Num_CI_ChangeTimes, P->CaseIsolationHouseEffectiveness, P);
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time after quarantine", P->HQ_HouseholdEffects_OverTime, P->Num_HQ_ChangeTimes, P->HQuarantineHouseEffect, P);
		Params::get_double_vec_ff(v, params, pre_params, "Relative household contact rates over time after place closure", P->PC_HouseholdEffects_OverTime, P->Num_PC_ChangeTimes, P->PlaceCloseHouseholdRelContact, P);
		Params::get_double_vec_ff(v, params, pre_params, "Residual household contacts after digital contact tracing isolation over time", P->DCT_HouseholdEffects_OverTime, P->Num_DCT_ChangeTimes, P->DCTCaseIsolationHouseEffectiveness, P);
	}

	///// place contact rates over time
	if (P->DoPlaces != 0)
	{
		//// soc dist
		Params::get_double_matrix(params, pre_params, "Relative place contact rates over time given social distancing by place type", P->SD_PlaceEffects_OverTime, P->Num_SD_ChangeTimes, P->PlaceTypeNum, 0, P);
		if (v || !Params::param_found(params, pre_params, "Relative place contact rates over time given social distancing by place type"))
			for (int ChangeTime = 0; ChangeTime < P->Num_SD_ChangeTimes; ChangeTime++) //// by default populate to values of P->SocDistPlaceEffect
				for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
					P->SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = P->SocDistPlaceEffect[PlaceType];

		//// enhanced soc dist
		Params::get_double_matrix(params, pre_params, "Relative place contact rates over time given enhanced social distancing by place type", P->Enhanced_SD_PlaceEffects_OverTime, P->Num_SD_ChangeTimes, P->PlaceTypeNum, 0, P);
		if (v || !Params::param_found(params, pre_params, "Relative place contact rates over time given enhanced social distancing by place type"))
			for (int ChangeTime = 0; ChangeTime < P->Num_SD_ChangeTimes; ChangeTime++) //// by default populate to values of P->EnhancedSocDistPlaceEffect
				for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
					P->Enhanced_SD_PlaceEffects_OverTime[ChangeTime][PlaceType] = P->EnhancedSocDistPlaceEffect[PlaceType];

		//// household quarantine
		Params::get_double_matrix(params, pre_params, "Residual place contacts over time after household quarantine by place type", P->HQ_PlaceEffects_OverTime, P->Num_HQ_ChangeTimes, P->PlaceTypeNum, 0, P);
		if (v || !Params::param_found(params, pre_params, "Residual place contacts over time after household quarantine by place type"))
			for (int ChangeTime = 0; ChangeTime < P->Num_HQ_ChangeTimes; ChangeTime++) //// by default populate to values of P->HQuarantinePlaceEffect
				for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
					P->HQ_PlaceEffects_OverTime[ChangeTime][PlaceType] = P->HQuarantinePlaceEffect[PlaceType];

		//// place closure
		Params::get_double_matrix(params, pre_params, "Proportion of places remaining open after closure by place type over time", P->PC_PlaceEffects_OverTime, P->Num_PC_ChangeTimes, P->PlaceTypeNum, 0, P);
		if (v || !Params::param_found(params, pre_params, "Proportion of places remaining open after closure by place type over time"))
			for (int ChangeTime = 0; ChangeTime < P->Num_PC_ChangeTimes; ChangeTime++) //// by default populate to values of P->PlaceCloseEffect
				for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
					P->PC_PlaceEffects_OverTime[ChangeTime][PlaceType] = P->PlaceCloseEffect[PlaceType];

		Params::get_double_matrix(params, pre_params, "Proportional attendance after closure by place type over time", P->PC_PropAttending_OverTime, P->Num_PC_ChangeTimes, P->PlaceTypeNum, 0, P);
		if (v || !Params::param_found(params, pre_params, "Proportional attendance after closure by place type over time"))
			for (int ChangeTime = 0; ChangeTime < P->Num_PC_ChangeTimes; ChangeTime++) //// by default populate to values of P->PlaceClosePropAttending
				for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
					P->PC_PropAttending_OverTime[ChangeTime][PlaceType] = P->PlaceClosePropAttending[PlaceType];
	}


	//// ****  compliance
	Params::get_double_vec_ff(v, params, pre_params, "Proportion of detected cases isolated over time", P->CI_Prop_OverTime, P->Num_CI_ChangeTimes, P->CaseIsolationProp, P);
	Params::get_double_vec_ff(v, params, pre_params, "Individual level compliance with quarantine over time", P->HQ_Individual_PropComply_OverTime, P->Num_HQ_ChangeTimes, P->HQuarantinePropIndivCompliant, P);
	Params::get_double_vec_ff(v, params, pre_params, "Household level compliance with quarantine over time", P->HQ_Household_PropComply_OverTime, P->Num_HQ_ChangeTimes, P->HQuarantinePropHouseCompliant, P);
	Params::get_double_vec_ff(v, params, pre_params, "Proportion of digital contacts who self-isolate over time", P->DCT_Prop_OverTime, P->Num_DCT_ChangeTimes, P->ProportionDigitalContactsIsolate, P);
	Params::get_int_vec_ff(v, params, pre_params, "Maximum number of contacts to trace per index case over time", P->DCT_MaxToTrace_OverTime, P->Num_DCT_ChangeTimes, P->MaxDigitalContactsToTrace, P);
	if (P->DoPlaces != 0)
	{
		//// ****  thresholds
		//// place closure (global threshold)
		Params::get_int_vec_ff(v, params, pre_params, "Place closure incidence threshold over time", P->PC_IncThresh_OverTime, P->Num_PC_ChangeTimes, P->PlaceCloseIncTrig1, P);
		Params::get_double_vec_ff(v, params, pre_params, "Place closure fractional incidence threshold over time", P->PC_FracIncThresh_OverTime, P->Num_PC_ChangeTimes, P->PlaceCloseFracIncTrig, P);
		Params::get_int_vec_ff(v, params, pre_params, "Trigger incidence per cell for place closure over time", P->PC_CellIncThresh_OverTime, P->Num_PC_ChangeTimes, P->PlaceCloseCellIncThresh1, P);
		for (int ChangeTime = 0; ChangeTime < P->Num_PC_ChangeTimes; ChangeTime++)
			if (P->PC_CellIncThresh_OverTime[ChangeTime] < 0) P->PC_CellIncThresh_OverTime[ChangeTime] = 1000000000; // allows -1 to be used as a proxy for no cell-based triggering
	}
	//// household quarantine
	Params::get_double_vec_ff(v, params, pre_params, "Household quarantine trigger incidence per cell over time", P->HQ_CellIncThresh_OverTime, P->Num_HQ_ChangeTimes, P->HHQuar_CellIncThresh, P);
	Params::get_double_vec_ff(v, params, pre_params, "Case isolation trigger incidence per cell over time", P->CI_CellIncThresh_OverTime, P->Num_CI_ChangeTimes, P->CaseIsolation_CellIncThresh, P);
	Params::get_int_vec_ff(v, params, pre_params, "Trigger incidence per cell for social distancing over time", P->SD_CellIncThresh_OverTime, P->Num_SD_ChangeTimes, P->SocDistCellIncThresh, P);

	//// **** Durations (later add Case isolation and Household quarantine)
	// place closure
	Params::get_double_vec_ff(v, params, pre_params, "Duration of place closure over time", P->PC_Durs_OverTime, P->Num_PC_ChangeTimes, P->PlaceCloseDurationBase, P);

	//// Guards: make unused change values in array equal to final used value
	if (P->VaryEfficaciesOverTime)
	{
		//// soc dist
		for (int SD_ChangeTime = P->Num_SD_ChangeTimes; SD_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; SD_ChangeTime++)
		{
			//// non-enhanced
			P->SD_SpatialEffects_OverTime[SD_ChangeTime] = P->SD_SpatialEffects_OverTime[P->Num_SD_ChangeTimes - 1];
			P->SD_HouseholdEffects_OverTime[SD_ChangeTime] = P->SD_HouseholdEffects_OverTime[P->Num_SD_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
				P->SD_PlaceEffects_OverTime[SD_ChangeTime][PlaceType] = P->SD_PlaceEffects_OverTime[P->Num_SD_ChangeTimes - 1][PlaceType];
			//// enhanced
			P->Enhanced_SD_SpatialEffects_OverTime[SD_ChangeTime] = P->Enhanced_SD_SpatialEffects_OverTime[P->Num_SD_ChangeTimes - 1];
			P->Enhanced_SD_HouseholdEffects_OverTime[SD_ChangeTime] = P->Enhanced_SD_HouseholdEffects_OverTime[P->Num_SD_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
				P->Enhanced_SD_PlaceEffects_OverTime[SD_ChangeTime][PlaceType] = P->Enhanced_SD_PlaceEffects_OverTime[P->Num_SD_ChangeTimes - 1][PlaceType];

			P->SD_CellIncThresh_OverTime[SD_ChangeTime] = P->SD_CellIncThresh_OverTime[P->Num_SD_ChangeTimes - 1];
		}

		//// case isolation
		for (int CI_ChangeTime = P->Num_CI_ChangeTimes; CI_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; CI_ChangeTime++)
		{
			P->CI_SpatialAndPlaceEffects_OverTime[CI_ChangeTime] = P->CI_SpatialAndPlaceEffects_OverTime[P->Num_CI_ChangeTimes - 1];
			P->CI_HouseholdEffects_OverTime[CI_ChangeTime] = P->CI_HouseholdEffects_OverTime[P->Num_CI_ChangeTimes - 1];
			P->CI_Prop_OverTime[CI_ChangeTime] = P->CI_Prop_OverTime[P->Num_CI_ChangeTimes - 1];
			P->CI_CellIncThresh_OverTime[CI_ChangeTime] = P->CI_CellIncThresh_OverTime[P->Num_CI_ChangeTimes - 1];
		}

		//// household quarantine
		for (int HQ_ChangeTime = P->Num_HQ_ChangeTimes; HQ_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; HQ_ChangeTime++)
		{
			P->HQ_SpatialEffects_OverTime[HQ_ChangeTime] = P->HQ_SpatialEffects_OverTime[P->Num_HQ_ChangeTimes - 1];
			P->HQ_HouseholdEffects_OverTime[HQ_ChangeTime] = P->HQ_HouseholdEffects_OverTime[P->Num_HQ_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
				P->HQ_PlaceEffects_OverTime[HQ_ChangeTime][PlaceType] = P->HQ_PlaceEffects_OverTime[P->Num_HQ_ChangeTimes - 1][PlaceType];

			P->HQ_Individual_PropComply_OverTime[HQ_ChangeTime] = P->HQ_Individual_PropComply_OverTime[P->Num_HQ_ChangeTimes - 1];
			P->HQ_Household_PropComply_OverTime[HQ_ChangeTime] = P->HQ_Household_PropComply_OverTime[P->Num_HQ_ChangeTimes - 1];

			P->HQ_CellIncThresh_OverTime[HQ_ChangeTime] = P->HQ_CellIncThresh_OverTime[P->Num_HQ_ChangeTimes - 1];
		}

		//// place closure
		for (int PC_ChangeTime = P->Num_PC_ChangeTimes; PC_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; PC_ChangeTime++)
		{
			P->PC_SpatialEffects_OverTime[PC_ChangeTime] = P->PC_SpatialEffects_OverTime[P->Num_PC_ChangeTimes - 1];
			P->PC_HouseholdEffects_OverTime[PC_ChangeTime] = P->PC_HouseholdEffects_OverTime[P->Num_PC_ChangeTimes - 1];
			for (int PlaceType = 0; PlaceType < P->PlaceTypeNum; PlaceType++)
			{
				P->PC_PlaceEffects_OverTime[PC_ChangeTime][PlaceType] = P->PC_PlaceEffects_OverTime[P->Num_PC_ChangeTimes - 1][PlaceType];
				P->PC_PropAttending_OverTime[PC_ChangeTime][PlaceType] = P->PC_PropAttending_OverTime[P->Num_PC_ChangeTimes - 1][PlaceType];
			}

			P->PC_IncThresh_OverTime[PC_ChangeTime] = P->PC_IncThresh_OverTime[P->Num_PC_ChangeTimes - 1];
			P->PC_FracIncThresh_OverTime[PC_ChangeTime] = P->PC_FracIncThresh_OverTime[P->Num_PC_ChangeTimes - 1];
			P->PC_CellIncThresh_OverTime[PC_ChangeTime] = P->PC_CellIncThresh_OverTime[P->Num_PC_ChangeTimes - 1];
		}

		//// digital contact tracing
		for (int DCT_ChangeTime = P->Num_DCT_ChangeTimes; DCT_ChangeTime < MAX_NUM_INTERVENTION_CHANGE_TIMES - 1; DCT_ChangeTime++)
		{
			P->DCT_SpatialAndPlaceEffects_OverTime[DCT_ChangeTime] = P->DCT_SpatialAndPlaceEffects_OverTime[P->Num_DCT_ChangeTimes - 1];
			P->DCT_HouseholdEffects_OverTime[DCT_ChangeTime] = P->DCT_HouseholdEffects_OverTime[P->Num_DCT_ChangeTimes - 1];
			P->DCT_Prop_OverTime[DCT_ChangeTime] = P->DCT_Prop_OverTime[P->Num_DCT_ChangeTimes - 1];
			P->DCT_MaxToTrace_OverTime[DCT_ChangeTime] = P->DCT_MaxToTrace_OverTime[P->Num_DCT_ChangeTimes - 1];
		}
	}

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	///// **** CFR SCALINGS OVER TIME (logic to set this up is the same as for VARIABLE EFFICACIES OVER TIME)
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// ****
	P->Num_CFR_ChangeTimes = Params::get_int(params, pre_params, "Num_CFR_ChangeTimes", 1, P);

	//// By default, initialize first change time to zero and all subsequent change times to occur after simulation time, i.e. single value e.g. of Critical CFR.
	P->CFR_ChangeTimes_CalTime[0] = 0;

	for (int ChangeTime = 1; ChangeTime < MAX_NUM_CFR_CHANGE_TIMES; ChangeTime++) P->CFR_ChangeTimes_CalTime[ChangeTime] = INT32_MAX;
	Params::get_int_vec(params, pre_params, "CFR_ChangeTimes_CalTime", P->CFR_ChangeTimes_CalTime, P->Num_CFR_ChangeTimes, INT32_MAX, P->Num_CFR_ChangeTimes, P);

	// Get various CFR scalings. 
	Params::get_double_vec(params, pre_params, "CFR_TimeScaling_Critical", P->CFR_TimeScaling_Critical, P->Num_CFR_ChangeTimes, 1, P->Num_CFR_ChangeTimes, P);
	Params::get_double_vec(params, pre_params, "CFR_TimeScaling_SARI", P->CFR_TimeScaling_SARI, P->Num_CFR_ChangeTimes, 1, P->Num_CFR_ChangeTimes, P);
	Params::get_double_vec(params, pre_params, "CFR_TimeScaling_ILI", P->CFR_TimeScaling_ILI, P->Num_CFR_ChangeTimes, 1, P->Num_CFR_ChangeTimes, P);

	if (P->FitIter == 0)
	{
		if (P->DoPlaces != 0)
		{
			P->KeyWorkerPopNum = Params::get_int(params, pre_params, "Number of key workers randomly distributed in the population", 0, P);
			Params::get_int_vec(params, pre_params, "Number of key workers in different places by place type", P->KeyWorkerPlaceNum, P->PlaceTypeNum, 0, NUM_PLACE_TYPES, P);
			Params::get_double_vec(params, pre_params, "Proportion of staff who are key workers per chosen place by place type", P->KeyWorkerPropInKeyPlaces, P->PlaceTypeNum, 1, NUM_PLACE_TYPES, P);
			P->KeyWorkerProphCellIncThresh = Params::get_int(params, pre_params, "Trigger incidence per cell for key worker prophylaxis", 1000000000, P);
			P->KeyWorkerProphTimeStartBase = Params::get_double(params, pre_params, "Key worker prophylaxis start time", USHRT_MAX / P->TimeStepsPerDay, P);
			P->KeyWorkerProphDuration = Params::get_double(params, pre_params, "Duration of key worker prophylaxis", 0, P);
			P->KeyWorkerProphRenewalDuration = Params::get_double(params, pre_params, "Time interval from start of key worker prophylaxis before policy restarted", P->KeyWorkerProphDuration, P);
			if (P->DoHouseholds != 0)
			{
				P->KeyWorkerHouseProp = Params::get_double(params, pre_params, "Proportion of key workers whose households are also treated as key workers", 0, P);
			}
			P->KeyWorkerProphRadius = Params::get_double(params, pre_params, "Minimum radius for key worker prophylaxis", 0, P);
		}
		else
		{
			P->KeyWorkerPopNum = 0;
			P->KeyWorkerProphTimeStartBase = 1e10;
		}

		//Added this to parameter list so that recording infection events (and the number to record) can easily be turned off and on: ggilani - 10/10/2014
		P->DoRecordInfEvents = Params::get_int(params, pre_params, "Record infection events", 0, P);
		if (P->DoRecordInfEvents != 0)
		{
			P->MaxInfEvents = Params::get_int(params, pre_params, "Max number of infection events to record", 1000, P);
			P->RecordInfEventsPerRun = Params::get_int(params, pre_params, "Record infection events per run", 0, P);
		}
		else
		{
			P->MaxInfEvents = 0;
		}
		//Include a limit to the number of infections to simulate, if this happens before time runs out
		P->LimitNumInfections = Params::get_int(params, pre_params, "Limit number of infections", 0, P);
		if (P->LimitNumInfections != 0)
		{
			P->MaxNumInfections = Params::get_int(params, pre_params, "Max number of infections", 60000, P);
		}
		//Add origin-destination matrix parameter
		P->DoOriginDestinationMatrix = Params::get_int(params, pre_params, "Output origin destination matrix", 0, P);

		P->MeanChildAgeGap = Params::get_int(params, pre_params, adm_params, "Mean child age gap", 2, P);
		P->MinAdultAge = Params::get_int(params, pre_params, adm_params, "Min adult age", 19, P);
		P->MaxMFPartnerAgeGap = Params::get_int(params, pre_params, adm_params, "Max MF partner age gap", 5, P);
		P->MaxFMPartnerAgeGap = Params::get_int(params, pre_params, adm_params, "Max FM partner age gap", 5, P);
		P->MinParentAgeGap = Params::get_int(params, pre_params, adm_params, "Min parent age gap", 19, P);
		P->MaxParentAgeGap = Params::get_int(params, pre_params, adm_params, "Max parent age gap", 44, P);
		P->MaxChildAge = Params::get_int(params, pre_params, adm_params, "Max child age", 20, P);
		P->OneChildTwoPersProb = Params::get_double(params, pre_params, adm_params, "One Child Two Pers Prob", 0.08, P);
		P->TwoChildThreePersProb = Params::get_double(params, pre_params, adm_params, "Two Child Three Pers Prob", 0.11, P);
		P->OnePersHouseProbOld = Params::get_double(params, pre_params, adm_params, "One Pers House Prob Old", 0.5, P);
		P->TwoPersHouseProbOld = Params::get_double(params, pre_params, adm_params, "Two Pers House Prob Old", 0.5, P);
		P->OnePersHouseProbYoung = Params::get_double(params, pre_params, adm_params, "One Pers House Prob Young", 0.23, P);
		P->TwoPersHouseProbYoung = Params::get_double(params, pre_params, adm_params, "Two Pers House Prob Young", 0.23, P);
		P->OneChildProbYoungestChildUnderFive = Params::get_double(params, pre_params, adm_params, "One Child Prob Youngest Child Under Five", 0.5, P);
		P->TwoChildrenProbYoungestUnderFive = Params::get_double(params, pre_params, adm_params, "Two Children Prob Youngest Under Five", 0.0, P);
		P->TwoChildrenProbYoungestUnderFive = Params::get_double(params, pre_params, adm_params, "Prob Youngest Child Under Five", 0, P);
		P->ZeroChildThreePersProb = Params::get_double(params, pre_params, adm_params, "Zero Child Three Pers Prob", 0.25, P);
		P->OneChildFourPersProb = Params::get_double(params, pre_params, adm_params, "One Child Four Pers Prob", 0.2, P);
		P->YoungAndSingleSlope = Params::get_double(params, pre_params, adm_params, "Young And Single Slope", 0.7, P);
		P->YoungAndSingle = Params::get_int(params, pre_params, adm_params, "Young And Single", 36, P);
		P->NoChildPersAge = Params::get_int(params, pre_params, adm_params, "No Child Pers Age", 44, P);
		P->OldPersAge = Params::get_int(params, pre_params, adm_params, "Old Pers Age", 60, P);
		P->ThreeChildFivePersProb = Params::get_double(params, pre_params, adm_params, "Three Child Five Pers Prob", 0.5, P);
		P->OlderGenGap = Params::get_int(params, pre_params, adm_params, "Older Gen Gap", 19, P);
	}

	if (P->DoOneGen != 0) P->DoOneGen = 1;
	P->ColourPeriod = 2000;
	P->MoveRestrRadius2 = P->MoveRestrRadius * P->MoveRestrRadius;
	P->SocDistRadius2 = P->SocDistRadius * P->SocDistRadius;
	P->VaccRadius2 = P->VaccRadius * P->VaccRadius;
	P->VaccMinRadius2 = P->VaccMinRadius * P->VaccMinRadius;
	P->TreatRadius2 = P->TreatRadius * P->TreatRadius;
	P->PlaceCloseRadius2 = P->PlaceCloseRadius * P->PlaceCloseRadius;
	P->KeyWorkerProphRadius2 = P->KeyWorkerProphRadius * P->KeyWorkerProphRadius;
	if (P->TreatRadius2 == 0) P->TreatRadius2 = -1;
	if (P->VaccRadius2 == 0) P->VaccRadius2 = -1;
	if (P->PlaceCloseRadius2 == 0) P->PlaceCloseRadius2 = -1;
	if (P->MoveRestrRadius2 == 0) P->MoveRestrRadius2 = -1;
	if (P->SocDistRadius2 == 0) P->SocDistRadius2 = -1;
	if (P->KeyWorkerProphRadius2 == 0) P->KeyWorkerProphRadius2 = -1;
	/*	if (P->TreatCellIncThresh < 1) P->TreatCellIncThresh = 1;
		if (P->CaseIsolation_CellIncThresh < 1) P->CaseIsolation_CellIncThresh = 1;
		if (P->DigitalContactTracing_CellIncThresh < 1) P->DigitalContactTracing_CellIncThresh = 1;
		if (P->HHQuar_CellIncThresh < 1) P->HHQuar_CellIncThresh = 1;
		if (P->MoveRestrCellIncThresh < 1) P->MoveRestrCellIncThresh = 1;
		if (P->PlaceCloseCellIncThresh < 1) P->PlaceCloseCellIncThresh = 1;
		if (P->KeyWorkerProphCellIncThresh < 1) P->KeyWorkerProphCellIncThresh = 1;
	*/

	//// Make unsigned short versions of various intervention variables. And scaled them by number of timesteps per day
	P->usHQuarantineHouseDuration = ((unsigned short int) (P->HQuarantineHouseDuration * P->TimeStepsPerDay));
	P->usVaccTimeToEfficacy = ((unsigned short int) (P->VaccTimeToEfficacy * P->TimeStepsPerDay));
	P->usVaccTimeEfficacySwitch = ((unsigned short int) (P->VaccTimeEfficacySwitch * P->TimeStepsPerDay));
	P->usCaseIsolationDelay = ((unsigned short int) (P->CaseIsolationDelay * P->TimeStepsPerDay));
	P->usCaseIsolationDuration = ((unsigned short int) (P->CaseIsolationDuration * P->TimeStepsPerDay));
	P->usCaseAbsenteeismDuration = ((unsigned short int) (P->CaseAbsenteeismDuration * P->TimeStepsPerDay));
	P->usCaseAbsenteeismDelay = ((unsigned short int) (P->CaseAbsenteeismDelay * P->TimeStepsPerDay));
	if (P->DoUTM_coords)
	{
		for (i = 0; i <= 1000; i++)
		{
			asin2sqx[i] = asin(sqrt(i / 1000.0));
			asin2sqx[i] = asin2sqx[i] * asin2sqx[i];
		}
		for (i = 0; i <= DEGREES_PER_TURN; i++)
		{
			t = PI * i / 180;
			sinx[i] = sin(t);
			cosx[i] = cos(t);
		}
	}
	if (P->R0scale != 1.0)
	{
		P->HouseholdTrans *= P->R0scale;
		if (P->FixLocalBeta) P->LocalBeta *= P->R0scale;
		P->R0 *= P->R0scale;
		for (int place_type = 0; place_type < P->PlaceTypeNum; place_type++) {
			P->PlaceTypeTrans[place_type] *= P->R0scale;
		}
		Files::xfprintf_stderr("Rescaled transmission coefficients by factor of %lg\n", P->R0scale);
	}
	Files::xfprintf_stderr("Parameters read\n");
	delete[] CountryNameBuf;
	delete[] AdunitListNamesBuf;
	delete[] CountryNames;
	delete[] AdunitListNames;
}

