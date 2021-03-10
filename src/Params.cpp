/** \file  Params.cpp
 *  \brief Support for parsing parameter files
 */

#include "Params.h"

ParamMap Params::read_params_map(const char* file)
{
  std::ifstream fin(file);  // The file to be read
  std::string buf;       // Buffer for a line of text

  ParamMap param_map;   // Resulting map of (key, value)
  std::string key;      // Current key being processed
  std::string value;    // Current value being built
  bool skipping = true;
  while (std::getline(fin, buf)) {    // Each line. (buffer includes new-line)
    buf.erase(0, buf.find_first_not_of(" \t\n\r"));     // Remove lead space
    buf.erase(buf.find_last_not_of(" \t\n\r") + 1);     // ... and trailing.

	if (skipping) {   // Currently, we're not doing anything interesting.
		if ((buf.length() > 1) && (buf.compare(0, 1, "[") == 0)) {   // We found a key
			skipping = false;
			key = buf.substr(1, buf.length() - 2);
		}
	}
	else
	{			// Not skipping...
		if ((buf.length() == 0) || (buf.compare(0, 1, "*") == 0) ||
			(buf.compare(0, 1, "=") == 0) || (buf.compare(0, 1, "[") == 0)) {
			value.erase(0, value.find_first_not_of(" \t\n\r"));
			value.erase(value.find_last_not_of(" \t\n\r") + 1);
			param_map.insert(ParamPair(key, value));
			value.clear();
			key.clear();
			if (buf.compare(0, 1, "[") == 0) {
				key = buf.substr(1, buf.length() - 2);
			}
			else
			{
				skipping = true;
			}
		}
		else {
			value.append(buf);
			value.append("\n");
		}
	}
  }

  if (!key.empty()) {
	  value.erase(0, value.find_first_not_of(" \t\n\r"));
	  value.erase(value.find_last_not_of(" \t\n\r") + 1);
      param_map.insert(ParamPair(key, value));
  }

  fin.close();
  return param_map;
}

/*****************************************************************************/

std::string Params::clp_overwrite(std::string value, Param* P) {
	if (value.at(0) != '#') {
	  return value;
  } else {
	  int clp_no = std::stoi(value.substr(1, std::string::npos));
		if ((clp_no < 0) || (clp_no > 99)) {
			ERR_CRITICAL_FMT("CLP %d is out of bounds reading parameters\n", clp_no);
		}
		return std::to_string(P->clP[clp_no]);
	}
}

/*****************************************************************************/

std::string Params::lookup_param(ParamMap &base, ParamMap &fallback,
                                 ParamMap &params, 
	                               std::string param_name, Param* P)
{

  ParamIter iter = params.find(param_name);
	if (iter != params.end()) {
		return Params::clp_overwrite(iter->second, P);
	}
	else {
		iter = fallback.find(param_name);
		if (iter != fallback.end()) {
			return Params::clp_overwrite(iter->second, P);
		}
		else if (base == fallback) {
		  return "NULL";
		} else {
			iter = base.find(param_name);
			if (iter != base.end()) {
				return Params::clp_overwrite(iter->second, P);
			}
			else
				return "NULL";
		}
	}
}

std::string Params::lookup_param(ParamMap &fallback, ParamMap &params,
                         std::string param_name, Param* P) {
  return Params::lookup_param(fallback, fallback, params, param_name, P);
}

/*****************************************************************************/

bool Params::param_found(ParamMap &base, ParamMap &fallback, ParamMap &params,
                         std::string param_name) {

	ParamIter iter = params.find(param_name);
	if (iter != params.end()) {
		return true;
	}
	else {
		iter = fallback.find(param_name);
		if (iter != fallback.end()) {
			return true;
		}
		else if (base == fallback) {
			return false;
		}
		else {
			iter = base.find(param_name);
			if (iter != base.end()) {
				return true;
			}
			else
				return false;
		}
	}
}

bool Params::param_found(ParamMap &fallback, ParamMap &params, std::string param_name) {
  return Params::param_found(fallback, fallback, params, param_name);
}

/*****************************************************************************************/

double Params::get_double(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                std::string param_name, double default_value,
                  bool err_on_missing, Param* P) {

  std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0) {
		std::string::size_type idx;
		return std::stod(str_value, &idx);
	}
	if (err_on_missing) {
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

int Params::get_int(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                     std::string param_name, int default_value,
	                     bool err_on_missing, Param* P) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);
	
	if (str_value.compare("NULL") != 0) {
		std::string::size_type idx;
		double dval = std::stod(str_value, &idx);
		int result;
		if (dval > (double) INT32_MAX) {
			Files::xfprintf_stderr("Warning: reducing int param %s (%s) to MAX_INT\n", param_name.c_str(), str_value.c_str());
			result = INT32_MAX;
		}
		else if (dval < (double)INT32_MIN) {
			Files::xfprintf_stderr("Warning: increasing int param %s (%s) to MIN_INT\n", param_name.c_str(), str_value.c_str());
			result = INT32_MIN;
		}
		else {
			result = std::stoi(str_value, &idx);
		}
		return result;
	}
	if (err_on_missing) {
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
	}
	else {
		return default_value;
	}
}

int Params::get_int(ParamMap &fallback, ParamMap &params, std::string param_name,
	                  int default_value, Param* P) {
	return Params::get_int(fallback, fallback, params, param_name, default_value,
		false, P);
}

int Params::get_int(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name,
	int default_value, Param* P) {
	return Params::get_int(base, fallback, params, param_name, default_value, false, P);
}

int Params::req_int(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                  std::string param_name, Param* P) {
	return Params::get_int(base, fallback, params, param_name, 0, true, P);
}

int Params::req_int(ParamMap &fallback, ParamMap &params,
	                  std::string param_name, Param* P) {
	return Params::req_int(fallback, fallback, params, param_name, P);
}

/***********************************************************************************/

void Params::get_double_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
                            std::string param_name, double* array, int expected,
                            double default_value, int default_size,
	                        bool err_on_missing, Param* P, bool force_fail) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if ((str_value.compare("NULL") != 0) && (!force_fail)) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int index = 0;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
			  if (index < expected) array[index] = std::stod(buffer, &idx);
			  index++;
			}
		}
		if (index != expected) {
			Files::xfprintf_stderr("Warning - Extra elements for %s (%d - only needed %d)\n",
				param_name.c_str(), index, expected);
		}
	}
	else {
		for (int i = 0; i < default_size; i++) array[i] = default_value;
		if (err_on_missing) {
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
		}
	}
}

void Params::get_double_vec(ParamMap &fallback, ParamMap &params,
	                          std::string param_name, double* array,
														int expected, double default_value,
														int default_size,
														Param* P) {

  Params::get_double_vec(fallback, fallback, params, param_name, array, expected,
	                       default_value, default_size, false, P, false);
}

void Params::req_double_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                          std::string param_name, double* array, int expected,
	                          Param* P) {
	Params::get_double_vec(base, fallback, params, param_name, array, expected,
	                       0, 0, true, P, false);
}

void Params::req_double_vec(ParamMap &fallback, ParamMap &params,
	                          std::string param_name, double* array, int expected,
	                          Param* P) {
	Params::get_double_vec(fallback, fallback, params, param_name, array, expected,
	                       0, 0, true, P, false);
}

void Params::get_double_vec_ff(bool force_fail, ParamMap &fallback, ParamMap &params,
                               std::string param_name, double* array, int expected,
															 double default_value, Param* P) {

  Params::get_double_vec(fallback, fallback, params, param_name, array, expected,
	                       default_value, expected, false, P, force_fail);
}

/***********************************************************************************/

void Params::get_int_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                       std::string param_name, int* array, int expected,
	                       int default_value, int default_size,
	                       bool err_on_missing, Param* P, bool force_fail) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if ((str_value.compare("NULL") != 0) && (!force_fail)) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int index = 0;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				double dval = std::stod(buffer, &idx);
				int result;
				if (dval > (double) INT32_MAX) {
					Files::xfprintf_stderr("Warning: reducing int param %s (%s) to MAX_INT\n", param_name.c_str(), buffer.c_str());
					result = INT32_MAX;
				}
				else if (dval < (double) INT32_MIN) {
					Files::xfprintf_stderr("Warning: increasing int param %s (%s) to MIN_INT\n", param_name.c_str(), buffer.c_str());
					result = INT32_MIN;
				}
				else {
					result = std::stoi(str_value, &idx);
				}

				array[index] = result;
				index++;
			}
		}
		if (index != expected) {
			ERR_CRITICAL_FMT("Expected %d elements, found %d for param %s\n",
				expected, index, param_name.c_str());
		}
	}
	else {
		for (int i = 0; i < default_size; i++) array[i] = default_value;
		if (err_on_missing) {
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
		}
	}
}

void Params::get_int_vec(ParamMap &fallback, ParamMap &params,
	                       std::string param_name, int* array,
	                       int expected, int default_value,
	                       int default_size, bool err_on_missing,
	                       Param* P, bool force_fail) {

	Params::get_int_vec(fallback, fallback, params, param_name, array, expected,
		                  default_value, default_size, err_on_missing, P, force_fail);
}

void Params::get_int_vec(ParamMap &fallback, ParamMap &params,
	                       std::string param_name, int* array,
	                       int expected, int default_value,
	                       int default_size, Param* P) {

	Params::get_int_vec(fallback, fallback, params, param_name, array, expected,
		default_value, default_size, false, P, false);
}

void Params::req_int_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                       std::string param_name, int* array, int expected,
	                       Param* P) {
	Params::get_int_vec(base, fallback, params, param_name, array, expected, 0, 0, false, P, false);
}

void Params::req_int_vec(ParamMap &fallback, ParamMap &params,
	                       std::string param_name, int* array, int expected,
	                       Param* P) {
	Params::get_int_vec(fallback, fallback, params, param_name, array, expected, 0, 0, false, P, false);
}

void Params::get_int_vec_ff(bool force_fail, ParamMap &fallback, ParamMap &params,
	                          std::string param_name, int* array, int expected,
	                          int default_value, Param* P) {

	Params::get_int_vec(fallback, fallback, params, param_name, array, expected,
		                  default_value, expected, false, P, force_fail);
}

/***********************************************************************************/

int Params::req_string_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                          std::string param_name, char** array, int expected,
	                          Param* P) {

	std::string str_value = Params::lookup_param(fallback, params, param_name, P);
	if (str_value.compare("NULL") != 0) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		int index = 0;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				array[index] = new char[buffer.length() + 1];
				strcpy(array[index], buffer.c_str());
				index++;
			}
		}
		return index;
	}
	else {
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
	}
}

int Params::req_string_vec(ParamMap &fallback, ParamMap &params,
	                          std::string param_name, char** array, int expected,
	                          Param* P) {
	return Params::req_string_vec(fallback, fallback, params, param_name, array, expected, P);
}

/***********************************************************************************/

void Params::get_double_matrix(ParamMap &base, ParamMap &fallback, ParamMap &params,
	                             std::string param_name, double** array, int sizex, int sizey,
	                             double default_value, bool err_on_missing, Param* P) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int x = 0;
		int y = 0;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				if ((y < sizey) && (x < sizex)) array[x][y] = std::stod(buffer, &idx);
				x++;
				if (x == sizex) {
					y++;
					x = 0;
				}
			}
		}
		if ((x >= sizex) || (y >= sizey)) {
			Files::xfprintf_stderr("Warning: more data available for %s than was read.\n", param_name.c_str());
		}
	}
	else {
		if (err_on_missing) {
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
		} else {
	    for (int x = 0; x < sizex; x++)
		    for (int y = 0; y < sizey; y++)
			    array[x][y] = default_value;
		}
	}
}

void Params::get_double_matrix(ParamMap &fallback, ParamMap &params,
	                             std::string param_name, double** array, int sizex,
	                             int sizey, double default_value, Param* P) {
	Params::get_double_matrix(fallback, fallback, params, param_name, array,
	                          sizex, sizey, default_value, false, P);
}

double** create_2d_double(int sizex, int sizey) {
	double** arr = new double* [sizex];
	for (int i = 0; i < sizex; i++)
		arr[i] = new double[sizey];
	return arr;
}

void Params::alloc_params(Param* P) {
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
