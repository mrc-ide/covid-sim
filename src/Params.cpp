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
  
  while (std::getline(fin, buf)) {    // Each line. (buffer includes new-line)

    buf.erase(0, buf.find_first_not_of(" \t\n\r"));     // Remove lead space
    buf.erase(buf.find_last_not_of(" \t\n\r") + 1);     // ... and trailing.

    if ((buf.length() == 0) || (buf.compare(0, 1, "[") == 0) ||
        (buf.compare(0, 1, "=") == 0) || (buf.compare(0, 1, "(") == 0)) {

      // Detected an empty line, or comment, or start of new key.
      // In any of these cases, if we were building a key-value before,
      // it is now finished and can be stored. 
      
      if (!key.empty()) {
        param_map.insert(ParamPair(key, value));
        key.clear();
        value.clear();
      }

      // If this is a new key, set it up.

      if (buf.compare(0, 1, "[") == 0)
        key = buf.substr(1, buf.length() - 2);

    } else {                  // It's a value (or part of one)
      if (!value.empty()) {   // If value already contains something, 
        value.append("\n");   // then separate with new-line in the map.
      }
      value.append(buf);
    }
  }

  // Deal with any straggler at the end of the file that didn't get terminated

  if (!key.empty()) {
    param_map.insert(ParamPair(key, value));
  }

  fin.close();
  return param_map;
}

/*****************************************************************************/

std::string Params::clp_overwrite(std::string value, Param P) {
	if (value.at(0) != '#') {
	  return value;
  } else {
	  int clp_no = std::stoi(value.substr(1, std::string::npos));
		if ((clp_no < 0) || (clp_no > 99)) {
			ERR_CRITICAL_FMT("CLP %d is out of bounds reading parameters\n", clp_no);
		}
		return std::to_string(P.clP[clp_no]);
	}
}

/*****************************************************************************/

std::string Params::lookup_param(ParamMap base, ParamMap fallback,
                                 ParamMap params, 
	                               std::string param_name, Param P)
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
		else if (base.empty()) {
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

std::string Params::lookup_param(ParamMap fallback, ParamMap params,
                         std::string param_name, Param P) {
  return Params::lookup_param(ParamMap(), fallback, params, param_name, P);
}

/*****************************************************************************/

bool Params::param_found(ParamMap base, ParamMap fallback, ParamMap params,
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
		else if (base.empty()) {
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

bool Params::param_found(ParamMap fallback, ParamMap params, std::string param_name) {
  return Params::param_found(ParamMap(), fallback, params, param_name);
}

/*****************************************************************************************/

double Params::get_double(ParamMap base, ParamMap fallback, ParamMap params,
	                std::string param_name, double default_value,
                  bool err_on_missing, Param P) {

  std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0) {
		std::string::size_type idx;
		return std::stod(str_value, &idx);
	}
	if (err_on_missing) {
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name);
	} else {
	  return default_value;
	}
}

double Params::get_double(ParamMap fallback, ParamMap params, std::string param_name,
                          double default_value, bool err_on_missing, Param P) {
	return Params::get_double(ParamMap(), fallback, params, param_name, default_value,
	                          err_on_missing, P);
}

double Params::req_double(ParamMap base, ParamMap fallback, ParamMap params,
                           std::string param_name, Param P) {
	return Params::get_double(base, fallback, params, param_name, 0, true, P);
}

double Params::req_double(ParamMap fallback, ParamMap params,
	std::string param_name, Param P) {
	return Params::req_double(ParamMap(), fallback, params, param_name, P);
}


/*****************************************************************************************/

int Params::get_int(ParamMap base, ParamMap fallback, ParamMap params,
	                     std::string param_name, int default_value,
	                     bool err_on_missing, Param P) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0) {
		std::string::size_type idx;
		return std::stoi(str_value, &idx);
	}
	if (err_on_missing) {
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name);
	}
	else {
		return default_value;
	}
}

int Params::get_int(ParamMap fallback, ParamMap params, std::string param_name,
	                  int default_value, bool err_on_missing, Param P) {
	return Params::get_int(ParamMap(), fallback, params, param_name, default_value,
		err_on_missing, P);
}

int Params::req_int(ParamMap base, ParamMap fallback, ParamMap params,
	                  std::string param_name, Param P) {
	return Params::get_int(base, fallback, params, param_name, 0, true, P);
}

int Params::req_int(ParamMap fallback, ParamMap params,
	                  std::string param_name, Param P) {
	return Params::req_int(ParamMap(), fallback, params, param_name, P);
}

/***********************************************************************************/

void Params::get_double_vec(ParamMap base, ParamMap fallback, ParamMap params,
                            std::string param_name, double* array, int expected,
														double default_value, int default_size,
	                          bool err_on_missing, Param P, bool force_fail = false) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if ((str_value.compare("NULL") != 0) && (!force_fail)) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int index = -1;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
			  index++;
			  array[index] = std::stod(buffer, &idx);
			}
		}
		if (index != expected) {
			ERR_CRITICAL_FMT("Expected %d elements, found %d for param %s\n",
			expected, index, param_name);
		}
	}
	else {
		for (int i = 0; i < default_size; i++) array[i] = default_value;
		if (err_on_missing) {
			ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name);
		}
	}
}

void Params::get_double_vec(ParamMap fallback, ParamMap params,
	                          std::string param_name, double* array,
														int expected, double default_value,
														int default_size, bool err_on_missing,
														Param P, bool force_fail = false) {

  Params::get_double_vec(ParamMap(), fallback, params, param_name, array, expected,
	                       default_value, default_size, err_on_missing, P, force_fail);
}

void Params::req_double_vec(ParamMap base, ParamMap fallback, ParamMap params,
	                          std::string param_name, double* array, int expected,
	                          Param P) {
	Params::get_double_vec(base, fallback, params, param_name, array, expected,
	                       0, 0, false, P);
}

void Params::req_double_vec(ParamMap fallback, ParamMap params,
	                          std::string param_name, double* array, int expected,
	                          Param P) {
	Params::get_double_vec(ParamMap(), fallback, params, param_name, array, expected,
	                       0, 0, false, P);
}

/***********************************************************************************/

void Params::get_int_vec(ParamMap base, ParamMap fallback, ParamMap params,
	                       std::string param_name, int* array, int expected,
	                       int default_value, int default_size,
	                       bool err_on_missing, Param P, bool force_fail = false) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if ((str_value.compare("NULL") != 0) && (!force_fail)) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int index = -1;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				index++;
				array[index] = std::stoi(buffer, &idx);
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

void Params::get_int_vec(ParamMap fallback, ParamMap params,
	                       std::string param_name, int* array,
	                       int expected, int default_value,
	                       int default_size, bool err_on_missing,
	                       Param P, bool force_fail = false) {

	Params::get_int_vec(ParamMap(), fallback, params, param_name, array, expected,
		                  default_value, default_size, err_on_missing, P, force_fail);
}

void Params::req_int_vec(ParamMap base, ParamMap fallback, ParamMap params,
	                       std::string param_name, int* array, int expected,
	                       Param P) {
	Params::get_int_vec(base, fallback, params, param_name, array, expected, 0, 0, false, P);
}

void Params::req_int_vec(ParamMap fallback, ParamMap params,
	                       std::string param_name, int* array, int expected,
	                       Param P) {
	Params::get_int_vec(ParamMap(), fallback, params, param_name, array, expected, 0, 0, false, P);
}

/***********************************************************************************/

void Params::req_string_vec(ParamMap base, ParamMap fallback, ParamMap params,
	                          std::string param_name, char** array, int expected,
	                          Param P) {

	std::string str_value = Params::lookup_param(fallback, params, param_name, P);
	if (str_value.compare("NULL") != 0) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		int index = -1;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				index++;
				array[index] = new char[buffer.length() + 1];
				strcpy(array[index], buffer.c_str());
			}
		}
		if (index != expected) {
			ERR_CRITICAL_FMT("Expected %d elements, found %d for param %s\n",
				expected, index, param_name.c_str());
		}
	}
}

void Params::req_string_vec(ParamMap fallback, ParamMap params,
	                          std::string param_name, char** array, int expected,
	                          Param P) {
	Params::req_string_vec(ParamMap(), fallback, params, param_name, array, expected, P);
}

/***********************************************************************************/

void Params::get_int_matrix(ParamMap base, ParamMap fallback, ParamMap params,
	                          std::string param_name, int** array, int sizex, int sizey,
										        Param P) {
	
	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int x = -1;
		int y = -1;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				x++;
				if (x == sizex) x = 0;
				if (x == 0) y++;
				array[x][y] = std::stoi(buffer, &idx);
			}
		}
		if ((x != (sizex - 1)) || (y != (sizey - 1))) {
			ERR_CRITICAL_FMT("Matrix read of %s expected to end on (%d,%d) - actually (%d,%d)\n", param_name.c_str(), sizex - 1, sizey - 1, x, y);
		}
	}
	else {
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
  }
}

void Params::get_int_matrix(ParamMap fallback, ParamMap params,
	                          std::string param_name, int** array, int sizex,
	                          int sizey, Param P) {

  Params::get_int_matrix(ParamMap(), fallback, params, param_name, array,  sizex, sizey, P);
}

/***********************************************************************************/

void Params::get_double_matrix(ParamMap base, ParamMap fallback, ParamMap params,
	                             std::string param_name, double** array, int sizex, int sizey,
	Param P) {

	std::string str_value = Params::lookup_param(base, fallback, params, param_name, P);

	if (str_value.compare("NULL") != 0) {
		std::stringstream str_stream(str_value);
		std::string buffer;
		std::string::size_type idx;
		int x = -1;
		int y = -1;
		while (str_stream >> buffer) {
			if (!buffer.empty()) {
				x++;
				if (x == sizex) x = 0;
				if (x == 0) y++;
				array[x][y] = std::stod(buffer, &idx);
			}
		}
		if ((x != (sizex - 1)) || (y != (sizey - 1))) {
			ERR_CRITICAL_FMT("Matrix read of %s expected to end on (%d,%d) - actually (%d,%d)\n", param_name.c_str(), sizex - 1, sizey - 1, x, y);
		}
	}
	else {
		ERR_CRITICAL_FMT("Required Parameter %s not found\n", param_name.c_str());
	}
}

void Params::get_double_matrix(ParamMap fallback, ParamMap params,
	std::string param_name, double** array, int sizex,
	int sizey, Param P) {
	Params::get_double_matrix(ParamMap(), fallback, params, param_name, array, sizex, sizey, P);
}
