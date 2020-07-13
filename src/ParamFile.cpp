#include <algorithm>
#include <cctype>
#include <cstdint>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

#include "Error.h"
#include "InverseCdf.h"
#include "ParamFile.h"
#include "Parsers.h"

// ltrim(), rtrim(), and trim() are from a highly up voted answer on StackOverflow
// https://stackoverflow.com/a/217605

// trim from start (in place)
static void ltrim(std::string& s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
			[](int ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
static void rtrim(std::string& s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(),
			[](int ch) { return !std::isspace(ch); }).base(), s.end());
}

// trim from both ends (in place)
static void trim(std::string& s)
{
	ltrim(s);
	rtrim(s);
}

ParamReader::ParamReader(std::string const& param_file, std::string const& preparam_file, std::string const& admin_file)
{
	// duplicate params from different files take priority over one another:
	//	 admin_file >> preparam_file >> param_file
	parse_param_file(param_file);
	parse_param_file(preparam_file);
	parse_param_file(admin_file);
}

bool ParamReader::extract_int(std::string const& param, int32_t& output, int32_t default_value)
{
	if (!exists(param))
	{
		std::cout << "Using default value: " << default_value << std::endl;
		output = default_value;
		return false;
	}

	parse_integer(param_value_map_[param], output, default_value);
	return true;
}

bool ParamReader::extract_bool(std::string const& param, bool& output, bool default_value)
{
	if (!exists(param))
	{
		std::cout << "Using default value: " << default_value << std::endl;
		output = default_value;
		return false;
	}

	parse_bool(param_value_map_[param], output, default_value);
	return true;
}

bool ParamReader::extract_double(std::string const& param, double& output, double default_value)
{
	if (!exists(param))
	{
		std::cout << "Using default value: " << default_value << std::endl;
		output = default_value;
		return false;
	}

	parse_double(param_value_map_[param], output, default_value);
	return true;
}

void ParamReader::extract_int_or_exit(std::string const& param, int& output)
{
	if (!exists(param))
		std::exit(1);

	parse_integer_or_exit(param_value_map_[param], output);
}

void ParamReader::extract_double_or_exit(std::string const& param, double& output)
{
	if (!exists(param))
		std::exit(1);

	parse_double_or_exit(param_value_map_[param], output);
}

bool ParamReader::extract_ints_no_default(std::string const& param, int32_t* output, std::size_t N)
{
	if (!exists(param))
		return false;

	std::istringstream iss(param_value_map_[param]);
	for (auto i = 0; i < N; i++)
	{
		std::string token;
		iss >> token;
		if (iss.fail())
		{
			std::cerr << "ERROR: Only got " << i << " out of expected " << N << " integers in line ("
					<< iss.str() << ") for parameter \"" << param << "\"" << std::endl;
			return false;
		}
		if (!parse_integer_no_default(token, output[i]))
		{
			std::cerr << "Error parsing integer for param \"" << param << "\"" << std::endl;
			return false;
		}
			
	}
	return true;
}

bool ParamReader::extract_doubles_no_default(std::string const& param, double* output, std::size_t N)
{
	if (!exists(param))
		return false;

	std::istringstream iss(param_value_map_[param]);
	for (auto i = 0; i < N; i++)
	{
		std::string token;
		iss >> token;
		if (iss.fail())
		{
			std::cerr << "ERROR: Only got " << i << " out of expected " << N << " doubles in line ("
					<< iss.str() << ") for parameter \"" << param << "\"" << std::endl;
			return false;
		}
		if (!parse_double_no_default(token, output[i]))
		{
			std::cerr << "Error parsing double for param \"" << param << "\"" << std::endl;
			return false;
		}
	}
	return true;
}

void ParamReader::extract_ints_or_exit(std::string const& param, int32_t* output, std::size_t N)
{
	if (extract_ints_no_default(param, output, N) == false)
		std::exit(1);
}

void ParamReader::extract_doubles_or_exit(std::string const& param, double* output, std::size_t N)
{
	if (extract_doubles_no_default(param, output, N) == false)
		std::exit(1);
}

bool ParamReader::extract_ints(std::string const& param, int32_t* output, std::size_t N, int32_t default_value)
{
	if (extract_ints_no_default(param, output, N) == false)
	{
		std::cerr << "Using default value: " << default_value << std::endl;
		std::fill_n(output, N, default_value);
		return false;
	}
	return true;
}

bool ParamReader::extract_doubles(std::string const& param, double* output, std::size_t N, double default_value)
{
	if (extract_doubles_no_default(param, output, N) == false)
	{
		std::cerr << "Using default value: " << default_value << std::endl;
		std::fill_n(output, N, default_value);
		return false;
	}
	return true;
}

bool ParamReader::cond_extract_ints(bool conditional, std::string const& param, int32_t* output, std::size_t N, int32_t default_value)
{
	if (conditional)
		return extract_ints(param, output, N, default_value);
	std::fill_n(output, N, default_value);
	return true;
}

bool ParamReader::cond_extract_doubles(bool conditional, std::string const& param, double* output, std::size_t N, double default_value)
{
	if (conditional)
		return extract_doubles(param, output, N, default_value);
	std::fill_n(output, N, default_value);
	return true;
}

bool ParamReader::extract_strings_no_default(std::string const& param, std::vector<std::string>& output, std::size_t N)
{
	if (!exists(param))
		return false;

	std::istringstream iss(param_value_map_[param]);
	std::string token;
	for (auto i = 0; i < N; i++)
	{
		iss >> token;
		if (iss.fail())
		{
			std::cerr << "ERROR: Only got " << i << " out of expected " << N << " strings in line ("
					<< iss.str() << ") for " << param << std::endl;
			return false;
		}
		output.push_back(std::move(token));
	}
	return true;
}

void ParamReader::extract_strings_or_exit(std::string const& param, std::vector<std::string>& output, std::size_t N)
{
	if (!extract_strings_no_default(param, output, N))
		std::exit(1);
}

void ParamReader::extract_string_matrix_or_exit(std::string const& param, std::vector<std::vector<std::string>>& output, std::size_t num_cols)
{
	if (!exists(param))
		std::exit(1);

	// iterate over every row of the matrix
	std::istringstream matrix(param_value_map_[param]);
	std::istringstream line_stream;
	for (std::string line; std::getline(matrix, line);)
	{
		// iterate over every token in each row and add it in a vector
		line_stream.str(line);
		line_stream.clear();
		std::vector<std::string> new_row;
		for (auto i = 0; i < num_cols; i++)
		{
			std::string token;
			line_stream >> token;
			if (line_stream.fail())
			{
				std::cerr << "ERROR: Only got " << i << " out of " << num_cols << " columns in line (" << line
						  << ") for parameter " << param << std::endl;
				std::exit(1);
			}
			new_row.push_back(std::move(token));
		}
		output.push_back(std::move(new_row));
	}
}

void ParamReader::extract_inverse_cdf(std::string const& param, InverseCdf& inverse_cdf, double start_value)
{
	if (!extract_doubles_no_default(param, inverse_cdf.get_values(), CDF_RES + 1))
	{
		inverse_cdf.set_neg_log(start_value);
	}
	inverse_cdf.assign_exponent();
}

bool ParamReader::exists(std::string const& param)
{
	auto it = param_value_map_.find(param);
	if (it == param_value_map_.cend())
	{
		std::cout << "ERROR: Parameter \"" << param << "\" not found." << std::endl;
		return false;
	}
	return true;
}

void ParamReader::parse_param_file(std::string const& param_file)
{
	std::ifstream param_stream(param_file);
	if (!param_stream)
	{
		std::cerr << "Could not open the specified param file: " << param_file << std::endl;
		return;
	}

	std::string param_line;
	while (std::getline(param_stream, param_line))
	{
		// strip all whitespace characters from the line
		trim(param_line);

		// skip empty lines
		if (param_line.empty())
			continue;

		// check if this line is a parameter
		if (param_line.front() == '[' && param_line.back() == ']') {
			// remove the brackets
			param_line.erase(param_line.begin());
			param_line.pop_back();
			// try to read the value of this parameter into a string
			std::string value_line;
			std::string tmp_line;
			// read every line after the parameter name until there's no more or a break
			while (std::getline(param_stream, tmp_line))
			{
				trim(tmp_line);
				// break if the next line is empty or doesn't start with a letter, number, or '-' character
				if (tmp_line.empty() || !(std::isalnum(tmp_line[0]) || tmp_line[0] == '-'))
					break;
				if (!value_line.empty())
					value_line.push_back('\n');
				value_line.append(tmp_line);
			}
			// if there is no value found for the parameter, print an error and continue
			if (value_line.empty())
			{
				std::cerr << "ERROR: Value is missing for parameter: " << param_line << std::endl;
				continue;
			}

			// store the parameter name and value for later use by ReadParams()
			auto it = param_value_map_.find(param_line);
			if (it != param_value_map_.cend())
			{
				std::cerr << "OVERRIDE: Param file (\"" << param_file << "\") is updating Key:["
						  << param_line << "]" << " Value:[" << value_line << "]" << std::endl;
				it->second = value_line;
			}
			else
			{
				param_value_map_.emplace(param_line, value_line);
			}
		}
	}
}
