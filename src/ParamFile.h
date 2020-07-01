#ifndef COVIDSIM_PARAM_FILE_HPP_INCLUDED_
#define COVIDSIM_PARAM_FILE_HPP_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * Class to handle reading and extracting values out of input parameter files
 * (e.g. pre-param.txt, param.txt, admin-params.txt) in a clean, type-safe, and
 * concise way.
 */
class ParamReader
{
public:
	/**
	 * Constructor which will open up the parameter files specified and read
	 * each parameter-value pair into a map.
	 *
	 * Duplicate parameter resolution is determined by the file they came from:
	 *
	 * admin_file >> preparam_file >> param_file
	 *
	 * Each parameter-value pair can be extracted out by one of the extract_*()
	 * functions by passing a storage variable of the type expected.
	 *
	 * @note: The parameter file names passed into this function can be empty.
	 */
	ParamReader(std::string const& param_file, std::string const& preparam_file, std::string const& admin_file);

	/**
	 * Extract a single 32-bit integer value or assign the default to `output`.
	 */
	bool extract_int(std::string const& param, int32_t& output, int32_t default_value);

	/**
	 * Extract a single double value or assign the default to `output`.
	 */
	bool extract_double(std::string const& param, double& output, double default_value);

	/**
	 * Extract and assign a single 32-bit integer value to `output` or exit the program.
	 */
	void extract_int_or_exit(std::string const& param, int32_t& output);

	/**
	 * Extract and assign a single double value to `output` or exit the program.
	 */
	void extract_double_or_exit(std::string const& param, double& output);

	/**
	 * Extract multiple 32-bit integer values to `N` numbers in `output`. 
	 */
	bool extract_ints_no_default(std::string const& param, int* output, std::size_t N);

	/**
	 * Extract multiple double values to `N` numbers in `output`. 
	 */
	bool extract_doubles_no_default(std::string const& param, double* output, std::size_t N);

	/**
	 * Extract and assign multiple 32-bit integer values to `N` numbers in `output` or exit the program.
	 */
	void extract_ints_or_exit(std::string const& param, int32_t* output, std::size_t N);

	/**
	 * Extract and assign multiple double values to `N` numbers in `output` or exit the program.
	 */
	void extract_doubles_or_exit(std::string const& param, double* output, std::size_t N);

	/**
	 * Extract multiple 32-bit integer values or assign the default to `N` numbers in `output`.
	 */
	bool extract_ints(std::string const& param, int32_t* output, std::size_t N, int32_t default_value);

	/**
	 * Extract multiple double values or assign the default to `N` numbers in `output`.
	 */
	bool extract_doubles(std::string const& param, double* output, std::size_t N, double default_value);

	/**
	 * Conditionally extract multiple 32-bit integer values or assign the default to `N` numbers in `output`.
	 */
	bool cond_extract_ints(bool conditional, std::string const& param, int32_t* output, std::size_t N, int32_t default_value);

	/**
	 * Conditionally extract multiple double values or assign the default to `N` numbers in `output`.
	 */
	bool cond_extract_doubles(bool conditional, std::string const& param, double* output, std::size_t N, double default_value);

	/**
	 * Extract exactly `N` strings into `output`, return false otherwise.
	 */
	bool extract_strings_no_default(std::string const& param, std::vector<std::string>& output, std::size_t N);
	
	/**
	 * Extract exactly `N` strings and add them to `output` or exit the program. 
	 */
	void extract_strings_or_exit(std::string const& param, std::vector<std::string>& output, std::size_t N);

	/**
	 * Extract a matrix of strings with `num_cols` columns per row and add them to `output` or exit the program. 
	 */
	void extract_string_matrix_or_exit(std::string const& param, std::vector<std::vector<std::string>>& output, std::size_t num_cols);

private:
	/**
	 * Checks if the `param` specified exists in the `m_param_value_map`.
	 */
	bool exists(std::string const& param);

	/**
	 * Open and reads in argument-value pairs from `param_file` into `m_param_value_map`.
	 */
	void parse_param_file(std::string const& param_file);

	std::unordered_map<std::string, std::string> m_param_value_map;
};

#endif // COVIDSIM_PARAM_FILE_HPP_INCLUDED_
