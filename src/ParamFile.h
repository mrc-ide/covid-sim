#ifndef COVIDSIM_PARAM_FILE_HPP_INCLUDED_
#define COVIDSIM_PARAM_FILE_HPP_INCLUDED_

#include <cstddef>
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
	 * Extract a single numeric value or assign the default to `output`.
	 */
	template<typename T>
	bool extract(std::string const& param, T& output, T default_value);

	/**
	 * Extract and assign a single numeric value to `output` or exit the program.
	 */
	template<typename T>
	void extract_or_exit(std::string const& param, T& output);

	/**
	 * Extract multiple numeric values to `N` numbers in `output`, return false otherwise.
	 */
	template<typename T>
	bool extract_multiple_no_default(std::string const& param, T* output, std::size_t N);

	/**
	 * Extract and assign multiple numeric values to `N` numbers in `output` or exit the program.
	 */
	template<typename T>
	void extract_multiple_or_exit(std::string const& param, T* output, std::size_t N);

	/**
	 * Extract multiple numeric values or assign the default to `N` numbers in `output`.
	 */
	template<typename T>
	bool extract_multiple(std::string const& param, T* output, std::size_t N, T default_value);

	/**
	 * Conditionally extract multiple numeric values or assign the default to `N` numbers in `output`.
	 */
	template<typename T>
	bool cond_extract_multiple(bool conditional, std::string const& param, T* output, std::size_t N, T default_value);

	/**
	 * Extract exactly `N` strings into `output`, return false otherwise.
	 */
	bool extract_multiple_strings_no_default(std::string const& param, std::vector<std::string>& output, std::size_t N);
	
	/**
	 * Extract exactly `N` strings and add them to `output` or exit the program. 
	 */
	void extract_multiple_strings_or_exit(std::string const& param, std::vector<std::string>& output, std::size_t N);

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

	/**
	 * Extracts the next value from `stream` into the output variable using the
	 * numeric limits of the `output` variable's type.
	 */
	template<typename T>
	bool raw_extract(std::istringstream& stream, T& output);

	std::unordered_map<std::string, std::string> m_param_value_map;
};

#endif // COVIDSIM_PARAM_FILE_HPP_INCLUDED_
