#ifndef COVIDSIM_PARSERS_H_INCLUDED_
#define COVIDSIM_PARSERS_H_INCLUDED_

#include <string>

/**
 * Parses and checks if the input string is a readable file on-disk.
 */
void parse_read_file(std::string const& input, std::string& output);

/**
 *  Parses and checks if the input string is a writable directory.
 */
void parse_write_dir(std::string const& input, std::string& output);

/**
 *  handles general string.
 */
void parse_string(std::string const& input, std::string& output);

/**
 * Parses and checks if the input string is an integer.
 */
bool parse_integer_no_default(std::string const& input, int& output);
bool parse_integer(std::string const& input, int& output, int default_value);
void parse_integer_or_exit(std::string const& input, int& output);

/**
 * Parses and checks if the input string is an double.
 */
bool parse_double_no_default(std::string const& input, double& output);
bool parse_double(std::string const& input, double& output, double default_value);
void parse_double_or_exit(std::string const& input, double& output);

#endif