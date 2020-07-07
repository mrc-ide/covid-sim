#ifndef COVIDSIM_CLI_H_INCLUDED_
#define COVIDSIM_CLI_H_INCLUDED_

#include <functional>
#include <map>
#include <string>

// only a forward-declartion, no need to pull in all of Param.h in this header
struct Param;

// The parse_* functions are kept outside of the CmdLineArgs class because they
// can be used to parse any string input and are not specific to CLI arguments.
// e.g. they may eventually be useful for parsing from the (pre)param files

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
void parse_integer(std::string const& input, int& output);

/**
 * Parses and checks if the input string is an double.
 */
void parse_double(std::string const& input, double& output);

class CmdLineArgs {
public:
	// Function prototype for a generic parser function
	using ParserFn = std::function<void(std::string const&)>;
	// Function prototype for a string parser function
	using StringParserFn = std::function<void(std::string const&, std::string&)>;

	/**
	 * Use this function when adding a new option to the CLI that needs a custom
	 * callback function (e.g. it does more than read a string or integer).
	 *
	 * This will insert the custom callback function into a map of other options.
	 */
	void add_custom_option(std::string const& option, ParserFn func, std::string const& doc);

	/**
	 * Use this function when adding a new double option to the CLI.
	 *
	 * This will bind the output variable to the parse_double() function and
	 * gets inserted into a map of other options. This provides a strong cohesion
	 * between an option name (i.e. 'R') with its variable (i.e. 'P.R0scale')
	 */
	void add_double_option(std::string const& option, double& output, std::string const& doc);

	/**
	 * Use this function when adding a new integral option to the CLI.
	 *
	 * This will bind the output variable to the parse_integer() function and
	 * gets inserted into a map of other options. This provides a strong cohesion
	 * between an option name (i.e. 'c') with its variable (i.e. 'P.MaxNumThreads')
	 */
	void add_integer_option(std::string const& option, int& output, std::string const& doc);

	/**
	 * Use this function when adding a new string option to the CLI.
	 *
	 * This will bind the output variable specified to a parser function which
	 * then gets inserted into a map of other options. This provides a strong
	 * cohesion between an option name (i.e. 'P') with a C++ variable (i.e. 'ParamFile')
	 */
	void add_string_option(std::string const& option, StringParserFn func, std::string& output, std::string const& doc);

	/**
	 * Call this function once all add_option() calls have been made to process
	 * the arguments passed in from the command-line.
	 */
	void parse(int argc, char* argv[], Param& P);

	void print_help();
	void print_help_and_exit();

	void print_detailed_help();
	void print_detailed_help_and_exit();

private:
	std::map<std::string, ParserFn> option_map_;
	std::map<std::string, std::string> doc_map_;
};

#endif
