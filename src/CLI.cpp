#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

#include "CLI.h"
#include "Error.h"
#include "Param.h"
#include "Parsers.h"

void CmdLineArgs::add_custom_option(std::string const& option, ParserFn func, std::string const& doc)
{
	if (option_map_.find(option) != option_map_.cend())
	{
		ERR_CRITICAL_FMT("Duplicate option specified %s\n", option.c_str());
	}
	option_map_.emplace(option, func);
	doc_map_.emplace(option, doc);
}

void CmdLineArgs::add_double_option(std::string const& option, double& output, std::string const& doc)
{
	add_custom_option(option, std::bind(parse_double, std::placeholders::_1, std::ref(output)), doc);
}

void CmdLineArgs::add_integer_option(std::string const& option, int& output, std::string const& doc)
{
	add_custom_option(option, std::bind(parse_integer, std::placeholders::_1, std::ref(output)), doc);
}

void CmdLineArgs::add_string_option(std::string const& option, StringParserFn func, std::string& output, std::string const& doc)
{
	add_custom_option(option, std::bind(func, std::placeholders::_1, std::ref(output)), doc);
}

void CmdLineArgs::parse(int argc, char* argv[], Param& P)
{
	// Detect if the user wants to print out the full help output
	if (argc >= 2)
	{
		std::string first(argv[1]);
		if (first.length() == 2 && first.compare("/H") == 0)
		{
			print_detailed_help_and_exit();
		}
	}
	if (argc < 7)
	{
		std::cerr << "Minimum number of arguments not met. Expected 6 got "
				  << (argc - 1) << "\n" << std::endl;
		print_help_and_exit();
	}

	// Get seeds.
	int i = argc - 4;
	parse_integer(argv[i], P.setupSeed1);
	parse_integer(argv[i+1], P.setupSeed2);
	parse_integer(argv[i+2], P.runSeed1);
	parse_integer(argv[i+3], P.runSeed2);

	for (i = 1; i < argc - 4; i++)
	{
		std::string argument(argv[i]);
		if (argument[0] != '/')
		{
			std::cerr << "Argument \"" << argument << "\" does not start with '/'" << std::endl;
			print_help_and_exit();
		}
		auto pos = argument.find_first_of(':');
		if (pos == std::string::npos)
		{
			std::cerr << "Argument \"" << argument << "\" did not have a ':' character" << std::endl;
			print_help_and_exit();
		}
		std::string option = argument.substr(1, pos - 1); // everything after / and up to :
		std::string value = argument.substr(pos + 1); // everything after :

		auto it = option_map_.find(option);
		if (it == option_map_.cend())
		{
			std::cout << "Unknown option \"" << option << "\". Skipping" << std::endl;
			continue;
		}

		it->second(value);
	}
}

static const std::string USAGE =
R"(SetupSeed1 SetupSeed2 RunSeed1 RunSeed2

All '/' arguments are followed by a colon ':' and then the value.
More detailed information can be found in docs/inputs-and-outputs.md.
)";

static const std::string DETAILED_USAGE =
R"(Named arguments:

		SetupSeed*	  RNG seeds when initializing
		RunSeed*		RNG seeds when running

Optional arguments:
)";

void CmdLineArgs::print_help()
{
	// Wes: Temporary hack to reduce the number of /CLP: on the help text.
	int clp_flag = 0;
	const std::string CLP = "CLP";
	std::stringstream ss;

	ss << "CovidSim";
	for (auto const& it : option_map_)
	{
		std::string option = it.first;
		std::string option_first_3 = option.substr(0, 3);

		if (option_first_3.compare(CLP) == 0) {
			if (clp_flag == 0) {
				ss << " [/CLP:00] .. [/CLP:99]";
				clp_flag = 1;
			}
		}
		else {
			ss << " [/" << option << ']';
		}
	}
	std::cerr << ss.str() << ' ' << USAGE << std::endl;
}

void CmdLineArgs::print_help_and_exit()
{
	print_help();
	std::cerr << "Use the '/H' argument to get a detailed listing of arguments." << std::endl;
	std::exit(1);
}

void CmdLineArgs::print_detailed_help()
{
	print_help();
	std::stringstream ss;
	for (auto const& it : doc_map_)
	{
		ss << '\t' << it.first << '\t' << it.second << '\n';
	}
	std::cerr << DETAILED_USAGE << '\n' << ss.str() << std::endl;
}

void CmdLineArgs::print_detailed_help_and_exit()
{
	print_detailed_help();
	std::exit(1);
}
