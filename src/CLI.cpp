#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <sstream>

#include "CLI.hpp"
#include "Param.h"

void parse_read_file(std::string const& input, std::string& output) {
    // check to see if the file exists and error out if it doesn't
    if (static_cast<bool>(std::ifstream(input)) == false) {
        std::cerr << input << " is not a file" << std::endl;
        PrintHelpAndExit();
    }
    output = input;
}

void parse_integer(std::string const& input, int& output) {
	std::istringstream iss(input);
	iss >> output;
    if (iss.fail()) {
		if (output == std::numeric_limits<int>::max()) {
			std::cerr << "OVERFLOW: detected a number larger than " << std::numeric_limits<int>::max()
					  << " when parsing an int: " << iss.str() << std::endl;
		}
        else if (output == std::numeric_limits<int>::min()) {
            std::cerr << "UNDERFLOW: detected a number smaller than " << std::numeric_limits<int>::min()
                      << " when parsing an int: " << iss.str() << std::endl;
        }
		else {
			std::cerr << "ERROR: Expected int, got " << iss.str() << std::endl;
		}
		PrintHelpAndExit();
	}
}

void parse_long(std::string const& input, long& output) {
	std::istringstream iss(input);
	iss >> output;
	if (iss.fail()) {
		if (output == std::numeric_limits<long>::max()) {
			std::cerr << "OVERFLOW: detected a number larger than " << std::numeric_limits<long>::max()
					  << " when parsing a long: " << iss.str() << std::endl;
		}
        else if (output == std::numeric_limits<long>::min()) {
            std::cerr << "UNDERFLOW: detected a number smaller than " << std::numeric_limits<long>::min()
                      << " when parsing a long: " << iss.str() << std::endl;
        }
		else {
			std::cerr << "ERROR: Expected long, got " << iss.str() << std::endl;
		}
		PrintHelpAndExit();
	}
}

template<class T>
void CmdLineArgs::add_option(std::string const&& option, std::function<void(std::string const&, T&)> func, T& output) {
    if (m_option_map.find(option) != m_option_map.cend()) {
        std::cerr << "Duplicate option specified " << option << ", ignoring..." << std::endl;
        return;
    }

    ParserFn parser_func = std::bind(func, std::placeholders::_1, std::ref(output));
    m_option_map.emplace(option, parser_func);
}

// Explicit template instantiations for the linker
// https://stackoverflow.com/questions/2351148/explicit-template-instantiation-when-is-it-used
template void CmdLineArgs::add_option<int>(std::string const&&, std::function<void(std::string const&, int&)>, int&);
template void CmdLineArgs::add_option<long>(std::string const&&, std::function<void(std::string const&, long&)>, long&);
template void CmdLineArgs::add_option<std::string>(std::string const&&, std::function<void(std::string const&, std::string&)>, std::string&);

int CmdLineArgs::parse(int argc, char* argv[], Param& P) {
    // Detect if the user wants to print out the full help output
    if (argc >= 2) {
        std::string first(argv[1]);
        if (first.length() == 2 && first.compare("/H") == 0) {
            PrintDetailedHelpAndExit();
        }
    }
    if (argc < 7) {
        std::cerr << "Minimum number of arguments not met. Expected 6 got "
                << (argc - 1) << "\n" << std::endl;
        PrintHelpAndExit();
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
		if (argument[0] != '/') {
			std::cerr << "Argument \"" << argument << "\" does not start with '/'" << std::endl;
			PrintHelpAndExit();
		}
        auto pos = argument.find_first_of(":");
        if (pos == std::string::npos) {
            std::cerr << "Argument \"" << argument << "\" did not have a ':' character" << std::endl;
            PrintHelpAndExit();
        }
        std::string option = argument.substr(1, pos - 1); // everything after / and up to :
        std::string value = argument.substr(pos + 1); // everything after :

        auto it = m_option_map.find(option);
        if (it == m_option_map.cend()) {
            std::cout << "Unknown option \"" << option << "\". Skipping" << std::endl;
            continue;
        }

        it->second(value);
    }

    return argc;
}

static const char* USAGE =
R"(CovidSim /O /P [/A] [/AP] [/c] [/C] [/CLP[1-6]] [/d] [/D] [/H] [/I] [/KO] [/KP] [/L] [/LS] [/M] [/PP] [/R] [/s] [/S] [/SS] [/T] SetupSeed1 SetupSeed2 RunSeed1 RunSeed2

All '/' arguments are followed by a colon ':' and then the value.
More detailed information can be found in docs/inputs-and-outputs.md.
)";

static const char* DETAILED_USAGE =
R"(Argument type:

    DIR     Directory path
    DOUBLE  Floating point number with double precision
    INT     Integer number
    FILE    File path

Required arguments:

    /O      DIR     Output file path prefix
    /P      FILE    Parameter file

    SetupSeed*      RNG seeds when initializing
    RunSeed*        RNG seeds when running

Optional arguments:

    /A      FILE    Administrative division
    /AP     FILE    Air travel data file
    /c      INT     Number of threads to use
    /C      INT     Sets the P.PlaceCloseIndepThresh parameter
    /CLP*   DOUBLE  Overwrites wildcards in parameter file
    /d      FILE    Regional demography file
    /D      FILE    Population density file
    /H              Display this help page
    /I      FILE    Intervention file
    /KO     DOUBLE  Scales the P.KernelOffsetScale parameter
    /KP     DOUBLE  Scales the P.KernelPowerScale parameter
    /L      FILE    Network file to load
    /LS     FILE    Snapshot file to load
    /M      FILE    Output density file
    /PP     FILE    Pre-parameter file
    /R      DOUBLE  R0 scaling
    /s      FILE    School file
    /S      FILE    Network file to save
    /SS     DOUBLE,FILE	Interval and file to save snapshots
    /T      FILE    Sets the P.PreControlClusterIdCaseThreshold
)";

void PrintHelpAndExit()
{
	std::cerr << USAGE << std::endl;
	std::cerr << "Use the '/H' argument to get a detailed listing of "
				 "required and optional arguments." << std::endl;
	std::exit(1);
}

void PrintDetailedHelpAndExit()
{
	std::cerr << USAGE << "\n" << DETAILED_USAGE << std::endl;
	std::exit(1);
}
