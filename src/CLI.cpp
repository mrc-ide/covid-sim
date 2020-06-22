#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>

#include "CLI.h"
#include "Param.h"

void parse_read_file(std::string const& input, std::string& output) {
    // check to see if the file exists and error out if it doesn't
    if (static_cast<bool>(std::ifstream(input)) == false) {
        std::cerr << input << " is not a file" << std::endl;
        std::exit(1);
    }
    output = input;
}

void parse_write_dir(std::string const& input, std::string& output) {
    // check to see if this prefix already exists as a file and error out
    if (static_cast<bool>(std::ifstream(input)) == true) {
        std::cerr << "Cannot use this prefix, this path already exists"
                     " as a file: " << input << std::endl;
        std::exit(1);
    }
    // TODO: add a platform-independent check to see if the prefix could
    // be added as a directory or file
    output = input;
}

template<typename T>
void parse_number(std::string const& input, T& output) {
    static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                  "Integral or floating point required.");

    std::istringstream iss(input);
    iss >> output;
    if (iss.fail()) {
        if (output == std::numeric_limits<T>::max()) {
            std::cerr << "OVERFLOW: detected a number larger than " << std::numeric_limits<T>::max()
                      << " when parsing " << iss.str() << std::endl;
        }
        else if (output == std::numeric_limits<T>::min()) {
            std::cerr << "UNDERFLOW: detected a number smaller than " << std::numeric_limits<T>::min()
                      << " when parsing " << iss.str() << std::endl;
        }
        else {
            std::cerr << "ERROR: Expected integral or floating point, got " << iss.str() << std::endl;
        }
        std::exit(1);
    }
}

void CmdLineArgs::add_custom_option(std::string const& option, ParserFn func, std::string const& doc)
{
    if (m_option_map.find(option) != m_option_map.cend()) {
        std::cerr << "Duplicate option specified " << option << ", ignoring..." << std::endl;
        return;
    }
    m_option_map.emplace(option, func);
    m_doc_map.emplace(option, doc);
}


template<typename T>
void CmdLineArgs::add_number_option(std::string const& option, T& output, std::string const& doc)
{
    add_custom_option(option, std::bind(parse_number<T>, std::placeholders::_1, std::ref(output)), doc);
}

// Explicit template instantiations for the linker
// https://stackoverflow.com/questions/2351148/explicit-template-instantiation-when-is-it-used
template void CmdLineArgs::add_number_option<double>(std::string const&, double&, std::string const&);
template void CmdLineArgs::add_number_option<int>(std::string const&, int&, std::string const&);

void CmdLineArgs::add_string_option(std::string const& option, StringParserFn func, std::string& output, std::string const& doc)
{
    add_custom_option(option, std::bind(func, std::placeholders::_1, std::ref(output)), doc);
}

void CmdLineArgs::parse(int argc, char* argv[], Param& P) {
    // Detect if the user wants to print out the full help output
    if (argc >= 2) {
        std::string first(argv[1]);
        if (first.length() == 2 && first.compare("/H") == 0) {
            print_detailed_help_and_exit();
        }
    }
    if (argc < 7) {
        std::cerr << "Minimum number of arguments not met. Expected 6 got "
                  << (argc - 1) << "\n" << std::endl;
        print_help_and_exit();
    }

    // Get seeds.
    int i = argc - 4;
    parse_number(argv[i], P.setupSeed1);
    parse_number(argv[i+1], P.setupSeed2);
    parse_number(argv[i+2], P.runSeed1);
    parse_number(argv[i+3], P.runSeed2);

    for (i = 1; i < argc - 4; i++)
    {
        std::string argument(argv[i]);
        if (argument[0] != '/') {
            std::cerr << "Argument \"" << argument << "\" does not start with '/'" << std::endl;
            print_help_and_exit();
        }
        auto pos = argument.find_first_of(':');
        if (pos == std::string::npos) {
            std::cerr << "Argument \"" << argument << "\" did not have a ':' character" << std::endl;
            print_help_and_exit();
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
}

static const std::string USAGE =
R"(SetupSeed1 SetupSeed2 RunSeed1 RunSeed2

All '/' arguments are followed by a colon ':' and then the value.
More detailed information can be found in docs/inputs-and-outputs.md.
)";

static const std::string DETAILED_USAGE =
R"(Named arguments:

        SetupSeed*      RNG seeds when initializing
        RunSeed*        RNG seeds when running

Optional arguments:
)";

void CmdLineArgs::print_help()
{
    std::stringstream ss;
    ss << "CovidSim";
    for (auto const& it : m_option_map) {
        ss << " [/" << it.first << ']';
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
    for (auto const& it : m_doc_map) {
        ss << '\t' << it.first << '\t' << it.second << '\n';
    }
    std::cerr << DETAILED_USAGE << '\n' << ss.str() << std::endl;
}

void CmdLineArgs::print_detailed_help_and_exit()
{
    print_detailed_help();
    std::exit(1);
}
