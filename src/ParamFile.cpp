#include <algorithm>
#include <cctype>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "ParamFile.hpp"

// ltrim(), rtrim(), and trim() are from a highly upvoted answer on StackOverflow
// https://stackoverflow.com/a/217605

// trim from start (in place)
static inline void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string& s) {
    ltrim(s);
    rtrim(s);
}

ParamReader::ParamReader(std::string const& param_file, std::string const& preparam_file, std::string const& admin_file)
{
    // params from certain CLI arguments take priority over one another
    // Admin File >> PreParam File >> Param File
    parse_param_file(param_file);
    parse_param_file(preparam_file);
    parse_param_file(admin_file);

/*
    for (auto const& it : m_param_value_map) {
        std::cout << "Key:[" << it.first << "] Value:[" << it.second << "]" << std::endl;
    }
*/
}

template<typename T>
bool ParamReader::extract(std::string const& param, T& output, T default)
{
    static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                  "Integral or floating point required.");

    auto it = m_param_value_map.find(param);
    if (it == m_param_value_map.cend())
    {
        std::cout << "Parameter \"" << param << "\" not found. Using default " << default << std::endl;
        output = default;
        return false;
    }

    std::istringstream iss(it->second);
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
            std::cerr << "ERROR: Expected integral type, got " << iss.str() << std::endl;
        }
        return false;
    }
    return true;
}

// Explicit template instantiations for the linker
// https://stackoverflow.com/questions/2351148/explicit-template-instantiation-when-is-it-used
template bool ParamReader::extract<double>(std::string const&, double&, double);
template bool ParamReader::extract<int>(std::string const&, int&, int);

template<typename T>
void ParamReader::extract_or_exit(std::string const& param, T& output)
{
    if (extract<T>(param, output, 0) == false)
    {
        std::exit(1);
    }
}

// Explicit template instantiations for the linker
// https://stackoverflow.com/questions/2351148/explicit-template-instantiation-when-is-it-used
template void ParamReader::extract_or_exit<double>(std::string const&, double&);
template void ParamReader::extract_or_exit<int>(std::string const&, int&);


void ParamReader::parse_param_file(std::string param_file)
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
            param_line.erase(param_line.cbegin());
            param_line.pop_back();
            // try to read the value of this parameter into a string
            std::string value_line;
            if (!std::getline(param_stream, value_line)) {
                std::cerr << "Error reading value for parameter: " << param_line << std::endl;
                break;
            }
            // store the parameter name and value for later use by ReadParams()
            std::cerr << "Key:[" << param_line << "] Value:[" << value_line << "]" << std::endl;
            m_param_value_map.emplace(param_line, value_line);
        }
    }
}
