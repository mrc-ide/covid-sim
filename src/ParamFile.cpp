#include <algorithm>
#include <cctype>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "ParamFile.hpp"

// ltrim(), rtrim(), and trim() are from a highly up voted answer on StackOverflow
// https://stackoverflow.com/a/217605

// trim from start (in place)
static inline void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            [](int ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
static inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            [](int ch) { return !std::isspace(ch); }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string& s) {
    ltrim(s);
    rtrim(s);
}

ParamReader::ParamReader(std::string const& param_file, std::string const& preparam_file, std::string const& admin_file)
{
    // duplicate params from different files take priority over one another:
    //     admin_file >> preparam_file >> param_file
    parse_param_file(param_file);
    parse_param_file(preparam_file);
    parse_param_file(admin_file);
}

template<typename T>
bool ParamReader::raw_extract(std::istringstream& stream, T& output)
{
    static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                  "Integral or floating point required.");

    stream >> output;
    if (stream.fail()) {
        if (output == std::numeric_limits<T>::max()) {
            std::cerr << "OVERFLOW: detected a number larger than " << std::numeric_limits<T>::max()
                      << " when parsing " << stream.str() << std::endl;
        }
        else if (output == std::numeric_limits<T>::min()) {
            std::cerr << "UNDERFLOW: detected a number smaller than " << std::numeric_limits<T>::min()
                      << " when parsing " << stream.str() << std::endl;
        }
        else {
            std::cerr << "ERROR: Expected integral type, got " << stream.str() << std::endl;
        }
        return false;
    }
    return true;
}

template<typename T>
bool ParamReader::extract(std::string const& param, T& output, T default_value)
{
    if (!exists(param))
    {
        std::cout << "Using default value: " << default_value << std::endl;
        output = default_value;
        return false;
    }

    std::istringstream iss(m_param_value_map[param]);
    return raw_extract(iss, output);
}

template<typename T>
void ParamReader::extract_or_exit(std::string const& param, T& output)
{
    if (!exists(param))
        std::exit(1);

    std::istringstream iss(m_param_value_map[param]);
    if (!raw_extract(iss, output))
        std::exit(1);
}

template<typename T>
bool ParamReader::extract_multiple(std::string const& param, T* output, std::size_t N, T default_value)
{
    if (!exists(param))
    {
        std::cerr << "Using default value: " << default_value << std::endl;
        std::fill_n(output, N, default_value);
        return false;
    }

    std::istringstream iss(m_param_value_map[param]);
    for (auto i = 0; i < N; i++)
    {
        if (!raw_extract(iss, output[i]))
        {
            std::cerr << "ERROR: Got " << i << " out of " << N << " parameters for "
                      << param << ".\nUsing default value: " << default_value << std::endl;
            std::fill_n(output, N, default_value);
            return false;
        }
    }
    return true;
}

template<typename T>
bool ParamReader::cond_extract_multiple(bool conditional, std::string const& param, T* output, std::size_t N, T default_value)
{
    if (conditional)
        return extract_multiple(param, output, N, default_value);
    std::fill_n(output, N, default_value);
    return true;
}

template<typename T>
bool ParamReader::extract_multiple_no_default(std::string const& param, T* output, std::size_t N)
{
    if (!exists(param))
    {
        return false;
    }

    std::istringstream iss(m_param_value_map[param]);
    for (auto i = 0; i < N; i++)
    {
        if (!raw_extract(iss, output[i]))
        {
            std::cerr << "ERROR: Got " << i << " out of " << N << " parameters for "
                      << param << std::endl;
            return false;
        }
    }
    return true;
}

template<typename T>
void ParamReader::extract_multiple_or_exit(std::string const& param, T* output, std::size_t N)
{
    if (!exists(param))
        std::exit(1);

    std::istringstream iss(m_param_value_map[param]);
    for (auto i = 0; i < N; i++)
    {
        if (!raw_extract(iss, output[i]))
        {
            std::cerr << "ERROR: Got " << i << " out of " << N << " parameters for " << param << std::endl;
            std::exit(1);
        }
    }
}

bool ParamReader::extract_multiple_strings_no_default(std::string const& param, std::vector<std::string>& output, std::size_t N)
{
    if (!exists(param))
        return false;

    std::istringstream iss(m_param_value_map[param]);
    std::string token;
    for (auto i = 0; i < N; i++)
    {
        iss >> token;
        if (iss.fail())
        {
            std::cerr << "ERROR: Only got " << i << " out of " << N << " parameters for " << param << std::endl;
            return false;
        }
        output.emplace_back(token);
    }
    return true;
}

void ParamReader::extract_multiple_strings_or_exit(std::string const& param, std::vector<std::string>& output, std::size_t N)
{
    if (!extract_multiple_strings_no_default(param, output, N))
        std::exit(1);
}

void ParamReader::extract_string_matrix_or_exit(std::string const& param, std::vector<std::vector<std::string>>& output, std::size_t num_cols)
{
    if (!exists(param))
        std::exit(1);

    // iterate over every row of the matrix
    std::istringstream matrix(m_param_value_map[param]);
    std::istringstream line_stream;
    for (std::string line; std::getline(matrix, line);)
    {
        // iterate over every token in each row and add it in a vector
        line_stream.str(line);
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

// Explicit template instantiations for the linker
// https://stackoverflow.com/questions/2351148/explicit-template-instantiation-when-is-it-used
template bool ParamReader::extract<double>(std::string const&, double&, double);
template bool ParamReader::extract<int>(std::string const&, int&, int);
template void ParamReader::extract_or_exit<double>(std::string const&, double&);
template void ParamReader::extract_or_exit<int>(std::string const&, int&);
template bool ParamReader::extract_multiple<double>(std::string const&, double*, std::size_t, double);
template bool ParamReader::extract_multiple<int>(std::string const&, int*, std::size_t, int);
template bool ParamReader::cond_extract_multiple<double>(bool, std::string const&, double*, std::size_t, double);
template bool ParamReader::cond_extract_multiple<int>(bool, std::string const&, int*, std::size_t, int);
template bool ParamReader::extract_multiple_no_default<double>(std::string const&, double*, std::size_t);
template bool ParamReader::extract_multiple_no_default<int>(std::string const&, int*, std::size_t);
template void ParamReader::extract_multiple_or_exit<double>(std::string const&, double*, std::size_t);
template void ParamReader::extract_multiple_or_exit<int>(std::string const&, int*, std::size_t);

bool ParamReader::exists(std::string const& param)
{
    auto it = m_param_value_map.find(param);
    if (it == m_param_value_map.cend())
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
            param_line.erase(param_line.cbegin());
            param_line.pop_back();
            // try to read the value of this parameter into a string
            std::string value_line;
            if (!std::getline(param_stream, value_line)) {
                std::cerr << "Error reading value for parameter: " << param_line << std::endl;
                break;
            }
            // store the parameter name and value for later use by ReadParams()
            auto it = m_param_value_map.find(param_line);
            if (it != m_param_value_map.cend())
            {
                std::cerr << "OVERRIDE: Param file (\"" << param_file << "\") is updating Key:["
                          << param_line << "]" << " Value:[" << value_line << "]" << std::endl;
                it->second = value_line;
            }
            else
            {
                m_param_value_map.emplace(param_line, value_line);
            }
        }
    }
}
