#ifndef COVIDSIM_PARAM_FILE_HPP_INCLUDED_
#define COVIDSIM_PARAM_FILE_HPP_INCLUDED_

#include <sstream>
#include <string>
#include <unordered_map>

class ParamReader {
public:
    ParamReader(std::string const& param_file, std::string const& preparam_file, std::string const& admin_file);

    template<typename T>
    bool extract(std::string const& param, T& output, T default_value);

    template<typename T>
    void extract_or_exit(std::string const& param, T& output);

    template<typename T>
    bool extract_multiple(std::string const& param, T* output, std::size_t N, T default_value);

    template<typename T>
    void extract_multiple_or_exit(std::string const& param, T* output, std::size_t N);

private:
    void parse_param_file(std::string param_file);

    template<typename T>
    bool raw_extract(std::istringstream& stream, T& output);

    std::unordered_map<std::string, std::string> m_param_value_map;
};

#endif // COVIDSIM_PARAM_FILE_HPP_INCLUDED_
