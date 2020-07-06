#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "Error.h"

void parse_read_file(std::string const& input, std::string& output)
{
	// check to see if the file exists and error out if it doesn't
	if (static_cast<bool>(std::ifstream(input)) == false)
	{
		ERR_CRITICAL_FMT("%s is not a file\n", input.c_str());
	}
	output = input;
}

void parse_write_dir(std::string const& input, std::string& output)
{
	// check to see if this prefix already exists as a file and error out
	if (static_cast<bool>(std::ifstream(input)) == true)
	{
		ERR_CRITICAL_FMT("Cannot use this prefix, this path already exists"
						 " as a file: %s\n", input.c_str());
	}
	// TODO: add a platform-independent check to see if the prefix could
	// be added as a directory or file
	output = input;
}

void parse_string(std::string const& input, std::string& output)
{
	output = input;
}

bool parse_integer_no_default(std::string const& input, int32_t& output)
{
	try
	{
		std::size_t pos;
		output = std::stoi(input, &pos);
		if (pos != input.size())
		{
			std::cerr << "Detected invalid characters after parsed integer: " << input << std::endl;
            return false;
		}
	}
	catch (const std::invalid_argument& e)
	{
		std::cerr << "EINVAL: Expected integer got %s" << input << std::endl;
        return false;
	}
	catch (const std::out_of_range& e)
	{
		std::cerr << "ERANGE: Input integer is out of range. Expected " << std::numeric_limits<int>::min()
                    << " to " << std::numeric_limits<int>::max() << ". Got " << input << std::endl;
        return false;
	}
    return true;
}

bool parse_integer(std::string const& input, int32_t& output, int32_t default_value)
{
    if (!parse_integer_no_default(input, output))
    {
        std::cerr << "Using default value: " << default_value << std::endl;
        output = default_value;
        return false;
    }
    return true;
}

void parse_integer_or_exit(std::string const& input, int32_t& output)
{
    if (!parse_integer_no_default(input, output))
        std::exit(1);
}

bool parse_double_no_default(std::string const& input, double& output)
{
	try
	{
		std::size_t pos;
		output = std::stod(input, &pos);
		if (pos != input.size())
		{
			std::cerr << "Detected invalid characters after parsed double: " << input << std::endl;
            return false;
		}
	}
	catch (const std::invalid_argument& e)
	{
		std::cerr << "EINVAL: Expected double got %s" << input << std::endl;
        return false;
	}
	catch (const std::out_of_range& e)
	{
		std::cerr << "ERANGE: Input integer is out of range. Expected " << std::numeric_limits<int>::min()
                    << " to " << std::numeric_limits<int>::max() << ". Got " << input << std::endl;
        return false;
	}
    return true;
}

bool parse_double(std::string const& input, double& output, double default_value)
{
    if (!parse_double_no_default(input, output))
    {
        std::cerr << "Using default value: " << default_value << std::endl;
        output = default_value;
        return false;
    }
    return true;
}

void parse_double_or_exit(std::string const& input, double& output)
{
    if (!parse_double_no_default(input, output))
        std::exit(1);
}
