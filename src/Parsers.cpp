#include <fstream>
#include <limits>
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

void parse_integer(std::string const& input, int& output)
{
	try
	{
		std::size_t pos;
		output = std::stoi(input, &pos);
		if (pos != input.size())
		{
			ERR_CRITICAL_FMT("Detected invalid characters after parsed integer: %s\n", input.c_str());
		}
	}
	catch (const std::invalid_argument& e)
	{
		ERR_CRITICAL_FMT("EINVAL: Expected integer got %s\n", input.c_str());
	}
	catch (const std::out_of_range& e)
	{
		ERR_CRITICAL_FMT("ERANGE: Input integer is out of range. Expected %d to %d. Got %s\n",
						std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), input.c_str());
	}
}

void parse_double(std::string const& input, double& output)
{
	try
	{
		std::size_t pos;
		output = std::stod(input, &pos);
		if (pos != input.size())
		{
			ERR_CRITICAL_FMT("Detected invalid characters after parsed double: %s\n", input.c_str());
		}
	}
	catch (const std::invalid_argument& e)
	{
		ERR_CRITICAL_FMT("EINVAL: Expected double got %s\n", input.c_str());
	}
	catch (const std::out_of_range& e)
	{
		ERR_CRITICAL_FMT("ERANGE: Input integer is out of range. Expected %.4e to %.4e. Got %s\n",
						std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), input.c_str());
	}
}