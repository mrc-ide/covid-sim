/** \file  ReadParams.h
 *  \brief Provide support for parsing param files.
 */

#ifndef READ_PARAMS_H_INCLUDED_
#define READ_PARAMS_H_INCLUDED_

#include <cstdint>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

#include "Dist.h"
#include "Error.h"
#include "Files.h"
#include "Memory.h"
#include "Model.h"
#include "Param.h"


typedef std::map<std::string, std::string> ParamMap;
typedef std::pair<std::string, std::string> ParamPair;
typedef std::map<std::string, std::string>::iterator ParamIter;

namespace Params
{

/** \brief             Read a parameter file (ini-like format) into a map
 *  \param  file       The filename to be read
 *  \return            A map of ParamPair (key, string) map
 */

  ParamMap read_params_map(const char* file);



/** \brief             Overwrite a value in the form "#00" with a command-line parameter
 *  \param  value      A string, which if it starts with "#" is interpreted as a CMD-line parameter.
 *  \param  P          The Param structure, containing the clP array.
 *  \return            The original string (if it's not a CLP), or otherwise, the matching CLP value.
 */

  std::string clp_overwrite(std::string value, Param* P);



/** \brief             Lookup a parameter string from one of three possible parameter maps
 *  \param params      The preferred map to find the parameter in
 *  \param fallback    A fallback if the parameter is not found in params
 *  \param base        A final map, should the parameter not be found in fallback.
 *  \param param_name  The name of the parameter to look up.
 *  \param P           The params (for looking up /CLP)
 *  \return            A string (possibly including '\n' for that parameter. If the result starts
 *                     with '#', the matching command-line parameter (/CP01:) will be returned.
 */

  std::string lookup_param(ParamMap &base, ParamMap &fallback, ParamMap &params,
                           std::string param_name, Param* P);



/** \brief             Overide of lookup_param with only one fallback option to find the parameter
 *  \param params      The preferred map to find the parameter in
 *  \param fallback    A fallback if the parameter is not found in params
 *  \param param_name  The name of the parameter to look up.
 *  \param P           The params (for looking up /CLP)
 *  \return            A string (possibly including '\n' for that parameter. If the result starts
 *                     with '#', the matching command-line parameter (/CP01:) will be returned.
 */

  std::string lookup_param(ParamMap &fallback, ParamMap &params, std::string param_name, Param* P);



/** \brief             Query whether a named parameter is found in one of three maps/
 *  \param params      The preferred map to find the parameter in
 *  \param fallback    A fallback if the parameter is not found in params
 *  \param base        A final map, should the parameter not be found in fallback.
 *  \param param_name  The name of the parameter to look up.
 *  \return            A boolean, TRUE only if the parameter is found.
 */

  bool param_found(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name);



/** \brief             Overide of param_found with only one fallback option to find the parameter
 *  \param params      The preferred map to find the parameter in
 *  \param fallback    A fallback if the parameter is not found in params
 *  \param param_name  The name of the parameter to look up.
 *  \return            A boolean, TRUE only if the parameter is found.
 */

  bool param_found(ParamMap &fallback, ParamMap &params, std::string param_name);



/** \brief                Retrieve parameter of type double from one of three parameter maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param default_value  Value to return if the parameter is missing
 *  \param err_on_missing If true, and the parameter is not found, then stop.
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               A boolean, TRUE only if the parameter is found.
 */

  double get_double(ParamMap &base, ParamMap &fallback, ParamMap &params,
                    std::string param_name, double default_value,
                    bool err_on_missing, Param* P);



/** \brief                Retrieve parameter of type double from one of three parameter maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param default_value  Value to return if the parameter is missing
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The parameter value (double)
 */

  double get_double(ParamMap &base, ParamMap &fallback, ParamMap &params,
                    std::string param_name, double default_value, Param* P);
 


/** \brief                A wrapper for get_double, taking only two parameter maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param default_value  Value to return if the parameter is missing
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The parameter value (double)
 */

  double get_double(ParamMap &fallback, ParamMap &params, std::string param_name,
    double default_value, Param* P);



/** \brief                A wrapper for get_double, which requires the parameter to be found
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The parameter value (double)
 */

  double req_double(ParamMap &base, ParamMap &fallback, ParamMap &params,
                    std::string param_name, Param* P);



/** \brief                A wrapper for req_double, where only two maps are needed.
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The parameter value (double)
 */


  double req_double(ParamMap &fallback, ParamMap &params,
                    std::string param_name, Param* P);



/** \brief                Retrieve parameter of type double from one of three parameter maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param default_value  Value to return if the parameter is missing
 *  \param err_on_missing If true, and the parameter is not found, then stop.
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The integer result
 */

  int get_int(ParamMap &base, ParamMap &fallback, ParamMap &params,
              std::string param_name, int default_value,
              bool err_on_missing, Param* P);



/** \brief                A wrapper for get_int, taking only two parameter maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param default_value  Value to return if the parameter is missing
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The integer result
 */

  int get_int(ParamMap &fallback, ParamMap &params, std::string param_name,
    int default_value, Param* P);


  /** \brief                Wrapper for get_int, three params, and default value
   *  \param params         The preferred map to find the parameter in
   *  \param fallback       A fallback if the parameter is not found in params
   *  \param param_name     The name of the parameter to look up.
   *  \param default_value  Value to return if the parameter is missing
   *  \param P              The Param struct, to loop up clP if necessary.
   *  \return               The integer result
   */

  int get_int(ParamMap &base, ParamMap &fallback, ParamMap &params, std::string param_name,
              int default_value, Param* P);



/** \brief                A wrapper for get_int, which requires the parameter exist
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The integer result
 */

  int req_int(ParamMap &base, ParamMap &fallback, ParamMap &params,
                 std::string param_name, Param* P);



/** \brief                A wrapper for req_int, where only two maps are needed.
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param P              The Param struct, to loop up clP if necessary.
 *  \return               The integer result.
 */

  int req_int(ParamMap &fallback, ParamMap &params,
              std::string param_name, Param* P);



/** \brief                Parse a parameter which is a vector of doubles 
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The double-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param default_value  The value to write if the param is not found
 *  \param default_size   The number of default_values to write if the param is not found
 *  \param err_on_missing Stop if the parameter is not found.
 *  \param force_fail     If true, just write the default value, don't look for params.
 *  \param P              The param struct, for looking up command-line args
 */

  void get_double_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
                      std::string param_name, double* array, int expected,
                      double default_value, int default_size,
                      bool err_on_missing, Param* P, bool force_fail);



/** \brief                Wrapper for get_double_vec only taking two param maps.
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The double-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param default_value  The value to write if the param is not found
 *  \param default_size   The number of default_values to write if the param is not found
 *  \param err_on_missing Stop if the parameter is not found.
 *  \param P              The param struct, for looking up command-line args
 */

  void get_double_vec(ParamMap &fallback, ParamMap &params,
                      std::string param_name, double* array, int expected,
                      double default_value, int default_size,
                      Param* P);




/** \brief                Require a vector of doubles of the given size.
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The double-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param P              The param struct, for looking up command-line args
 */

  void req_double_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
                      std::string param_name, double* array, int expected,
                      Param* P);



/** \brief                Wrapper for req_double_vec with only two maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The double-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param P              The param struct, for looking up command-line args
 */

 void req_double_vec(ParamMap &fallback, ParamMap &params,
                     std::string param_name, double* array, int expected,
                     Param* P);



/** \brief                Wrapper for req_double_vec allowing "forced failure"
 *                        in which case default values are written
 *  \param force_fail     If true, write defaults even if param exists
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The double-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param default_value  Default value to write if forced (or if params absent)
 *  \param P              The param struct, for looking up command-line args
 */

 void get_double_vec_ff(bool force_fail, ParamMap &fallback, ParamMap &params,
   std::string param_name, double* array, int expected,
   double default_value, Param* P);



 /** \brief              A wrapper for req_int, where only two maps are needed.
*  \param params         The preferred map to find the parameter in
*  \param fallback       A fallback if the parameter is not found in params
*  \param param_name     The name of the parameter to look up.
*  \param P              The Param struct, to loop up clP if necessary.
*  \return               The integer result.
*/

 int req_int(ParamMap &fallback, ParamMap &params,
   std::string param_name, Param* P);



/** \brief                Parse a parameter which is a vector of doubles
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The int-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param default_value  The value to write if the param is not found
 *  \param default_size   The number of default_values to write if the param is not found
 *  \param err_on_missing Stop if the parameter is not found.
 *  \param force_fail     If true, just write the default value, don't look for params.
 *  \param P              The param struct, for looking up command-line args
 */

 void get_int_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
                  std::string param_name, int* array, int expected,
                  int default_value, int default_size,
                  bool err_on_missing, Param* P, bool force_fail);



/** \brief                Wrapper for get_int_vec only taking two param maps.
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The int-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param default_value  The value to write if the param is not found
 *  \param default_size   The number of default_values to write if the param is not found
 *  \param err_on_missing Stop if the parameter is not found.
 *  \param force_fail     Optional - if true, just write the default value, don't look for params.
 *  \param P              The param struct, for looking up command-line args
 */

 void get_int_vec(ParamMap &fallback, ParamMap &params,
                  std::string param_name, int* array, int expected,
                  int default_value, int default_size,
                  bool err_on_missing, Param* P, bool force_fail);



 /** \brief                Wrapper for get_int_vec only taking two param maps, with default values
  *  \param params         The preferred map to find the parameter in
  *  \param fallback       A fallback if the parameter is not found in params
  *  \param param_name     The name of the parameter to look up.
  *  \param array          The int-array into which to put the results
  *  \param expected       The number of elements we expect to find
  *  \param default_value  The value to write if the param is not found
  *  \param default_size   The number of default_values to write if the param is not found
  *  \param P              The param struct, for looking up command-line args
  */

 void get_int_vec(ParamMap &fallback, ParamMap &params,
   std::string param_name, int* array, int expected,
   int default_value, int default_size,
   Param* P);





/** \brief                Require a vector of ints of the given size.
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The int-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param P              The param struct, for looking up command-line args
 */

 void req_int_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
                  std::string param_name, int* array, int expected,
                  Param* P);



/** \brief                Wrapper for req_int_vec with only two maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The int-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param P              The param struct, for looking up command-line args
 */

 void req_int_vec(ParamMap &fallback, ParamMap &params,
                  std::string param_name, int* array, int expected,
                  Param* P);


/** \brief                Wrapper for req_int_vec allowing "forced failure"
 *                        in which case default values are written
 *  \param force_fail     If true, write defaults even if param exists
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The int-array into which to put the results
 *  \param expected       The number of elements we expect to find
 *  \param default_value  Default value to write if forced (or if params absent)
 *  \param P              The param struct, for looking up command-line args
 */

 void get_int_vec_ff(bool force_fail, ParamMap &fallback, ParamMap &params,
                     std::string param_name, int* array, int expected,
                     int default_value, Param* P);


 /** \brief                Require a vector of strings as a parameter
  *  \param params         The preferred map to find the parameter in
  *  \param fallback       A fallback if the parameter is not found in params
  *  \param base           A final map, should the parameter not be found in fallback.
  *  \param param_name     The name of the parameter to look up.
  *  \param array          An array of char* to write the results into
  *  \param expected       The number of elements we expect to find
  *  \param P              The param struct, for looking up command-line args
  *  \return               Number of strings parsed*
  */

  int req_string_vec(ParamMap &base, ParamMap &fallback, ParamMap &params,
                      std::string param_name, char** array, int expected, Param* P);


/** \brief                Wrapper for req_string_vec with only two maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          An array of char* to write the results into
 *  \param expected       The number of elements we expect to find
 *  \param P              The param struct, for looking up command-line args
 *  \return               Number of strings parsed
 */

  int req_string_vec(ParamMap &fallback, ParamMap &params,
                      std::string param_name, char** array, int expected, Param* P);



/** \brief                Read a matrix of doubles from a parameter
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param base           A final map, should the parameter not be found in fallback.
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The array to write into
 *  \param sizex          Size of first dimesion in A[x][y]
 *  \param sizey          Size of second dimension in A[x][y]
 *  \param default_value  Double to fill array with if no param is found
 *  \param err_on_missing If true, stop if no matching param is found
 *  \param P              The param struct, for looking up command-line args
 */

  void get_double_matrix(ParamMap &base, ParamMap &fallback, ParamMap &params,
                         std::string param_name, double** array, int sizex, int sizey,
                         double default_value, bool err_on_missing, Param* P);

/** \brief                Wrapper for get_double_matrix with only two maps
 *  \param params         The preferred map to find the parameter in
 *  \param fallback       A fallback if the parameter is not found in params
 *  \param param_name     The name of the parameter to look up.
 *  \param array          The array to write into
 *  \param sizex          Size of first dimesion in A[x][y]
 *  \param sizey          Size of second dimension in A[x][y]
 *  \param default_value  Double to fill array with if no param is found
 *  \param P              The param struct, for looking up command-line args
 */

  void get_double_matrix(ParamMap &fallback, ParamMap &params,
                         std::string param_name, double** array, int sizex,
                         int sizey, double default_value, Param* P);

/** \brief                Parse an inverse CDF
 */
   void get_inverse_cdf(ParamMap fallback, ParamMap params, const char* icdf_name, InverseCdf* inverseCdf, Param* P, double start_value);

/** \brief                Allocate memory for 2-D (or higher) objects in params
 */
  void alloc_params(Param* P);


  /** \brief                Top-level call for ReadParams.
 */
  void airport_params(ParamMap adm_params, ParamMap pre_params, ParamMap params, Param* P);


/** \brief                Top-level call for ReadParams.
 */
  void ReadParams(std::string const& ParamFile, std::string const& PreParamFile, std::string const& AdUnitFile, Param* P, AdminUnit* AdUnits);

} // namespace Params

#endif // READ_PARAMS_H_INCLUDED_
