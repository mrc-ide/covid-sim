#include <gtest/gtest.h>
#include "Files.h"
#include "ReadParams.h"
#include "Param.h"
#include <string>


ParamMap test_map(std::string contents) {
  Files::xfopen("test_param_tmp.txt", "w");
  Files::fprintf("%s", contents);
  Files::xfclose();
  ParamMap P = Params::read_params_map("test_param_tmp.txt");
  Files::xremove("test_param_tmp.txt");
  return P;
}

TEST(ReadParams, scalar_int) {
  Param P;
  ParamMap PM = test_map("[A]\n42[C]99");
  EXPECT_EQ(Params::get_int(PM, PM, "A", 0, &P), 42);    // Basic usage
  EXPECT_EQ(Params::get_int(PM, PM, "C", 0, &P), 99);    // End of file, no last line
  EXPECT_EQ(Params::get_int(PM, PM, "B", 1, &P), 1);\
  EXPECT_EQ(Params::get_int(PM, PM, PM, "A", 0, &P), 42);
  EXPECT_EQ(Params::get_int(PM, PM, PM, "B", 1, &P), 1);

  EXPECT_EQ(Params::get_int_ff(false, PM, PM, "A", 0, &P), 42);
  EXPECT_EQ(Params::get_int_ff(true, PM, PM, "A", 37, &P), 37);
  EXPECT_EQ(Params::req_int(PM, PM, PM, "A", 0, &P), 42);
  EXPECT_EQ(Params::req_int(PM, PM, "A", 0, &P), 42);
}
