#include <gtest/gtest.h>
#include "Dist.h"
#include "Error.h"
#include "Files.h"
#include "Memory.h"
#include "Model.h"
#include "ReadParams.h"
#include "Param.h"

#include <string>



ParamMap test_map(std::string contents) {
  FILE* f = Files::xfopen("test_param_tmp.txt", "w");
  Files::xfprintf(f, "%s", contents);
  Files::xfclose(f);
  ParamMap P = Params::read_params_map("test_param_tmp.txt");
  Files::xremove("test_param_tmp.txt");
  return P;
}

TEST(ReadParams, scalar_int) {
  Param P;
  ParamMap PM = test_map("[A]\n42[C]99");
  ParamMap PM2 = test_map("");
  EXPECT_EQ(Params::get_int(PM2, PM, "A", 0, &P), 42);      // Found in first (2)
  EXPECT_EQ(Params::get_int(PM, PM2, "A", 0, &P), 42);      // Found in fallback (2)
  EXPECT_EQ(Params::get_int(PM2, PM2, PM, "A", 0, &P), 42); // Found in first (3)
  EXPECT_EQ(Params::get_int(PM2, PM, PM2, "A", 0, &P), 42); // Found in fallback (3)
  EXPECT_EQ(Params::get_int(PM, PM2, PM2, "A", 0, &P), 42); // Found in base (3)
  EXPECT_EQ(Params::get_int(PM, PM, "B", 7, &P), 7);        // Not found in (2)
  EXPECT_EQ(Params::get_int(PM, PM, PM, "B", 7, &P), 7);    // Not found in (3)

  EXPECT_EQ(Params::get_int(PM, PM, "C", 0, &P), 99);       // End of file, no lineberak at end

  EXPECT_EQ(Params::get_int_ff(false, PM, PM, "A", 0, &P), 42); // Force failure (false)
  EXPECT_EQ(Params::get_int_ff(true, PM, PM, "A", 37, &P), 37); // Force failure (true)

  EXPECT_EQ(Params::req_int(PM, PM, PM, "A", &P), 42);     // Require in first (3)
  EXPECT_EQ(Params::req_int(PM, PM, PM2, "A", &P), 42);    // Require in fallback (3)
  EXPECT_EQ(Params::req_int(PM, PM2, PM2, "A", &P), 42);   // Require in base (3)
  EXPECT_EQ(Params::req_int(PM, PM, "A", &P), 42);         // Require in first (2)
  EXPECT_EQ(Params::req_int(PM, PM2, "A", &P), 42);        // Require in fallback (3)
}
