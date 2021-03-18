#include <gtest/gtest.h>
#include "Dist.h"
#include "Error.h"
#include "Files.h"
#include "Memory.h"
#include "Model.h"
#include "ReadParams.h"
#include "Param.h"

#include <string>

Param P;

ParamMap test_map(std::string contents) {
  FILE* f = Files::xfopen("test_param_tmp.txt", "w");
  Files::xfprintf(f, "%s", contents.c_str());
  Files::xfclose(f);
  ParamMap P = Params::read_params_map("test_param_tmp.txt");
  Files::xremove("test_param_tmp.txt");
  return P;
}

TEST(ReadParams, scalar_int) {
  ParamMap PM = test_map("[A]\n42\n[C]\n99");
  ParamMap PM2 = test_map("");
  EXPECT_EQ(PM.size(), 2);
  EXPECT_EQ(PM2.size(), 0);
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


TEST(ReadParams, scalar_double) {
  ParamMap PM = test_map("[A]\n4.2\n[C]\n9.9");
  ParamMap PM2 = test_map("");
  EXPECT_EQ(PM.size(), 2);
  EXPECT_EQ(PM2.size(), 0);
  EXPECT_EQ(Params::get_double(PM2, PM, "A", 0, &P), 4.2);      // Found in first (2)
  EXPECT_EQ(Params::get_double(PM, PM2, "A", 0, &P), 4.2);      // Found in fallback (2)
  EXPECT_EQ(Params::get_double(PM2, PM2, PM, "A", 0, &P), 4.2); // Found in first (3)
  EXPECT_EQ(Params::get_double(PM2, PM, PM2, "A", 0, &P), 4.2); // Found in fallback (3)
  EXPECT_EQ(Params::get_double(PM, PM2, PM2, "A", 0, &P), 4.2); // Found in base (3)
  EXPECT_EQ(Params::get_double(PM, PM, "B", 7, &P), 7);        // Not found in (2)
  EXPECT_EQ(Params::get_double(PM, PM, PM, "B", 7, &P), 7);    // Not found in (3)
  EXPECT_EQ(Params::get_double(PM, PM, "C", 0, &P), 9.9);       // End of file, no lineberak at end
  EXPECT_EQ(Params::req_double(PM, PM, PM, "A", &P), 4.2);     // Require in first (3)
  EXPECT_EQ(Params::req_double(PM, PM, PM2, "A", &P), 4.2);    // Require in fallback (3)
  EXPECT_EQ(Params::req_double(PM, PM2, PM2, "A", &P), 4.2);   // Require in base (3)
  EXPECT_EQ(Params::req_double(PM, PM, "A", &P), 4.2);         // Require in first (2)
  EXPECT_EQ(Params::req_double(PM, PM2, "A", &P), 4.2);        // Require in fallback (3)
}

TEST(ReadParams, scalar_wildcard) {
  ParamMap PM = test_map("[A]\n#1\n[C]\n#2");
  P.clP[1] = 21;
  P.clP[2] = 4.2;
  EXPECT_EQ(Params::get_int(PM, PM, "A", 0, &P), 21);
  EXPECT_EQ(Params::get_double(PM, PM, "C", 0, &P), 4.2);
}

TEST(ReadParams, int_vec) {
  ParamMap PM = test_map("[A]\n-1 2 4 6 8\n");
