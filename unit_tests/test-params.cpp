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

bool match_int(int* res, int a, int b, int c, int d, int e) {
  return (res[0] == a) && (res[1] == b) && (res[2] == c) && (res[3] == d) && (res[4] == e);
}

TEST(ReadParams, int_vec) {

  // Basic lookups

  ParamMap PM = test_map("[A]\n-1 2 4 6 8\n[E]\n42\n");
  ParamMap PM2 = test_map("[B]\n3 5 7 0 42\n");
  ParamMap PM3 = test_map("[C]\n9 8 7 6 5\n[B]\n3 4 5 6 7\n");
  int* res = new int[5];
  Params::get_int_vec(PM2, PM, "A", res, 5, 0, 5, &P);
  EXPECT_TRUE(match_int(res, -1, 2, 4, 6, 8));
  Params::get_int_vec(PM2, PM, "B", res, 5, 0, 5, &P);
  EXPECT_TRUE(match_int(res, 3, 5, 7, 0, 42));

  // Lookups with defaults

  Params::get_int_vec(PM2, PM, "D", res, 5, 77, 5, &P);
  EXPECT_TRUE(match_int(res, 77, 77, 77, 77, 77));

  // Required (no default)

  Params::req_int_vec(PM2, PM, "B", res, 5, &P);
  EXPECT_TRUE(match_int(res, 3, 5, 7, 0, 42));
  Params::req_int_vec(PM3, PM2, PM, "B", res, 5, &P);
  EXPECT_TRUE(match_int(res, 3, 5, 7, 0, 42));
  Params::req_int_vec(PM3, PM, PM, "B", res, 5, &P);
  EXPECT_TRUE(match_int(res, 3, 4, 5, 6, 7));

  // Length 1

  Params::req_int_vec(PM, PM2, "E", res, 1, &P);
  EXPECT_EQ(res[0], 42);

  // Force-fail

  Params::get_int_vec_ff(false, PM, PM2, "A", res, 5,17, &P);
  EXPECT_TRUE(match_int(res, -1, 2, 4, 6 ,8));
  Params::get_int_vec_ff(false, PM, PM2, "Z", res, 5, 17, &P);
  EXPECT_TRUE(match_int(res, 17, 17, 17, 17, 17));
  Params::get_int_vec_ff(true, PM, PM2, "A", res, 5, 18, &P);
  EXPECT_TRUE(match_int(res, 18, 18, 18, 18, 18));
  Params::get_int_vec_ff(true, PM, PM2, "Z", res, 5, 19, &P);
  EXPECT_TRUE(match_int(res, 19, 19, 19, 19, 19));
}

bool ish(double x, double y) {
  return (fabs(x - y) < 1E-10);
}

bool match_dbl(double* res, double a, double b, double c, double d, double e) {
  return ish(res[0], a) && ish(res[1], b) && ish(res[2], c) && ish(res[3], d) && ish(res[4], e);
}


TEST(ReadParams, double_vec) {

  // Basic lookups

  ParamMap PM = test_map("[A]\n-1.1 2.1 4.1 6.1 8.1\n[E]\n42\n");
  ParamMap PM2 = test_map("[B]\n3.2 5.3 7.4 0 42.6\n");
  ParamMap PM3 = test_map("[C]\n9.1 8.2 7.3 6.4 5.5\n[B]\n3.6 4.7 5.8 6.9 7.0\n");
  double* res = new double[5];
  Params::get_double_vec(PM2, PM, "A", res, 5, 0, 5, &P);
  EXPECT_TRUE(match_dbl(res, -1.1, 2.1, 4.1, 6.1, 8.1));
  Params::get_double_vec(PM2, PM, "B", res, 5, 0, 5, &P);
  EXPECT_TRUE(match_dbl(res, 3.2, 5.3, 7.4, 0, 42.6));

  // Lookups with defaults

  Params::get_double_vec(PM2, PM, "D", res, 5, 7.7, 5, &P);
  EXPECT_TRUE(match_dbl(res, 7.7, 7.7, 7.7, 7.7, 7.7));

  // Required (no default)

  Params::req_double_vec(PM2, PM, "B", res, 5, &P);
  EXPECT_TRUE(match_dbl(res, 3.2, 5.3, 7.4, 0, 42.6));
  Params::req_double_vec(PM3, PM2, PM, "B", res, 5, &P);
  EXPECT_TRUE(match_dbl(res, 3.2, 5.3, 7.4, 0, 42.6));
  Params::req_double_vec(PM3, PM, PM, "B", res, 5, &P);
  EXPECT_TRUE(match_dbl(res, 3.6, 4.7, 5.8, 6.9, 7.0));

  // Length 1

  Params::req_double_vec(PM, PM2, "E", res, 1, &P);
  EXPECT_NEAR(res[0], 42, 1e-10);

  // Force-fail

  Params::get_double_vec_ff(false, PM, PM2, "A", res, 5, 17, &P);
  EXPECT_TRUE(match_dbl(res, -1.1, 2.1, 4.1, 6.1 ,8.1));
  Params::get_double_vec_ff(false, PM, PM2, "Z", res, 5, 17, &P);
  EXPECT_TRUE(match_dbl(res, 17, 17, 17, 17, 17));
}

TEST(ReadParams, string_vec) {
  ParamMap PM = test_map("[A]\nAbc Bcd\nCde Def Efg\n");
  ParamMap PM2 = test_map("[F]\nGhi\n[A]");
  ParamMap PM3 = test_map("[H]\nIjk Jkl\nKlm\nNop");
  char** res = new char*[5];
  for (int i=0; i<5; i++) res[i] = new char[4];
  Params::req_string_vec(PM3, PM2, PM, "A", res, 5, &P);
  EXPECT_STREQ(res[0], "Abc");
  EXPECT_STREQ(res[1], "Bcd");
  EXPECT_STREQ(res[2], "Cde");
  EXPECT_STREQ(res[3], "Def");
  EXPECT_STREQ(res[4], "Efg");
  Params::req_string_vec(PM3, PM2, "H", res, 4, &P);
  EXPECT_STREQ(res[0], "Ijk");
  EXPECT_STREQ(res[1], "Jkl");
  EXPECT_STREQ(res[2], "Klm");
  EXPECT_STREQ(res[3], "Nop");
}

//  void get_double_matrix(ParamMap &base, ParamMap &fallback, ParamMap &params,
//                         std::string param_name, double** array, int sizex, int sizey,
//                         double default_value, bool err_on_missing, Param* P);

//  void get_double_matrix(ParamMap &fallback, ParamMap &params,
//                         std::string param_name, double** array, int sizex,
//                         int sizey, double default_value, Param* P);
