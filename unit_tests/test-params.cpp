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

ParamMap test_map(std::string contents, int test_no) {
  std::string fname = "test_param_tmp";
  fname.append(std::to_string(test_no));
  fname.append(".txt");
  FILE* f = Files::xfopen(fname.c_str(), "w");
  Files::xfprintf(f, "%s", contents.c_str());
  Files::xfclose(f);
  ParamMap P = Params::read_params_map(fname.c_str());
  Files::xremove(fname.c_str());
  return P;
}

TEST(ReadParams, scalar_int) {
  ParamMap PM = test_map("[A]\n42\n[C]\n99", 1);
  ParamMap PM2 = test_map("", 1);
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
  ParamMap PM = test_map("[A]\n4.2\n[C]\n9.9", 2);
  ParamMap PM2 = test_map("", 2);
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
  ParamMap PM = test_map("[A]\n#1\n[C]\n#2", 3);
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

  ParamMap PM = test_map("[A]\n-1 2 4 6 8\n[E]\n42\n", 4);
  ParamMap PM2 = test_map("[B]\n3 5 7 0 42\n", 4);
  ParamMap PM3 = test_map("[C]\n9 8 7 6 5\n[B]\n3 4 5 6 7\n[F]\n#4 #2 7 #3 #1", 4);
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

    // Wildcards

  P.clP[1] = 2;
  P.clP[2] = 4;
  P.clP[3] = 8;
  P.clP[4] = 16;

  Params::req_int_vec(PM, PM2, PM3, "F", res, 5, &P);
  EXPECT_TRUE(match_int(res, 16, 4, 7, 8, 2));
}

bool ish(double x, double y) {
  return (fabs(x - y) < 1E-10);
}

bool match_dbl(double* res, double a, double b, double c, double d, double e) {
  return ish(res[0], a) && ish(res[1], b) && ish(res[2], c) && ish(res[3], d) && ish(res[4], e);
}


TEST(ReadParams, double_vec) {

  // Basic lookups

  ParamMap PM = test_map("[A]\n-1.1 2.1 4.1 6.1 8.1\n[E]\n42\n", 5);
  ParamMap PM2 = test_map("[B]\n3.2 5.3 7.4 0 42.6\n[F]\n#1 #2 #3 #4 17.1\n", 5);
  ParamMap PM3 = test_map("[C]\n9.1 8.2 7.3 6.4 5.5\n[B]\n3.6 4.7 5.8 6.9 7.0\n", 5);
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

  // Wildcards

  P.clP[1] = 3.14;
  P.clP[2] = 6.28;
  P.clP[3] = 0.707;
  P.clP[4] = 1.414;

  Params::req_double_vec(PM, PM2, PM3, "F", res, 5, &P);
  EXPECT_TRUE(match_dbl(res, 3.14, 6.28, 0.707, 1.414, 17.1));

}

TEST(ReadParams, string_vec) {
  ParamMap PM = test_map("[A]\nAbc Bcd\nCde Def Efg\n", 6);
  ParamMap PM2 = test_map("[F]\nGhi\n[A]", 6);
  ParamMap PM3 = test_map("[H]\nIjk Jkl\nKlm\nNop", 6);
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

TEST(ReadParams, double_matrix) {
  ParamMap PM = test_map("[M]\n1.1 2.2 3.3 4.4 5.5\n2.2 3.3 4.4 5.5 6.6\n3.3 4.4 5.5 6.6 7.7\n4.4 5.5 6.6 7.7 8.8\n5.5 6.6 7.7 8.8 9.9\n[MW]\n#1 0 0 0 0\n0 #2 0 0 0\n0 0 #3 0 0\n0 0 0 #4 0\n0 0 0 10.5 #5\n", 7);
  double** M = new double*[5];
  for (int i=0; i<5; i++) M[i] = new double[5];

  Params::get_double_matrix(PM, PM, "M", M, 5, 5, 0, &P);
  EXPECT_TRUE(match_dbl(M[0], 1.1, 2.2, 3.3, 4.4, 5.5));
  EXPECT_TRUE(match_dbl(M[1], 2.2, 3.3, 4.4, 5.5, 6.6));
  EXPECT_TRUE(match_dbl(M[2], 3.3, 4.4, 5.5, 6.6, 7.7));
  EXPECT_TRUE(match_dbl(M[3], 4.4, 5.5, 6.6, 7.7, 8.8));
  EXPECT_TRUE(match_dbl(M[4], 5.5, 6.6, 7.7, 8.8, 9.9));

  Params::get_double_matrix(PM, PM, "Q", M, 5, 5, 3.14, &P);
  EXPECT_TRUE(match_dbl(M[0], 3.14, 3.14, 3.14, 3.14, 3.14));
  EXPECT_TRUE(match_dbl(M[1], 3.14, 3.14, 3.14, 3.14, 3.14));
  EXPECT_TRUE(match_dbl(M[2], 3.14, 3.14, 3.14, 3.14, 3.14));
  EXPECT_TRUE(match_dbl(M[3], 3.14, 3.14, 3.14, 3.14, 3.14));
  EXPECT_TRUE(match_dbl(M[4], 3.14, 3.14, 3.14, 3.14, 3.14));

// Wildcards

  P.clP[1] = 1.5;
  P.clP[2] = 3.5;
  P.clP[3] = 7.5;
  P.clP[4] = 9.5;
  P.clP[5] = 11.5;

  Params::get_double_matrix(PM, PM, "MW", M, 5, 5, 0, &P);
  EXPECT_TRUE(match_dbl(M[0], 1.5, 0, 0, 0, 0));
  EXPECT_TRUE(match_dbl(M[1], 0, 3.5, 0, 0, 0));
  EXPECT_TRUE(match_dbl(M[2], 0, 0, 7.5, 0, 0));
  EXPECT_TRUE(match_dbl(M[3], 0, 0, 0, 9.5, 10.5));
  EXPECT_TRUE(match_dbl(M[4], 0, 0, 0, 0, 11.5));

  for (int i=0; i<5; i++) delete[] M[i];
  delete[] M;
}




//  void get_double_matrix(ParamMap &fallback, ParamMap &params,
//                         std::string param_name, double** array, int sizex,
//                         int sizey, double default_value, Param* P);
