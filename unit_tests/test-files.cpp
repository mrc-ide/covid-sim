#include <gtest/gtest.h>
#include "Error.h"
#include "Files.h"


TEST(Files, BasicTests) {

  // Tests normal use of fread_big, fwrite_big, xfopen, xfclose, xremove, xrename.

  uint64_t* buf = new uint64_t[10];
  uint64_t* compare = new uint64_t[10];
  for (uint64_t i=0; i<10; i++) {
    buf[i] = 180000000000000000L + i;
  }

  FILE* f = Files::xfopen("test_basictest.txt", "wb");
  Files::fwrite_big(buf, sizeof(uint64_t), 10, f);
  Files::xfclose(f);

  Files::xrename("test_basictest.txt", "renamed_test_basictest.txt");
  FILE* f2 = Files::xfopen("renamed_test_basictest.txt", "rb");
  Files::fread_big(compare, sizeof(uint64_t), 10, f2);
  Files::xfclose(f2);
  Files::xremove("renamed_test_basictest.txt");

  int matches = 0;
  for (int i = 0; i < 10; i++) {
    matches += (buf[i] == compare[i]) ? 1 : 0;
  }

  delete [] compare;
  delete [] buf;

  EXPECT_EQ(10, matches);
}

TEST(FilesDeathTests, file_open_fail) {
  ASSERT_DEATH({
    FILE* f = Files::xfopen("a_non_existent_file.txt", "rb");
  }, "Error .* opening file .*");
}

TEST(FilesDeathTests, file_close_fail) {
  ASSERT_DEATH({
    Files::xfclose(NULL);
  }, "Error .* closing filestream .*");
}

TEST(Files, file_open_fallback_exists) {
  FILE* f = Files::xfopen("prefer_this.txt", "w");
  Files::xfprintf(f, "1\n");
  Files::xfclose(f);
  FILE* f2 = Files::xfopen("to_this.txt", "w");
  Files::xfprintf(f2, "2\n");
  Files::xfclose(f);

  f2 = Files::xfopen("to_this.txt", "w");
  int result;

  ////////////////////////////////////////////////////////////////////////

  FILE* fav = Files::xfopen_fallback("prefer_this.txt", "r", f2);
  Files::xfscanf(fav, 1, "%d", &result);
  Files::xfclose(fav);
  EXPECT_EQ(result, 1);

  Files::xremove("prefer_this.txt");

  FILE* fail = Files::xfopen_fallback("prefer_this.txt", "r", f2);
  Files::xfscanf(fail, 1, "%d", &result);
  Files::xfclose(fail);
  EXPECT_EQ(result, 2);

  Files::xfclose(f2);

  FILE* nullfile = Files::xfopen_fallback("prefer_this.txt", "r", NULL);
  EXPECT_EQ(nullfile, null);


  //////////////////////////////////////////////////////////////////////

  FILE* f_if_exists = Files::xfopen_if_exists("to_this.txt", "r");
  Files::xfscanf(f_if_exists, 1, "%d", &result);
  Files::xfclose(f_if_exists);
  EXPECT_EQ(result, 2);
  
  Files::xremove("to_this.txt");
  f_if_exists = Files::xfopen_if_exists("to_this.txt", "r");
  EXPECT_EQ(f_if_exists, 0);
}

////////////////////////////////////////////////////////////////////////

TEST(FilesDeathTests, printf_fail) {
  ASSERT_DEATH({
    Files::xfprintf(NULL, "%d", 123);
  }, "Error .* doing fprintf .*");
}

TEST(FilesDeathTests, scanf_fail) {
  int result;
  ASSERT_DEATH({
    Files::xfscanf(NULL, 1, "%d", &result);
  }, "Error fsscanf looking for .*, expected .* matches, got .*");
}

TEST(FilesDeathTests, rename_fail) {
  int result;
  ASSERT_DEATH({
    Files::xrename("this_doesn_exist.txt", "neither_does_this.txt");
  }, "Error .* renaming file .* to fsscanf looking for .*, expected .* matches, got .*");
}

TEST(FilesDeathTests, remove_fail) {
  ASSERT_DEATH({
    Files::xremove("this_doesn_exist.txt");
  }, "Error .* removing file .*");
}

TEST(Files, sprintf) {
  char* buf = new char[10];
  Files::xsprintf(buf, "%d", 123);
  EXPECT_EQ(atoi(buf), 123);
  delete[] buf;
}


TEST(FilesDeathTests, sprintf_fail) {
  char* buf = new char[10];
  ASSERT_DEATH({
    Files::xsprintf(buf, "%d");
  }, "Error .* doing sprintf .*");
  delete[] buf;
}

TEST(Files, sscanf) {
  char* buf = new char[10];
  buf[0] = '1';
  buf[1] = '2';
  buf[2] = '3';
  buf[3] = 0;
  int result = 0;
  Files::xsscanf(buf, 1, "%d", &result);
  EXPECT_EQ(results, 123);
  delete[] buf;
}

TEST(FilesDeathTests, sscanf_fail) {
  char* buf = new char[10];
  ASSERT_DEATH({
    Files::xsscanf(buf, 1, "%d");
  }, "Error, xsscanf looking for .* in .*, expected .* matches, got .*");
  delete[] buf;
}
