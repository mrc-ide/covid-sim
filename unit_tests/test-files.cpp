#include <gtest/gtest.h>
#include "Error.h"
#include "Files.h"


TEST(Files, BasicTests) {
  uint64_t* buf = new uint64_t[10];
  uint64_t* compare = new uint64_t[10];
  for (int i=0; i<10; i++) {
    buf[i] = 2e30 + i;
  }

  FILE* f = Files::xfopen("test_basictest.txt", "wb");
  Files::fwrite_big(buf, sizeof(uint64_t), 10, f);
  Files::xfclose(f);

  Files::xrename("test_basictest.txt", "renamed_test_basictest.txt");
  FILE* f2 = Files::xfopen("renamed_test_basictest.txt", "wb");
  Files::fread_big(compare, sizeof(uint64_t), 10, f2);
  Files::xfclose(f2);
  Files::xremove("renamed_test_basictest.txt");

  int matches = 0;
  for (int i=0; i<10; i++) {
    matches += (buf[i] == compare[i]) ? 1 : 0;
  }

  delete [] compare;
  delete [] buf;

  EXPECT_EQ(10, matches);
}

