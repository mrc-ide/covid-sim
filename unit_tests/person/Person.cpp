#include "Models/Person.h"

// Third party headers
#include "gtest/gtest.h"

namespace CovidSim
{
  namespace Test
  {
    namespace Person
    {
      TEST(Person, SimpleTests)
      {
        Person  P;
        EXPECT_TRUE(P.is_alive());
        P.set_susceptible();
        EXPECT_TRUE(P.is_suscesptible());
        P.set_case();
        EXPECT_TRUE(P.is_case());
        P.set_dead();
        EXPECT_TRUE(P.is_dead();
      }
    }
  }
}
