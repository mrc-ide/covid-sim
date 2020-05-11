#include "geometry/Geometry.h"

// Third party headers
#include "gtest/gtest.h"

using namespace CovidSim::Geometry;

namespace CovidSim
{
  namespace Test
  {
    namespace Geometry
    {
      TEST(Geometry_Vector, OnCopy)
      {
        Vector2d v1(9.0, 8.5), v2;
        v2 = v1;
        EXPECT_EQ(v1.x(), v2.x()) << "X value not copied";
        EXPECT_EQ(v1.y(), v2.y()) << "Y value not copied";
      }

      TEST(Geometry_Vector, SetMax)
      {
        Maximum2<double> v1(9.0, 6.5);
        Vector2d v2(7.0, 8.0);
        v1.add(v2, 0.0);
        EXPECT_EQ(v1.x(), 9.0) << "X value not max";
        EXPECT_EQ(v1.y(), 8.0) << "Y value not max";
      }

      TEST(Geometry_Vector, SetMin)
      {
        Minimum2<double> v1(9.0, 6.5);
        Vector2d v2(7.0, 8.0);
        v1.add(v2);
        EXPECT_EQ(v1.x(), 7.0) << "X value not min";
        EXPECT_EQ(v1.y(), 6.5) << "Y value not min";
      }
    }
  }
}
