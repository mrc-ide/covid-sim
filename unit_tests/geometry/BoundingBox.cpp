#include "BoundingBox.h"

// Third party headers
#include "gtest/gtest.h"

using namespace Geometry;

namespace CovidSim
{
  namespace Test
  {
    namespace Geometry
    {
      TEST(Geometry_Vector, SetMax)
      {
        Maximum2<double> v1(9.0, 6.5);
        Vector2d v2(7.0, 8.0);
        v1.expand(v2, 0.0);
        EXPECT_EQ(v1.x, 9.0) << "X value not max";
        EXPECT_EQ(v1.y, 8.0) << "Y value not max";
      }

      TEST(Geometry_Vector, SetMin)
      {
        Minimum2<double> v1(9.0, 6.5);
        Vector2d v2(7.0, 8.0);
        v1.expand(v2);
        EXPECT_EQ(v1.x, 7.0) << "X value not min";
        EXPECT_EQ(v1.y, 6.5) << "Y value not min";
      }

      TEST(Geometry_BoundingBox, Inside)
      {
        BoundingBox2<double> box;
        box.top_right() = Vector2d(9.0, 8.5);
        box.bottom_left() = Vector2d(5.0, 6.0);
        EXPECT_TRUE( box.inside(Vector2d(7.0, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(3.0, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(9.5, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(7.0, 5.0)));
        EXPECT_FALSE(box.inside(Vector2d(7.0, 9.0)));
      }

      TEST(Geometry_BoundingBox, Width)
      {
        BoundingBox2<double> box;
        box.top_right() = Vector2d(9.0, 8.5);
        box.bottom_left() = Vector2d(5.0, 6.0);
        EXPECT_EQ(4.0, box.width());
      }

      TEST(Geometry_BoundingBox, Height)
      {
        BoundingBox2<double> box;
        box.top_right() = Vector2d(9.0, 8.5);
        box.bottom_left() = Vector2d(5.0, 6.0);
        EXPECT_EQ(2.5, box.height());
      }

      TEST(Geometry_BoundingBox2d, Reset)
      {
        BoundingBox2d box;
        box.reset();
        EXPECT_LT(box.top_right().x, box.bottom_left().x);
        EXPECT_LT(box.top_right().y, box.bottom_left().y);
      }

      TEST(Geometry_BoundingBox2d, Add)
      {
        BoundingBox2d box;
        box.reset();
        box.expand(Vector2d(9.0, 8.5));
        box.expand(Vector2d(5.0, 6.0));
        EXPECT_TRUE( box.inside(Vector2d(7.0, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(3.0, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(9.5, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(7.0, 5.0)));
        EXPECT_FALSE(box.inside(Vector2d(7.0, 9.0)));
      }
    }
  }
}
