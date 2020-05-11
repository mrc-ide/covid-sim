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

      TEST(Geometry_BoundingBox2d, Outside)
      {
        BoundingBox2d box;
        box.top_right() = Vector2d(9.0, 8.5);
        box.bottom_left() = Vector2d(5.0, 6.0);
        EXPECT_FALSE(box.outside(Vector2f(7.0f, 7.0f)));
        EXPECT_TRUE( box.outside(Vector2f(3.0f, 7.0f)));
        EXPECT_TRUE( box.outside(Vector2f(9.5f, 7.0f)));
        EXPECT_TRUE( box.outside(Vector2f(7.0f, 5.0f)));
        EXPECT_TRUE( box.outside(Vector2f(7.0f, 9.0f)));
      }

      TEST(Geometry_BoundingBox2d, Reset)
      {
        BoundingBox2d box;
        box.reset();
        EXPECT_LT(box.top_right().x(), box.bottom_left().x());
        EXPECT_LT(box.top_right().y(), box.bottom_left().y());
      }

      TEST(Geometry_BoundingBox2d, Add)
      {
        BoundingBox2d box;
        box.reset();
        box.add(Vector2d(9.0, 8.5));
        box.add(Vector2d(5.0, 6.0));
        EXPECT_TRUE( box.inside(Vector2d(7.0, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(3.0, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(9.5, 7.0)));
        EXPECT_FALSE(box.inside(Vector2d(7.0, 5.0)));
        EXPECT_FALSE(box.inside(Vector2d(7.0, 9.0)));
      }

      TEST(Geometry_BitMapBounds, Inside)
      {
        BitMapBounds box;
        box.width() = 8;
        box.height() = 4;
        EXPECT_TRUE( box.inside(Vector2i( 6,  2)));
        EXPECT_FALSE(box.inside(Vector2i(-1,  2)));
        EXPECT_FALSE(box.inside(Vector2i( 6, -1)));
        EXPECT_FALSE(box.inside(Vector2i( 9,  2)));
        EXPECT_FALSE(box.inside(Vector2i( 6,  6)));
      }
    }
  }
}
