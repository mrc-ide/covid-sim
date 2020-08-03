#ifndef COVIDSIM_GEOMETRY_BOUNDING_BOX_H_INCLUDED_
#define COVIDSIM_GEOMETRY_BOUNDING_BOX_H_INCLUDED_

// Local headers
#include "Vector2.h"
#include "Size.h"

// Standard library headers
#include <cmath>

namespace CovidSim { namespace Geometry {
  /// General class for representing a 2D maximum point
  template<typename T>
  class Maximum2 : public Vector2<T>
  {
  public:
    /// Default constructor initialises nothing
    Maximum2()
      : Vector2<T>()
    {
    }

    /// Construct with the values
    /// \param x X coordinate
    /// \param y Y coordinate
    Maximum2(T x, T y)
      : Vector2<T>(x, y)
    {
    }

    /// Expand so it is the maximum of this and the other point
    /// \param p Point
    /// \param offset Value to add to component of p
    void expand(const Vector2<T>& p, T offset)
    {
      if (p.x >= this->x) this->x = p.x + offset;
      if (p.y >= this->y) this->y = p.y + offset;
    }

    /// Is the point inside the box?
    /// \param p Point
    /// \return true if inside
    bool inside(const Vector2<T>& p) const
    {
      return p.x < this->x && p.y < this->y;
    }

    /// Expand to the grid
    /// \param width width
    /// \param height height
    void to_grid(const double& width, const double& height)
    {
      this->x = ceil(this->x / width ) * width;
      this->y = ceil(this->y / height) * height;
    }
  };

  /// General class for representing a 2D minimum point
  template<typename T>
  class Minimum2 : public Vector2<T>
  {
  public:
    /// Default constructor initialises nothing
    Minimum2()
      : Vector2<T>()
    {
    }

    /// Construct with the values
    /// \param x X coordinate
    /// \param y Y coordinate
    Minimum2(T x, T y)
      : Vector2<T>(x, y)
    {
    }

    /// Expand so it is the minimum of this and the other point
    /// \param p Point
    void expand(const Vector2<T>& p)
    {
      if (p.x < this->x) this->x = p.x;
      if (p.y < this->y) this->y = p.y;
    }

    /// Is the point inside the box?
    /// \param p Point
    /// \return true if inside
    bool inside(const Vector2<T>& p) const
    {
      return p.x >= this->x && p.y >= this->y;
    }

    /// Expand to the grid
    /// \param width width
    /// \param height height
    void to_grid(const double& width, const double& height)
    {
      this->x = floor(this->x / width ) * width;
      this->y = floor(this->y / height) * height;
    }
  };

  /// General class for representing an axis aligned bounding box
  template<typename T>
  class BoundingBox2
  {
  protected:
    /// The bottom left corner of the bounding box
    Minimum2<T> bottom_left_;

    /// The top right corner of the bounding box
    Maximum2<T> top_right_;

  public:
    /// Getter
    /// \return The bottom left corner of the bounding box
    Vector2<T>& bottom_left()
    {
      return bottom_left_;
    }

    /// Getter
    /// \return The bottom left corner of the bounding box
    const Vector2<T>& bottom_left() const
    {
      return bottom_left_;
    }

    /// Getter
    /// \return The top right corner of the bounding box
    Vector2<T>& top_right()
    {
      return top_right_;
    }

    /// Getter
    /// \return The top right corner of the bounding box
    const Vector2<T>& top_right() const
    {
      return top_right_;
    }

    /// Is the point inside the box?
    /// \param p Point
    /// \return true if inside
    bool inside(const Vector2<T>& p) const
    {
      return bottom_left_.inside(p) && top_right_.inside(p);
    }

    /// Expand to the grid
    /// \param width width
    /// \param height height
    void to_grid(const double& width, const double& height)
    {
      top_right_.to_grid(width, height);
      bottom_left_.to_grid(width, height);
    }

    /// Calculate the width
    /// \return The width
    double width() const
    {
      return top_right_.x - bottom_left_.x;
    }

    /// Calculate the height
    /// \return The height
    double height() const
    {
      return top_right_.y - bottom_left_.y;
    }
  };

  /// Specialisation of BoundingBox for double
  class BoundingBox2d : public BoundingBox2<double>
  {
  public:
    /// Set the minimum to working infinity and maximum to minus working infinity
    void reset()
    {
      bottom_left_.x = bottom_left_.y = 1e10;
      top_right_.x = top_right_.y = -1e10;
    }

    /// Expand the bounds to fit the given point
    /// \param p Point
    void expand(const Vector2d& p)
    {
      bottom_left_.expand(p);
      top_right_.expand(p, 1e-6);
    }
  };
}}

#endif
