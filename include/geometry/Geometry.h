#pragma once

#include <Eigen/Dense>

namespace CovidSim
{
  namespace Geometry
  {
    typedef Eigen::Vector2i Vector2i;
    typedef Eigen::Vector2d Vector2d;
    typedef Eigen::Vector2f Vector2f;

    typedef Eigen::DiagonalMatrix<double, 2> DiagonalMatrix2d;

    /// General class for representing a 2D maximum point
    template<typename T>
    class Maximum2 : public Eigen::Matrix<T, 2, 1>
    {
    public:
      /// Default constructor initialises nothing
      Maximum2() = default;

      /// Construct with the values
      /// \param x X coordinate
      /// \param y Y coordinate
      Maximum2(T x, T y)
        : Eigen::Matrix<T, 2, 1>(x, y)
      {
      }

      /// Adjust so it is the maximum of this and the other point
      /// \param p Point
      /// \param offset Value to add to component of p
      void add(const Eigen::Matrix<T, 2, 1>& p, T offset)
      {
        if (p.x() >= x()) x() = p.x() + offset;
        if (p.y() >= y()) y() = p.y() + offset;
      }

      /// Is the point inside the box?
      /// \param p Point
      /// \return true if inside
      bool inside(const Eigen::Matrix<T, 2, 1>& p) const
      {
        return p.x() < x() && p.y() < y();
      }

      /// Expand to the grid
      /// \param width width
      /// \param height height
      void to_grid(const double& width, const double& height)
      {
        x() = ceil(x() / width ) * width;
        y() = ceil(y() / height) * height;
      }
    };

    /// General class for representing a 2D minimum point
    template<typename T>
    class Minimum2 : public Eigen::Matrix<T, 2, 1>
    {
    public:
      /// Default constructor initialises nothing
      Minimum2() = default;

      /// Construct with the values
      /// \param x X coordinate
      /// \param y Y coordinate
      Minimum2(T x, T y)
        : Eigen::Matrix<T, 2, 1>(x, y)
      {
      }

      /// Adjust so it is the minimum of this and the other point
      /// \param p Point
      void add(const Eigen::Matrix<T, 2, 1>& p)
      {
        if (p.x() < x()) x() = p.x();
        if (p.y() < y()) y() = p.y();
      }

      /// Is the point inside the box?
      /// \param p Point
      /// \return true if inside
      bool inside(const Eigen::Matrix<T, 2, 1>& p) const
      {
        return p.x() >= x() && p.y() >= y();
      }

      /// Expand to the grid
      /// \param width width
      /// \param height height
      void to_grid(const double& width, const double& height)
      {
        x() = floor(x() / width ) * width;
        y() = floor(y() / height) * height;
      }
    };

    /// General class for representing an axis aligned bounding box
    /// \todo Replace with a third party GIS library?
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
      Eigen::Matrix<T, 2, 1>& bottom_left()
      {
        return bottom_left_;
      }

      /// Getter
      /// \return The bottom left corner of the bounding box
      const Eigen::Matrix<T, 2, 1>& bottom_left() const
      {
        return bottom_left_;
      }

      /// Getter
      /// \return The top right corner of the bounding box
      Eigen::Matrix<T, 2, 1>& top_right()
      {
        return top_right_;
      }

      /// Getter
      /// \return The top right corner of the bounding box
      const Eigen::Matrix<T, 2, 1>& top_right() const
      {
        return top_right_;
      }

      /// Is the point inside the box?
      /// \param p Point
      /// \return true if inside
      bool inside(const Eigen::Matrix<T, 2, 1>& p) const
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
        return top_right_.x() - bottom_left_.x();
      }

      /// Calculate the height
      /// \return The height
      double height() const
      {
        return top_right_.y() - bottom_left_.y();
      }
    };

    /// Specialisation of BoundingBox for double
    class BoundingBox2d : public BoundingBox2<double>
    {
    public:
      /// Is the point outside the box?
      /// \param p Point
      /// \return true if outside
      /// \todo This is not !inside() so is there a bug?
      bool outside(const Eigen::Vector2f& p) const
      {
        return p.x() < bottom_left_.x() || p.x() > top_right_.x() || p.y() < bottom_left_.y() || p.y() > top_right_.y();
      }

      /// Set the minimum to working infinity and maximum to minus working infinity
      void reset()
      {
        bottom_left_.x() = bottom_left_.y() = 1e10;
        top_right_.x() = top_right_.y() = -1e10;
      }

      /// Add a point to the bounds
      /// \param p Point
      void add(const Eigen::Vector2d& p)
      {
        bottom_left_.add(p);
        top_right_.add(p, 1e-6);
      }
    };

    /// As inheritance is evil use composition to change the behaviour of BoundingBox<int> for a BitMap
    /// Also add in the height in pixels of the entire bitmap output
    class BitMapBounds
    {
    private:
      /// The wrapped bounding box class
      BoundingBox2<int> box_;

      /// Height in pixels of the entire bitmap output, including both the spectrum at the top and the map area
      int total_height_; 

    public:
      /// Constructor
      BitMapBounds()
      {
        box_.bottom_left().x() = 0;
        box_.bottom_left().y() = 0;
      }

      /// Getter
      /// \return Width
      int width() const
      {
        return box_.top_right().x();
      }

      /// Setter
      /// \return Width
      int& width()
      {
        return box_.top_right().x();
      }

      /// Getter
      /// \return Height
      int height() const
      {
        return box_.top_right().y();
      }

      /// Setter
      /// \return Height
      int& height()
      {
        return box_.top_right().y();
      }

      /// Getter
      /// \return Height of the entire bitmap output
      int total_height() const
      {
        return total_height_;
      }

      /// Setter
      /// \return Height of the entire bitmap output
      int& total_height()
      {
        return total_height_;
      }

      /// Is the point inside the box?
      /// \param p Point
      /// \return true if inside
      bool inside(const Eigen::Vector2i& p) const
      {
        return box_.inside(p);
      }
    };
  }
}
