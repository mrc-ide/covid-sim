#ifndef COVIDSIM_GEOMETRY_DISTANCE_H_INCLUDED_
#define COVIDSIM_GEOMETRY_DISTANCE_H_INCLUDED_

#include "BoundingBox.h"
#include "Vector2.h"
#include "Size.h"

#include <memory>

namespace CovidSim
{
  namespace Geometry
  {
    /// ABC for distance calculations
    class Distance
    {
    public:
      /// Calculate the distance squared between two rectangular grid indices
      /// \param l index
      /// \param m index
      /// \param stride number of indices in a row
      /// \param scales conversion to global co-ordinates
      /// \param minimum whether to calculate the minimum between cells
      /// \return distance squared
      virtual double distance_squared(int l, int m, int stride, const Size<double>& scales, bool minimum = false) = 0;

      /// Calculate the distance squared between two points
      /// \param a point in global co-ordinates
      /// \param b point in global co-ordinates
      /// \return distance squared
      virtual double distance_squared(const Vector2d& a, const Vector2d& b) = 0;

      /// Calculate the distance squared between two points
      /// \param a point in global co-ordinates
      /// \param b point in global co-ordinates
      /// \return distance squared
      double distance_squared(const Vector2f& a, const Vector2f& b)
      {
        return distance_squared(Vector2d(a.x, a.y), Vector2d(b.x, b.y));
      }
    };

    /// Provide a factory method for creating a concrete Distance class
    class DistanceFactory
    {
    public:
      /// Create a concrete distance
      /// \param utm whether to perform UTM calculations
      /// \param periodic_boundaries if not UTM whether the geometry is periodic or Eulerian
      /// \param spatial_bounding_box for UTM the bounds to which any lat long co-ordinates are relative
      /// \param global for periodic the horizontal and vertical periods in global units
      /// \return instance of a concrete distance
      static std::shared_ptr<Distance> create(bool utm, bool periodic_boundaries, const BoundingBox2d& spatial_bounding_box, const Size<double>& global);
    };
  }
}

#endif
