#include "geometry/Distance.h"

namespace CovidSim
{
  namespace Geometry
  {
    /// Concrete class for performing UTM distance calculations
    class UTMDistance : public Distance
    {
      /// Math constant defined as the ratio of a circle's circumference to its diameter.
      ///
      /// TODO: since all calculations using this constant are being automatically
      /// type-casted to double, should the precision be extended for more accuracy in
      /// the simulations?
      ///
      /// Eventually could be replaced with C++20's std::numbers::pi.
      /// https://en.cppreference.com/w/cpp/header/numbers
      static constexpr double pi_ = 3.1415926535; // full double precision: 3.14159265358979323846

      /// An arc minute of latitude along any line of longitude in meters.
      ///
      /// Also known as the International Nautical Mile.
      ///
      /// @see https://en.wikipedia.org/wiki/Nautical_mile
      static constexpr int nmi_ = 1852;

      /// The number of arc minutes in one degree.
      ///
      /// @see https://en.wikipedia.org/wiki/Minute_and_second_of_arc
      static constexpr int arcminutes_per_degree_ = 60;

      /// The number of degrees in a complete rotation.
      ///
      /// @see https://en.wikipedia.org/wiki/Turn_(angle)
      static constexpr int degrees_per_turn_ = 360;

      /// The earth's circumference in meters.
      ///
      /// The units of cancellation:
      ///    meters/minute * minutes/degree * degrees = meters
      static constexpr int earth_circumference_ = nmi_ * arcminutes_per_degree_ * degrees_per_turn_;

      /// The earth's diameter in meters.
      static constexpr double earth_diameter_ = earth_circumference_ / pi_;

      /// The Earth's radius in meters.
      static constexpr double earth_radius_ = earth_diameter_ / 2;

      /// The bounds to which any lat long co-ordinates are relative
      const BoundingBox2d& spatial_bounding_box_;

      /// The lookup table for sin(x)
      double sin_x_[degrees_per_turn_ + 1];

      /// The lookup table for cos(x)
      double cos_x_[degrees_per_turn_ + 1];

      /// The lookup table for asin(sqrt(x)) * asin(sqrt(x))
      double asin2_sq_x_[1001];

    public:
      /// Construct a class for UTM distance calculations
      /// \param spatial_bounding_box the bounds to which any lat long co-ordinates are relative
      UTMDistance(const BoundingBox2d& spatial_bounding_box)
        : spatial_bounding_box_(spatial_bounding_box)
      {
      }

      /// Create the lookup tables
      void initialise()
      {
        for (int i = 0; i <= 1000; i++)
        {
          asin2_sq_x_[i] = asin(sqrt(i / 1000.0));
          asin2_sq_x_[i] = asin2_sq_x_[i] * asin2_sq_x_[i];
        }

        for (int i = 0; i <= degrees_per_turn_; i++)
        {
          double t = pi_ * i / 180;
          sin_x_[i] = sin(t);
          cos_x_[i] = cos(t);
        }
      }

      /// Calculate the distance squared between two rectangular grid indices
      /// \param l index
      /// \param m index
      /// \param stride number of indices in a row
      /// \param scales conversion to global co-ordinates
      /// \param minimum whether to calculate the minimum between cells
      /// \return distance squared
      double distance_squared(int l, int m, int stride, const Size<double>& scales, bool minimum) override
      {
        int i = l;
        int j = m;
        if (minimum)
        {
          if (scales.width * ((double)abs(m / stride - l / stride)) > pi_)
          {
            if (m / stride > l / stride)
              j += stride;
            else if (m / stride < l / stride)
              i += stride;
          }
          else
          {
            if (m / stride > l / stride)
              i += stride;
            else if (m / stride < l / stride)
              j += stride;
          }

          if (m % stride > l % stride)
            i++;
          else if (m % stride < l % stride)
            j++;
        }

        CovidSim::Geometry::Vector2d u(scales.width * fabs((double)(i / stride)), scales.height * fabs((double)(i % stride)));
        CovidSim::Geometry::Vector2d v(scales.width * fabs((double)(j / stride)), scales.height * fabs((double)(j % stride)));
        return distance_squared(u, v);
      }

      /// Calculate the distance squared between two points
      /// \param a lat long point in degrees relative to the spatial bounding box
      /// \param b lat long point in degrees relative to the spatial bounding box
      /// \return distance squared
      double distance_squared(const Vector2d& a, const Vector2d& b) override
      {
        double x1 = a.x;
        double y1 = a.y;
        double x2 = b.x;
        double y2 = b.y;

        double x, y, cy1, cy2, yt, xi, yi;

        x = fabs(x1 - x2) / 2;
        y = fabs(y1 - y2) / 2;
        xi = floor(x);
        yi = floor(y);
        x -= xi;
        y -= yi;
        x = (1 - x) * sin_x_[(int)xi] + x * sin_x_[((int)xi) + 1];
        y = (1 - y) * sin_x_[(int)yi] + y * sin_x_[((int)yi) + 1];
        yt = fabs(y1 + spatial_bounding_box_.bottom_left().y);
        yi = floor(yt);
        cy1 = yt - yi;
        cy1 = (1 - cy1) * cos_x_[((int)yi)] + cy1 * cos_x_[((int)yi) + 1];
        yt = fabs(y2 + spatial_bounding_box_.bottom_left().y);
        yi = floor(yt);
        cy2 = yt - yi;
        cy2 = (1 - cy2) * cos_x_[((int)yi)] + cy2 * cos_x_[((int)yi) + 1];
        x = fabs(1000 * (y * y + x * x * cy1 * cy2));
        xi = floor(x);
        x -= xi;
        y = (1 - x) * asin2_sq_x_[((int)xi)] + x * asin2_sq_x_[((int)xi) + 1];
        return 4 * earth_radius_ * earth_radius_ * y;
      }
    };

    /// Concrete class for performing periodic distance calculations
    class PeriodicDistance : public Distance
    {
      /// The horizontal and vertical periods in global units
      const Size<double>& global_;

    public:
      /// Construct a class for periodic distance calculations
      /// \param global the horizontal and vertical periods in global units
      PeriodicDistance(const Size<double>& global)
        : global_(global)
      {
      }

      /// Calculate the distance squared between two rectangular grid indices
      /// \param l index
      /// \param m index
      /// \param stride number of indices in a row
      /// \param scales conversion to global co-ordinates
      /// \param minimum whether to calculate the minimum between cells
      /// \return distance squared
      double distance_squared(int l, int m, int stride, const Size<double>& scales, bool minimum) override
      {
        int i = l;
        int j = m;
        if (minimum)
        {
          if (scales.width * ((double)abs(m / stride - l / stride)) > global_.width * 0.5)
          {
            if (m / stride > l / stride)
              j += stride;
            else if (m / stride < l / stride)
              i += stride;
          }
          else
          {
            if (m / stride > l / stride)
              i += stride;
            else if (m / stride < l / stride)
              j += stride;
          }

          if (scales.height * ((double)abs(m % stride - l % stride)) > global_.height * 0.5)
          {
            if (m % stride > l % stride)
              j++;
            else if (m % stride < l % stride)
              i++;
          }
          else
          {
            if (m % stride > l % stride)
              i++;
            else if (m % stride < l % stride)
              j++;
          }
        }

        CovidSim::Geometry::Vector2d u(scales.width * fabs((double)(i / stride)), scales.height * fabs((double)(i % stride)));
        CovidSim::Geometry::Vector2d v(scales.width * fabs((double)(j / stride)), scales.height * fabs((double)(j % stride)));
        return distance_squared(u, v);
      }

      /// Calculate the distance squared between two points
      /// \param a point in global co-ordinates
      /// \param b point in global co-ordinates
      /// \return distance squared
      double distance_squared(const Vector2d& a, const Vector2d& b) override
      {
        auto diff = a - b;
        if (diff.x > global_.width * 0.5)
          diff.x = global_.width - diff.x;

        if (diff.y > global_.height * 0.5)
          diff.y = global_.height - diff.y;

        return diff.length_squared();
      }
    };

    /// Concrete class for performing Eulerian distance calculations
    class EulerDistance : public Distance
    {
    public:
      /// Calculate the distance squared between two rectangular grid indices
      /// \param l index
      /// \param m index
      /// \param stride number of indices in a row
      /// \param scales conversion to global co-ordinates
      /// \param minimum whether to calculate the minimum between cells
      /// \return distance squared
      double distance_squared(int l, int m, int stride, const Size<double>& scales, bool minimum) override
      {
        int i = l;
        int j = m;
        if (minimum)
        {
          if (m / stride > l / stride)
            i += stride;
          else if (m / stride < l / stride)
            j += stride;

          if (m % stride > l % stride)
            i++;
          else if (m % stride < l % stride)
            j++;
        }

        CovidSim::Geometry::Vector2d u(scales.width * fabs((double)(i / stride)), scales.height * fabs((double)(i % stride)));
        CovidSim::Geometry::Vector2d v(scales.width * fabs((double)(j / stride)), scales.height * fabs((double)(j % stride)));
        return distance_squared(u, v);
      }

      /// Calculate the distance squared between two points
      /// \param a point in global co-ordinates
      /// \param b point in global co-ordinates
      /// \return distance squared
      double distance_squared(const Vector2d& a, const Vector2d& b) override
      {
        auto diff = a - b;
        return diff.length_squared();
      }
    };

    std::shared_ptr<Distance> DistanceFactory::create(bool utm, bool periodic_boundaries, const BoundingBox2d& spatial_bounding_box, const Size<double>& global)
    {
      if (utm)
      {
        auto dist = std::make_shared<UTMDistance>(spatial_bounding_box);
        dist->initialise();
        return dist;
      }

      if (periodic_boundaries)
        return std::make_shared<PeriodicDistance>(global);

      return std::make_shared<EulerDistance>();
    }
  }
}
