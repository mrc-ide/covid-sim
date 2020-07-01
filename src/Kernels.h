#ifndef COVIDSIM_KERNELS_H_INCLUDED_
#define COVIDSIM_KERNELS_H_INCLUDED_

#include <vector>

struct Cell;

namespace CovidSim
{
  /// \todo TBD1 is specific to this file and is to be changed to something meaningful during a later refactoring
  namespace TBD1
  {
    /// A probability distribution
    struct KernelStruct
    {
      /// Which distribution to use
      int type_;

      /// Distribution parameter
      double shape_;

      /// Representative distance
      double scale_;

      /// Distribution parameter
      double p3_;

      /// Distribution parameter
      double p4_;

      /// \param r2 The distance squared
      /// \return Probability
      double exponential(double r2) const;

      /// \param r2 The distance squared
      /// \return Probability
      double power(double r2) const;

      /// \param r2 The distance squared
      /// \return Probability
      double power_b(double r2) const;

      /// \param r2 The distance squared
      /// \return Probability
      double power_us(double r2) const;

      /// \param r2 The distance squared
      /// \return Probability
      double power_exp(double r2) const;

      /// Gaussian distribution a.k.a. normal distribution and bell curve
      /// \param r2 The distance squared
      /// \return Probability
      double gaussian(double r2) const;

      /// Step function
      /// \param r2 The distance squared
      /// \return Probability
      double step(double r2) const;
    };

    /// \brief To speed up calculation of kernel values we provide a couple of lookup
    /// tables.
    ///
    /// lookup_ is a table of lookups. lookup_[0] is the kernel
    /// value at a distance of 0, and lookup_.back() is the kernel value at the
    /// largest possible distance.
    ///
    /// hi_res_ is a higher-resolution table of lookups, also of lookup_.size()
    /// elements.  hi_res_[n * expansion_factor_] corresponds to lookup_[n] for
    /// n in [0, lookup_.size() / expansion_factor_]
    ///
    /// Graphically:
    ///
    /// \code
    /// Distance 0 ...                                                 Bound Box diagonal
    /// lookup_[0] ... lookup_[lookup_.size() / expansion_factor_] ... lookup_.back()
    /// hi_res_[0] ... hi_res_.back()
    /// \endcode
    struct KernelLookup
    {
    private:
      /// Kernel lookup table
      std::vector<double> lookup_;

      /// Hi-res kernel lookup table for closer distances
      std::vector<double> hi_res_;

      /// Longest distance / lookup_.size()
      double delta_;

    public:
      /// Size of kernel lookup table
      int size_ = 4000000;

      /// The hi-res kernel extends only to the longest distance / expansion_factor_
      int expansion_factor_;

      /// Resize the vectors and calculate delta_
      /// \param longest_distance The longest distance to lookup
      void setup(double longest_distance);

      /// Set the values in the lookup table
      /// \param norm Value to divide the probability by
      /// \param kernel The kernel to use when calculating the lookup tables
      void init(double norm, KernelStruct& kernel);

      /// Set the values in the cell lookup table
      /// \param lookup The kernel lookup table
      /// \param cell_lookup The cell lookup table
      /// \param cell_lookup_size Number of Cell*
      static void init(const KernelLookup& lookup, Cell **cell_lookup, int cell_lookup_size);

      /// Perform a lookup
      /// \param r2 The distance squared
      /// \return Probability
      double num(double r2) const;
    };
  }
}

#endif // COVIDSIM_KERNELS_H_INCLUDED_
