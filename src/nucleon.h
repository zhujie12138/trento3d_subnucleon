// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef NUCLEON_H
#define NUCLEON_H

#include <array>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <iostream>
#include "fast_exp.h"
#include "fwd_decl.h"
#include "random.h"

namespace trento {

class Nucleon;

/// \rst
/// Encapsulates properties shared by all nucleons: transverse thickness
/// profile, cross section, fluctuations.  Responsible for sampling
/// nucleon-nucleon participation with given `\sigma_{NN}`.
/// \endrst
class NucleonProfile {
 public:
  /// Instantiate from the configuration.
  explicit NucleonProfile(const VarMap& var_map);

  /// The radius at which the nucleon profile is truncated.
  double radius() const;

  /// The maximum impact parameter for participation.
  double max_impact() const;

  /// Randomly fluctuate the profile.  Should be called prior to evaluating the
  /// thickness function for a new nucleon.
  void fluctuate();

  /// Corners of the tile enclosing the nucleon thickness./zzzj
  std::array<double, 4> boundary(const Nucleon& nucleon) const;

  /// Compute the thickness function at a (squared) distance from the profile
  /// center.
  double thickness(const Nucleon& nucleon, double x, double y) const;

  /// WK: same as above, but without the Gamma fluctuation, 
  /// used in the calculation of binary collision density 
  double deterministic_thickness(double distance_sqr) const;

  /// WK: return Tpp given bpp^2
  double norm_Tpp(double bpp_sqr) const;

  /// Randomly determine if a pair of nucleons participates.
  bool participate(Nucleon& A, Nucleon& B) const;

 private:
  /// Sample constituent positions inside the nucleon.
  void sample_constituent_positions(Nucleon& nucleon) const;

  /// Fast exponential for calculating the thickness profile.
  const FastExp<double> fast_exp_;

  /// Width of Gaussian thickness function.
  const double width_sqr_;

  /// Gaussian nucleon width.
  const double nucleon_width_;

  /// Gaussian constituent width.
  const double constituent_width_;

  /// Number of constituents inside the nucleon.
  const std::size_t constituent_number_;

  /// Gaussian width of constituent position sampling distribution.
  const double sampling_width_;

  /// Truncate the Gaussian at this radius.
  const double trunc_radius_sqr_;

  /// Maximum impact parameter for participants.
  const double max_impact_sqr_;

  /// Constituent width squared.
  const double constituent_width_sq_;

  /// Calculate thickness out to this distance from constituent center.
  const double constituent_radius_sq_;

  /// Nuclear opacity parameter
  double sigma_partonic_;

  /// Cache (-1/2w^2) for use in the thickness function exponential.
  /// Yes, this actually makes a speed difference.../zj comment: consti width
  const double neg_one_div_two_width_sqr_;

  /// Thickness function prefactor = 1/(constituent_number*2*pi*w^2)
  double prefactor_;
 
  /// bool variable to calcualte Ncoll
  bool with_ncoll_;

  /// WK 1/4w^2
  const double neg_one_div_four_width_sqr_;

  /// WK 1/4pi
  const double one_div_four_pi_;

  /// Dimensionless parameter set to reproduce the inelastic nucleon-nucleon
  /// cross section \sigma_{NN}.  Calculated in constructor.
  double cross_sec_param_;

  /// Fluctuation distribution.
  mutable std::gamma_distribution<double> fluct_dist_;

  /// Gaussian distribution for sampling constituent positions.
  mutable std::normal_distribution<double> constituent_position_dist_;

};

/// \rst
/// Represents a single nucleon.  Stores its position and whether or not it's a
/// participant.  These properties are globally readable, but can only be set
/// through ``Nucleus`` and ``NucleonProfile``.
/// \endrst
class Nucleon {
 public:
  /// Only a default constructor is necessary\---the class is designed to be
  /// constructed once and repeatedly updated.
  Nucleon() = default;

  /// Whether or not this nucleon is a participant.
  bool is_participant() const;

  /// Whether or not its constituents have been sampled.
  bool constituents_exist() const;

  /// The transverse \em x position.
  double x() const;

  /// The transverse \em y position.
  double y() const;

  /// The longitudinal \em z position.
  double z() const;

 private:
  /// A Nucleus must be able to set its Nucleon positions.
  friend class Nucleus;
  friend class NucleonProfile;

  /// The NucleonProfile samples participants so must be able to set
  /// participation status.///correspond to the class NucleonCommon in TRENTo2d
  //friend bool NucleonProfile::participate(Nucleon&, Nucleon&) const;

  /// Set the position and reset participant status to false.
  void set_position(double x, double y, double z);

  /// Internal storage of the position.
  double x_, y_, z_;

  /// Mark as a participant.
  void set_participant(); // Trento2d doesn't have this func

  /// Internal storage of participant status.
  bool participant_ = false;

  /// Whether or not nucleon's constituents are sampled.
  bool constituents_exist_ = false;

  /// Constituent transverse position and fluctuation prefactor.
  struct Constituent {
    double x, y;
    double fluctuation;
  };

  /// Vector of constituent positions and fluctuation prefactors.
  std::vector<Constituent> constituents_;
};

// These functions are short, called very often, and account for a large
// fraction of the total computation time, so request inlining.

// Trivial helper function.
template <typename T>
inline constexpr T sqr(const T& value) {
  return value * value;
}

// Nucleon inline member functions

inline double Nucleon::x() const {
  return x_;
}

inline double Nucleon::y() const {
  return y_;
}

inline double Nucleon::z() const {
  return z_;
}

inline bool Nucleon::is_participant() const {
  return participant_;
}

inline bool Nucleon::constituents_exist() const {
  return constituents_exist_;
}

inline void Nucleon::set_position(double x, double y, double z) {
  x_ = x;
  y_ = y;
  z_ = z;
  participant_ = false;
  constituents_exist_ = false;
}

inline void Nucleon::set_participant() {//corresponding to trento2d's NucleonCommon::set_participant
  participant_ = true;
}

// NucleonProfile inline member functions

inline double NucleonProfile::radius() const {
  return std::sqrt(trunc_radius_sqr_);
}

inline double NucleonProfile::max_impact() const {
  return std::sqrt(max_impact_sqr_);
}

inline std::array<double, 4>
NucleonProfile::boundary(const Nucleon& nucleon) const {
  auto constituent = nucleon.constituents_.begin();

  // Initialize the boundary with the position of the first constituent.
  auto xmin = constituent->x, xmax = constituent->x;
  auto ymin = constituent->y, ymax = constituent->y;

  // Check the remaining constituents and update the boundary accordingly.
  // Using this instead of something from std::algorithm because it finds all
  // four quantities {xmin, xmax, ymin, ymax} in a single pass over the constituents.
  for (std::advance(constituent, 1); constituent != nucleon.constituents_.end(); ++constituent) {
    auto x = constituent->x;
    if (x < xmin)
      xmin = x;
    else if (x > xmax)
      xmax = x;

    auto y = constituent->y;
    if (y < ymin)
      ymin = y;
    else if (y > ymax)
      ymax = y;
  }

  // Remember to add and subtract the constituent radius.
  auto r = std::sqrt(constituent_radius_sq_);

  return {xmin - r, xmax + r, ymin - r, ymax + r};
}


inline void NucleonProfile::fluctuate() {//this fluctuation function is not used in constituent mode
  prefactor_ = fluct_dist_(random::engine) *
     math::double_constants::one_div_two_pi / width_sqr_;
}



inline double NucleonProfile::thickness(
  const Nucleon& nucleon, double x, double y) const {

  auto t = 0.;

  for (const auto& constituent : nucleon.constituents_) {
    auto fluct = constituent.fluctuation;
    auto distance_sq = sqr(x - constituent.x) + sqr(y - constituent.y);
    if (distance_sq < constituent_radius_sq_)
      t += fluct * fast_exp_(distance_sq*neg_one_div_two_width_sqr_);
  }

  return prefactor_ * t;
}

// WK
inline double NucleonProfile::deterministic_thickness(double distance_sqr) const {
  if (distance_sqr > trunc_radius_sqr_)
    return 0.;
  return math::double_constants::one_div_two_pi / width_sqr_ 
		* fast_exp_(neg_one_div_two_width_sqr_*distance_sqr);
}

// WK/zzzj need to figure out
inline double NucleonProfile::norm_Tpp(double bpp_sqr) const  {
  return one_div_four_pi_ / width_sqr_ 
		* fast_exp_(neg_one_div_four_width_sqr_*bpp_sqr);
}

inline bool NucleonProfile::participate(Nucleon& A, Nucleon& B) const {
  // If both nucleons are already participants, there's nothing to do, unless in Ncoll mode
  if (A.is_participant() && B.is_participant() && (! with_ncoll_))
    {//std::cout << "skip one" << std::endl;
        return true;}

  double dx = A.x() - B.x();
  double dy = A.y() - B.y();
  //std::cout<<"("<<A.x()<<","<<A.y()<<");("<<B.x()<<","<<B.y()<<")"<<std::endl;
  double distance_sqr = dx*dx + dy*dy;

  // Check if nucleons are out of range.
  if (distance_sqr > max_impact_sqr_)
    {//std::cout << "false" << std::endl;
    return false;}

  sample_constituent_positions(A);
  sample_constituent_positions(B);

  auto overlap = 0.;
  for (const auto& qA : A.constituents_) {
    for (const auto& qB : B.constituents_) {
      auto distance_sq = sqr(qA.x - qB.x) + sqr(qA.y - qB.y);
      overlap += std::exp(-.25*distance_sq/constituent_width_sq_);
    }
  }
  auto one_minus_prob = std::exp(
          -sigma_partonic_ * prefactor_/(2.*constituent_number_) * overlap
      );
  
  if (one_minus_prob < random::canonical<double>()) {
    A.set_participant();
    B.set_participant();
    return true;
  } 

  return false;
}

inline void NucleonProfile::sample_constituent_positions(Nucleon& nucleon) const {
  if (nucleon.constituents_exist())
    return;

  nucleon.constituents_.resize(static_cast<std::size_t>(constituent_number_));

  double xcom = 0.0;
  double ycom = 0.0;

  // Sample nucleon constituent positions//zjtest
  for (auto&& constituent : nucleon.constituents_) {
    auto xloc = constituent_position_dist_(random::engine);
    auto yloc = constituent_position_dist_(random::engine);

    constituent.x = xloc;
    constituent.y = yloc;
    constituent.fluctuation = fluct_dist_(random::engine);

    xcom += xloc;
    ycom += yloc;
  }

  xcom /= constituent_number_;
  ycom /= constituent_number_;

  // Place nucleon at the desired position
  for (auto&& constituent : nucleon.constituents_) {
    constituent.x += (nucleon.x() - xcom);
    constituent.y += (nucleon.y() - ycom);
  }
  
  //for (auto&& constituent : nucleon.constituents_) {//zjtest
  //  constituent.x = (nucleon.x() - xcom);
  //  constituent.y = (nucleon.y() - ycom);
  //}
  nucleon.constituents_exist_ = true;
}

class MonteCarloCrossSection {
 public:
  MonteCarloCrossSection(const VarMap& var_map);

  double operator() (const double sigma_partonic) const;

 private:
  const double nucleon_width_;
  const double constituent_width_;
  const std::size_t constituent_number_;
  const double sampling_width_;
  const double max_impact_;
  const double constituent_width_sq_;
  const double prefactor_;

  const std::size_t n_max = 10000000;
  const std::size_t cache_size = 1000000;
  const std::size_t n_loops = 10;
  const int n_pass = 10000;
  const double tolerance = 0.001;
};

}  // namespace trento

#endif  // NUCLEON_H
