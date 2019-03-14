#include "parabolic_parametrization.hpp"

#include<cmath>

using namespace parametricbem2d;

ParabolicParametrization::ParabolicParametrization(
    const CoefficientsList &coeffs, double tmin, double tmax)
    : coeffs_(coeffs), tmin_(tmin), tmax_(tmax) {}

Eigen::Vector2d ParabolicParametrization::operator()(double t) const {
  t = t * (tmax_ - tmin_) / 2 +
      (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  Eigen::VectorXd output = Eigen::VectorXd::Zero(2);
  for (unsigned i = 0 ; i < coeffs_.cols() ; ++i) {
    Eigen::Vector2d contri = std::pow(t,i) * coeffs_.col(i);
    //std::cout << "contri for i " << i << " = " << contri << std::endl;
    output += std::pow(t,i) * coeffs_.col(i);
  }
  return output;
}

Eigen::Vector2d ParabolicParametrization::Derivative(double t) const {
  t = t * (tmax_ - tmin_) / 2 +
      (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  Eigen::VectorXd output = Eigen::VectorXd::Zero(2);
  for (unsigned i = 1 ; i < coeffs_.cols() ; ++i) {
    output += i * std::pow(t,i-1) * coeffs_.col(i);
  }
  return output * (tmax_ - tmin_) / 2;
}

Eigen::Vector2d ParabolicParametrization::DoubleDerivative(double t) const {
  t = t * (tmax_ - tmin_) / 2 +
      (tmax_ + tmin_) / 2; // converting to the range [tmin,tmax]
  Eigen::VectorXd output = Eigen::VectorXd::Zero(2);
  for (unsigned i = 2 ; i < coeffs_.cols() ; ++i) {
    output += i * (i-1) * std::pow(t,i-2) * coeffs_.col(i);
  }
  return output * (tmax_ - tmin_) / 2 * (tmax_ - tmin_) / 2;
}

PanelVector ParabolicParametrization::split(unsigned int) const {
  PanelVector nul;
  return nul;
}

double ParabolicParametrization::length() const { return 1; }
