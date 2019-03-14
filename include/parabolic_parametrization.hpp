#ifndef PARABOLICPARAMHPP
#define PARABOLICPARAMHPP

#include "abstract_parametrized_curve.hpp"
#include<Eigen/Dense>

using namespace parametricbem2d;

class ParabolicParametrization : public AbstractParametrizedCurve {
public:
  using CoefficientsList = typename Eigen::Matrix<double,2,3>;

  ParabolicParametrization(const CoefficientsList& coeffs, double tmin = -1, double tmax = 1);
  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d operator()(double) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d Derivative(double) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  Eigen::Vector2d DoubleDerivative(double) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  PanelVector split(unsigned int) const;

  /**
   * See documentation in AbstractParametrizedCurve
   */
  double length() const;

private:
  const double tmin_;
  const double tmax_;
  const CoefficientsList coeffs_;
};

#endif
