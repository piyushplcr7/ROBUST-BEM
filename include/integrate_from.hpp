#ifndef INTEGRATEFROMHPP
#define INTEGRATEFROMHPP

#include <vector>
#include <Eigen/Dense>

template <typename integrand>
double ComputeIntegralFrom(const std::vector<double> &points,
                           const std::vector<double> &weights,
                           const integrand &integ) {
  unsigned N = points.size();
  double integral = 0.;
  for (unsigned i = 0; i < N; ++i)
    integral += weights[i] * integ(points[i]);
  return integral;
}

template <typename integrand>
double ComputeIntegralFrom(const Eigen::VectorXd &points,
                           const Eigen::VectorXd &weights,
                           const integrand &integ) {
  unsigned N = points.rows();
  double integral = 0.;
  for (unsigned i = 0; i < N; ++i)
    integral += weights(i) * integ(points(i));
  return integral;
}

#endif
