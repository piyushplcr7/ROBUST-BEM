#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"
#include "integral_gauss.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "shifted_logweight_quadrature.hpp"
#include "single_layer.hpp"
#include <Eigen/Dense>
#include <limits>
#include <vector>

using namespace parametricbem2d;
double eps = 0.0000000000001;

double f(double x) { return x; }

double integrand(double r) { return log(r * r + eps * eps) * f(r); }

double integrand1(double r) { return log(r * r / eps / eps + 1) * f(r); }

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

int main() {
  typedef std::numeric_limits<double> dbl;
  std::string filename = "gen_weight_test.txt";
  std::cout.precision(dbl::max_digits10);
  std::ofstream out(filename);
  // out << std::setw(30)<< "# order" << std::setw(30) << "gauss" <<
  // std::setw(30) << "gen" <<std::endl; out << std::setw(30)<< "# eps" <<
  // std::setw(30) << "error" <<std::endl;

  /*double exact_value = 0.5 * ( (4+eps*eps)*(log(4+eps*eps)-1) -
  (eps*eps)*(log(eps*eps)-1) ); Eigen::VectorXd pts,wts; std::tie(pts,wts) =
  getQuadrature(eps); std::vector<double> points =
  {0.02367475777864181,0.1394787677804392,0.3463462452121395,0.6247025163185129,0.946163825859407,1.276624941168976,1.579500517231165,1.819101092762548,1.964528016323514};
  std::vector<double> weights =
  {-0.6525916039749501,-1.021896788127096,-1.091252499429839,-0.9899987532944808,-0.8011133386283669,-0.585617593186547,-0.3844881537602865,-0.2178876355301909,-0.08599147153172651};

  for (unsigned order = 2 ; order < 50 ; order = order+2) {
    double gauss_error = fabs(exact_value-ComputeIntegral(integrand,0,2,order));
    double gen_error =
  fabs(exact_value-log(10)*ComputeIntegral(f,0,2,order)-ComputeIntegralFrom(pts,wts,f));
    //double gen_error =
  fabs(exact_value-2*log(eps)*ComputeIntegral(f,0,2,order)-ComputeIntegral(integrand1,0,2,order));
    std::cout << std::setw(3)<< order << std::setw(30) << gauss_error <<
  std::setw(30) << gen_error << std::endl;

  }*/

  Eigen::VectorXd e1 = Eigen::VectorXd::LinSpaced(10000, 0, 0.1);
  Eigen::VectorXd e2 = Eigen::VectorXd::LinSpaced(100, 0.1, 1);
  Eigen::VectorXd ess(e1.rows() + e2.rows() - 1);
  ess << e1.segment(1, e1.rows() - 1), e2;
  /*Eigen::VectorXd ess(14);
  for (unsigned i = 0 ; i < 14 ; ++i) {
     ess(i) = std::pow(0.1,(i+1));
     std::cout << "ess(i) = " << ess(i) << std::endl;
  }*/
  for (unsigned i = 0; i < ess.rows(); ++i) {
    double exact_value =
        0.5 * ((4 + ess(i) * ess(i)) * (log(4 + ess(i) * ess(i)) - 1) -
               (ess(i) * ess(i)) * (log(ess(i) * ess(i)) - 1));
    Eigen::VectorXd pts, wts;
    std::tie(pts, wts) = getQuadrature(ess(i));
    double gen_error =
        fabs(exact_value - log(10) * ComputeIntegral(f, 0, 2, 9) -
             ComputeIntegralFrom(pts, wts, f));
    out << std::setw(30) << ess(i) << std::setw(30) << gen_error << std::endl;
  }
  // QuadRule logweightQR = getLogWeightQR(2,9);
  // std::cout << "weights:" << std::endl << logweightQR.w << std::endl;
  // std::cout << "nodeess:" << std::endl << logweightQR.x << std::endl;
  return 0;
}
