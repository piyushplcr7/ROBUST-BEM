#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"
#include "integral_gauss.hpp"
#include "line_near_singularity.hpp"
#include "logweight_quadrature.hpp"
#include "parabolic_parametrization.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "shifted_logweight_quadrature.hpp"
#include "singleLayerPotential.hpp"
#include "single_layer.hpp"
#include <Eigen/Dense>

using namespace parametricbem2d;

int main() {
  typedef std::numeric_limits<double> dbl;
  std::cout.precision(dbl::max_digits10);
  double tol = 1e-9;
  double eps = 0.071612800000000004;
  DiscontinuousSpace<0> space;
  std::string file_name = "results_solution_3.txt";
  std::ofstream output(file_name);
  output << "#order" << std::setw(15) << "error" << std::endl;
  int orderstart = 2;
  Eigen::Vector2d point1, point2;
  point1 << 0, 0;
  point2 << -eps, 1;
  parametricbem2d::ParametrizedLine panel1(point1, point2);
  // Panel 2
  Eigen::Vector2d point3, point4;
  point3 << eps, 0;
  point4 << 2 * eps, 1;
  parametricbem2d::ParametrizedLine panel2(point3, point4);
  double exact_value = computeVij(point1, point2, point3, point4, 0.);
  exact_value = 0.23802816757460257;
  double old, neww;
  int convergence_order;
  Eigen::MatrixXd coeffs1(2, 3);
  coeffs1 << 0, 1, 0, 1, 0, 1;
  std::cout << "coeffs1 \n" << coeffs1 << std::endl;
  Eigen::MatrixXd coeffs2(2, 3);
  coeffs2 << 1, 1, 0, 1.4, 0, -1;
  std::cout << "coeffs1 \n" << coeffs2 << std::endl;
  ParabolicParametrization parabola1(coeffs1, 0, 1);
  ParabolicParametrization parabola2(coeffs2, -1, 0);
  parametricbem2d::ParametrizedCircularArc c1(Eigen::Vector2d(0, 0), 1.,
                                              -M_PI / 2, M_PI / 2);
  parametricbem2d::ParametrizedCircularArc c2(Eigen::Vector2d(0, 0), 1. + eps,
                                              -M_PI / 2, M_PI / 2);
  // Finding out the quadrature order for convergence
  for (int order = 2; order < 200; order += 2) {
    double i1, i2;
    neww = evaluateL(parabola1, parabola2, space, order, eps)(0, 0);
    output << order << std::setw(15) << fabs(exact_value - neww) << std::endl;
    std::cout << order << "          " << fabs(exact_value - neww) << std::endl;
  }
  output.close();
  return 0;
}

/*// integral of type int_{0}^{a} log(z+eps) F(z) dz
template<typename T>
 double shifted_logweighted(const T& integrand,double a, double epsilon,unsigned
order) { std::function<double(double)> shifted_integrand = [&] (double x) {
     return integrand(x-std::sqrt(epsilon));
   };
   return ComputeLogIntegral(shifted_integrand,a+epsilon,order) -
ComputeLogIntegral(shifted_integrand,epsilon,order);
 }*/
