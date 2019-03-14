#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "singleLayerPotential.hpp"
#include "single_layer.hpp"
#include "parabolic_parametrization.hpp"
#include <Eigen/Dense>

using namespace parametricbem2d;

// Function to evaluate the lowest order single layer entry for two panels using
// quadrature of given order
double evaluate(AbstractParametrizedCurve &panel1,
                AbstractParametrizedCurve &panel2, AbstractBEMSpace &space,
                unsigned int order) {
  //
  QuadRule GaussQR = getGaussQR(order);
  Eigen::MatrixXd interaction_matrix = single_layer::ComputeIntegralGeneral(
      panel1, panel2, space, GaussQR);
  return interaction_matrix(0, 0);
}

int main() {
  typedef std::numeric_limits<double> dbl;
  std::cout.precision(dbl::max_digits10);
  DiscontinuousSpace<0> space;
  std::string file_name = "mp3_results.txt";
  std::ofstream output(file_name);
  output << "#order" << std::setw(15) << "error" << std::endl;
  double eps = .001;
  // Panel 1
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
  parametricbem2d::ParametrizedCircularArc c1(Eigen::Vector2d(0,0),1.,-M_PI/2,M_PI/2);
  parametricbem2d::ParametrizedCircularArc c2(Eigen::Vector2d(0,0),1.+eps,-M_PI/2,M_PI/2);
  double old, neww;
  int convergence_order;
  Eigen::MatrixXd coeffs1(2,3);
  coeffs1 << 0,1,0,1,0,1;
  std::cout << "coeffs1 \n" << coeffs1 << std::endl;
  Eigen::MatrixXd coeffs2(2,3);
  coeffs2 << 1,1,0,1.4,0,-1;
  std::cout << "coeffs1 \n" << coeffs2 << std::endl;
  ParabolicParametrization parabola1(coeffs1,0,1);
  ParabolicParametrization parabola2(coeffs2,-1,0);
  for (int order = 2; order < 500; order += 2) {
    neww = evaluate(parabola1, parabola2, space, order);
    output << order << std::setw(15) << fabs(exact_value - neww) << std::endl;
    //std::cout << std::setw(5) << order << std::setw(25) << neww << std::endl;
    std::cout << std::setw(5) << order << std::setw(25) << fabs(exact_value - neww) << std::endl;
  }
  output.close();

  return 0;
}
