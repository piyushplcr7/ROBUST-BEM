#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "single_layer.hpp"
#include <Eigen/Dense>

using namespace parametricbem2d;

// Function to evaluate the lowest order single layer entry for two panels using
// quadrature of given order
double evaluate(AbstractParametrizedCurve &panel1,
                AbstractParametrizedCurve &panel2, AbstractBEMSpace &space,
                unsigned int order) {
  QuadRule GaussQR = getGaussQR(order);
  Eigen::MatrixXd interaction_matrix =
      single_layer::InteractionMatrix(panel1, panel2, space, GaussQR);
  return interaction_matrix(0, 0);
}

int main() {
  // Preparing the output
  std::string file_name = "mp2_results.txt";
  std::ofstream output(file_name);
  output << "#order" << std::setw(15) << "error" << std::endl;
  // Setting the output precision
  typedef std::numeric_limits<double> dbl;
  std::cout.precision(dbl::max_digits10);
  double tol = 1e-9;
  // Defining the Point Near Singularity problem
  // Defining the circular arc panel
  Eigen::VectorXd center(2);
  center << 0, 0;
  // Quarter circle of radius 1 in first quadrant
  parametricbem2d::ParametrizedCircularArc panel1(center, 1., 0, M_PI / 2);
  // BEM Space used
  DiscontinuousSpace<0> space;
  // Exact integral value (to be calculated from overkill and stored manually)
  double exact_value = 0;
  // The distance between the two panels
  double eps = 0.001;
  // Defining the linear panels
  Eigen::Vector2d point1, point2;
  point1 << 2 + eps, 1;
  point2 << 1 + eps, 0;
  parametricbem2d::ParametrizedLine panel2(point1, point2);
  double old, neww;
  // Looping over quadrature order
  for (int order = orderstart; order < 500; order += 2) {
    neww = evaluate(panel1, panel2, space, order);
    output << order << std::setw(15) << fabs(exact_value - neww) << std::endl;
    std::cout << std::setw(5) << order << std::setw(25)
              << fabs(exact_value - neww) << std::endl;
  }
  output.close();
  return 0;
}
