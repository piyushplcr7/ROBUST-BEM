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
#include "logweight_quadrature.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "point_near_singularity.hpp"
#include "shifted_logweight_quadrature.hpp"
#include "single_layer.hpp"
#include <Eigen/Dense>

using namespace parametricbem2d;

int main() {
  std::string file_name = "results_solution_22.txt";
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
  old = evaluate(panel1, panel2, space, 4);
  // Finding out the quadrature order for convergence
  for (int order = orderstart; order < 500; order += 2) {
    neww = evaluateP(panel1, panel2, space, order, eps);
    output << order << std::setw(15) << fabs(exact_value - neww) << std::endl;
    std::cout << std::setw(15) << order << std::setw(15)
              << fabs(exact_value - neww) << std::endl;
  }
  output.close();

  return 0;
}
