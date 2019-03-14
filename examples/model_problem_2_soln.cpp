#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <limits>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "single_layer.hpp"
#include "point_near_singularity.hpp"
#include <Eigen/Dense>

using namespace parametricbem2d;
double epsilon = 0.001;
double eps = epsilon;
// Function to evaluate the lowest order single layer entry for two panels using
// quadrature of given order
double evaluate(AbstractParametrizedCurve &pi, AbstractParametrizedCurve &pi_p,
                AbstractBEMSpace &space, unsigned int N) {
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space

  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // when transforming the parametrizations from [-1,1]->\Pi to local
      // arclength parametrizations [0,|\Pi|] -> \Pi, swap is used to ensure
      // that the common point between the panels corresponds to the parameter 0
      // in both arclength parametrizations
      // bool swap = (pi(1) != pi_p(-1));
      bool swap = true;
      // Panel lengths for local arclength parametrization. Actual values are
      // not required so a length of 1 is used for both the panels
      double length_pi = 1.; // Length for panel pi
      double length_pi_p = 1.;
      /*swap ? pi_p.Derivative(1).norm() / pi.Derivative(-1).norm()
           : pi_p.Derivative(-1).norm() /
                 pi.Derivative(1).norm(); // Length for panel pi_p*/

      // Lambda expressions for the functions F,G and D(r,phi) in eq. 1.4.172
      auto F = [&](double t_pr) { // Function associated with panel pi_p
        // Transforming the local arclength parameter to standard parameter
        // range [-1,1] using swap
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s_pr) { // Function associated with panel pi
        // Transforming the local arclength parameter to standard parameter
        // range [-1,1] using swap
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      auto D_r_phi = [&](double r, double phi) { // eq. 1.4.172
        double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        // Transforming to local arclength parameter range
        double s_pr = r * cos(phi);
        // Transforming to standard parameter range [-1,1] using swap
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        // Transforming to local arclength parameter range
        double t_pr = r * sin(phi);
        // Transforming to standard parameter range [-1,1] using swap
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        /*if (r > tol) // Away from singularity, simply use the formula
          return (pi(s) - pi_p(t)).squaredNorm() / r / r;
        else // Near singularity, use analytically evaluated limit for r -> 0
          //return 1 - sin(2 * phi) * pi.Derivative(s).dot(pi_p.Derivative(t));
          return 4*pi.Derivative(s).squaredNorm() - sin(2 * phi) *
        pi.Derivative(s).dot(pi_p.Derivative(t));*/
        if (r > sqrt_epsilon) // Away from singularity, simply use the formula
          return (pi(s) - pi_p(t)).squaredNorm() / (r * r + epsilon * epsilon);
        else // Near singularity, use analytically evaluated limit for r -> 0
          // return 1 - sin(2 * phi) * pi.Derivative(s).dot(pi_p.Derivative(t));
          return 1.;
      };

      // Getting Gauss Legendre Quadrature weights and nodes
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) =
          gauleg(-1, 1, N, std::numeric_limits<double>::epsilon());

      // The two integrals in eq. 1.4.172 have to be further split into two
      // parts part 1 is where phi goes from 0 to alpha part 2 is where phi goes
      // from alpha to pi/2
      double alpha = atan(length_pi_p / length_pi); // the split point

      // i_IJ -> Integral I, part J
      double i11 = 0., i21 = 0., i12 = 0., i22 = 0.;
      // part 1 (phi from 0 to alpha)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature node into phi
        double phi = alpha / 2 * (1 + points(i));
        // Computing inner integral with fixed phi
        // Inner integral for double integral 1, evaluated with Gauss Legendre
        // quadrature
        double inner1 = 0.;
        // Inner integral for double integral 2, evaluated with Log weighted
        // Gauss quadrature
        double inner2 = 0.;
        // Upper limit for inner 'r' integral
        double rmax = length_pi / cos(phi);
        // Evaluating the inner 'r' integral
        Eigen::RowVectorXd weights_i, points_i;
        std::tie(points_i, weights_i) =
            gauleg(2 * log(epsilon), log(rmax * rmax + epsilon * epsilon), N,
                   std::numeric_limits<double>::epsilon());
        for (unsigned int j = 0; j < N; ++j) {
          // Getting Quadrature weights and nodes for Log weighted Gauss
          // quadrature
          // QuadRule logweightQR = getLogWeightQR(rmax, N);
          // Evaluating inner2 using Log weighted Gauss quadrature
          double x = points_i(j);
          double r = std::sqrt(exp(x) - epsilon * epsilon);
          inner2 +=
              weights_i(j) * x * exp(x) * F(r * sin(phi)) * G(r * cos(phi));

          // Evaluating inner1 using Gauss Legendre quadrature
          r = rmax / 2 * (1 + points(j));
          inner1 += weights(j) * r * log(D_r_phi(r, phi)) * F(r * sin(phi)) *
                    G(r * cos(phi));
        }
        // Multiplying the integral with appropriate constants for
        // transformation to r from Gauss Legendre nodes
        inner1 *= rmax / 2;
        // Multiplying the integrals with appropriate constants for
        // transformation to phi from Gauss Legendre nodes
        i11 += weights(i) * inner1 * alpha / 2;
        i21 += weights(i) * inner2 * alpha / 2;
      }

      // part 2 (phi from alpha to pi/2)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature node into phi (alpha to pi/2)
        double phi =
            points(i) * (M_PI / 2. - alpha) / 2. + (M_PI / 2. + alpha) / 2.;
        // Computing inner integral with fixed phi
        // Inner integral for double integral 1, evaluated with Gauss Legendre
        // quadrature
        double inner1 = 0.;
        // Inner integral for double integral 2, evaluated with Log weighted
        // Gauss quadrature
        double inner2 = 0.;
        // Upper limit for inner 'r' integral
        double rmax = length_pi_p / sin(phi);
        // Evaluating the inner 'r' integral
        Eigen::RowVectorXd weights_i, points_i;
        std::tie(points_i, weights_i) =
            gauleg(2 * log(epsilon), log(rmax * rmax + epsilon * epsilon), N,
                   std::numeric_limits<double>::epsilon());
        for (unsigned int j = 0; j < N; ++j) {
          // Getting Quadrature weights and nodes for Log weighted Gauss
          // quadrature
          // QuadRule logweightQR = getLogWeightQR(rmax, N);
          // Evaluating inner2 using Log weighted Gauss quadrature
          double x = points_i(j);
          double r = std::sqrt(exp(x) - epsilon * epsilon);
          inner2 +=
              weights_i(j) * x * exp(x) * F(r * sin(phi)) * G(r * cos(phi));

          // Evaluating inner1 using Gauss Legendre quadrature
          r = rmax / 2 * (1 + points(j));
          inner1 += weights(j) * r * log(D_r_phi(r, phi)) * F(r * sin(phi)) *
                    G(r * cos(phi));
        }
        // Multiplying the integral with appropriate constants for
        // transformation to r from Gauss Legendre quadrature nodes
        inner1 *= rmax / 2;
        // Multiplying the integrals with appropriate constants for
        // transformation to phi from Gauss Legendre quadrature nodes
        i12 += weights(i) * inner1 * (M_PI / 2. - alpha) / 2.;
        i22 += weights(i) * inner2 * (M_PI / 2. - alpha) / 2.;
      }
      // Summing up the parts to get the final integral
      double integral = 0.5 * (i11 + i12) + 0.25 * (i21 + i22);
      // Multiplying the integral with appropriate constants for transformation
      // to local arclength variables
      integral *= 4 / length_pi / length_pi_p;
      // Filling up the matrix entry
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix(0, 0);
}

int main() {
  std::string file_name = "results_solution_2.txt";
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
    neww = evaluate(panel1, panel2, space, order);
    output << order << std::setw(15) << fabs(exact_value - neww) << std::endl;
    std::cout << "value " << fabs(exact_value - neww) << " for order " << order
              << std::endl;
  }
  output.close();
  return 0;
}
