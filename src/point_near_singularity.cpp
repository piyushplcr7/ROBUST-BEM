#include "point_near_singularity.hpp"

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "integrate_from.hpp"
#include <Eigen/Dense>

Eigen::MatrixXd evaluateP(AbstractParametrizedCurve &pi,
                          AbstractParametrizedCurve &pi_p,
                          AbstractBEMSpace &space, unsigned int N,
                          const double &epsilon) {
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
      bool swap = true; // Assuming Pi(-1) and Pi'(1) are the closest points
      // Panel lengths for local arclength parametrization. Actual values are
      // not required so a length of 1 is used for both the panels
      double length_pi = 1.;   // Length for panel pi
      double length_pi_p = 1.; // Length for panel pi_p

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

        if (r > sqrt_epsilon) // Away from singularity, simply use the formula
          return (pi(s) - pi_p(t)).squaredNorm() / (r * r + epsilon * epsilon);
        else // Near singularity, use analytically evaluated limit for r -> 0
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
        // Evaluating the inner 'r' integral by transforming the integration
        // range to [0,2] in order to use the shifted log weighted quadrature
        double a = rmax / 2;
        double ess = eps / a;
        Eigen::VectorXd pts, wts;
        // Getting the shifted log weighted quadrature rule
        std::tie(pts, wts) = getQuadrature(ess);
        // Number of Quadrature points
        unsigned Ni = pts.rows();
        // Inner integrand without the weight function
        auto inner_integrand = [&](double x) {
          double r = a * x;
          return x * F(r * sin(phi)) * G(r * cos(phi));
        };
        // Residual integral due to normalizing the log argument by 10
        inner2 +=
            a * a * log(10 * a * a) * ComputeIntegral(inner_integrand, 0, 2, N);
        // The normalized integrand computed from the weighted quadrature
        inner2 += a * a * ComputeIntegralFrom(pts, wts, inner_integrand);
        for (unsigned int j = 0; j < N; ++j) {
          // Evaluating inner1 using Gauss Legendre quadrature
          double r = rmax / 2 * (1 + points(j));
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
        // Evaluating the inner 'r' integral by transforming the integration
        // range to [0,2] in order to use the shifted log weighted quadrature
        double a = rmax / 2;
        double ess = eps / a;
        Eigen::VectorXd pts, wts;
        // Getting the shifted log weighted quadrature rule
        std::tie(pts, wts) = getQuadrature(ess);
        // Number of Quadrature points
        unsigned Ni = pts.rows();
        // Inner integrand without the weight function
        auto inner_integrand = [&](double x) {
          double r = a * x;
          return x * F(r * sin(phi)) * G(r * cos(phi));
        };
        // Residual integral due to normalizing the log argument by 10
        inner2 +=
            a * a * log(10 * a * a) * ComputeIntegral(inner_integrand, 0, 2, N);
        // The normalized integrand computed from the weighted quadrature
        inner2 += a * a * ComputeIntegralFrom(pts, wts, inner_integrand);
        for (unsigned int j = 0; j < N; ++j) {
          // Evaluating inner1 using Gauss Legendre quadrature
          double r = rmax / 2 * (1 + points(j));
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
      double integral = 0.5 * (i11 + i12) + 0.5 * (i21 + i22);
      // Multiplying the integral with appropriate constants for transformation
      // to local arclength variables
      integral *= 4 / length_pi / length_pi_p;
      // Filling up the matrix entry
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}
