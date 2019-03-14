#include "line_near_singularity.hpp"

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
#include "parabolic_parametrization.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "shifted_logweight_quadrature.hpp"
#include "singleLayerPotential.hpp"
#include "single_layer.hpp"
#include "integrate_from.hpp"
#include <Eigen/Dense>

// Function to evaluate the lowest order single layer entry for two panels using
// quadrature of given order
Eigen::MatrixXd evaluateL(const parametricbem2d::AbstractParametrizedCurve &pi,
                         const parametricbem2d::AbstractParametrizedCurve &pi_p,
                         const parametricbem2d::AbstractBEMSpace &space,
                         const unsigned int &N, const double &eps) {
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Lambda expression for functions F and G in \f$\eqref{eq:Vidp}\f$
      auto F = [&](double t) { // Function associated with panel pi_p
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) { // Function associated with panel pi
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      // The usual integrand containing the near singularity is transformed by
      // splitting the logarithm into two parts. Lambda expression for the 1st
      // integrand in this split is given by:
      auto integrand1 = [&](double s, double t) {
        double s_st;
        // Transformed argument of the logarithm
        s_st =
            (pi(s) - pi_p(t)).squaredNorm() / ((s - t) * (s - t) + eps * eps);
        return -1. / 4. / M_PI * log(s_st) * F(t) * G(s); // 1/4 for squarednorm
      };

      // Lambda expression for integrand 1 in z,w coordinates
      auto integrand1_wz = [&](double z, double w) {
        // Including the jacobian of transformation to z,w coordinates = 0.5
        return 0.5 * integrand1(0.5 * (w - z), 0.5 * (w + z)) +
               0.5 * integrand1(0.5 * (w + z), 0.5 * (w - z));
      };

      // transformation of the integral1 to z_pr,w coordinates
      auto integrand1_wzpr = [&](double zpr, double w) {
        // Relation between z_pr and z
        double z = eps * zpr;
        return integrand1_wz(z, w) * eps; // Including the jacobian = eps
      };

      // Integrand 1 presents a Layer function which can be integrated by
      // splitting The integration domain [0,2] into [0,eps] & [eps,2] Function
      // to integrate integrand for z \f$\in\ [0,\epsilon]f$ using quadrature of
      // order N
      auto integrate1_0_e = [&](unsigned N) {
        // points and weights for outer integral
        Eigen::RowVectorXd w_o, p_o;
        std::tie(p_o, w_o) =
            gauleg(0, 1, N, std::numeric_limits<double>::epsilon());
        double integral = 0;
        for (unsigned i = 0; i < N; ++i) {
          double zpr = p_o(i);
          // points and weights for inner integral
          Eigen::RowVectorXd w_i, p_i;
          std::tie(p_i, w_i) = gauleg(-2 + zpr * eps, 2 - zpr * eps, N,
                                      std::numeric_limits<double>::epsilon());
          for (unsigned j = 0; j < N; ++j) {
            double w = p_i(j);
            integral += w_o(i) * w_i(j) * integrand1_wzpr(zpr, w);
          }
        }
        return integral;
      };

      // Function to integrate integrand 1 for z \f$\in\ [\epsilon, 2]f$ using
      // quadrature of order N
      auto integrate1_e_2 = [&](unsigned N) {
        // points and weights for outer integral
        Eigen::RowVectorXd w_o, p_o;
        std::tie(p_o, w_o) =
            gauleg(eps, 2, N, std::numeric_limits<double>::epsilon());
        double integral = 0;
        for (unsigned i = 0; i < N; ++i) {
          double z = p_o(i);
          // points and weights for inner integral
          Eigen::RowVectorXd w_i, p_i;
          std::tie(p_i, w_i) =
              gauleg(-2 + z, 2 - z, N, std::numeric_limits<double>::epsilon());
          for (unsigned j = 0; j < N; ++j) {
            double w = p_i(j);
            integral += w_o(i) * w_i(j) * integrand1_wz(z, w);
          }
        }
        return integral;
      };

      // Getting the quadrature rule for the second part of the transformed
      // integral
      Eigen::VectorXd pts, wts;
      std::tie(pts, wts) = getQuadrature(eps);

      // Evaluating the inner integral of the second part as a function of z
      auto f = [&](double z) {
        // Inner lambda function which keeps z fixed for integrating in w.
        auto inner = [&](double w) {
          return 0.5 * F(0.5 * (w - z)) * G(0.5 * (w + z)) +
                 0.5 * F(0.5 * (w + z)) * G(0.5 * (w - z));
          ;
        };
        // Computing the inner integral with proper limits
        return ComputeIntegral(inner, -2 + z, 2 - z, N);
      };

      // Second integral
      double integral2 = -1. / 4. / M_PI *
                         (log(10) * ComputeIntegral(f, 0, 2, N) +
                          ComputeIntegralFrom(pts, wts, f));

      // Filling up the matrix entry
      interaction_matrix(i, j) =
          integrate1_0_e(N) + integrate1_e_2(N) + integral2;
    }
  }
  return interaction_matrix;
}
