#ifndef SHIFTEDLOGWEIGHTQUAD
#define SHIFTEDLOGWEIGHTQUAD

#include <Eigen/Dense>
#include <utility>
#include <vector>

/**
 * This function is used to get the quadrature rule for the weighted integrals
 * of the following type: \f$\int_{0}^{2}log(z^{2}+\epsilon^{2})F(z)\f$ where F
 * is analytic and \f$\epsilon \in (0,1]\f$. The Quadrature rule returned by the
 * function is of a fixed order equal to 9. It is evaluated using Lagrangian
 * interpolation of the Quadrature rule evaluated at 100 Chebychev nodes in the
 * range [0,1].
 *
 * @param eps The value for epsilon in the integral mentioned above
 * @return std::pair type containing Eigen::VectorXd types containing the
 *         weights and nodes
 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> getQuadrature(const double &eps);

std::pair<Eigen::VectorXd, Eigen::VectorXd> getQuadrature2(const double &eps);

#endif // SHIFTEDLOGWEIGHTQUAD
