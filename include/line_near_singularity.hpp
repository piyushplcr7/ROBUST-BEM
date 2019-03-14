#ifndef LINENS
#define LINENS

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include <Eigen/Dense>

/**
 * This function is used to evaluate integrals induced by Single Layer BIO in
 * BEM formulations in the case of a line near-singularity, which is the case
 * when the point of minimum distance between two smooth non linear panels lies
 * strictly inside the parameter domain.
 *
 * This function assumes that the point of closest distance lies on the diagonal
 * of parameter domain, which can be achieved through mapping. Evaluation is
 * done using transformation techniques. This function uses the framework
 * provided by 2D-parametric BEM Library.
 *
 * @param pi The first parametrized panel, derived from
 * AbstractParametrizedCurve
 * @param pi_p The second parametrized panel, derived from
 * AbstractParametrizedCurve
 * @param space Trial/Test BEM space for evaluating the integral
 * @param i Index for the first shape function
 * @param j Index for the second shape function
 * @param N Quadrature order to be used for evaluation
 * @param eps Minimum distance between the two panels
 * @return double type with the integral value
 */

Eigen::MatrixXd
evaluateL(const parametricbem2d::AbstractParametrizedCurve &pi,
         const parametricbem2d::AbstractParametrizedCurve &pi_p,
         const parametricbem2d::AbstractBEMSpace &space,
         const unsigned int &N, const double &eps);

#endif
