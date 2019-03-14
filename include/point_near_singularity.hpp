#ifndef POINTNS
#define POINTNS

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include <Eigen/Dense>

/**
 * This function is used to evaluate integrals induced by Single Layer BIO in
 * BEM formulations in the case of a point near-singularity, which is the case
 * when the point of minimum distance between two panels is at their ends.
 * Evaluation is done using transformation techniques. This function uses the
 * framework provided by 2D-parametric BEM Library.
 *
 * IMPORTANT: This function assumes the following condition: \f$\Pi(-1)\f$ and
 * \f$\Pi'(1)\f$ are the closest points with distance equal to \f$\epsilon\f$
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
evaluateP(const parametricbem2d::AbstractParametrizedCurve &pi,
         const parametricbem2d::AbstractParametrizedCurve &pi_p,
         const parametricbem2d::AbstractBEMSpace &space,
         const unsigned int &N, const double &epsilon);

#endif
