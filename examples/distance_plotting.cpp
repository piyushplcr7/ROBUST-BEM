#include "parabolic_parametrization.hpp"

#include<iomanip>
#include<iostream>
#include<Eigen/Dense>
#include<fstream>
#include<string>

int main() {
  std::string filename = "distances.txt";
  std::ofstream output(filename);
  output << std::setw(15) << "#s" << std::setw(15) << "t" << std::setw(15) << "dist" << std::endl;
  double a = 0.3;
  unsigned N = 100;
  Eigen::MatrixXd coeffs1(2,3);
  coeffs1 << 0,1,0,1,0,1;
  std::cout << "coeffs1 \n" << coeffs1 << std::endl;
  Eigen::MatrixXd coeffs2(2,3);
  coeffs2 << 1,1,0,1.4,0,-1;
  std::cout << "coeffs1 \n" << coeffs2 << std::endl;
  ParabolicParametrization parabola1(coeffs1,0,1);
  ParabolicParametrization parabola2(coeffs2,-1,0);
  Eigen::VectorXd ts = Eigen::VectorXd::LinSpaced(N,-1,1);
  for (unsigned i = 0 ; i < N ; ++i) {
    for (unsigned j = 0 ; j < N ; ++j) {
      output << std::setw(15) << ts(i) << std::setw(15) << ts(j) << std::setw(15) << (parabola1(ts(i))-parabola2(ts(j))).norm() << std::endl;
    }
  }
  std::cout << "parabola1(1) = " << parabola1(1) << std::endl;
  std::cout << "parabola1(1) = " << parabola2(1) << std::endl;
  output.close();
  return 0;
}
