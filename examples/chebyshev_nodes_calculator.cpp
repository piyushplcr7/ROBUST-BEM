/**
 * This is a simple script to calculate Chebychev Nodes and save them to be
 * directly used in Mathematica
 */
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

int main() {
  unsigned num_nodes = 80;
  typedef std::numeric_limits<double> dbl;
  std::string filename = "chebyshev_nodes.txt";
  std::ofstream output(filename);
  output.precision(dbl::max_digits10);
  double a = 0;
  double b = 1;
  output << "{";
  for (unsigned i = 0; i < num_nodes; ++i) {
    output.setf(std::ios::fixed, std::ios::floatfield);
    output.setf(std::ios::showpoint);
    output << 0.5 * (a + b) +
                  0.5 * (b - a) * std::cos((2 * i + 1) * M_PI / 2 / num_nodes);
    if (i != num_nodes - 1)
      output << ", ";
  }
  output << "}";
  output.close();
  return 0;
}
