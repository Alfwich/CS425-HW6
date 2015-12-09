#include <iostream>
#include "utils.h"

void PrintMatrix(const TSPSolver::TMatrix& matrix)
{
    for (auto row : matrix) {
        for (auto col : row) {
            if(col > 1000) {
              std::cout << "-" << "  ";
            } else {
              if(col >= 10) {
                std::cout << col << " ";
              } else {
                std::cout << col << "  ";
              }
            }
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}
