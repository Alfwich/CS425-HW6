#include "tspsolver.h"
#include <vector>
#include <iostream>

TMatrix ReadCostMatrix(); // Read in a cost matrix from the user
void PrintMatrix(const TMatrix& matrix); // Print a cost matrix
void SolveTSP(const TMatrix& matrix); // Solve a cost matrix

int main(int argc, char** argv)
{
    TMatrix matrix = ReadCostMatrix();

    std::cout << "Read matrix: " << std::endl;
    PrintMatrix(matrix);

    SolveTSP(matrix);

    return 0;
}

TMatrix ReadCostMatrix()
{
    size_t numCities;

    std::cout << "Enter the number of cities: ";
    std::cin >> numCities;

    TMatrix matrix(numCities, std::vector<double>(numCities));

    for (size_t i = 0; i < numCities; ++i) {
        std::cout << "Enter elements of row " << i << " : ";

        for (size_t j = 0; j < numCities; ++j) {
            double tmp;
            std::cin >> tmp;
            matrix[i][j] = tmp;
        }
    }

    return matrix;
}

void PrintMatrix(const TMatrix& matrix)
{
    for (auto row : matrix) {
        for (auto col : row) {
            std::cout << col << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void SolveTSP(const TMatrix& matrix)
{
    CTSPSolver solver;

    SStep* solutionRoot = solver.solve(matrix.size(), matrix);
    std::string sortedPath = solver.getSortedPath("City");

    std::cout << "sorted path: " << sortedPath << std::endl;
}

