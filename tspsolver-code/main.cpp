#include "tspsolver.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace TSPSolver;

TMatrix ReadCostMatrix(); // Read in a cost matrix from the user
TMatrix ReadCostMatrixFromFile(std::ifstream& file); // Read in a cost matrix from a file
void ReadMatrixRow(TMatrix& matrix, size_t row, size_t numCities, std::istream& is = std::cin);

void PrintMatrix(const TMatrix& matrix); // Print a cost matrix
void SolveTSP(const TMatrix& matrix); // Solve a cost matrix

int main(int argc, char** argv)
{
    char yn;
    std::cout << "Read const matrix from a file? (Y/N): ";
    std::cin >> yn;

    const bool readFromFile = (yn == 'Y' || yn == 'y');
    TMatrix matrix;

    if (readFromFile) {
        std::string filename;
        std::cout << "Enter filename: ";
        std::cin >> filename;

        std::ifstream ifs;
        ifs.open(filename);
        if (!ifs.good()) {
            std::cerr << "Error opening filename: " << filename;
            exit(EXIT_FAILURE);
        }

        matrix = ReadCostMatrixFromFile(ifs);
        ifs.close();
    } else {
        matrix = ReadCostMatrix();
    }

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
        ReadMatrixRow(matrix, i, numCities);
    }

    return matrix;
}

TMatrix ReadCostMatrixFromFile(std::ifstream& file)
{
    size_t numCities;
    file >> numCities;

    TMatrix matrix(numCities, std::vector<double>(numCities));
    for (size_t i = 0; i < numCities; ++i) {
        ReadMatrixRow(matrix, i, numCities, file);
    }

    return matrix;
}

void ReadMatrixRow(TMatrix& matrix, size_t row, size_t numCities, std::istream& is)
{
    for (size_t i = 0; i < numCities; ++i) {
        std::string tmp;
        is >> tmp;

        if (tmp == "-") {
            matrix[row][i] == INFINITY;
        } else {
            matrix[row][i] = atof(tmp.c_str());
        }
    }
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

int getTotalPath(SStep* solutionRoot) {
  return 0;
}

void SolveTSP(const TMatrix& matrix)
{
    CTSPSolver solver;

    SStep* solutionRoot = solver.solve(matrix.size(), matrix);
    std::string sortedPath = solver.getSortedPath("City");

    std::cout << "sorted path: " << sortedPath << std::endl;

    int totalPathLength = solver.getTotalSteps();
    std::cout << "Total Path Length: " << totalPathLength << std::endl;
}
