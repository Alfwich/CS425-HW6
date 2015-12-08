/*
 *  TSPSG: TSP Solver and Generator
 *  Copyright (C) 2007-2014 Oleksii Serdiuk <contacts[at]oleksii[dot]name>
 *
 *  $Id: $Format:%h %ai %an$ $
 *  $URL: http://tspsg.info/ $
 *
 *  This file is part of TSPSG.
 *
 *  TSPSG is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  TSPSG is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TSPSG.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include "tspsolver.h"

//! \internal \brief A short for maximum double, used internally in the solution algorithm.
#define MAX_DOUBLE std::numeric_limits<double>::max()

const int numThreads = 4;

// Join a list of strings by separator
std::string ListJoin(std::list<std::string>& list, std::string separator=" ") {
  std::ostringstream oss;
  std::copy(list.begin(), std::prev(list.end()), std::ostream_iterator<std::string>(oss, separator.c_str()));
  oss << *list.rbegin();
  return oss.str();
}

std::string getCityName(std::string prefix, int index) {
  std::stringstream ss;
  ss << prefix << " " << index;
  return ss.str();
}

namespace TSPSolver {

/*!
 * \brief Constructs CTSPSolver object.
 * \param parent A parent object.
 */
CTSPSolver::CTSPSolver() {

}

/*!
 * \brief Returns the sorted optimal path, starting from City 1.
 * \param city A string that represents city elements in the path.
 * \param separator A string that represents separators between cities in the path.
 * \return A string, containing sorted optimal path.
 *
 *  The resulting path will be in the form \a city+\a separator+\a city+...+\a separator+\a city.
 *  \c \%1 in \a city will be replaced by the city number.
 */
std::string CTSPSolver::getSortedPath(const std::string &separator) const {
    if (!root || route.empty() || (route.size() != nCities))
        return std::string();

    int i = 0; // We start from City 1
    std::list<std::string> path;
    path.push_back("City 1");

    while ((i = route.at(i)) != 0) {
        path.push_back(getCityName("City", i+1));
    }

    // And finish in City 1, too
    path.push_back("City 1");

    return ListJoin(path, separator);
}

/*!
 * \brief Returns a total number of steps in the current solution.
 * \return A total number of steps or \c 0 if no solution.
 * \note This is not always the same as the number of cities.
 */
int CTSPSolver::getTotalSteps() const {
    return total;
}

int CTSPSolver::getTotalCost() const {
  return totalCost;
}

/*!
 * \brief Indicates whether or not the solution is definitely optimal.
 * \return \c true if the solution is definitely optimal, otherwise \c false.
 *
 *  The solution may need some further iterations to determine whether or not it is optimal.
 *  In such cases this function returns \c false.
 */
bool CTSPSolver::isOptimal() const {
    return !mayNotBeOptimal;
}

/*!
 * \brief Solves the given task.
 * \param numCities Number of cities in the task.
 * \param task The matrix of city-to-city travel costs.
 * \return Pointer to the root of the solution tree.
 *
 * \todo TODO: Comment the algorithm.
 */
SStep* CTSPSolver::solve(int numCities, const TMatrix &task) {
    if (numCities < 3) {
        return NULL;
    }

    nCities = numCities;
    SStep *step = new SStep(); // Initial node for the solution
    step->matrix = task; // Copy the initial cost matrix

    normalize(step->matrix); // Replace the infinity references in the cost matrix
                             // by the maximum available double value.

    step->price = align(step->matrix); // align the cost matrix and assign a lower bound
    root = step; // Sets the root node for the tree

    SStep *left, *right; // Some config
    int nRow, nCol;
    bool firstStep = true;
    double check = INFINITY;
    total = 0;
    while (route.size() < nCities) {
        // For this step setup the alternative paths while also setting the nRow nCol for the next transition
        step->alts = findCandidate(step->matrix, nRow, nCol);

        // Continuously refine the step's cost matrix while we have subcycles
        // Updates the price and alternatives
        while (hasSubCycles(nRow,nCol)) {
            step->matrix[nRow][nCol] = INFINITY;
            step->price += align(step->matrix);
            step->alts = findCandidate(step->matrix,nRow,nCol);
        }

        // A path could not be generated
        if ((nRow == -1) || (nCol == -1)) {
            return NULL;
        }

        // Create the right child for the current step. This will invalidate the nRow and nCol
        // for the right child and update the right child's cost. We also invalidate the nRow, nCol
        // from further consideration
        right = new SStep();
        right->pNode = step;
        right->matrix = step->matrix;
        for (int k = 0; k < nCities; k++) {
            if (k != nCol) {
                right->matrix[nRow][k] = INFINITY;
            }

            if (k != nRow) {
                right->matrix[k][nCol] = INFINITY;
            }
        }
        right->price = step->price + align(right->matrix);

        // Forbid the selected route to exclude its reuse in next steps.
        right->matrix[nCol][nRow] = INFINITY;
        right->matrix[nRow][nCol] = INFINITY;

        // Create the left child for the current step and invalidate nRow and nCol
        left = new SStep();
        left->pNode = step;
        left->matrix = step->matrix;
        left->matrix[nRow][nCol] = INFINITY;
        left->price = step->price + align(left->matrix);

        step->candidate.nRow = nRow;
        step->candidate.nCol = nCol;
        step->plNode = left;
        step->prNode = right;

        // This matrix is not used anymore. Restoring infinities back.
        denormalize(step->matrix);

        // Route with (nRow, nCol) path is cheaper
        if (right->price <= left->price) {
            step->next = SStep::RightBranch;
            step = right;
            route[nRow] = nCol;

            if (firstStep) {
                check = left->price;
                firstStep = false;
            }

        // Route without (nRow,nCol) path is cheaper
        } else {
            step->next = SStep::LeftBranch;
            step = left;

            if (firstStep) {
                check = right->price;
                firstStep = false;
            }
        }

        total++;
    }

    mayNotBeOptimal = (check < step->price);
    totalCost = step->price;

    return root;
}

CTSPSolver::~CTSPSolver()
{
    if (root != NULL) {
        deleteTree(root);
    }
}

/* Privates **********************************************************/

// matrix: The cost matrix to align
// Returns a lower bound for the cost matrix
double CTSPSolver::align(TMatrix &matrix) {
    double r = 0;
    double min;

    // Do row subtraction from the matrix with the min row value per row
    for (int k = 0; k < nCities; k++) {
        min = findMinInRow(k, matrix, numThreads);
        if (min > 0) {
            r += min;
            if (min < MAX_DOUBLE) {
                subRow(matrix, k, min);
            }
        }
    }

    // Do col subtraction from the matrix with the min col value per col
    for (int k = 0; k < nCities; k++) {
        min = findMinInCol(k, matrix, numThreads);
        if (min > 0) {
            r += min;
            if (min < MAX_DOUBLE) {

                subCol(matrix, k, min);
            }
        }
    }

    return (r != MAX_DOUBLE) ? r : INFINITY;
}

// Cleanup Process
void CTSPSolver::deleteTree(SStep *&root) {
    if (root == NULL)
        return;
    SStep *step = root;
    SStep *parent;
    while(true) {
        if (step->plNode != NULL) {
            // We have left child node - going inside it
            step = step->plNode;
            step->pNode->plNode = NULL;
            continue;
        } else if (step->prNode != NULL) {
            // We have right child node - going inside it
            step = step->prNode;
            step->pNode->prNode = NULL;
            continue;
        } else {
            // We have no child nodes. Deleting the current one.
            parent = step->pNode;
            delete step;
            if (parent != NULL) {
                // Going back to the parent node.
                step = parent;
            } else {
                // We came back to the root node. Finishing.
                root = NULL;
                break;
            }
        }
    }
}

// Replaces every instance of MAX_DOUBLE with INFINITY in the cost matrix
// matrix: The cost matrix
void CTSPSolver::denormalize(TMatrix &matrix) const {
    for (int r = 0; r < nCities; r++) {
        for (int c = 0; c < nCities; c++) {
            if (matrix.at(r).at(c) == MAX_DOUBLE) {
                matrix[r][c] = INFINITY;
            }
        }
    }
}

// Return an vector of canidate path selections
std::vector<SStep::SCandidate> CTSPSolver::findCandidate(const TMatrix &matrix, int &nRow, int &nCol) const {
    nRow = -1;
    nCol = -1;
    std::vector<SStep::SCandidate> alts;
    SStep::SCandidate cand;
    double h = -1;
    double sum;

    // For each row and col check for canidate paths
    for (int r = 0; r < nCities; r++) {
        for (int c = 0; c < nCities; c++) {
            // Choose edges that are the min for the row and column
            if (matrix.at(r).at(c) == 0) {
                sum = findMinInRow(r, matrix, numThreads, c) + findMinInCol(c, matrix, numThreads, r);

                // If we found another min value for the row or col then reset the alternatives
                // vector and set the next row and col to this node
                if (sum > h) {
                    h = sum;
                    nRow = r;
                    nCol = c;
                    alts.clear();
                // If there is an equivalent path without subcycles push it into the alternative vector
                } else if ((sum == h) && !hasSubCycles(r ,c)) {
                    cand.nRow = r;
                    cand.nCol = c;
                    alts.push_back(cand);
                }
            }
        }
    }

    return alts;
}

// nCol: The col number
// matrix: The cost matrix
// exr: row to exclude from the min calculation
double CTSPSolver::findMinInCol(int nCol, const TMatrix &matrix, int numThreads, int exr) const {
    double total_min = INFINITY;

    #pragma omp parallel num_threads(numThreads)
    {
      double min = INFINITY;
      #pragma omp for
      for (int k = 0; k < nCities; k++) {
          if ((k != exr) && (min > matrix.at(k).at(nCol))) {
              min = matrix.at(k).at(nCol);
          }
      }

      #pragma omp critical
      total_min = fmin(min, total_min);
    }

    return (total_min == INFINITY) ? 0 : total_min;
}

// nRow: The row number
// matrix: The cost matrix
// exc: column to exclude from the min calculation
double CTSPSolver::findMinInRow(int nRow, const TMatrix &matrix, int numThreads, int exc) const {
    double total_min = INFINITY;

    #pragma omp parallel num_threads(numThreads)
    {
        double min = INFINITY;
        #pragma omp for
        for (int k = 0; k < nCities; k++) {
            if (((k != exc)) && (min > matrix.at(nRow).at(k))) {
                min = matrix.at(nRow).at(k);
            }
        }

        #pragma omp critical
        total_min = fmin(min, total_min);
    }

    return (total_min == INFINITY) ? 0 : total_min;
}

// Detects if a route has sub cycles
bool CTSPSolver::hasSubCycles(int nRow, int nCol) const {
    if ((nRow < 0) || (nCol < 0) || route.empty() || !(route.size() < nCities - 1) || !route.count(nCol)) {
        return false;
    }

    int i = nCol;
    while(true) {
        if ((i = route.at(i)) == nRow) {
            return true;
        }

        if (!route.count(i)) {
            return false;
        }
    }

    return false;
}

// Replaces every instance of INFINITY with MAX_DOUBLE in the cost matrix
// matrix: The cost matrix
void CTSPSolver::normalize(TMatrix &matrix) const {
    for (int r = 0; r < nCities; r++) {
        for (int c = 0; c < nCities; c++) {
            if (matrix.at(r).at(c) == INFINITY) {
                matrix[r][c] = MAX_DOUBLE;
            }
        }
    }
}

// matrix: The cost matrix
// nCol: The col to operate on
// val: The amount to remove from each element in the col
void CTSPSolver::subCol(TMatrix &matrix, int nCol, double val) {
    for (int k = 0; k < nCities; k++) {
        if (k != nCol) {
            matrix[k][nCol] -= val;
        }
    }
}

// matrix: The cost matrix
// nRow: The row to operate on
// val: The amount to remove from each element in the row
void CTSPSolver::subRow(TMatrix &matrix, int nRow, double val) {
    for (int k = 0; k < nCities; k++) {
        if (k != nRow) {
            matrix[nRow][k] -= val;
        }
    }
}

} // end namespace
