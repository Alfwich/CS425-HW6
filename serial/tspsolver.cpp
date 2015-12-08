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
#include "utils.h"

//! \internal \brief A short for maximum double, used internally in the solution algorithm.
#define MAX_DOUBLE std::numeric_limits<double>::max()

// Join a list of strings by separator
namespace {
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
}

namespace TSPSolver {

/*!
 * \brief Constructs CTSPSolver object.
 * \param parent A parent object.
 */
CTSPSolver::CTSPSolver() {
  numThreads = 1;
  total = 0;
  totalCost = 0;
}

CTSPSolver::CTSPSolver(int numThreads) {
  this->numThreads = 1;
  total = 0;
  totalCost = 0;
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
    step->matrix = task;       // Copy the initial cost matrix

    // Align the matrix and set the price of the first step to a lower bound for the
    // entire algorithm.
    step->price = align(step->matrix);
    root = step;

    SStep *left, *right;
    int nRow, nCol;
    bool firstStep = true;
    double check = INFINITY;
    total = 0;
    while (route.size() < nCities) {
        // For this step setup the alternative paths while also setting the nRow nCol for the next transition
        step->alts = findCandidate(step->matrix, nRow, nCol);

        // Continually align the matrix while the path has any subcycles. We also
        // eliminate the currently considered path as we can determine this is not within
        // the solution. For each align we need to update the currenly lower bound.
        while (hasSubCycles(nRow,nCol)) {
            step->matrix[nRow][nCol] = INFINITY; // Eliminate path from consideration
            step->price += align(step->matrix);  // Update lower bound
            step->alts = findCandidate(step->matrix,nRow,nCol); // Get new best path
        }

        // A path could not be generated without subcycles; our algorithm has failed
        if ((nRow == -1) || (nCol == -1)) {
            return NULL;
        }

        // Create the inclusion transition. This will compute a new cost matrix with the
        // the selected path as being within the soulution.
        right = new SStep();
        right->pNode = step;
        right->matrix = step->matrix;

        // Remove the to and from paths for the respective cities as we are selecting a path.
        // This effectivly reduces the matrix from a NxN => (N-1)x(N-1)
        for (int k = 0; k < nCities; k++) {
            if (k != nCol) {
                right->matrix[nRow][k] = INFINITY;
            }

            if (k != nRow) {
                right->matrix[k][nCol] = INFINITY;
            }
        }

        // By considering the path compute the new lower bound for the matrix
        right->price = step->price + align(right->matrix);

        // Remove the path to and symetrical path from from the cost matrix
        right->matrix[nCol][nRow] = INFINITY;
        right->matrix[nRow][nCol] = INFINITY;

        // Create the left child for the current step and invalidate nRow and nCol
        left = new SStep();
        left->pNode = step;
        left->matrix = step->matrix;

        // Exclude this path from the solution
        left->matrix[nRow][nCol] = INFINITY;

        // Update the lower bound for this cost matrix with the path excluded
        left->price = step->price + align(left->matrix);

        // Book-keeping for the old GUI program
        step->candidate.nRow = nRow;
        step->candidate.nCol = nCol;
        step->plNode = left;
        step->prNode = right;

        // If the right price is cheaper or equal then we will include the path into the solution
        // and add this transition to the route. Then we will repeat the algorithm from the right child.
        if (right->price <= left->price) {
            step->next = SStep::RightBranch;
            step = right;
            route[nRow] = nCol;

            if (firstStep) {
                check = left->price;
                firstStep = false;
            }

        // The exclusion path is cheaper and we will repeat the algorithm using the left child.
        // IMPORTANT: We do not modify the route
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

    // On out first transition the step price is an estimated lower bound for the entire algorithm.
    // If our final steps price is greater or equal then we cannot guarantee the solution is optimum.
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
        min = findMinInRow(k, matrix);
        if (min > 0) {
            r += min;
            if (min < MAX_DOUBLE) {
                subRow(matrix, k, min);
            }
        }
    }

    // Do col subtraction from the matrix with the min col value per col
    for (int k = 0; k < nCities; k++) {
        min = findMinInCol(k, matrix);
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

// Using the cost matrix provided configure nRow and nCol to locate the canidate that maximises
// the difference between the inclusion and exclusion branch. We return a vector of alternative
// canidate paths
std::vector<SStep::SCandidate> CTSPSolver::findCandidate(const TMatrix &matrix, int &nRow, int &nCol) const {
    nRow = -1;
    nCol = -1;
    std::vector<SStep::SCandidate> alts;
    SStep::SCandidate cand;

    // Our best difference so far
    double h = -1;
    double sum;

    // For each row and col check for canidate paths. Because we have performed aligns to the
    // cost matrix we are only concerned with 0 values
    for (int r = 0; r < nCities; r++) {
        for (int c = 0; c < nCities; c++) {
            if (matrix.at(r).at(c) == 0) {
                // Find the cost change for the exclusion branch by selecting this nRow, nCol
                sum = findMinInRow(r, matrix, c) + findMinInCol(c, matrix, r);

                // If the difference is greater then we want to choose this path.
                if (sum > h) {
                    h = sum;
                    nRow = r;
                    nCol = c;
                    alts.clear();
                // Alternativly we are finding paths of equal difference and will store these. The
                // optimum soulution could use one of these paths rathen than the one we have selected.
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
double CTSPSolver::findMinInCol(int nCol, const TMatrix &matrix, int exr) const {
    double min = INFINITY;
    for (int k = 0; k < nCities; k++) {
        if ((k != exr) && (min > matrix.at(k).at(nCol))) {
            min = matrix.at(k).at(nCol);
        }
    }

    return (min == INFINITY) ? 0 : min;
}

// nRow: The row number
// matrix: The cost matrix
// exc: column to exclude from the min calculation
double CTSPSolver::findMinInRow(int nRow, const TMatrix &matrix, int exc) const {
    double min = INFINITY;

    for (int k = 0; k < nCities; k++) {
        if (((k != exc)) && (min > matrix.at(nRow).at(k))) {
            min = matrix.at(nRow).at(k);
        }
    }

    return (min == INFINITY) ? 0 : min;
}

// Detects if a route has sub cycles by making sure that all destination
// paths do not point to the origin
bool CTSPSolver::hasSubCycles(int nRow, int nCol) const {
    if ((nRow < 0) || (nCol < 0) || route.empty() || !(route.size() < nCities - 1) || !route.count(nCol)) {
        return false;
    }

    int destination = nCol;
    while(true) {
        destination = route.at(destination);
        if (destination == nRow) {
            return true;
        }

        if (!route.count(destination)) {
            return false;
        }
    }

    return false;
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
