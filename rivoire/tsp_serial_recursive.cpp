#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include "cities.h"
using namespace std;

/* Recursive, depth-first implementation of tree search. */
/* ----------------------------------------------------- */


// Class to keep track of tour so far
//
// Right now, we only have 2 Tour objects (the current Tour,
// and the best Tour we've seen).
//
// Non-recursive implementations will require us to save more tours.
class Tour {
    public:
        // Used to initialize new tours. 
        Tour() { cities.push_back(0); mileage = 0; }

        // Used to initialize best tour. Set to a high
        // cost so that the first complete tour will become
        // the new "best."
        Tour(unsigned int cost) { mileage = cost; } 

        // Copy constructor -- in recursive algorithm,
        // only used to update the "best" tour
	Tour(const Tour& t) {
            cities = t.cities;
            mileage = t.mileage;
        }

        // ----- Accessors -----

        unsigned int get_num_cities() { return cities.size(); }

        // Returns True if City c is in the tour
        bool is_present(unsigned int c) {
            return (find(cities.begin(), cities.end(), c) != cities.end());
        }

        // Returns True if this tour is better than Tour t
        bool better_than(const Tour& t) { return (mileage < t.mileage); }

        // ----- Adding and removing cities -----
        
        // Compute change in mileage if City c is added to the list
        unsigned int mileage_delta(unsigned int c) {
            unsigned int last = cities.back(); // current last city
            unsigned int test_mileage = mileage_table[last][c];
         
            // If this would be our last city, add the mileage
            // to return home
            if (cities.size() + 1 == num_cities) {
                test_mileage += mileage_table[c][0];
            }
            return test_mileage;
         }

 
         // Add city unconditionally (without checking anything) 
         void add_city(unsigned int c) {
             mileage += mileage_delta(c);
             cities.push_back(c);
         }

         // Add city conditionally
         // Check to make sure that (a) it is not already present
         // and (b) it would not lead to a tour that is worse
         // than the best tour so far
         // Return True if city was added
         bool add_city(unsigned int c, const Tour& best_so_far) {
             // If city is already in the list, return false
             if (is_present(c)) {
                 return false;
             }
             // If this would be worse than our best completed tour,
             // return false
             if (mileage + mileage_delta(c) > best_so_far.mileage) {
#if DEBUG
                 cout << "Not adding " << city_names[c]
                      << " -- mileage would be " << mileage + mileage_delta(c)
                      << " -- best is " << best_so_far.mileage << endl;
#endif
                 return false;
             }
             add_city(c);
#if DEBUG
             cout << "Added " << city_names[c] << endl;
             cout << "Tour is now " << *this;
#endif
             return true;      
         }

         // Delete the last city
         void remove_last() {
             unsigned int c = cities.back();
             cities.pop_back();

             // Subtract the mileage that c added to the cost
             mileage -= mileage_delta(c);
#if DEBUG
             cout << "Removed " << city_names[c] << endl;
             cout << "Tour is now " << *this;
#endif
         }

         // Output operator: print cost and list of names
         friend ostream& operator << (ostream &outs, Tour t) {
             cout << "\nCost = " << t.mileage << endl;
             for (vector<unsigned int>::iterator it = t.cities.begin();
                  it != t.cities.end(); it++) {
                 cout << city_names[*it] << endl;
             }
             cout << endl;
             return outs;
         }
	
    private:
        vector<unsigned int> cities;  // Visited cities, in order
        unsigned int mileage; // Total mileage
} best(UINT_MAX); // global for best tour so far - initialized high


/* Do a depth-first search to find the best tour 
 * starting with the partially complete tour t
 */
void DFS(Tour& t) {
    if (!t.better_than(best)) return; // if we've already seen something better

    // If we've added all the cities and we're still here, this must be
    // the best tour
    if (t.get_num_cities() == num_cities) {
        best = t;
#ifdef DEBUG
        cout << "New best:" << best;
#endif

    }
    // Otherwise, try to add the next city and recursively call DFS
    // if successful
    else {
        for (unsigned int i=0; i < num_cities; i++) {
            if (t.add_city(i, best)) {
                DFS(t);
                t.remove_last();
            }
        }
    }
}


int main() { 
    Tour t;              // create an initial tour
    DFS(t);
    cout << best;
}
