#ifndef CITIES_H
#define CITIES_H

using namespace std;

/* debug flag */
#define DEBUG 1

/* City parameters to set up:
 * can do this automatically using Python script
 * 
 * num_cities = number of cities
 * city_names = names of all cities; the first name should be the "hometown"
 *     where tours start and end
 * mileage_table: each row should contain the mileages FROM that city
 *     TO the cities in each column. The distances from A to B and B to A
 *     do not have to be equal (see Python script).
 *  */

const unsigned int num_cities = 5;

const string city_names[num_cities] = {
	"Rohnert Park",
	"San Francisco",
	"Petaluma",
	"Oakland",
	"Santa Rosa"
};

const unsigned int mileage_table[][num_cities] = {
	{0, 48, 11, 52, 10},
	{48, 0, 39, 12, 55},
	{11, 39, 0, 42, 18},
	{52, 12, 42, 0, 58},
	{10, 55, 18, 58, 0}
};

#endif
