"""
Helper program to generate the mileage table array from a list of cities

Will prompt you for number of city, names of cities, and mileages
and will spit out C code for you to copy and paste

There is no error checking -- be sure that
the number of cities and the mileages are positive ints.

This assumes that the distance between two cities doesn't depend on the order of
the cities -- modify main to pass a False to get_mileages if that's not true.
"""


def get_mileages(city_list, symmetric):
    """
    Get the mileage list from the user
    Parameters:
       city_list: the list of city names
       symmetric: True if distances are the same in both directions
    Returns: 2D mileage list
    """
    N = len(city_list)

    # Initialize list
    mileage_list = [[0 for i in range(N)] for j in range(N)]
    for i in range(len(city_list)):
        for j in range(i + 1, len(city_list)):
            if i != j:
                mileage = int(input('Distance from ' + city_list[i] + ' to ' +
                                    city_list[j] + ': '))
                mileage_list[i][j] = mileage
                if not(symmetric):
                    mileage = int(input('Distance from ' + city_list[i] + ' to ' +
                                         city_list[j] + ': '))
                mileage_list[j][i] = mileage 
        
    return mileage_list


def write_c(city_names, mileages):
    print("\nconst unsigned int num_cities = ", len(city_names), ";", sep="")
    print()
    print("const string city_names[num_cities] = {")
    print(",\n".join(map(lambda x: "\t\"" + x + "\"", city_names)))
    print("};")
    print()
    print("const unsigned int mileage_table[][num_cities] = {")

    for i in range(len(city_names)):
        print("\t{", end="")
        print(", ".join([str(mileages[i][j]) for j in range(len(city_names))]),
              end="")
        if i < len(city_names) - 1:
            print("},")
        else:
            print("}")

    print("};")


def main():
    num_cities = int(input('Number of cities: '))

    # Get names
    city_names = []
    for i in range(num_cities):
        city_names.append(input('Name of city ' + str(i + 1) + ': '))

    mileages = get_mileages(city_names, True)

    write_c(city_names, mileages)


main()
