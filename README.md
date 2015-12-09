# Subject
Solution to the Traveling Salesman Problem using a greedy branch-and-bound inclusion-exclusion algorithm.

# Group Members
* Cameron Hall
* Arthur Wuterich

# Contents
There are two implementations included:
* Serial: Adaptation of the algorithm from [TSP Solver Generator](http://tspsg.info/) to execute command line.
* Parallel: Extension of the serial version to adopt OpenMP constructs to make useage of parallel resources.

# Useage
* Both the serial and parallel version had a python driver in the root of the folder to demonstrate their performance.
* To generate new example cost matrices use the gentable.py program in the python folder

# About
```
Sonoma State University: CS425-HW6
Fall 2015
Professor Riviore
```
#Resources
* [Algorithm Basis](https://github.com/leppa/tspsg)
* [TSP Wiki Link](https://simple.wikipedia.org/wiki/Travelling_salesman_problem)
* [B&B Wiki Link](https://en.wikipedia.org/wiki/Branch_and_bound)
* [CS Document covering B&B TSP](http://cs.indstate.edu/cpothineni/alg.pdf)
* [1HR YouTube Video covering B&B TSP(primary)](https://www.youtube.com/watch?v=-cLsEHP0qt0A)
* [1HR YouTube Video covering B&B TSP(alternative)](https://www.youtube.com/watch?v=nN4K8xA8ShM)
