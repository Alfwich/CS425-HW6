import os
import sys
import re

EXECUTABLE = "./executable.x"
USAGE = "usage: " + sys.argv[0] + " <num_threads>"

TESTAMOUNTS = [50, 250]
TESTFORMAT = "test.%d.txt"
TESTFILES = map(lambda x: TESTFORMAT % x, TESTAMOUNTS)
TESTFILEDIR = "./test/"

def parseTestOutput(output):
    output = output.splitlines()

    duration = [int(s) for s in output[0].split(' ') if s.isdigit()][0]
    pathlen = [int(s) for s in output[-1].split(' ') if s.isdigit()][0]

    return (duration, pathlen)

def runTest(testFilename, numThreads):
    tempfile = ".results.tmp"
    command = EXECUTABLE + " " + testFilename + " " + numThreads + " > " + tempfile
    os.system(command)
    
    with open(tempfile, 'r') as f:
        results = f.read()
        os.remove(tempfile)
        return parseTestOutput(results)

def main():
    if len(sys.argv) < 2:
        print(USAGE)
        sys.exit(1)
    
    numThreads = sys.argv[1]

    for filename in TESTFILES:
        results = [runTest(TESTFILEDIR + filename, numThreads) for i in range(len(TESTFILES))]
        durations = [item[0] for item in results]
        pathlens = [item[1] for item in results]

        avgDuration = sum(durations) / float(len(durations))
        print "average duration for test file {0}: {1} ms".format(filename, str(avgDuration))
        print "path length: {0}".format(pathlens[0])

if __name__ == "__main__":
    main()
