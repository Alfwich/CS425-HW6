import sys, os, random

def main():
    if len(sys.argv) < 2:
        print "Useage: 'gentable.py {size}'"
        return

    size = int(sys.argv[1])
    print size
    for row in genTable(size):
        print row

def genRow(size, pos):
    row = [ str(random.randint(1,10)) for i in range(size)]
    row[pos] = "-"
    return row

def genTable(size):
    return [" ".join(genRow(size, i)) for i in range(size)]


if __name__ == "__main__":
    main()
