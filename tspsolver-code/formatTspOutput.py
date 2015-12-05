import sys, os

def main():
    if len(sys.argv) < 3:
        print "Useage: 'formatTspOutput.py {inputfile} {size}'"
        return

    size = int(sys.argv[2])
    result = [str(size)]
    with open(sys.argv[1]) as f:
        lineBuffer = []
        for line in f.readlines():
            line = line.strip()
            if len(str(line)) > 0:
                if line == "---":
                    line = "-"
                lineBuffer.append(str(line))
                if len(lineBuffer) == size:
                    result.append( " ".join(lineBuffer))
                    lineBuffer = []

    if len(result) == size+1:
        stringResult = "\n".join(result)
        with open(sys.argv[1],"w") as f:
            f.write(stringResult)
    else:
        print "Size of the result list was not corretly. No changes made."



if __name__ == "__main__":
    main()
