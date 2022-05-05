import sys

fas = sys.argv[1]
out = sys.argv[2]

with open(fas, "r") as f:
    with open(out, "w") as o:
        for line in f:
            if line.startswith(">"):
                line = line.split("|")[0] + "\n"
            o.write(line)
