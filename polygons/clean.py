# Clean and validate the output of RPG.
#
# To clean and validate the polygon data stored in FILE, run
#
#  python3 clean.py FILE
#
# Note that the polygon must be stored in the "line" format by RPG.

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename", metavar = "FILE", help = "File with polygon data.")
args = parser.parse_args()

try:
    with open(args.filename) as f:
        lines = f.readlines()
except IOError:
    print("Error: Could not open \"" + args.filename + "\" for reading.")
    raise SystemExit(1)

if len(lines) <= 4:
    print("Error: Not enough lines in file.")
    raise SystemExit(1)

points = []

for (i, line) in enumerate(lines[1:-1]):
    parts = line.split()

    if len(parts) != 2:
        print("Error: Line " + str(i + 1) + " is not formatted properly.")
        raise SystemExit(1)

    try:
        x = float(parts[0])
        y = float(parts[1])
    except:
        print("Error: Line " + str(i + 1) + " is not formatted properly.")
        raise SystemExit(1)

    points.append((x, y))

try:
    with open(args.filename, "w") as f:
        for (i, point) in enumerate(points):
            line = "{:0.15f} {:0.15f}".format(point[0], point[1])
            f.write(line)

            if i != len(points) - 1:
                f.write("\n")
except IOError:
    print("Error: Could not open \"" + args.filename + "\" for writing.")
    raise SystemExit(1)
