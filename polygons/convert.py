# Convert a polygon data file into a compatible format for Julia.
#
# To convert the polygon data stored in FILE, run
#
#  python3 convert.py FILE
#
# Note that the polygon must be stored in the format specified in "display.py".

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

if len(lines) <= 2:
    print("Error: Not enough lines in file.")
    raise SystemExit(1)

for (i, line) in enumerate(lines):
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

    output = "[{:0.15f}, {:0.15f}]".format(x, y)

    if i != len(lines) - 1:
        output += ","
    else:
        output += "]"

    if i == 0:
        output = "[" + output

    print(output)
