# Displays a polygon stored in a file.
#
# To display the polygon stored in FILE, run the command
#
#   python3 display.py FILE
#
# The file should be formatted as follows, where (xi,yi) is the ith vertex of the
# polygon:
#
#   x1 y1
#   x2 y2
#   .. ..
#   xN yN
#
# There should be no other lines in the file, and no repeated vertices.

import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("filename", metavar = "FILE", help = "File with polygon data.")

args = parser.parse_args()

try:
    with open(args.filename) as f:
        lines = f.readlines()
except IOError:
    print("Error: Could not open \"" + args.filename + "\" for reading.")
    raise SystemExit(1)

if len(lines) <= 1:
    print("Error: Not enough lines in file.")
    raise SystemExit(1)

xPoints = []
yPoints = []

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

    xPoints.append(x)
    yPoints.append(y)

plt.axes().set_aspect("equal")
plt.axis("off")
plt.fill(xPoints, yPoints, facecolor = "violet", edgecolor = "purple", linewidth = 1)
plt.show()
