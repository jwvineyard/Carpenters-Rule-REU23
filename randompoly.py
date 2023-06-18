import alphashape
import matplotlib.pyplot as plt
import random
import shapely

points = [(random.random(), random.random()) for i in range(100)]
polygon = alphashape.alphashape(points, 7.5)

if type(polygon) == shapely.MultiPolygon:
    for count, poly in enumerate(polygon.geoms):
        print("Plotting polygon " + str(count) + "...")
        print("Coordinates: " + str(list(zip(*poly.normalize().exterior.xy))))
        plt.plot(*poly.exterior.xy)
        plt.show()
else:
    print("Coordinates: " + str(list(zip(*polygon.normalize().exterior.xy))))
    plt.plot(*polygon.exterior.xy)
    plt.show()
