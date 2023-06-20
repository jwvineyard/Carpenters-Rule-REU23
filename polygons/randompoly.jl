using LibGEOS
using Random

function find_intersection(poly)
	n = size(poly, 1)

	for i in range(1, n - 1, step = 1)
		p = LibGEOS.createLineString([poly[i], poly[i + 1]])

		for j in range(i + 2, n - 1, step = 1)
			q = LibGEOS.createLineString([poly[j], poly[j + 1]])

			if LibGEOS.intersects(p, q)
				return (i, j)
			end
		end

		q = LibGEOS.createLineString([poly[n], poly[1]])

		if i != 1 && i != n - 1 && LibGEOS.intersects(p, q)
			return (i, n)
		end
	end

	return (-1, -1)
end

function two_opt_heuristic(vertices)
	poly = shuffle(vertices)

	(i, j) = find_intersection(poly)

	while (i, j) != (-1, -1)
		tmp = poly[i]

		if j == size(poly)[1]
			poly[i] = poly[1]
			poly[1] = tmp
		else
			poly[i] = poly[j + 1]
			poly[j + 1] = tmp
		end

		(i, j) = find_intersection(poly)
	end

	return poly
end
