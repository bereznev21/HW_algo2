from shapely.geometry import Point, Polygon, MultiPoint
from main import load_projected_points, load_projected_district, closest_pair_bruteforce, SCALE

import math

points = load_projected_points("points.txt")
rings = load_projected_district("district.txt")

mp = MultiPoint([(p[0], p[1]) for p in points])
hull_poly = mp.convex_hull
hull_vertices_count = len(hull_poly.exterior.coords) - 1
hull_twice_area = int(round(hull_poly.area * 2))

print(f"Число вершин оболочки: {hull_vertices_count}")
print(f"Удвоенная площадь оболочки: {hull_twice_area}")

polys = [Polygon(r) for r in rings]
total_district_twice_area = sum(int(round(p.area * 2)) for p in polys)

inside_shapely = 0
for x, y, _ in points:
    pt = Point(x, y)
    if any(p.contains(pt) for p in polys):
        inside_shapely += 1

print(f"Удвоенная площадь района: {total_district_twice_area}")
print(f"Точек внутри района: {inside_shapely}")


bf_d2, bf_pair, bf_calcs = closest_pair_bruteforce(points)
bf_dist_m = math.sqrt(bf_d2) / SCALE
bf_ids = sorted([bf_pair[0][2], bf_pair[1][2]])
print(f"Полный перебор пар  : d² = {bf_d2} ({bf_dist_m:.4f} м), пара = {bf_ids}")
