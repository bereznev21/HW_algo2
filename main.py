import sys
import json
import math
import argparse
import time
import matplotlib.pyplot as plt

R_EARTH = 6371000.0
SCALE = 10 

def orient(a, b, c):

    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def polygon_area2(ring):
    #Формула шнурования
    n = len(ring)
    if n < 3:
        return 0
    s = 0
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return s


def convex_hull(points):
    #Алгоритм Эндрю (монотонная цепочка)
    unique_map = {}
    for x, y, obj_id in points:
        if (x, y) not in unique_map:
            unique_map[(x, y)] = (x, y, obj_id)

    pts = sorted(unique_map.values(), key=lambda p: (p[0], p[1]))
    n = len(pts)

    if n <= 1:
        return pts, True

    lower = []
    for p in pts:
        while len(lower) >= 2 and orient(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and orient(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    hull = lower[:-1] + upper[:-1]

    is_degenerate = len(hull) < 3 or abs(polygon_area2([(p[0], p[1]) for p in hull])) == 0
    return hull, is_degenerate


def sq_dist(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    return dx * dx + dy * dy


def closest_pair_bruteforce(points):
    #Поиск ближайшей пары полным перебором
    n = len(points)
    min_d2 = float('inf')
    best_pair = (None, None)
    calc_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            calc_count += 1
            d2 = sq_dist(points[i], points[j])
            if d2 < min_d2:
                min_d2 = d2
                best_pair = (points[i], points[j])
    return min_d2, best_pair, calc_count


def closest_pair_dc(points):
    #Поиск ближайшей пары методом «разделяй и властвуй»
    pts_x = sorted(points, key=lambda p: (p[0], p[1]))
    pts_y = sorted(points, key=lambda p: (p[1], p[0]))
    calc_count = 0

    def solve(px, py):
        nonlocal calc_count
        n = len(px)
        if n <= 3:
            min_d2 = float('inf')
            best_pair = (None, None)
            for i in range(n):
                for j in range(i + 1, n):
                    calc_count += 1
                    d2 = sq_dist(px[i], px[j])
                    if d2 < min_d2:
                        min_d2 = d2
                        best_pair = (px[i], px[j])
            return min_d2, best_pair

        mid = n // 2
        mid_point = px[mid]

        left_set = set(id(p) for p in px[:mid])
        py_left = [p for p in py if id(p) in left_set]
        py_right = [p for p in py if id(p) not in left_set]

        d2_left, pair_left = solve(px[:mid], py_left)
        d2_right, pair_right = solve(px[mid:], py_right)

        if d2_left < d2_right:
            min_d2 = d2_left
            best_pair = pair_left
        else:
            min_d2 = d2_right
            best_pair = pair_right

        strip = [p for p in py if (p[0] - mid_point[0]) ** 2 < min_d2]

        m = len(strip)
        for i in range(m):
            for j in range(i + 1, m):
                if (strip[j][1] - strip[i][1]) ** 2 >= min_d2:
                    break
                calc_count += 1
                d2 = sq_dist(strip[i], strip[j])
                if d2 < min_d2:
                    min_d2 = d2
                    best_pair = (strip[i], strip[j])

        return min_d2, best_pair

    min_d2, best_pair = solve(pts_x, pts_y)
    return min_d2, best_pair, calc_count


def point_in_district_raycast(px, py, rings):
    inside = False
    for ring in rings:
        n = len(ring)
        for i in range(n):
            x1, y1 = ring[i]
            x2, y2 = ring[(i + 1) % n]

            if (y1 > py) != (y2 > py):
                if y2 > y1:
                    hit = (px - x1) * (y2 - y1) < (x2 - x1) * (py - y1)
                else:
                    hit = (px - x1) * (y2 - y1) > (x2 - x1) * (py - y1)
                if hit:
                    inside = not inside
    return inside


def load_projected_points(filename):
    points = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            obj_id = parts[0]
            x, y = int(parts[1]), int(parts[2])
            points.append((x, y, obj_id))
    return points


def load_projected_district(filename):
    rings = []
    current_ring = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_ring:
                    rings.append(current_ring)
                    current_ring = []
            else:
                parts = line.split()
                current_ring.append((int(parts[0]), int(parts[1])))
        if current_ring:
            rings.append(current_ring)
    return rings


def cmd_project(args):
    with open(args.points_json, 'r', encoding='utf-8') as f:
        osm_data = json.load(f)

    elements = osm_data.get('elements', [])
    raw_points = []
    skipped = 0
    skip_reasons = {}

    for el in elements:
        el_type = el.get('type')
        el_id = el.get('id')
        obj_id = f"{el_type}/{el_id}"
        if 'lat' in el and 'lon' in el:
            raw_points.append((el['lat'], el['lon'], obj_id))
        elif 'center' in el and 'lat' in el['center'] and 'lon' in el['center']:
            raw_points.append((el['center']['lat'], el['center']['lon'], obj_id))
        else:
            skipped += 1
            reason = "отсутствуют координаты (нет ни lat/lon, ни center)"
            skip_reasons[reason] = skip_reasons.get(reason, 0) + 1

    lats = [p[0] for p in raw_points]
    lons = [p[1] for p in raw_points]
    lat0 = (min(lats) + max(lats)) / 2.0
    lon0 = (min(lons) + max(lons)) / 2.0

    phi0 = math.radians(lat0)
    lam0 = math.radians(lon0)
    cos_phi0 = math.cos(phi0)

    def proj(lat, lon):
        phi = math.radians(lat)
        lam = math.radians(lon)
        x = R_EARTH * (lam - lam0) * cos_phi0
        y = R_EARTH * (phi - phi0)
        return int(round(x * SCALE)), int(round(y * SCALE))

    max_abs_coord = 0
    with open(args.out_points, 'w', encoding='utf-8') as f:
        for lat, lon, obj_id in raw_points:
            ix, iy = proj(lat, lon)
            max_abs_coord = max(max_abs_coord, abs(ix), abs(iy))
            f.write(f"{obj_id} {ix} {iy}\n")

    with open(args.district_geojson, 'r', encoding='utf-8') as f:
        geo_data = json.load(f)

    rings_out = []
    features = geo_data.get('features', [])
    for feat in features:
        geom = feat.get('geometry', {})
        gtype = geom.get('type')
        coords = geom.get('coordinates', [])
        polygons = [coords] if gtype == 'Polygon' else (coords if gtype == 'MultiPolygon' else [])

        for poly in polygons:
            for ring_idx, ring in enumerate(poly):
                # В GeoJSON порядок [долгота, широта]
                pts = [proj(pt[1], pt[0]) for pt in ring]
                # Убираем повторение первой вершины в конце кольца
                if len(pts) > 1 and pts[0] == pts[-1]:
                    pts.pop()
                for ix, iy in pts:
                    max_abs_coord = max(max_abs_coord, abs(ix), abs(iy))

                area2 = polygon_area2(pts)
                # Внешнее кольцо (ring_idx == 0) должно быть против часовой стрелки (area2 > 0)
                # Кольца дыр (ring_idx > 0) — по часовой стрелке (area2 < 0)
                if ring_idx == 0 and area2 < 0:
                    pts.reverse()
                elif ring_idx > 0 and area2 > 0:
                    pts.reverse()
                rings_out.append(pts)

    with open(args.out_district, 'w', encoding='utf-8') as f:
        for idx, ring in enumerate(rings_out):
            if idx > 0:
                f.write("\n")
            for ix, iy in ring:
                f.write(f"{ix} {iy}\n")

    print(f"Всего объектов в выгрузке: {len(elements)}")
    print(f"Пропущено объектов: {skipped}" + (
        f" ({skip_reasons})" if skipped else " (у всех есть координаты или center)"))
    print(f"Центр проекции (lat0, lon0): ({lat0:.7f}, {lon0:.7f})")
    print(f"Единица округления: {'дециметр (0.1 м)' if SCALE == 10 else 'метр (1 м)'}")
    print(f"Наибольшая по модулю координата: {max_abs_coord}")


def cmd_hull(args):
    points = load_projected_points(args.points_file)
    unique_coords = set((p[0], p[1]) for p in points)
    hull, is_degenerate = convex_hull(points)

    print(f"Общее число точек: {len(points)}")
    print(f"Число различных точек: {len(unique_coords)}")

    if is_degenerate:
        print("ВНИМАНИЕ: Выпуклая оболочка вырождена (все точки на одной прямой или точек < 3). Площадь = 0.")

    print(f"Число вершин оболочки: {len(hull)}")
    print("Вершины оболочки (против часовой стрелки):")
    for idx, (x, y, obj_id) in enumerate(hull, 1):
        print(f"  {idx:2d}. {obj_id} ({x}, {y})")

    hull_ring = [(p[0], p[1]) for p in hull]
    twice_area = abs(polygon_area2(hull_ring)) if not is_degenerate else 0
    area_km2 = (twice_area / 2.0) / (SCALE ** 2) / 1e6

    perimeter_units = 0.0
    if len(hull) > 1:
        for i in range(len(hull)):
            p1 = hull[i]
            p2 = hull[(i + 1) % len(hull)]
            perimeter_units += math.hypot(p2[0] - p1[0], p2[1] - p1[1])
        if len(hull) == 2:
            perimeter_units /= 2.0
    perimeter_km = (perimeter_units / SCALE) / 1000.0

    print(f"Удвоенная площадь (в целых ед.): {twice_area}")
    print(f"Площадь оболочки: {area_km2:.4f} км²")
    print(f"Периметр оболочки: {perimeter_km:.4f} км")


def cmd_closest(args):
    points = load_projected_points(args.points_file)
    if len(points) < 2:
        print("Недостаточно точек для поиска пары.")
        return

    t0 = time.perf_counter()
    if args.bruteforce:
        min_d2, best_pair, calcs = closest_pair_bruteforce(points)
        method_name = "Полный перебор пар"
    else:
        min_d2, best_pair, calcs = closest_pair_dc(points)
        method_name = "Разделяй и властвуй"
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    dist_m = math.sqrt(min_d2) / SCALE
    p1, p2 = best_pair
    print(f"Метод: {method_name}")
    print(f"Расстояние: {dist_m:.2f} м (квадрат в целых ед.: {min_d2})")
    print(f"Объект 1: {p1[2]} ({p1[0]}, {p1[1]})")
    print(f"Объект 2: {p2[2]} ({p2[0]}, {p2[1]})")
    print(f"Число вычислений расстояния: {calcs}")
    print(f"Время выполнения: {elapsed_ms:.3f} мс")


def cmd_inside(args):
    points = load_projected_points(args.points_file)
    rings = load_projected_district(args.district_file)

    total_twice_area = 0
    print("Кольца района:")
    for idx, ring in enumerate(rings, 1):
        a2 = polygon_area2(ring)
        total_twice_area += a2
        print(f"  Кольцо {idx} (вершин: {len(ring)}): удвоенная ориентированная площадь = {a2}")

    area_km2 = (total_twice_area / 2.0) / (SCALE ** 2) / 1e6

    inside_count = 0
    for x, y, _ in points:
        if point_in_district_raycast(x, y, rings):
            inside_count += 1
    outside_count = len(points) - inside_count

    print(f"Суммарная удвоенная площадь района: {total_twice_area}")
    print(f"Площадь района: {area_km2:.4f} км²")
    print(f"Точек внутри района: {inside_count}")
    print(f"Точек снаружи района: {outside_count}")


def cmd_plot(args):

    points = load_projected_points(args.points_file)
    rings = load_projected_district(args.district_file)
    hull, _ = convex_hull(points)
    _, best_pair, _ = closest_pair_dc(points)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)

    for idx, ring in enumerate(rings):
        rx = [p[0] / (SCALE * 1000) for p in ring] + [ring[0][0] / (SCALE * 1000)]
        ry = [p[1] / (SCALE * 1000) for p in ring] + [ring[0][1] / (SCALE * 1000)]
        label = "Граница Советского района" if idx == 0 else None
        ax.plot(rx, ry, color='royalblue', linewidth=1.5, label=label)
        ax.fill(rx, ry, color='cornflowerblue', alpha=0.2)

    if len(hull) >= 2:
        hx = [p[0] / (SCALE * 1000) for p in hull] + [hull[0][0] / (SCALE * 1000)]
        hy = [p[1] / (SCALE * 1000) for p in hull] + [hull[0][1] / (SCALE * 1000)]
        ax.plot(hx, hy, color='darkorange', linestyle='--', linewidth=1.8, label="Выпуклая оболочка")

    px = [p[0] / (SCALE * 1000) for p in points]
    py = [p[1] / (SCALE * 1000) for p in points]
    ax.scatter(px, py, color='darkslategray', s=15, alpha=0.8, zorder=3, label="Аптеки (все точки)")

    if best_pair[0] and best_pair[1]:
        p1, p2 = best_pair
        cx = [p1[0] / (SCALE * 1000), p2[0] / (SCALE * 1000)]
        cy = [p1[1] / (SCALE * 1000), p2[1] / (SCALE * 1000)]
        ax.scatter(cx, cy, color='crimson', s=70, marker='*', zorder=5,
                   label=f"Ближайшая пара ({p1[2]}, {p2[2]})")

    ax.set_title("Аптеки Казани и граница Советского района")
    ax.set_xlabel("X (км от центра проекции)")
    ax.set_ylabel("Y (км от центра проекции)")
    ax.set_aspect('equal')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='lower left')

    plt.tight_layout()
    plt.savefig(args.out_image)
    print(f"График успешно сохранён в файл: {args.out_image}")


def main():
    parser = argparse.ArgumentParser(description="Геометрическая утилита анализа объектов OSM")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_proj = subparsers.add_parser("project", help="Проекция координат в целые числа")
    p_proj.add_argument("points_json", help="Входной JSON файл точек из Overpass")
    p_proj.add_argument("district_geojson", help="Входной GeoJSON файл границы района")
    p_proj.add_argument("-p", "--out-points", default="points.txt", help="Выходной файл точек")
    p_proj.add_argument("-d", "--out-district", default="district.txt", help="Выходной файл района")
    p_proj.set_defaults(func=cmd_project)

    p_hull = subparsers.add_parser("hull", help="Построение выпуклой оболочки")
    p_hull.add_argument("points_file", help="Промежуточный файл точек")
    p_hull.set_defaults(func=cmd_hull)

    p_clos = subparsers.add_parser("closest", help="Поиск ближайшей пары")
    p_clos.add_argument("points_file", help="Промежуточный файл точек")
    p_clos.add_argument("-b", "--bruteforce", action="store_true", help="Использовать полный перебор пар")
    p_clos.set_defaults(func=cmd_closest)

    p_ins = subparsers.add_parser("inside", help="Подсчет точек внутри района и его площади")
    p_ins.add_argument("points_file", help="Промежуточный файл точек")
    p_ins.add_argument("district_file", help="Промежуточный файл района")
    p_ins.set_defaults(func=cmd_inside)

    p_plot = subparsers.add_parser("plot", help="Отрисовка карты")
    p_plot.add_argument("points_file", help="Промежуточный файл точек")
    p_plot.add_argument("district_file", help="Промежуточный файл района")
    p_plot.add_argument("-o", "--out-image", default="map.png", help="Выходной файл картинки")
    p_plot.set_defaults(func=cmd_plot)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()