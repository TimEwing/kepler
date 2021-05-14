from math import sqrt, asin, pi

## Translated from https://stackoverflow.com/a/32698993/8239386
# returns the positive root of intersection of line y = h with circle centered at the origin and radius r
def section(h, r=1):
    assert(r >= 0) # assume r is positive, leads to some simplifications in the formula below (can factor out r from the square root)
    # http://www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+%3D+h
    if h < r:
        return sqrt(r * r - h * h)
    else:
        return 0

# indefinite integral of circle segment
def g(x, h, r=1): 
    # http://www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+-+h
    return .5 * (sqrt(1 - x * x / (r * r)) * x * r + r * r * asin(x / r) - 2 * h * x)

# area of intersection of an infinitely tall box with left edge at x0, right edge at x1, bottom edge at h and top edge at infinity, with circle centered at the origin with radius r
def area_infinite(x0, x1, h, r):
    if x0 > x1:
        x0, x1 = x1, x0 # this must be sorted otherwise we get negative area
    s = section(h, r)
    return g(max(-s, min(s, x1)), h, r) - g(max(-s, min(s, x0)), h, r) # integrate the area

# area of the intersection of a finite box with a circle centered at the origin with radius r
def area_finite(x0, x1, y0, y1, r):
    if y0 > y1:
        y0, y1 = y1, y0 # swap y0,y1 to simplify the reasoning
    if y0 < 0:
        if y1 < 0:
            # the box is completely under, just flip it above and try again
            return area_finite(x0, x1, -y0, -y1, r)
        else:
            # the box is both above and below, divide it to two boxes and go again
            return area_finite(x0, x1, 0, -y0, r) + area_finite(x0, x1, 0, y1, r)
    else:
        # y0 >= 0, which means that y1 >= 0 also (y1 >= y0) because of the swap at the beginning
        assert(y1 >= 0)
        # area of the lower box minus area of the higher box
        return area_infinite(x0, x1, y0, r) - area_infinite(x0, x1, y1, r)

# area of the intersection of a general box with a general circle
def circular_rect_intersect(x0, y0, x1, y1, cx, cy, r):
    # move circle to 0,0
    x0 -= cx
    x1 -= cx
    y0 -= cy
    y1 -= cy

    return area_finite(x0, x1, y0, y1, r)

## Translated from https://stackoverflow.com/a/14646734/8239386
def circular_intersect(x0, y0, r0, x1, y1, r1):
    rr0 = r0 ** 2
    rr1 = r1 ** 2
    d = sqrt((x1 - x0)**2 + (y1 - y0)**2)

    # Circles do not overlap
    if (d > r1 + r0):
        return 0

    # Circle1 is completely inside circle0
    if d <= Math.abs(r0 - r1) and r0 >= r1:
        # Return area of circle1
        return pi * rr1

    # Circle0 is completely inside circle1
    if d <= Math.abs(r0 - r1) and r0 < r1:
        # Return area of circle0
        return pi * rr0

    # Circles partially overlap
    phi = (Math.acos((rr0 + (d * d) - rr1) / (2 * r0 * d))) * 2
    theta = (Math.acos((rr1 + (d * d) - rr0) / (2 * r1 * d))) * 2
    area1 = 0.5 * theta * rr1 - 0.5 * rr1 * Math.sin(theta)
    area2 = 0.5 * phi * rr0 - 0.5 * rr0 * Math.sin(phi)

    # Return area of intersection
    return area1 + area2

## Translated from https://stackoverflow.com/a/42803692/8239386
def circular_intersect_points(x0, y0, r0, x1, y1, r1):
    d = sqrt((x1-x0)**2 + (y1-y0)**2)
    if d == 0:
        return None
    a = (r0**2 - r1**2 + d**2) / (2*d)
    if a**2 > r0**2:
        return None
    h = sqrt(r0**2 - a**2)

    x2 = x0 + a * (x1 - x0) / d
    y2 = y0 + a * (y1 - y0) / d

    p1 = x2 + h*(y1-y0)/d, y2 - h*(x1-x0)/d
    p2 = x2 - h*(y1-y0)/d, y2 + h*(x1-x0)/d
    return p1, p2


## Translated from nowhere because I actually wrote this one. Whoa.
def circular_rect_intersect_points(x0, y0, x1, y1, cx, cy, r):
    if x0 > x1:
        x0,x1 = x1,x0
    if y0 > y1:
        y0,y1 = y1,y0

    def grid_intersect_x(cx, cy, r, x):
        if abs(cx - x) > r:
            return None, None
        h = sqrt(r**2 - (cx - x)**2)
        return cy + h, cy - h

    def grid_intersect_y(cx, cy, r, y):
        if abs(cy - y) > r:
            return None, None
        w = sqrt(r**2 - (cy - y)**2)
        return cx + w, cx - w

    point_pair_list = [
        list(zip([x0, x0], grid_intersect_x(cx, cy, r, x0))),
        list(zip([x1, x1], grid_intersect_x(cx, cy, r, x1))),
        list(zip(grid_intersect_x(cx, cy, r, y0), [y0, y0])),
        list(zip(grid_intersect_x(cx, cy, r, y1), [y1, y1])),
    ]
    # Flatten lists
    points = [x for sublist in point_pair_list for x in sublist if not None in x]
    # Exclude intersects outside the box (since we're looking at lines, not line segments)
    points = [[x,y] for x,y in points if x0<=x<=x1 and y0<=y<=y1]
    return points

def three_way_intersect(x0, y0, x1, y1, cx0, cy0, r0, cx1, cy1, r1):
    x_list = []
    y_list = []

    # Split at circular intersects
    circle_intersects = circular_intersect_points(cx0, cy0, r0, cx1, cy1, r1)
    # Remove intersections outside of box
    if circle_intersects is not None:
        circle_intersects = [[x,y] for x,y in circle_intersects if x0<=x<=x1 and y0<=y<=y1]
    # Add remaining intersections to coord lists
    if circle_intersects:
        circle_x_list, circle_y_list = zip(*circle_intersects)
        x_list += list(circle_x_list)
        y_list += list(circle_y_list)

    # Split at box intersects
    # Calculate intersects and split to x/y lists for circle 0
    intersect_list = circular_rect_intersect_points(x0, y0, x1, y1, cx0, cy0, r0)
    if intersect_list:
        rect_x_list, rect_y_list = zip(*intersect_list)
        # Add to master lists
        x_list += rect_x_list
        y_list += rect_y_list
    # Calculate intersects and split to x/y lists for circle 1
    intersect_list = circular_rect_intersect_points(x0, y0, x1, y1, cx1, cy1, r1)
    if intersect_list:
        rect_x_list, rect_y_list = zip(*intersect_list)
        # Add to master lists
        x_list += rect_x_list
        y_list += rect_y_list

    # Add box coords
    x_list += [x0, x1]
    y_list += [y0, y1]

    # Sort, remove non-unique
    x_list = sorted(list(set(x_list)))
    y_list = sorted(list(set(y_list)))

    # Make new list of boxes
    boxes = []
    for new_x0, new_x1 in zip(x_list[:-1], x_list[1:]):
        for new_y0, new_y1 in zip(y_list[:-1], y_list[1:]):
            boxes.append([new_x0, new_y0, new_x1, new_y1])

    # Sum intersections
    total_area = 0
    for box_x0, box_y0, box_x1, box_y1 in boxes:
        c0_area = circular_rect_intersect(box_x0, box_y0, box_x1, box_y1, cx0, cy0, r0)
        c1_area = circular_rect_intersect(box_x0, box_y0, box_x1, box_y1, cx1, cy1, r1)
        total_area += min(c0_area, c1_area)
    return total_area

def pixel_area(pix_x0, pix_y0, pix_x1, pix_y1, sol_x, sol_y, sol_r, obj_x, obj_y, obj_r):
    solar_area = circular_rect_intersect(
        x0=pix_x0,
        y0=pix_y0,
        x1=pix_x1,
        y1=pix_y1,
        cx=sol_x,
        cy=sol_y,
        r=sol_r,
    )
    object_area = three_way_intersect(
        x0=pix_x0,
        y0=pix_y0,
        x1=pix_x1,
        y1=pix_y1,
        cx0=sol_x,
        cy0=sol_y,
        r0=sol_r,
        cx1=obj_x,
        cy1=obj_y,
        r1=obj_r,
    )
    return solar_area - object_area