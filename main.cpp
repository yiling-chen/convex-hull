// C++ implementation of Quick Hull algorithm
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

struct Point2D {
    int x, y;

    Point2D (int _x, int _y) : x(_x), y(_y) {}

    friend bool operator<(const Point2D& p1, const Point2D& p2) {
        return p1.x < p2.x ? p1.x < p2.x : p1.y < p2.y;
    }
    friend bool operator==(const Point2D& p1, const Point2D& p2) {
        return p1.x == p2.x && p1.y == p2.y;
    }
};

/**
 * Determine the orientation of p w.r.t the line segment joining p1, p2
 *
 * Slope of line segment L'(p1-p2) : delta = (p2.y - p1.y) / (p2.x - p1.x)
 * Slope of line segment L (p1-p)  : sigma = (p.y  - p1.y) / (p.x  - p1.x)
 *
 * If delta < sigma, p is left to L'
 * If delta = sigma, p, p1, p2 are collinear.
 * If delta < sigma, p is right to L'
 *
 * @param p1
 * @param p2
 * @param p
 * @return 0 (collinear), 1 (right), -1 (left)
 */
int orientation(Point2D p1, Point2D p2, Point2D p) {
    int val = (p2.y - p1.y) * (p.x - p1.x) - (p2.x - p1.x) * (p.y - p1.y);
    if (val == 0) return 0;
    return val > 0 ? -1 : 1;
}

/**
 * Calculate an approximate measure of the distance from p to the line segment joining p1, p2
 *
 * Slope of line segment L'(p1-p2) : delta = (p2.y - p1.y) / (p2.x - p1.x)
 * Slope of line segment L (p1-p)  : sigma = (p.y  - p1.y) / (p.x  - p1.x)
 *
 * The larger abs(delta - sigma) is, the further p is away from L'
 *
 * @param p1
 * @param p2
 * @param p
 * @return approximate distance
 */
int approxDistToLine(Point2D p1, Point2D p2, Point2D p) {
    return abs((p2.y - p1.y) * (p.x - p1.x) - (p2.x - p1.x) * (p.y - p1.y));
}

/**
 *  Square distance between p and q
 * @param p
 * @param q
 * @return square distance
 */
int distSquare(Point2D p, Point2D q)
{
    return (p.y - q.y) * (p.y - q.y) + (p.x - q.x) * (p.x - q.x);
}


/**
 * Given a point set P and a line L formed by p1 and p2.
 * Do convex check by finding all points in P lying on the
 * specified side of L
 *
 * @param P   : input point set
 * @param p1  : end point of L
 * @param p2  : end point of L
 * @param side : specified side of L
 */
void findHull(vector<Point2D> &P, set<Point2D> &hull, Point2D p1, Point2D p2, int side)
{
    int idx = -1;
    int max_dist = 0;

    // finding the point farthest from L and also on the specified side of L.
    for (int i = 0; i < P.size(); i++) {
        int temp = approxDistToLine(p1, p2, P[i]);
        if (orientation(p1, p2, P[i]) == side && temp > max_dist) {
            idx = i;
            max_dist = temp;
        }
    }

    // If no point is on the specified side, p1 and p2 belongs to convex hull.
    // Stop recursion.
    if (idx == -1) {
        hull.insert(p1);
        hull.insert(p2);
        return;
    }

    // Recursively partition P with the two lines joining P[idx] and p1, p2.
    // Note here we have to specify the "outside" part of the triangle formed
    // by P[idx], p1, p2 to perform convex check.
    // For example, the opposite side of p2 w.r.t. line of P[idx] -> p1
    findHull(P, hull, P[idx], p1, -orientation(P[idx], p1, p2));
    findHull(P, hull, P[idx], p2, -orientation(P[idx], p2, p1));
}

/**
 * Quick Hull algorithm
 * @param points : point set P
 * @param hull   : convext hull of P
 */
void quickHull(vector<Point2D> &points, set<Point2D> &hull)
{
    if (points.size() < 3) {
        cout << "Convex hull does not exist.\n";
        return;
    }

    // Finding the point with minimum and maximum x-coordinate
    int min_p = 0, max_p = 0;
    for (int i = 1; i < points.size(); i++) {
        if (points[i].x < points[min_p].x)
            min_p = i;
        if (points[i].x > points[max_p].x)
            max_p = i;
    }

    // Recursively find convex hull points on both sides of
    // the line joining points[min_p] and points[max_p]
    findHull(points, hull, points[min_p], points[max_p], 1);
    findHull(points, hull, points[min_p], points[max_p], -1);
}

/**
 * Given p, q, r are collinear, check if q is on line p-r.
 * @param p
 * @param q
 * @param r
 * @return true (q is on line p-r), false (otherwise)
 */
bool onSegment(Point2D p, Point2D q, Point2D r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
}

/**
 * Line intersection test of line p1-q1 and line p2-q2
 *
 * @param p1
 * @param q1
 * @param p2
 * @param q2
 * @return true (intersect), false (not intersect)
 */
bool isIntersect(Point2D p1, Point2D q1, Point2D p2, Point2D q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // p2 and q2 should be on opposite side of line p1-q1
    // p1 and q1 should be on opposite side of line p2-q2
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases: p1-q1 and p2-q2 partially overlap
    // p1, q1 and p2 are collinear and p2 lies on segment p1-q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1-q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2-q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2-q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false;
}

/**
 * Ray casting algorithm for point in polygon problem
 * https://en.wikipedia.org/wiki/Point_in_polygon
 *
 * For a given point p and polygon S, cast a horizontal line
 * from p to infinity and count how many times the ray intersects
 * with the boundary of S.
 *
 * If odd, p is inside of S.
 * If even, p is outside of S.
 *
 * @param polygon : vector of points forming a polygon
 * @param p : point to test
 * @return  true (inside), false (outside)
 */
bool point_in_polygon(vector<Point2D> &polygon, Point2D p) {
    const auto n = polygon.size();

    // Check polygon size
    if (n < 3)  return false;

    // Create a point in the infinity for ray casting
    Point2D extreme = {numeric_limits<int>::max(), p.y};

    // For each edge of the polygon, check intersection
    int count = 0;
    for (int i = 0; i < n; i++) {
        int j = (i == n-1) ? 0 : i+1;

        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if (isIntersect(polygon[i], polygon[j], p, extreme)) {
            // Boundary case:
            // If p is collinear with line segment p_i to p_j,
            // then check if it lies on segment.
            // If yes (p is on the edge of S), return true,
            // otherwise false
            if (orientation(polygon[i], polygon[j], p) == 0)
                return onSegment(polygon[i], p, polygon[j]);

            count++;
        }
    }

    return count % 2 == 1;  // odd => inside
}

void test_convex_hull() {
    vector<Point2D> outer = {{2, 0}, {-2, 0}, {0, 2}, {0, -2}};
    vector<Point2D> points = outer;

    points.emplace_back(1, 0);
    points.emplace_back(-1, 0);
    points.emplace_back(0, 1);
    points.emplace_back(0, -1);

    set<Point2D> hull;
    quickHull(points, hull);

    // convex hull should equal to the outer points
    assert(hull.size() == outer.size());
    for (auto p : hull) {
        assert(find(outer.begin(), outer.end(), p) != outer.end());
    }

    // newly added points should be in the convex hull
    hull.clear();
    points.emplace_back(0, 3);
    quickHull(points, hull);
    assert(hull.find(Point2D(0, 3)) != hull.end());
    assert(hull.find(Point2D(0, 2)) == hull.end());

    vector<Point2D> vhull;
    for (auto p : hull) {
        vhull.push_back(p);
    }

    assert(point_in_polygon(vhull, Point2D(0, 0)));
    assert(!point_in_polygon(vhull, Point2D(5, 0)));

    // (-1,1) lies on edge of convex hull
    assert(point_in_polygon(vhull, Point2D(-1,-1)));
    // (0,3) is a corner of convex hull
    assert(point_in_polygon(vhull, Point2D(0,3)));

    cout << "[INFO] test_convex_hull() passed!\n";
}

int main()
{
    vector<Point2D> points = {{97, 96}, {77, 195}, {111, 268}, {202, 317},
                             {338, 314}, {388, 209}, {272, 83}, {197, 162},
                              {275, 175}, {163, 210}, {255, 262}};

    cout << "Input points:" << endl;
    for (auto p : points) {
        cout << "(" << p.x << ", " << p.y << "), ";
    }
    cout << endl;

    set<Point2D> hull;
    quickHull(points, hull);

    cout << "Convex Hull:" << endl;
    for (auto p : hull) {
        cout << "(" << p.x << ", " << p.y << "), ";
    }
    cout << endl << endl;

    test_convex_hull();

    return 0;
}
