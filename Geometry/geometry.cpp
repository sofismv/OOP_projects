#define _USE_MATH_DEFINES
#include "geometry.h"
#include <cmath>
#include <vector>
#include <iostream>

///////////////////////////////////////// POINT /////////////////////////////////////////

bool areSame(const double &first, const double &second) {
  return std::abs(first - second) < epsilon;
}

bool Point::operator==(const Point &other) const {
  return areSame(x, other.x) && areSame(y, other.y);
}

bool Point::operator!=(const Point &other) const {
  return !(*this == other);
}
double findDistance(const Point &first, const Point &second) {
  return sqrt((first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
}

void rotatePoint(Point &point, const Point &center, double &angle) {
  double x = center.x + cos(toRadians(angle)) * (point.x - center.x) - sin(toRadians(angle)) * (point.y - center.y);
  double y = center.y + sin(toRadians(angle)) * (point.x - center.x) + cos(toRadians(angle)) * (point.y - center.y);
  point.x = x;
  point.y = y;
}
void reflexPoint(Point &point, const Point &center) {
  double add_x = center.x - point.x;
  double add_y = center.y - point.y;
  point.x += 2 * add_x;
  point.y += 2 * add_y;
}
void scalePoint(Point &point, const Point &center, double &coef) {
  Point vector(point.x - center.x, point.y - center.y);
  vector.x *= coef;
  vector.y *= coef;
  point.x = center.x + vector.x;
  point.y = center.y + vector.y;
}

///////////////////////////////////////// LINE /////////////////////////////////////////

Line::Line(const Point &first, const Point &second) {
  a = first.y - second.y;
  b = second.x - first.x;
  c = -second.x * first.y + second.y * first.x;
}
Line::Line(double slope, double shift) {
  a = -slope;
  b = 1;
  c = -shift;
}
Line::Line(const Point &point, double slope) {
  Point second(point.x + 1, point.y + slope);
  *this = Line(point, second);
}
bool Line::operator==(const Line &other) const {
  if (areSame(a / b, other.a / other.b) && areSame(b / c, other.b / other.c)) {
    return true;
  }
  return false;
}
bool Line::operator!=(const Line &other) const {
  return !(*this == other);
}
double Line::getSlope() const {
  return -(a / b);
}
Point findIntersectionPoint(const Line &first, const Line &second) {
  double determinant = first.a * second.b - first.b * second.a;
  double x = (first.b * second.c - first.c * second.b) / determinant;
  double y = (-first.a * x - first.c) / first.b;
  return Point(x, y);
}
void reflexPoint(Point &point, const Line &line) {
  Line perpendicular(point, -1 / line.getSlope());
  const Point center = findIntersectionPoint(line, perpendicular);
  reflexPoint(point, center);
}

///////////////////////////////////////// ELLIPSE /////////////////////////////////////////

double Ellipse::a() const {
  return dist / 2;
}
double Ellipse::b() const {
  return sqrt(a() * a() - c() * c());
}
double Ellipse::c() const {
  return findDistance(focus.first, focus.second) / 2;
}
Ellipse::Ellipse() {
  focus = std::make_pair(Point(0, 0), Point(0, 0));
  dist = 0;
}
Ellipse::Ellipse(Point focus1, Point focus2, double _dist) {
  focus = std::make_pair(focus1, focus2);
  dist = _dist;
}
std::pair<Point, Point> Ellipse::focuses() const {
  return focus;
}
std::pair<Line, Line> Ellipse::directrices() const {
  double max_from_a_and_b = std::max(a(), b());
  Line directrice_1(1, 0, max_from_a_and_b / eccentricity());
  Line directrice_2(1, 0, -max_from_a_and_b / eccentricity());
  return std::make_pair(directrice_1, directrice_2);
}
double Ellipse::eccentricity() const {
  return c() / a();
}
Point Ellipse::center() const {
  return Point((focus.first.x + focus.second.x) / 2, (focus.first.y + focus.second.y) / 2);
}
double Ellipse::perimeter() const {
  return M_PI * (3 * (a() + b()) - sqrt((3 * a() + b()) * (a() + 3 * b())));
}
double Ellipse::area() const {
  return M_PI * a() * b();
}
bool Ellipse::operator==(const Shape &another) const {
  Ellipse other;
  try {
    other = dynamic_cast<const Ellipse &>(another);
  } catch (...) {
    return false;
  }
  if (((focus.first == other.focus.first && focus.second == other.focus.second)
      || (focus.first == other.focus.second && focus.second == other.focus.first)) && areSame(dist, other.dist)) {
    return true;
  }
  return false;
}
bool Ellipse::operator!=(const Shape &another) const {
  return !(*this == another);
}
bool Ellipse::isSimilarTo(const Shape &another) const {
  Ellipse other;
  try {
    other = dynamic_cast<const Ellipse &>(another);
  } catch (...) {
    return false;
  }
  if (a() / b() == other.a() / other.b() || a() / b() == other.b() / other.a()) {
    return true;
  }
  return false;
}
bool Ellipse::isCongruentTo(const Shape &another) const {
  return isSimilarTo(another) && area() == another.area();
}
double Ellipse::focal_distance() const {
  return dist;
}
bool Ellipse::containsPoint(Point point) const {
  if (findDistance(point, focus.first) + findDistance(point, focus.second) < focal_distance() + epsilon) {
    return true;
  }
  return false;
}
void Ellipse::rotate(Point center, double angle) {
  rotatePoint(focus.first, center, angle);
  rotatePoint(focus.second, center, angle);
}
void Ellipse::reflex(Point center) {
  reflexPoint(focus.first, center);
  reflexPoint(focus.second, center);
}
void Ellipse::reflex(Line axis) {
  reflexPoint(focus.first, axis);
  reflexPoint(focus.second, axis);
}
void Ellipse::scale(Point center, double coefficient) {
  scalePoint(focus.first, center, coefficient);
  scalePoint(focus.first, center, coefficient);
  dist *= coefficient;
}

///////////////////////////////////////// CIRCLE /////////////////////////////////////////


double Circle::radius() const {
  return focal_distance() / 2;
}
void Circle::reflex(Point _center) {
  Point new_center = focuses().first;
  reflexPoint(new_center, _center);
  *this = Circle(new_center, focal_distance() / 2);
}
void Circle::rotate(Point _center, double angle) {
  Point new_center = focuses().first;
  rotatePoint(new_center, _center, angle);
  *this = Circle(new_center, focal_distance() / 2);
}
void Circle::reflex(Line axis) {
  Point new_center = focuses().first;
  reflexPoint(new_center, axis);
  *this = Circle(new_center, focal_distance() / 2);
}
void Circle::scale(Point _center, double coefficient) {
  Point new_center = focuses().first;
  scalePoint(new_center, _center, coefficient);
  *this = Circle(new_center, focal_distance() / 2);
}
double Circle::perimeter() const {
  return M_PI * focal_distance();
}
double Circle::area() const {
  return M_PI * (focal_distance() / 2) * (focal_distance() / 2);
}

///////////////////////////////////////// POLYGON /////////////////////////////////////////

Polygon::Polygon(const std::vector<Point> &points) {
  vertices = points;
}
size_t Polygon::verticesCount() const {
  return vertices.size();
}
std::vector<Point> Polygon::getVertices() const {
  return vertices;
}
double Polygon::perimeter() const {
  double perimeter = 0;
  for (size_t i = 0; i < verticesCount() - 1; ++i) {
    perimeter += findDistance(vertices[i], vertices[i + 1]);
  }
  perimeter += findDistance(vertices[0], vertices[verticesCount() - 1]);
  return perimeter;
}
double Polygon::area() const {
  double area = 0.0;
  for (size_t i = 0; i < vertices.size(); ++i) {
    area += (vertices[i].x + vertices[(i + 1) % verticesCount()].x)
        * (vertices[(i + 1) % verticesCount()].y - vertices[i].y);
  }
  area /= 2;
  return abs(area);
}
bool Polygon::isConvex() const {
  size_t n = verticesCount();
  if (n < 4) {
    return true;
  }
  bool sign = false;
  for (size_t i = 0; i < n; ++i) {
    double crossProduct =
        (vertices[(i + 2) % n].x - vertices[(i + 1) % n].x) * (vertices[i].y - vertices[(i + 1) % n].y)
            - (vertices[(i + 2) % n].y - vertices[(i + 1) % n].y) * (vertices[i].x - vertices[(i + 1) % n].x);
    if (i == 0) {
      sign = crossProduct > 0;
    } else {
      if (sign != (crossProduct > 0)) {
        return false;
      }
    }
  }
  return true;
}
void Polygon::rotate(Point center, double angle) {
  for (size_t i = 0; i < verticesCount(); ++i) {
    rotatePoint(vertices[i], center, angle);
  }
}
void Polygon::reflex(Point center) {
  for (size_t i = 0; i < verticesCount(); ++i) {
    reflexPoint(vertices[i], center);
  }
}
void Polygon::reflex(Line axis) {
  for (size_t i = 0; i < verticesCount(); ++i) {
    reflexPoint(vertices[i], axis);
  }
}
void Polygon::scale(Point center, double coefficient) {
  for (size_t i = 0; i < verticesCount(); ++i) {
    scalePoint(vertices[i], center, coefficient);
  }
}
bool Polygon::operator==(const Shape &another) const {
  Polygon shape;
  try {
    shape = dynamic_cast<const Polygon &>(another);
  } catch (...) {
    return false;
  }
  if (shape.verticesCount() != verticesCount()) {
    return false;
  }
  size_t n = verticesCount();
  std::vector<Point> vertices = getVertices();
  std::vector<Point> other_vertices = shape.getVertices();
  size_t first_vertex = n;
  for (size_t i = 0; i < n; ++i) {
    if (other_vertices[i] == vertices[0]) {
      first_vertex = i;
      break;
    }
  }
  if (first_vertex == n) {
    return false;
  }
  if (other_vertices[(first_vertex + 1) % n] == vertices[1]) {
    for (size_t i = 0; i < n; ++i) {
      if (other_vertices[(first_vertex + i) % n] != vertices[i]) {
        return false;
      }
    }
    return true;
  } else {
    for (size_t i = 0; i < n; ++i) {
      if (other_vertices[(n + first_vertex - i) % n] != vertices[i]) {
        return false;
      }
    }
    return true;
  }
}
bool Polygon::operator!=(const Shape &another) const {
  return !(*this == another);
}
bool Polygon::containsPoint(Point point) const {
  double sum = 0;
  for (size_t i = 0; i < verticesCount(); ++i) {
    Point A = vertices[i];
    Point B = vertices[(i + 1) % verticesCount()];
    Line A_point(A, point);
    Line B_point(B, point);
    double cos_angle = (A_point.a * B_point.a + A_point.b * B_point.b)
        / (sqrt(A_point.a * A_point.a + A_point.b * A_point.b) * sqrt(B_point.a * B_point.a + B_point.b * B_point.b));
    double angle = acos(cos_angle);
    sum += angle;
  }
  return areSame(sum, 2 * M_PI);
}
bool Polygon::isSimilarTo(const Shape &another) const {
  Polygon shape;
  try {
    shape = dynamic_cast<const Polygon &>(another);
  } catch (...) {
    return false;
  }
  if (shape.verticesCount() != verticesCount()) {
    return false;
  }
  size_t n = verticesCount();
  std::vector<double> edges;
  std::vector<double> another_edges;
  for (size_t i = 0; i < n; ++i) {
    edges.push_back(findDistance(vertices[i], vertices[(i + 1) % n]));
    another_edges.push_back(findDistance(shape.vertices[i], shape.vertices[(i + 1) % n]));
  }
  for (size_t i = 0; i < n; i++) {
    double k = edges[i] / another_edges[0];
    bool flag = true;
    for (size_t j = 0; j < n; ++j) {
      if (!areSame(edges[(j + i) % n] / another_edges[j], k)) {
        flag = false;
        break;
      }
    }
    if (flag) {
      return true;
    }
  }
  for (size_t i = 0; i < n; i++) {
    double k = edges[i] / another_edges[0];
    bool flag = true;
    for (size_t j = 0; j < n; ++j) {
      if (!areSame(edges[(n + i - j) % n] / another_edges[j], k)) {
        flag = false;
        break;
      }
    }
    if (flag) {
      return true;
    }
  }
  return false;
}
bool Polygon::isCongruentTo(const Shape &another) const {
  Polygon shape;
  try {
    shape = dynamic_cast<const Polygon &>(another);
  } catch (...) {
    return false;
  }
  if (shape.verticesCount() != verticesCount()) {
    return false;
  }
  return (isSimilarTo(another) && area() == another.area());
}
Polygon::Polygon(std::initializer_list<Point>
                 points) {
  for (auto point : points) {
    vertices.push_back(point);
  }
}

///////////////////////////////////////// TRIANGLE /////////////////////////////////////////

Circle Triangle::circumscribedCircle() const {
  std::vector<Point> vertices = getVertices();
  Point A = vertices[0];
  Point B = vertices[1];
  Point C = vertices[2];

  double D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
  double center_x = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y)
      + (C.x * C.x + C.y * C.y) * (A.y - B.y)) / D;
  double center_y = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x)
      + (C.x * C.x + C.y * C.y) * (B.x - A.x)) / D;

  Point center(center_x, center_y);
  return Circle(center, findDistance(center, A));
}
Circle Triangle::inscribedCircle() const {
  std::vector<Point> vertices = getVertices();
  Point A = vertices[0];
  Point B = vertices[1];
  Point C = vertices[2];

  double a = findDistance(B, C);
  double b = findDistance(A, C);
  double c = findDistance(A, B);

  double center_x = (a * A.x + b * B.x + c * C.x) / (a + b + c);
  double center_y = (a * A.y + b * B.y + c * C.y) / (a + b + c);
  Point center(center_x, center_y);
  double rad = 2 * area() / (a + b + c);
  return Circle(center, rad);
}
Point Triangle::centroid() const {
  std::vector<Point> vertices = getVertices();
  Point A = vertices[0];
  Point B = vertices[1];
  Point C = vertices[2];
  return Point((A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3);
}
double det3(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
  double result = 0;
  result += a * e * i + b * f * g + c * d * h;
  result -= g * e * c + h * f * a + d * b * i;
  return result;
}
Point Triangle::orthocenter() const {
  std::vector<Point> vertices = getVertices();
  Point A = vertices[0];
  Point B = vertices[1];
  Point C = vertices[2];
  Line AC(A, C);
  Line BH(B, -1 / AC.getSlope());
  Line AB(A, B);
  Line CH(C, -1 / AB.getSlope());
  Point H = Point(findIntersectionPoint(BH, CH));
  return H;
}
Line Triangle::EulerLine() const {
  return Line(orthocenter(), circumscribedCircle().center());
}
Circle Triangle::ninePointsCircle() const {
  Point H = orthocenter();
  Point O = circumscribedCircle().center();
  Point center((H.x + O.x) / 2, (H.y + O.y) / 2);
  return Circle(center, circumscribedCircle().radius() / 2);
}

///////////////////////////////////////// RECTANGLE /////////////////////////////////////////

Point Rectangle::center() const {
  return findIntersectionPoint(diagonals().first, diagonals().second);
}
std::pair<Line, Line> Rectangle::diagonals() const {
  std::vector<Point> vertices = getVertices();
  return std::pair<Line, Line>(Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3]));
}
Rectangle::Rectangle(Point
                     x, Point
                     y, double
                     coef) {
  Point cent((x.x + y.x) / 2, (x.y + y.y) / 2);
  coef = std::max(coef, 1.0 / coef);
  double diag_len = findDistance(x, y);
  double width = diag_len / sqrt(1.0 + coef * coef);

  Line diag(x, y);
  Point p(diag.a, diag.b);
  Point diag_vec(y.x - x.x, y.y - x.y);
  Point height_vec(diag_vec.x / (1.0 + coef * coef), diag_vec.y / (1.0 + coef * coef));
  Point h(x.x + height_vec.x, x.y + height_vec.y);

  p = Point(p.x / sqrt(p.x * p.x + p.y * p.y), p.y / sqrt(p.x * p.x + p.y * p.y));
  double case_true_x = h.x + p.x * (diag_len * coef / (coef * coef + 1.0));
  double case_true_y = h.y + p.y * (diag_len * coef / (coef * coef + 1.0));

  double case_false_x = h.x - p.x * (diag_len * coef / (coef * coef + 1.0));
  double case_false_y = h.y - p.y * (diag_len * coef / (coef * coef + 1.0));

  Point case_true(case_true_x, case_true_y);
  Point case_false(case_false_x, case_false_y);

  Point a;
  if (areSame(findDistance(case_true, x), width) && !areSame(p.x * diag_vec.y - p.y * diag_vec.x, 0.))
    a = case_true;
  Point b(a.x + (cent.x - a.x) * 2, a.y + (cent.y - a.y) * 2);
  *this = Rectangle(x, a, y, b);
}

///////////////////////////////////////// SQUARE /////////////////////////////////////////

Circle Square::circumscribedCircle() const {
  return Circle(center(), findDistance(center(), getVertices()[0]));
}
Circle Square::inscribedCircle() const {
  return Circle(center(), findDistance(getVertices()[1], getVertices()[0]) / 2);
}