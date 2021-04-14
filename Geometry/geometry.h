#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>

const double epsilon = 1e-5;

double toRadians(double angle) {
  return (angle / 180) * M_PI;
}

///////////////////////////////////////// POINT /////////////////////////////////////////

struct Point {
  double x;
  double y;
  Point() : x(0), y(0) {};
  Point(const double &x, const double &y) : x(x), y(y) {};
  bool operator==(const Point &other) const;
  bool operator!=(const Point &other) const;
};

///////////////////////////////////////// LINE /////////////////////////////////////////

class Line {
 public:
  double a, b, c;
  Line() : a(0), b(0), c(0) {};
  Line(double a, double b, double c) : a(a), b(b), c(c) {};
  Line(const Point &first, const Point &second);
  Line(double slope, double shift);
  Line(const Point &point, double slope);
  bool operator==(const Line &other) const;
  bool operator!=(const Line &other) const;

  double getSlope() const;
};

///////////////////////////////////////// SHAPE /////////////////////////////////////////

class Shape {
 public:
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool operator==(const Shape &another) const = 0;
  virtual bool operator!=(const Shape &another) const = 0;
  virtual bool isCongruentTo(const Shape &another) const = 0;
  virtual bool isSimilarTo(const Shape &another) const = 0;
  virtual bool containsPoint(Point point) const = 0;

  virtual void rotate(Point center, double angle) = 0;
  virtual void reflex(Point center) = 0;
  virtual void reflex(Line axis) = 0;
  virtual void scale(Point center, double coefficient) = 0;

  virtual ~Shape() {};
};

///////////////////////////////////////// ELLIPSE /////////////////////////////////////////

class Ellipse : public Shape {
 private:
  std::pair<Point, Point> focus;
  double dist;
  double a() const;
  double b() const;
  double c() const;

 public:
  Ellipse();
  Ellipse(Point focus1, Point focus2, double _dist);
  std::pair<Point, Point> focuses() const;
  double eccentricity() const;
  std::pair<Line, Line> directrices() const;
  Point center() const;

  double focal_distance() const;

  double perimeter() const override;
  double area() const override;
  bool operator==(const Shape &another) const override;
  bool operator!=(const Shape &another) const override;
  bool isCongruentTo(const Shape &another) const override;
  bool isSimilarTo(const Shape &another) const override;
  bool containsPoint(Point point) const override;

  void rotate(Point center, double angle) override;
  void reflex(Point center) override;
  void reflex(Line axis) override;
  void scale(Point center, double coefficient) override;

  ~Ellipse() override {};
};

///////////////////////////////////////// CIRCLE /////////////////////////////////////////

class Circle : public Ellipse {
 public:
  Circle(const Point &center, double radius) : Ellipse(center, center, 2 * radius) {};
  double radius() const;

  double perimeter() const override;
  double area() const override;

  void rotate(Point _center, double angle) override;
  void reflex(Point _center) override;
  void reflex(Line axis) override;
  void scale(Point _center, double coefficient) override;
};

///////////////////////////////////////// POLYGON /////////////////////////////////////////

class Polygon : public Shape {
 private:
  std::vector<Point> vertices;
 public:
  Polygon() = default;
  Polygon(const std::vector<Point> &points);
  Polygon(std::initializer_list<Point> points);
  size_t verticesCount() const;
  std::vector<Point> getVertices() const;
  bool isConvex() const;

  double perimeter() const override;
  double area() const override;
  bool operator==(const Shape &another) const override;
  bool operator!=(const Shape &another) const override;
  bool isCongruentTo(const Shape &another) const override;
  bool isSimilarTo(const Shape &another) const override;
  bool containsPoint(Point point) const override;

  void rotate(Point center, double angle) override;
  void reflex(Point center) override;
  void reflex(Line axis) override;
  void scale(Point center, double coefficient) override;
};

///////////////////////////////////////// TRIANGLE /////////////////////////////////////////

class Triangle : public Polygon {
 public:
  Triangle(Point A, Point B, Point C) : Polygon({A, B, C}) {};

  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;

  Point centroid() const;
  Point orthocenter() const;

  Line EulerLine() const;
  Circle ninePointsCircle() const;
};

///////////////////////////////////////// RECTANGLE /////////////////////////////////////////

class Rectangle : public Polygon {
 public:
  Rectangle(Point A, Point C, double coef);
  Rectangle(Point a, Point b, Point c, Point d) : Polygon({a, b, c, d}) {};
  Point center() const;
  std::pair<Line, Line> diagonals() const;

};

///////////////////////////////////////// SQUARE /////////////////////////////////////////

class Square : public Rectangle {
 public:
  Square(Point A, Point C) : Rectangle(A, C, 1) {};
  Circle circumscribedCircle() const;
  Circle inscribedCircle() const;
};