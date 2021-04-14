#include <iostream>
#include <vector>
#include <string>

class BigInteger {
 private:
  static const int BASE = 10000;
  static const int BASE_LENGTH = 4;
  std::vector<unsigned int> digits;
  bool sign = false;

 public:
  // Constructors
  BigInteger();
  BigInteger(const long long &number);
  BigInteger(const std::string &str_number);
  std::string toString() const;
  BigInteger &operator=(const int &num);

  // Arithmetics
  friend BigInteger operator+(const BigInteger &left, const BigInteger &right);
  friend BigInteger operator-(const BigInteger &left, const BigInteger &right);
  friend BigInteger operator*(const BigInteger &left, const BigInteger &right);
  friend BigInteger operator/(const BigInteger &left, const BigInteger &right);
  friend BigInteger operator%(const BigInteger &left, const BigInteger &right);

  // Arithmetics + assignment
  BigInteger &operator+=(const BigInteger &other);
  BigInteger &operator-=(const BigInteger &other);
  BigInteger &operator*=(const BigInteger &other);
  BigInteger &operator/=(const BigInteger &other);
  BigInteger &operator%=(const BigInteger &other);

  // Unary operators
  BigInteger &operator++();
  BigInteger operator++(int);
  BigInteger &operator--();
  BigInteger operator--(int);
  BigInteger operator-() const;

  // Comparison
  friend bool operator<(const BigInteger &left, const BigInteger &right);
  friend bool operator>(const BigInteger &left, const BigInteger &right);
  friend bool operator<=(const BigInteger &left, const BigInteger &right);
  friend bool operator>=(const BigInteger &left, const BigInteger &right);
  friend bool operator==(const BigInteger &left, const BigInteger &right);
  friend bool operator!=(const BigInteger &left, const BigInteger &right);

  explicit operator bool() const;
  friend std::ostream &operator<<(std::ostream &out, const BigInteger &num);
  friend std::istream &operator>>(std::istream &in, BigInteger &num);
};

class Rational {
 private:
  bool sign;
  BigInteger numerator;
  BigInteger denominator;

  BigInteger gcd(const BigInteger &first, const BigInteger &second) const;
  static const unsigned int DOUBLE_DECIMAL_DIGITS_NUM = 16;
  static const unsigned int BI_BASE = 10000;
  static const unsigned int BI_BASE_LENGTH = 4;

 public:
  // Constructors
  Rational();
  Rational(int num);
  Rational(const BigInteger &num);
  std::string toString() const;
  std::string asDecimal(unsigned int precision) const;
  explicit operator double() const;

  // Arithmetics
  friend Rational operator+(const Rational &left, const Rational &right);
  friend Rational operator-(const Rational &left, const Rational &right);
  friend Rational operator*(const Rational &left, const Rational &right);
  friend Rational operator/(const Rational &left, const Rational &right);

  // Arithmetics + assigment
  Rational &operator+=(const Rational &other);
  Rational &operator-=(const Rational &other);
  Rational &operator*=(const Rational &other);
  Rational &operator/=(const Rational &other);

  // Unary operator
  Rational operator-() const;

  // Comparison
  friend bool operator<(const Rational &left, const Rational &right);
  friend bool operator>(const Rational &left, const Rational &right);
  friend bool operator<=(const Rational &left, const Rational &right);
  friend bool operator>=(const Rational &left, const Rational &right);
  friend bool operator==(const Rational &left, const Rational &right);
  friend bool operator!=(const Rational &left, const Rational &right);
};

