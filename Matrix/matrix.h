#include <iostream>
#include <vector>
#include <string>

class BigInteger {
 private:
  static const int BASE = 1000000000;
  static const int BASE_LENGTH = 9;
  std::vector<long long> digits;
  bool sign = false;
  BigInteger divide(size_t i) const;
  BigInteger multiplicate(size_t i) const;

 public:

  friend int main();
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

  bool isEven() const;
  BigInteger &double_bi();
  BigInteger &split_in_half();
  void changeTo();

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

  static const unsigned int DOUBLE_DECIMAL_DIGITS_NUM = 16;
  static const unsigned int BI_BASE = 1000000000;
  static const unsigned int BI_BASE_LENGTH = 9;

 public:
  // Constructors
  BigInteger gcd(const BigInteger &first, const BigInteger &second) const;
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

template<int N>
class Finite {
 public:
  Finite<N>(int number_);
  Finite<N>();
  Finite<N> operator+=(const Finite<N> &other);
  Finite<N> operator-=(const Finite<N> &other);
  Finite<N> operator*=(const Finite<N> &other);
  Finite<N> operator/=(const Finite<N> &other);

  Finite<N> &operator++();
  Finite<N> operator++(int);
  Finite<N> &operator--();
  Finite<N> operator--(int);

  explicit operator int() const;

  friend bool operator==(const Finite<N> &left, const Finite<N> &right) {
    return left.number == right.number;
  }
  friend bool operator!=(const Finite<N> &left, const Finite<N> &right) {
    return !(left == right);
  }
  friend bool operator>(const Finite<N> &left, const Finite<N> &right) {
    return left.number > right.number;
  }
  friend bool operator<(const Finite<N> &left, const Finite<N> &right) {
    return left.number < right.number;
  }
  inline friend Finite<N> operator+(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> result(first.number);
    result += second;
    return result;
  }
  inline friend Finite<N> operator-(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> result(first.number);
    result -= second;
    return result;
  }
  inline friend Finite<N> operator*(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> result(first.number);
    result *= second;
    return result;
  }
  inline friend Finite<N> operator/(const Finite<N> &first, const Finite<N> &second) {
    Finite<N> result(first.number);
    result /= second;
    return result;
  }
 private:
  uint64_t number = 0;
};

template<size_t M, size_t N, typename Field = Rational>
class Matrix {
 public:
  Matrix<M, N, Field>();
  template<typename T>
  Matrix<M, N, Field>(const std::vector<std::vector<T>> &);

  Matrix<N, M, Field> transposed() const;
  Field det() const;
  Field trace() const;
  void invert();
  Matrix<N, M, Field> inverted();
  std::vector<Field> getRow(size_t i);
  std::vector<Field> getColumn(size_t j);
  size_t rank() const;

  Matrix<M, N, Field> operator+=(const Matrix<M, N, Field> &other);
  Matrix<M, N, Field> operator-=(const Matrix<M, N, Field> &other);
  Matrix<M, N, Field> operator*=(Field number);
  Matrix<M, N, Field> operator*=(const Matrix<N, N, Field> &other);
  std::vector<Field> &operator[](size_t i);
  const std::vector<Field> &operator[](size_t i) const;

  friend bool operator==(const Matrix<M, N, Field> &first, const Matrix<M, N, Field> &second) {
    if (first.array == second.array) {
      return true;
    }
    return false;
  }
  friend bool operator!=(const Matrix<M, N, Field> &first, const Matrix<M, N, Field> &second) {
    return !(first == second);
  }
 public:
  std::vector<std::vector<Field>> array;
  size_t height = M;
  size_t width = N;
};
