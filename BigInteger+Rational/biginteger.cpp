#include <iostream>
#include <vector>
#include <string>
#include "biginteger.h"

// BIG INTEGER

BigInteger::BigInteger() {
  *this = 0;
}

BigInteger operator ""_bi(const char *num) {
  return BigInteger(std::string(num));
}

BigInteger::BigInteger(const long long &number) {
  long long tmp_number = number;
  digits = std::vector<unsigned int>();
  if (number < 0) {
    sign = true;
    tmp_number *= -1;
  }
  while (tmp_number > 0) {
    digits.push_back(tmp_number % BASE);
    tmp_number /= BASE;
  }
}

BigInteger::BigInteger(const std::string &str_number) {
  std::string number;
  if (str_number[0] == '-' && str_number != "-0") {
    sign = true;
    number = str_number.substr(1);
  } else {
    number = str_number;
  }
  unsigned int number_end = number.length();
  while (number_end > BASE_LENGTH) {
    unsigned int number_begin = number_end - BASE_LENGTH;
    std::string digit = number.substr(number_begin, BASE_LENGTH);
    digits.push_back(stoi(digit));
    number_end -= BASE_LENGTH;
  }
  std::string digit = number.substr(0, number_end);
  digits.push_back(stoi(digit));
}

std::string BigInteger::toString() const {
  if (digits.empty()) {
    return "0";
  }
  std::string result = "";
  if (sign) {
    result += "-";
  }
  for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
    if (it == digits.rbegin()) {
      result += std::to_string(*it);
    } else {
      for (unsigned int i = 0; i < BASE_LENGTH - std::to_string(*it).length(); ++i) {
        result += "0";
      }
      result += std::to_string(*it);
    }
  }
  return result;
}
BigInteger &BigInteger::operator=(const int &num) {
  *this = BigInteger(num);
  return *this;
}
BigInteger operator+(const BigInteger &left, const BigInteger &right) {
  BigInteger result = left;
  result += right;
  return result;
}
BigInteger operator-(const BigInteger &left, const BigInteger &right) {
  BigInteger result = left;
  result -= right;
  return result;
}
BigInteger operator%(const BigInteger &left, const BigInteger &right) {
  BigInteger result = left;
  result %= right;
  return result;
}

BigInteger &BigInteger::operator+=(const BigInteger &other) {
  if (other == 0) {
    return *this;
  }
  if (sign && !other.sign) {
    sign = !sign;
    *this -= other;
    if (*this != 0) {
      sign = !sign;
    }
    return *this;
  } else {
    if (!sign && other.sign) {
      *this -= -other;
      return *this;
    }
  }
  BigInteger result;
  unsigned int addition = 0;
  unsigned int max_size = digits.size() > other.digits.size() ? digits.size() : other.digits.size();
  unsigned int first;
  unsigned int second;
  unsigned int cur_digit;

  for (unsigned int i = 0; i < max_size; ++i) {
    first = i < digits.size() ? digits[i] : 0;
    second = i < other.digits.size() ? other.digits[i] : 0;
    cur_digit = (first + second + addition) % BASE;
    if (i < digits.size()) {
      digits[i] = cur_digit;
    } else {
      digits.push_back(cur_digit);
    }
    addition = (first + second + addition) / BASE;
  }
  if (addition) {
    digits.push_back(addition);
  }
  return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &other) {
  if (other == 0) {
    return *this;
  }
  if (sign && !other.sign) {
    sign = !sign;
    *this += other;
    sign = !sign;
    return *this;
  }
  if (!sign && other.sign) {
    *this += -other;
    return *this;
  }
  if (sign && other.sign && *this > other) {
    BigInteger result = -other;
    result -= -*this;
    *this = result;
    return *this;
  }
  if (!sign && !other.sign && *this < other) {
    BigInteger result = other;
    result -= *this;
    *this = result;
    sign = !sign;
    return *this;
  }
  unsigned int loan = 0;
  unsigned int max_size = digits.size();
  unsigned int first;
  unsigned int second;
  unsigned int cur_digit;
  for (unsigned int i = 0; i < max_size; ++i) {
    first = i < digits.size() ? digits[i] : 0;
    second = i < other.digits.size() ? other.digits[i] : 0;
    cur_digit = (2 * BASE + first - second - loan) % BASE;
    if (i < digits.size()) {
      digits[i] = cur_digit;
    } else {
      digits.push_back(cur_digit);
    }
    loan = 2 - (2 * BASE + first - second - loan) / BASE;
  }
  while (!digits.empty() && digits[digits.size() - 1] == 0) {
    digits.pop_back();
  }
  if (digits.empty()) {
    sign = false;
  }
  return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &other) {
  BigInteger divider = other;
  bool tmp_sign = sign;
  sign = false;
  divider.sign = false;
  *this -= (*this / divider) * divider;
  sign = tmp_sign;
  return *this;
}

BigInteger &BigInteger::operator++() {
  if (*this == 0) {
    *this = 1;
    return *this;
  }
  if (*this == -1) {
    *this = 0;
    return *this;
  }
  if (sign) {
    sign = !sign;
    --*this;
    sign = !sign;
    return *this;
  }
  bool need_addition = true;
  unsigned int cur_digit = 0;
  while (need_addition) {
    if (digits.size() == cur_digit) {
      digits.push_back(0);
    }
    ++digits[cur_digit];
    if (digits[cur_digit] < BASE) {
      need_addition = false;
    } else {
      digits[cur_digit] = 0;
    }
    ++cur_digit;
  }
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger num = *this;
  ++*this;
  return num;
}
BigInteger &BigInteger::operator--() {
  if (*this == 0) {
    *this = -1;
    return *this;
  }
  if (*this == 1) {
    *this = 0;
    return *this;
  }
  if (sign) {
    sign = !sign;
    ++*this;
    sign = !sign;
    return *this;
  }
  unsigned int cur_digit = 0;
  bool need_loan = true;
  while (need_loan) {
    if (digits[cur_digit] > 0) {
      --digits[cur_digit];
      need_loan = false;
    } else {
      digits[cur_digit] = BASE - 1;
    }
    ++cur_digit;
  }
  if (digits[digits.size() - 1] == 0) {
    digits.pop_back();
  }
  return *this;
}
BigInteger BigInteger::operator--(int) {
  BigInteger cp = *this;
  --*this;
  return *this;
}
BigInteger BigInteger::operator-() const {
  BigInteger result = *this;
  result.sign = !sign;
  return result;
}
bool operator<=(const BigInteger &left, const BigInteger &right) {
  if (left.digits.empty() && (!right.sign || right.digits.empty())) {
    return true;
  }
  unsigned int max_size = left.digits.size() > right.digits.size() ? left.digits.size() : right.digits.size();
  if (left.sign && !right.sign) {
    return true;
  }
  if (!left.sign && right.sign) {
    return false;
  }
  if (left.sign && right.sign) {
    return -left >= -right;
  }
  unsigned int first;
  unsigned int second;
  for (unsigned int i = max_size; i > 0; --i) {
    first = i - 1 < left.digits.size() ? left.digits[i - 1] : 0;
    second = i - 1 < right.digits.size() ? right.digits[i - 1] : 0;
    if (first == second) {
      continue;
    } else {
      return first <= second;
    }
  }
  return true;
}
bool operator>=(const BigInteger &left, const BigInteger &right) {
  return right <= left;
}
bool operator==(const BigInteger &left, const BigInteger &right) {
  return left <= right && right <= left;
}
bool operator!=(const BigInteger &left, const BigInteger &right) {
  return !(left == right);
}
bool operator<(const BigInteger &left, const BigInteger &right) {
  return left <= right && left != right;
}
bool operator>(const BigInteger &left, const BigInteger &right) {
  return right < left;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &num) {
  out << num.toString();
  return out;
}
std::istream &operator>>(std::istream &in, BigInteger &num) {
  std::string str;
  in >> str;
  num = BigInteger(str);
  return in;
}
BigInteger::operator bool() const {
  return !digits.empty();
}
BigInteger &BigInteger::operator*=(const BigInteger &other) {
  BigInteger result;
  bool res_sign = sign != other.sign;
  for (unsigned int i = 0; i < digits.size() + other.digits.size() + 1; ++i) {
    result.digits.push_back(0);

  }
  for (unsigned int i = 0; i < digits.size(); ++i) {
    for (unsigned int j = 0; j < other.digits.size(); ++j) {
      result.digits[i + j] += digits[i] * other.digits[j];
    }
  }
  for (unsigned int i = 0; i < result.digits.size(); ++i) {
    result.digits[i + 1] += (result.digits[i] / BASE);
    result.digits[i] %= BASE;
  }
  while (!result.digits.empty() && result.digits[result.digits.size() - 1] == 0) {
    result.digits.pop_back();
  }
  *this = result;
  sign = !(*this == 0) && res_sign;
  return *this;
}
BigInteger operator*(const BigInteger &left, const BigInteger &right) {
  BigInteger result = left;
  result *= right;
  return result;
}

BigInteger &BigInteger::operator/=(const BigInteger &other) {
  if (other == 0) {
    std::cerr << "Division by zero" << std::endl;
  }
  if (other == 1) {
    return *this;
  }
  BigInteger result;
  bool result_sign = sign != other.sign;
  BigInteger dividend = *this;
  dividend.sign = false;
  BigInteger divisor = other;
  divisor.sign = false;
  BigInteger BASE_SHIFT;
  while (dividend >= divisor) {
    BASE_SHIFT = 1;
    while (dividend >= divisor * BASE) {
      divisor *= BASE;
      BASE_SHIFT *= BASE;
    }
    while (dividend >= divisor) {
      dividend -= divisor;
      result += BASE_SHIFT;
    }
    divisor = other;
    divisor.sign = false;
  }
  while (!result.digits.empty() && result.digits[result.digits.size() - 1] == 0) {
    result.digits.pop_back();
  }
  result.sign = result_sign;
  if (result.digits.empty()) {
    result.sign = false;
  }
  *this = result;
  return *this;
}
BigInteger operator/(const BigInteger &left, const BigInteger &right) {
  BigInteger result = left;
  result /= right;
  return result;
}

// RATIONAL

Rational::Rational() {
  numerator = 0;
  denominator = 1;
  sign = false;
}

Rational::Rational(int num) {
  numerator = num >= 0 ? num : -num;
  denominator = 1;
  sign = num < 0;
}
Rational::Rational(const BigInteger &num) {
  numerator = num >= 0 ? num : -num;
  denominator = 1;
  sign = num < 0;
}
std::string Rational::toString() const {
  std::string result = "";
  if (sign) {
    result += "-";
  }
  result += numerator.toString();
  if (denominator != 1) {
    result += "/";
    result += denominator.toString();
  }
  return result;
}
BigInteger Rational::gcd(const BigInteger &first, const BigInteger &second) const {
  if (first == 0 || second == 0) {
    return first + second;
  }
  BigInteger left = first;
  BigInteger right = second;

  left = left >= 0 ? left : -left;
  right = right >= 0 ? right : -right;

  while (left != 0 && right != 0) {
    if (left > right) {
      left %= right;
    } else {
      right %= left;
    }
  }
  return left + right;
}
Rational &Rational::operator+=(const Rational &other) {
  if (sign && !other.sign) {
    sign = !sign;
    *this -= other;
    sign = !(numerator == 0) && !sign;
    return *this;
  } else {
    if (!sign && other.sign) {
      *this -= -other;
      return *this;
    }
  }
  numerator = numerator * other.denominator + denominator * other.numerator;
  denominator *= other.denominator;
  BigInteger gcd_ = gcd(numerator, denominator);
  numerator /= gcd_;
  denominator /= gcd_;
  if (numerator == 0 || -numerator == 0) {
    numerator = 0;
    sign = false;
  }
  return *this;
}
Rational operator+(const Rational &left, const Rational &right) {
  Rational result = left;
  result += right;
  return result;
}
Rational &Rational::operator-=(const Rational &other) {
  if (other == 0) {
    return *this;
  }
  if (sign && !other.sign) {
    sign = !sign;
    *this += other;
    sign = !(numerator == 0) && !sign;
    return *this;
  }
  if (!sign && other.sign) {
    *this += -other;
    return *this;
  }
  if (sign && other.sign && *this > other) {
    Rational result = -other;
    result -= -*this;
    *this = result;
    return *this;
  }
  if (!sign && !other.sign && *this < other) {
    Rational result = other;
    result -= *this;
    *this = result;
    sign = !(numerator == 0) && !sign;
    return *this;
  }
  numerator = numerator * other.denominator - denominator * other.numerator;
  denominator *= other.denominator;
  BigInteger gcd_ = gcd(numerator, denominator);
  numerator /= gcd_;
  denominator /= gcd_;
  if (numerator == 0 || -numerator == 0) {
    numerator = 0;
    sign = false;
  }
  return *this;
}
Rational operator-(const Rational &left, const Rational &right) {
  Rational result = left;
  result -= right;
  return result;
}

Rational &Rational::operator*=(const Rational &other) {
  if (numerator == 0) {
    *this = 0;
    return *this;
  }
  sign = sign != other.sign;
  numerator *= other.numerator;
  denominator *= other.denominator;
  BigInteger gcd_ = gcd(numerator, denominator);
  numerator /= gcd_;
  denominator /= gcd_;
  return *this;
}
Rational operator*(const Rational &left, const Rational &right) {
  Rational result = left;
  result *= right;
  return result;
}
Rational &Rational::operator/=(const Rational &other) {
  if (numerator == 0) {
    *this = 0;
    return *this;
  }
  sign = sign != other.sign;
  numerator *= other.denominator;
  denominator *= other.numerator;
  BigInteger gcd_ = gcd(numerator, denominator);
  numerator /= gcd_;
  denominator /= gcd_;
  return *this;
}
Rational operator/(const Rational &left, const Rational &right) {
  Rational result = left;
  result /= right;
  return result;
}
Rational Rational::operator-() const {
  Rational res = *this;
  if (*this != 0) {
    res.sign = !sign;
  }
  return res;
}
bool operator<=(const Rational &left, const Rational &right) {
  if (left.sign && !right.sign) {
    return true;
  }
  if (!left.sign && right.sign) {
    return false;
  }
  if (left.sign && right.sign) {
    return -left >= -right;
  }
  return left.numerator * right.denominator <= left.denominator * right.numerator;
}
bool operator>=(const Rational &left, const Rational &right) {
  return right <= left;
}
bool operator<(const Rational &left, const Rational &right) {
  return left <= right && left != right;
}
bool operator>(const Rational &left, const Rational &right) {
  return right < left;
}
bool operator==(const Rational &left, const Rational &right) {
  return left >= right && left <= right;
}
bool operator!=(const Rational &left, const Rational &right) {
  return !(left == right);
}
std::string Rational::asDecimal(unsigned int precision = 0) const {
  BigInteger num = numerator;
  BigInteger den = denominator;
  std::string int_part = (sign ? "-" : "") + (num / den).toString();
  std::string frac_part = "";
  while (frac_part.length() < precision) {
    num %= den;
    num *= BI_BASE;
    std::string addition = (num / den).toString();
    while (addition.length() < BI_BASE_LENGTH) {
      addition = "0" + addition;
    }
    frac_part += addition;
  }
  while (frac_part.length() > precision) {
    frac_part.pop_back();
  }
  if (precision == 0) {
    return int_part;
  }
  return int_part + "." + frac_part;
}
Rational::operator double() const {
  return stod(asDecimal(DOUBLE_DECIMAL_DIGITS_NUM));
}
