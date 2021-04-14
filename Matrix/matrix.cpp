#include "matrix.h"
#include <vector>
#include <string>
#include <iostream>

///////////////////////////////////FINITE///////////////////////////////////////////////////

template<int N>
Finite<N>::Finite(int number_) {
  number = (N + number_ % N) % N;
}
template<int N>
Finite<N> Finite<N>::operator+=(const Finite<N> &other) {
  number = (number + other.number) % N;
  return *this;
}
template<int N>
Finite<N> Finite<N>::operator-=(const Finite<N> &other) {
  number = (N + number - other.number) % N;
  return *this;
}
template<int N>
Finite<N> Finite<N>::operator*=(const Finite<N> &other) {
  number = (number % N * other.number % N) % N;
  return *this;
}
template<int N>
Finite<N> Finite<N>::operator/=(const Finite<N> &other) {
  for (size_t i = 0; i < N - 2; i++) {
    number = (number * other.number) % N;
  }
  return *this;
}
template<int N>
Finite<N> &Finite<N>::operator++() {
  *this += 1;
  return *this;
}
template<int N>
Finite<N> Finite<N>::operator++(int) {
  Finite<N> result = *this;
  ++*this;
  return result;
}
template<int N>
Finite<N> &Finite<N>::operator--() {
  *this -= 1;
  return *this;
}
template<int N>
Finite<N> Finite<N>::operator--(int) {
  Finite<N> result = *this;
  --*this;
  return result;
}
template<int N>
Finite<N>::operator int() const {
  return number;
}
template<int N>
Finite<N>::Finite() {
  number = 0;
}


///////////////////////////////////MATRIX///////////////////////////////////////////////////

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix() {
  array = std::vector<std::vector<Field>>(height, std::vector<Field>(width));
  for (size_t i = 0; i < height; ++i) {
    for (size_t j = 0; j < width; ++j) {
      if (i == j) {
        array[i][j] = 1;
      } else {
        array[i][j] = 0;
      }
    }
  }
}
template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> &first, const Matrix<M, N, Field> &second);
template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> &first, const Matrix<M, N, Field> &second);
template<size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field> &first, const Matrix<N, K, Field> &second);
template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(Field number, const Matrix<M, N, Field> &matrix);
template<size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field> &matrix, Field number);

template<size_t M, size_t N, typename Field>
template<typename T>
Matrix<M, N, Field>::Matrix(const std::vector<std::vector<T>> &vector) {
  array = std::vector<std::vector<Field>>(height, std::vector<Field>(width));
  for (size_t i = 0; i < height; ++i) {
    for (size_t j = 0; j < width; ++j) {
      Field element(vector[i][j]);
      array[i][j] = element;
    }
  }
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field> &other) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      array[i][j] += other.array[i][j];
    }
  }
  return *this;
}
template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field> &other) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      array[i][j] -= other.array[i][j];
    }
  }
  return *this;
}
template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> &first, const Matrix<M, N, Field> &second) {
  Matrix<M, N, Field> result = first;
  result += second;
  return result;
}
template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> &first, const Matrix<M, N, Field> &second) {
  Matrix<M, N, Field> result = first;
  result -= second;
  return result;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator*=(Field number) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      array[i][j] *= number;
    }
  }
  return *this;
}
template<size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
  Matrix<N, M, Field> result;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      result.array[j][i] = array[i][j];
    }
  }
  return result;
}
template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator*=(const Matrix<N, N, Field> &other) {
  *this = *this * other;
  return *this;
}
template<size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const {
  static_assert(M == N, "Not square matrix");
  std::vector<std::vector<Field>> copy_array = array;
  Field det = 1;
  for (size_t i = 0; i < M; ++i) {
    size_t k = i;
    for (size_t j = i + 1; j < M; ++j)
      if (copy_array[j][i] > copy_array[k][i]) {
        k = j;
      }
    if (copy_array[k][i] == 0) {
      det = 0;
      break;
    }
    std::swap(copy_array[i], copy_array[k]);
    if (i != k) {
      det *= -1;
    }
    det *= copy_array[i][i];
    for (size_t j = i + 1; j < M; ++j) {
      copy_array[i][j] /= copy_array[i][i];
    }
    for (size_t j = 0; j < M; ++j) {
      if (j != i && copy_array[j][i] > 0) {
        for (size_t t = i + 1; t < M; ++t) {
          copy_array[j][t] -= copy_array[i][t] * copy_array[j][i];
        }
      }
    }
  }
  return det;
}
template<size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::invert() {
  static_assert(M == N, "Not square matrix");
  Matrix<M, 2 * M, Field> Gauss_Jordan_matrix;
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < M; ++j) {
      Gauss_Jordan_matrix[i][j] = array[i][j];
    }
  }
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < M; ++j) {
      if (i == j) {
        Gauss_Jordan_matrix[i][j + M] = 1;
      } else {
        Gauss_Jordan_matrix[i][j + M] = 0;
      }
    }
  }
  for (size_t i = M - 1; i > 0; i--) {
    if (Gauss_Jordan_matrix[i - 1][0] < Gauss_Jordan_matrix[i][0]) {
      std::swap(Gauss_Jordan_matrix[i - 1], Gauss_Jordan_matrix[i]);
    }
  }
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < M; j++) {
      if (i != j) {
        Field d = Gauss_Jordan_matrix[j][i] / Gauss_Jordan_matrix[i][i];
        for (size_t k = 0; k < M * 2; k++) {
          Gauss_Jordan_matrix[j][k] -= d * Gauss_Jordan_matrix[i][k];
        }
      }
    }
  }
  for (size_t i = 0; i < M; i++) {
    for (size_t j = M; j < 2 * M; ++j) {
      Gauss_Jordan_matrix[i][j] /= Gauss_Jordan_matrix[i][i];
    }
  }
  for (size_t i = 0; i < M; i++) {
    for (size_t j = 0; j < M; ++j) {
      array[i][j] = Gauss_Jordan_matrix[i][j + M];
    }
  }
}

template<size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::inverted() {
  Matrix<N, M, Field> copy_matrix = *this;
  copy_matrix.invert();
  return copy_matrix;
}
template<size_t M, size_t N, typename Field>
std::vector<Field> &Matrix<M, N, Field>::operator[](size_t i) {
  return array[i];
}
template<size_t M, size_t N, typename Field>
const std::vector<Field> &Matrix<M, N, Field>::operator[](size_t i) const {
  return array[i];
}
template<size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getRow(size_t i) {
  return array[i];
}
template<size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getColumn(size_t j) {
  std::vector<Field> column;
  for (size_t i = 0; i < M; ++i) {
    column.push_back(array[i][j]);
  }
  return column;
}
template<size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const {
  static_assert(M == N, "Not square matrix");
  Field trace = 0;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < M; ++j) {
      if (i == j) {
        trace += array[i][j];
      }
    }
  }
  return trace;
}
template<size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
  Matrix<M, N, Field> copy_matrix = *this;
  size_t rank = 0;
  std::vector<bool> rowSelected(M, false);
  for (size_t j = 0; j < N; ++j) {
    size_t i;
    for (i = 0; i < M; ++i) {
      if (!rowSelected[i] && copy_matrix[i][j] != 0) {
        break;
      }
    }
    if (i != M) {
      ++rank;
      rowSelected[i] = true;
      for (size_t k = j + 1; k < N; ++k) {
        copy_matrix[i][k] /= copy_matrix[i][j];
      }
      for (size_t k = 0; k < M; ++k) {
        if (k != i && copy_matrix[k][i] != 0) {
          for (size_t p = j + 1; p < N; ++p) {
            copy_matrix[k][p] -= copy_matrix[i][p] * copy_matrix[k][j];
          }
        }
      }
    }
  }
  return rank;
}

template<size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field> &first, const Matrix<N, K, Field> &second) {
  Matrix<M, K, Field> result;
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      result.array[i][j] = 0;
    }
  }
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < N; ++k) {
        result.array[i][j] += first.array[i][k] * second.array[k][j];
      }
    }
  }
  return result;
}
template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(Field number, const Matrix<M, N, Field> &matrix) {
  Matrix<M, N, Field> result = matrix;
  result *= number;
  return result;
}
template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field> &matrix, Field number) {
  return number * matrix;
}

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

///////////////////////////////////BIGINTEGER///////////////////////////////////////////////

BigInteger::BigInteger() {
  *this = 0;
}

BigInteger operator ""_bi(const char *num) {
  return BigInteger(std::string(num));
}

BigInteger::BigInteger(const long long &number) {
  long long tmp_number = number;
  digits = std::vector<long long>();
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
  long long addition = 0;
  long long max_size = digits.size() > other.digits.size() ? digits.size() : other.digits.size();
  long long first;
  long long second;
  long long cur_digit;

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
  long long loan = 0;
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

void BigInteger::changeTo() {
  long long addition;
  for (size_t i = 0; i < digits.size(); ++i) {
    if (digits[i] < 0) {
      if (i == digits.size() - 1) {
        sign = !sign;
        for (long long &digit : digits) {
          digit *= -1;
        }
        changeTo();
        return;
      }
      digits[i] = digits[i] % BASE;
      addition = (-digits[i]) / BASE;
      if (digits[i] < 0) {
        digits[i] += BASE;
        ++addition;
      }
      digits[i + 1] -= addition;
    } else {
      addition = digits[i] / BASE;
      digits[i] = digits[i] % BASE;
      if (i == digits.size() - 1) {
        if (addition != 0) {
          digits.push_back(addition);
        }
      } else {
        digits[i + 1] += addition;
      }
    }
  }
}

bool BigInteger::isEven() const {
  return (digits[0] % 2 == 0);
}
BigInteger &BigInteger::double_bi() {
  int carry = 0;
  for (size_t i = 0; i < digits.size() || carry; ++i) {
    if (i == digits.size()) {
      digits.push_back(0);
    }
    long long current = carry + digits[i] * 2;
    digits[i] = current % BASE;
    carry = current / BASE;
  }
  while (!digits.empty() && digits.back() == 0) {
    digits.pop_back();
  }
  return *this;
}
BigInteger &BigInteger::split_in_half() {
  int carry = 0;
  for (int i = digits.size() - 1; i >= 0; --i) {
    long long current = digits[i] + carry * BASE;
    digits[i] = current / 2;
    carry = current % 2;
  }
  while (!digits.empty() && digits.back() == 0) {
    digits.pop_back();
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

BigInteger BigInteger::divide(size_t i) const {
  BigInteger res = 0;
  for (size_t j = i; j < digits.size(); ++j) {
    res.digits.push_back(digits[j]);
  }
  return res;
}

BigInteger BigInteger::multiplicate(size_t i) const {
  BigInteger res = 0;
  for (size_t j = 0; j < i; ++j) {
    res.digits.push_back(0);
  }
  for (long long digit : digits) {
    res.digits.push_back(digit);
  }
  return res;
}
BigInteger& BigInteger::operator/=(const BigInteger& other) {
  BigInteger numerator = *this > 0 ? *this : -*this;
  BigInteger denominator = other > 0 ? other : -other;
  if (numerator.digits.size() < denominator.digits.size()) {
    *this = 0;
    return *this;
  }
  std::vector<long long> nums(0);
  digits.clear();
  BigInteger div, sub;
  long long digit, digit_, pivot;
  for (size_t i = numerator.digits.size() - denominator.digits.size() + 1; i > 0; --i) {
    div = numerator.divide(i - 1);
    digit = 0;
    digit_ = BASE;
    while (digit_ - 1 > digit) {
      pivot = (digit + digit_) / 2;
      if (div >= pivot * denominator) {
        digit = pivot;
      } else {
        digit_ = pivot;
      }
    }
    sub = digit * denominator;
    numerator -= sub.multiplicate(i - 1);
    nums.push_back(digit);
  }
  size_t size = nums.size();
  for (size_t i = 1; i <= size; ++i) {
    digits.push_back(nums[size - i]);
  }
  if (nums.empty()) {
    digits.push_back(0);
  }
  while (digits.size() > 1 && digits.back() == 0) {
    digits.pop_back();
  }
  if (digits.size() == 1 && digits.back() == 0) {
    sign = false;
  }
  sign = sign != other.sign;
  return *this;
}
BigInteger operator/(const BigInteger &left, const BigInteger &right) {
  BigInteger result = left;
  result /= right;
  return result;
}

///////////////////////////////////////RATIONAL///////////////////////////////////////////

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
  if (first == 0) {
    return second;
  }
  if (second == 0) {
    return first;
  }
  if (first == second) {
    return first;
  }
  if (first == 1 || second == 1) {
    return 1;
  }
  BigInteger ans = 1;
  BigInteger first_ = first, second_ = second;
  while (first_.isEven() && second_.isEven()) {
    ans.double_bi();
    first_.split_in_half();
    second_.split_in_half();
  }
  while (first_ != 0) {
    while (first_.isEven()) {
      first_.split_in_half();
    }
    while (second_.isEven()) {
      second_.split_in_half();
    }
    if (first_ >= second_) {
      first_ -= second_;
    } else {
      second_ -= first_;
    }
  }
  ans *= second_;
  return ans;
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
std::istream &operator>>(std::istream &in, Rational &num) {
  BigInteger str;
  in >> str;
  num = Rational(str);
  return in;
}
std::ostream &operator<<(std::ostream &out, const Rational &num) {
  out << num.toString();
  return out;
}