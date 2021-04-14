#include "string.h"
#include <iostream>

void String::shrink() {
  if (capacity == size + 1) {
    capacity *= 2;
    char *tmp_str = new char[capacity];
    strcpy(tmp_str, str);
    delete[] str;
    str = tmp_str;
  }
}

void String::cut_capacity() {
  if (8 * (size + 1) < 3 * capacity) {
    capacity /= 2;
    char *tmp_str = new char[capacity];
    strcpy(tmp_str, str);
    delete[] str;
    str = tmp_str;
  }
}

String::String() : size(0), capacity(1) {
  str = new char[1];
}

String::String(const char *c_str) : size(strlen(c_str)), capacity(strlen(c_str) + 1), str(new char[strlen(c_str) + 1]) {
  strcpy(str, c_str);
}

String::String(size_t char_num, char symbol) {
  size = char_num;
  capacity = size + 1;
  str = new char[capacity];
  memset(str, symbol, char_num);
  str[size] = '\0';
}

String::String(const String &other) {
  size = other.size;
  capacity = size + 1;
  str = new char[capacity];
  strcpy(str, other.str);
}

template<class T>
void String::swap(T &first, T &second) {
  if (&first != &second) {
    T tmp = first;
    first = second;
    second = tmp;
  }
}

void String::swap(String &s) {
  swap(size, s.size);
  swap(str, s.str);
  swap(capacity, s.capacity);
}

String &String::operator=(String other) {
  swap(other);
  return *this;
}

size_t String::length() const {
  return size;
}

void String::push_back(char symbol) {
  shrink();
  str[size++] = symbol;
  str[size] = '\0';
}

void String::pop_back() {
  cut_capacity();
  str[--size] = '\0';
}

char String::front() const {
  return str[0];
}

char String::back() const {
  return str[size - 1];
}

char &String::front() {
  return str[0];
}

char &String::back() {
  return str[size - 1];
}

String::~String() {
  delete[] str;
}

size_t String::find(const String &substring) const {
  char *ptr = strstr(str, substring.str);
  if (!ptr) {
    return size;
  }
  return ptr - str;
}

size_t String::rfind(const String &substring) const {
  for (size_t i = size - substring.size + 1; i > 0; --i) {
    if (strncmp(str + i - 1, substring.str, substring.size) == 0) {
      return i - 1;
    }
  }
  return size;
}
String String::substr(size_t start, size_t count) const {
  String res;
  for (size_t i = start; i < start + count; ++i) {
    res += str[i];
  }
  return res;
}
bool String::empty() const {
  return size == 0;
}

void String::clear() {
  delete[] str;
  size = 0;
  capacity = 1;
  str = new char[1];
  str[size] = '\0';
}

char &String::operator[](size_t pos) {
  return str[pos];
}

char String::operator[](size_t pos) const {
  return str[pos];
}

String &String::operator+=(char symbol) {
  push_back(symbol);
  return *this;
}

String &String::operator+=(const String &other) {
  size_t old_size = size;
  size += other.size;
  if (size + 1 > capacity) {
    shrink();
  }
  memcpy(str + old_size, other.str, other.size);
  str[size] = '\0';
  return *this;
}

String operator+(char symbol, const String &str) {
  String res(1, symbol);
  res += str;
  return res;
}

String operator+(const String &str, char symbol) {
  String res = str;
  res += symbol;
  return res;
}

String operator+(const String &str1, const String &str2) {
  String res = str1;
  res += str2;
  return res;
}

bool operator==(const String &str1, const String &str2) {
  if (str1.size != str2.size) {
    return false;
  }
  return strcmp(str1.str, str2.str) == 0;
}

std::ostream &operator<<(std::ostream &out, const String &str) {
  out << str.str;
  return out;
}

std::istream &operator>>(std::istream &in, String &str) {
  str.clear();
  char symbol;
  while (in.get(symbol) && !isspace(symbol) && symbol != '\0') {
    str.push_back(symbol);
  }
  return in;
}