#include <cstring>
#include <iostream>

class String {
 public:
  String();
  explicit String(const char *);
  String(size_t, char);
  String(const String &);
  ~String();
  String &operator=(String);

  template<class T> inline
  void swap(T &first, T &second);
  void swap(String& s);
  size_t length() const;
  void push_back(char);
  void pop_back();
  char front() const;
  char back() const;
  char &front();
  char &back();
  size_t find(const String &) const;
  size_t rfind(const String &) const;
  String substr(size_t, size_t) const;
  bool empty() const;
  void clear();

  char &operator[](size_t pos);
  char operator[](size_t pos) const;
  String &operator+=(char);
  String &operator+=(const String &);
  friend bool operator==(const String &, const String &);

  friend std::ostream &operator<<(std::ostream &, const String &);
  friend std::istream &operator>>(std::istream &in, String &str);

 private:
  char *str;
  size_t size;
  size_t capacity;
  void shrink();
  void cut_capacity();

};