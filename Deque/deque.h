#include <vector>
#include <algorithm>

template<typename T>
class Deque {
 private:

  std::vector<T *> buffer;
  const static size_t chunk_size = 8;
  size_t size_ = 0;
  size_t buffer_size = 0;
  size_t head_chunk_ind = 0;
  size_t head_buf_ind = 0;

 private:
  void resize(size_t new_buffer_size) {
    if (new_buffer_size < buffer_size) {
      return;
    }
    std::vector<T *> temp_buffer;

    size_t new_head_buf_ind = (new_buffer_size - buffer_size) / 2;
    for (size_t i = 0; i < new_head_buf_ind; ++i) {
      temp_buffer.push_back(reinterpret_cast<T *>(new int8_t[chunk_size * sizeof(T)]));
    }
    for (size_t i = 0; i < buffer_size; ++i) {
      temp_buffer.push_back(buffer[i]);
    }
    for (size_t i = new_head_buf_ind + buffer_size; i < new_buffer_size; ++i) {
      temp_buffer.push_back(reinterpret_cast<T *>(new int8_t[chunk_size * sizeof(T)]));
    }
    std::swap(buffer, temp_buffer);
    buffer_size = new_buffer_size;
    head_chunk_ind += new_head_buf_ind * chunk_size;
    head_chunk_ind %= chunk_size;
    head_buf_ind = new_head_buf_ind;
  }

 public:

  // constructors

  Deque() {
    resize(1);
  }

  Deque(size_t size) {
    resize(size / chunk_size + 1);
    for (size_t i = 0; i < size; ++i) {
      push_back(T());
    }
  }

  Deque(size_t size, const T &value) : Deque() {
    for (size_t i = 0; i < size; ++i) {
      push_back(value);
    }
  }

  Deque(const Deque &other) : Deque() {
    try {
      for (size_t i = 0; i < other.size_; ++i) {
        push_back(other[i]);
      }
    } catch (...) {
      throw;
    }
  }

  void swap(Deque &other) {
    std::swap(size_, other.size_);
    std::swap(buffer_size, other.buffer_size);
    std::swap(head_chunk_ind, other.head_chunk_ind);
    std::swap(head_buf_ind, other.head_buf_ind);
    std::swap(buffer, other.buffer);
  }

  Deque &operator=(const Deque &other) {
    Deque deque_copy = other;
    swap(deque_copy);
    return *this;
  }

  size_t size() const {
    return size_;
  }

  T &operator[](size_t pos) {
    return buffer[head_buf_ind + (head_chunk_ind + pos) / chunk_size][(head_chunk_ind + pos) % chunk_size];
  }

  const T &operator[](size_t pos) const {
    return buffer[head_buf_ind + (head_chunk_ind + pos) / chunk_size][(head_chunk_ind + pos) % chunk_size];
  }

  T &at(size_t pos) {
    if (pos >= size_) {
      throw std::out_of_range("...");
    }
    return buffer[head_buf_ind + (head_chunk_ind + pos) / chunk_size][(head_chunk_ind + pos) % chunk_size];
  }

  const T &at(size_t pos) const {
    if (pos >= size_) {
      throw std::out_of_range("...");
    }
    return buffer[head_buf_ind + (head_chunk_ind + pos) / chunk_size][(head_chunk_ind + pos) % chunk_size];
  }

  // iterators

  template<bool is_const>
  struct common_iterator {
   private:
    friend class Deque;
    friend class iterator;
    friend class const_iterator;

    size_t cur_buff_index;
    size_t cur_chunk_index;
    std::conditional_t<is_const, const T *, T *> ptr;

    const std::vector<T *> *vector_pointer;

   public:

    common_iterator(size_t cur_buff_index,
                    size_t cur_chunk_index,
                    std::conditional_t<is_const, const T *, T *> ptr,
                    const std::vector<T *> *vector_pointer)
        : cur_buff_index(cur_buff_index), cur_chunk_index(cur_chunk_index), ptr(ptr), vector_pointer(vector_pointer) {}

    common_iterator(size_t index, std::conditional_t<is_const, const Deque &, Deque &> deque_) {
      cur_chunk_index = (deque_.head_chunk_ind + index) % chunk_size;
      cur_buff_index = deque_.head_buf_ind + (deque_.head_chunk_ind + index) / chunk_size;
      ptr = deque_.buffer[cur_buff_index] + cur_chunk_index;
      vector_pointer = &deque_.buffer;
    }

    common_iterator &operator++() {
      ++cur_chunk_index;
      if (cur_chunk_index == chunk_size) {
        ++cur_buff_index;
        ptr = (*vector_pointer)[cur_buff_index];
        cur_chunk_index = 0;
      } else {
        ++ptr;
      }
      return *this;
    }

    common_iterator operator++(int) {
      common_iterator copy = *this;
      ++*this;
      return copy;
    }

    common_iterator &operator--() {
      if (cur_chunk_index == 0u) {
        cur_buff_index--;
        ptr = (*vector_pointer)[cur_buff_index] + chunk_size - 1;
        cur_chunk_index = chunk_size - 1;
      } else {
        --cur_chunk_index;
        --ptr;
      }
      return *this;
    }

    common_iterator operator--(int) {
      common_iterator copy = *this;
      --*this;
      return copy;
    }

    common_iterator &operator+=(int x) {
      size_t index = cur_buff_index * chunk_size + cur_chunk_index + x;
      cur_chunk_index = index % chunk_size;
      cur_buff_index = index / chunk_size;
      ptr = (*vector_pointer)[cur_buff_index] + cur_chunk_index;
      return *this;
    }

    common_iterator &operator-=(int x) {
      size_t index = cur_buff_index * chunk_size + cur_chunk_index - x;
      cur_chunk_index = index % chunk_size;
      cur_buff_index = index / chunk_size;
      ptr = (*vector_pointer)[cur_buff_index] + cur_chunk_index;
      return *this;
    }

    common_iterator operator+(int x) {
      common_iterator copy = *this;
      copy += x;
      return copy;
    }

    common_iterator operator-(int x) {
      common_iterator copy = *this;
      copy -= x;
      return copy;
    }

    bool operator==(const common_iterator &other) const {
      return ptr == other.ptr;
    }

    bool operator!=(const common_iterator &other) const {
      return ptr != other.ptr;
    }

    bool operator<(const common_iterator &other) const {
      return cur_buff_index < other.cur_buff_index
          || (cur_buff_index == other.cur_buff_index && cur_chunk_index < other.cur_chunk_index);
    }

    bool operator<=(const common_iterator &other) const {
      return *this < other || *this == other;
    }

    bool operator>(const common_iterator &other) const {
      return other < *this;
    }

    bool operator>=(const common_iterator &other) const {
      return other <= *this;
    }

    int operator-(const common_iterator &other) const {
      size_t index = cur_buff_index * chunk_size + cur_chunk_index;
      size_t other_index = other.cur_buff_index * chunk_size + other.cur_chunk_index;
      return index - other_index;
    }

    std::conditional_t<is_const, const T &, T &> operator*() const {
      return *ptr;
    }

    std::conditional_t<is_const, const T *, T *> operator->() const {
      return ptr;
    }

    operator common_iterator<true>() {
      return common_iterator<true>(cur_buff_index, cur_chunk_index, ptr, vector_pointer);
    }
  };

  using const_iterator = common_iterator<true>;
  using iterator = common_iterator<false>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() {
    return iterator(0, *this);
  }

  const_iterator begin() const {
    return const_iterator(0, *this);
  }

  iterator end() {
    return iterator(size_, *this);
  }

  const_iterator end() const {
    return const_iterator(size_, *this);
  }

  const_iterator cbegin() const {
    return const_iterator(0, *this);
  }

  const_iterator cend() const {
    return const_iterator(size_, *this);
  }

  reverse_iterator rbegin() {
    return std::make_reverse_iterator(end());
  }

  reverse_iterator rend() {
    return std::make_reverse_iterator(begin());
  }

  const_reverse_iterator crbegin() const {
    return std::make_reverse_iterator(cend());
  }

  const_reverse_iterator crend() const {
    return std::make_reverse_iterator(cbegin());
  }

  void push_back(T value) {
    if (head_buf_ind + (head_chunk_ind + size_) / chunk_size == buffer_size) {
      resize(3 * buffer_size);
    }
    try {
      new(buffer[head_buf_ind + (head_chunk_ind + size_) / chunk_size]
              + (head_chunk_ind + size_) % chunk_size) T(value);
      ++size_;
    } catch (...) {
      throw;
    }
  }

  void push_front(T value) {
    if (head_chunk_ind == 0 && head_buf_ind == 0) {
      resize(3 * buffer_size);
    }
    try {
      new(buffer[head_chunk_ind == 0 ? head_buf_ind - 1 : head_buf_ind]
              + (head_chunk_ind - 1 + chunk_size) % chunk_size) T(value);
      ++size_;
      head_buf_ind = head_chunk_ind == 0 ? head_buf_ind - 1 : head_buf_ind;
      head_chunk_ind = (head_chunk_ind - 1 + chunk_size) % chunk_size;
    } catch (...) {
      throw;
    }
  }

  void pop_back() {
    size_t tail_buf_ind = head_buf_ind + (head_chunk_ind + size_) / chunk_size;
    size_t tail_chunk_ind = (head_chunk_ind + size_) % chunk_size;
    (buffer[tail_buf_ind] + tail_chunk_ind)->~T();
    size_--;
  }

  void pop_front() {
    (buffer[head_buf_ind] + head_chunk_ind)->~T();
    head_buf_ind += (head_chunk_ind + 1) / chunk_size;
    head_chunk_ind = (head_chunk_ind + 1) % chunk_size;
    size_--;
  }

  void erase(iterator it) {
    iterator it_ = it;
    iterator next_it = it_;
    ++next_it;
    for (; it_ < end(); ++it_, ++next_it) {
      std::swap(*it_, *(next_it));
    }
    pop_back();
  }

  void insert(iterator it, const T &value) {
    push_back(value);
    iterator it_ = end();
    --it_;
    iterator prev_it = it_;
    --prev_it;
    for (; it_ > it; --it_, --prev_it) {
      std::swap(*it_, *(prev_it));
    }
  }
};
