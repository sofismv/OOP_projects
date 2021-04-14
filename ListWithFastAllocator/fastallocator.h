#include <iostream>
#include <vector>

template<size_t chunkSize>
class FixedAllocator {
  struct Block {
    char data[chunkSize];
  };
  std::vector<Block *> buffer;
  size_t size_ = 1;
  const size_t buffer_size = 1024;
  std::vector<Block *> deallocated_pointers;
  size_t current_index = 0;

  void double_buffer() {
    if (size_ < buffer_size) {
      size_ *= 2;
    }
    auto *empty_block = new Block[size_];
    current_index = 0;
    buffer.push_back(empty_block);
  }

  FixedAllocator() {
    double_buffer();
  }

 public:

  void *allocate() {
    if (!deallocated_pointers.empty()) {
      Block *ptr = deallocated_pointers.back();
      deallocated_pointers.pop_back();
      return ptr;
    }
    if (current_index == size_ - 1) {
      double_buffer();
    }
    return buffer.back() + (current_index++);
  }

  void deallocate(void *ptr) {
    if (deallocated_pointers.size() < 10000000) {
      deallocated_pointers.push_back(static_cast<Block *>(ptr));
    }
  }

  FixedAllocator(const FixedAllocator &) = delete;

  static FixedAllocator &getInstance() {
    static FixedAllocator instance;
    return instance;
  }

  ~FixedAllocator() {
    for (auto block : buffer) {
      delete[] block;
    }
  }

};

template<typename T>
class FastAllocator {
  static FixedAllocator<sizeof(T)> *fixed_allocator_;

 public:
  using value_type = T;
  using pointer = T *;
  using const_pointer = const T *;
  using reference = T &;
  using const_reference = const T &;

  FastAllocator<T>() = default;
  FastAllocator<T>(const FastAllocator<T> &) {};

  template<typename U>
  FastAllocator<T>(const FastAllocator<U> &) {}

  template<class U>
  struct rebind {
    using other = FastAllocator<U>;
  };

  T *allocate(size_t n) {
    if (n == 1 && sizeof(T) <= 24) {
      return static_cast<T *>(fixed_allocator_->getInstance().allocate());
    } else {
      return std::allocator<T>().allocate(n);
    }
  }

  void deallocate(T *ptr, size_t n) {
    if (n == 1 && sizeof(T) <= 24) {
      fixed_allocator_->getInstance().deallocate(static_cast<void *>(ptr));
    } else {
      std::allocator<T>().deallocate(ptr, n);
    }
  }
};

template<typename A, typename B>
bool operator==(const FastAllocator<A>& a, const FastAllocator<B>& b) {
  return &a == &b;
}

template<typename A, typename B>
bool operator!=(const FastAllocator<A>& a, const FastAllocator<B>& b) {
  return !(a == b);
}

template<typename T, typename Allocator = std::allocator<T> >
class List {
  struct Node {
    T value;
    Node *prev = nullptr;
    Node *next = nullptr;

    Node() = default;
    Node(const T &value_) : value(value_) {}
  };
  using node_alloc = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;

  Node *fake;
  size_t size_ = 0;
  node_alloc allocator;

 public:

  explicit List(const Allocator &alloc = Allocator()) : allocator(alloc){
    fake = allocator.allocate(1);
    fake->prev = fake;
    fake->next = fake;
  }

  List(size_t count, const T &value, const Allocator &alloc = Allocator()) : List(alloc) {
    for (size_t i = 0; i < count; ++i) {
      push_back(value);
    }
  }

  List(size_t count, const Allocator &alloc = Allocator()) : List(alloc) {
    for (size_t i = 0; i < count; ++i) {
      Node *new_node = std::allocator_traits<node_alloc>::allocate(allocator, 1);
      std::allocator_traits<node_alloc>::construct(allocator, new_node);
      Node *new_node_prev = fake->prev;
      Node *new_node_next = fake;
      new_node->next = new_node_next;
      new_node->prev = new_node_prev;
      fake->prev->next = new_node;
      fake->prev = new_node;
    }
    size_ = count;
  }

  List(const List &other)
      : List(std::allocator_traits<Allocator>::select_on_container_copy_construction(other.get_allocator())) {
    for (auto it = other.cbegin(); it != other.cend(); ++it) {
      push_back(*it);
    }
  }

  List &operator=(const List &other) {
    bool alloc_copy = std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value;
    if (this == &other) {
      return *this;
    }
    while (size_ != 0) {
      pop_back();
    }
    allocator.deallocate(fake, 1);
    if (alloc_copy) {
      allocator = other.get_allocator();
    }
    fake = allocator.allocate(1);
    fake->next = fake;
    fake->prev = fake;
    for (auto it = other.cbegin(); it != other.cend(); ++it) {
      push_back(*it);
    }
    return *this;
  }

  size_t size() const {
    return size_;
  }

  void push_back(const T &value) {
    Node *new_node = std::allocator_traits<node_alloc>::allocate(allocator, 1);
    std::allocator_traits<node_alloc>::construct(allocator, new_node, value);
    Node *new_node_prev = fake->prev;
    Node *new_node_next = fake;
    new_node->next = new_node_next;
    new_node->prev = new_node_prev;
    fake->prev->next = new_node;
    fake->prev = new_node;
    size_++;
  }
  void push_front(const T &value) {
    Node *new_node = std::allocator_traits<node_alloc>::allocate(allocator, 1);
    std::allocator_traits<node_alloc>::construct(allocator, new_node, value);
    Node *new_node_prev = fake;
    Node *new_node_next = fake->next;
    new_node->next = new_node_next;
    new_node->prev = new_node_prev;
    fake->next->prev = new_node;
    fake->next = new_node;
    size_++;
  }
  void pop_back() {
    erase(--end());
  }
  void pop_front() {
    erase(begin());
  }

  Allocator get_allocator() const {
    return allocator;
  }

  template<bool is_const>
  struct common_iterator {
   private:
    friend class List;
    friend class iterator;
    friend class const_iterator;

    Node *ptr;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = typename std::conditional<is_const, const T, T>::type;
    using pointer = typename std::conditional<is_const, const T *, T *>::type;
    using reference = typename std::conditional<is_const, const T &, T &>::type;
    using iterator_category = std::bidirectional_iterator_tag;

    common_iterator() = default;
    common_iterator(Node *ptr_) : ptr(ptr_) {}

    common_iterator &operator++() {
      ptr = ptr->next;
      return *this;
    }
    common_iterator &operator--() {
      ptr = ptr->prev;
      return *this;
    }
    common_iterator operator++(int) {
      common_iterator copy = *this;
      ptr = ptr->next;
      return copy;
    }
    common_iterator operator--(int) {
      common_iterator copy = *this;
      ptr = ptr->prev;
      return copy;
    }
    bool operator==(const common_iterator &other) const {
      return ptr == other.ptr;
    }
    bool operator!=(const common_iterator &other) const {
      return ptr != other.ptr;
    }
    reference operator*() {
      return ptr->value;
    }
    pointer operator->() {
      return &ptr->value;
    }
    operator common_iterator<true>() {
      return common_iterator<true>(ptr);
    }
  };
  using const_iterator = common_iterator<true>;
  using iterator = common_iterator<false>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() {
    return iterator(fake->next);
  }

  const_iterator begin() const {
    return const_iterator(fake->next);
  }

  iterator end() {
    return iterator(fake);
  }

  const_iterator end() const {
    return const_iterator(fake);
  }

  const_iterator cbegin() const {
    return const_iterator(fake->next);
  }

  const_iterator cend() const {
    return const_iterator(fake);
  }

  const_reverse_iterator rbegin() const {
    return std::make_reverse_iterator(end());
  }

  const_reverse_iterator rend() const {
    return std::make_reverse_iterator(begin());
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

  void erase(const_iterator it) {
    Node *next = it.ptr->next;
    Node *prev = it.ptr->prev;
    prev->next = next;
    next->prev = prev;
    std::allocator_traits<node_alloc>::destroy(allocator, it.ptr);
    allocator.deallocate(it.ptr, 1);
    --size_;
  }

  void insert(const_iterator it, const T &value) {
    Node *new_node = std::allocator_traits<node_alloc>::allocate(allocator, 1);
    std::allocator_traits<node_alloc>::construct(allocator, new_node, value);
    Node *new_node_next = it.ptr;
    Node *new_node_prev = it.ptr->prev;
    new_node_next->prev = new_node;
    new_node->prev = new_node_prev;
    new_node->next = new_node_next;
    new_node_prev->next = new_node;
    ++size_;
  }

  ~List() {
    auto it = begin();
    for (size_t i = 0; i < size_; ++i) {
      auto prev = it++;
      std::allocator_traits<node_alloc>::destroy(allocator, prev.ptr);
      allocator.deallocate(prev.ptr, 1);
    }
    allocator.deallocate(fake, 1);
  }
};

