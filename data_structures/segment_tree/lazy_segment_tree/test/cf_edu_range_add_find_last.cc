// Problem:
// https://codeforces.com/edu/course/2/lesson/5/3/practice/contest/280799/problem/C
// Submission:
// https://codeforces.com/edu/course/2/lesson/5/3/practice/contest/280799/submission/240235172

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun)  // NOLINT(*forwarding-reference*)
      : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

namespace ds {
namespace lazy_segment_tree {

namespace internal {

template <typename T>
constexpr bool has_binary_plus = requires(const T& lhs, const T& rhs) {
  { lhs + rhs } -> std::same_as<T>;
};

template <typename T>
concept monoid = std::is_default_constructible_v<T> && has_binary_plus<T>;

template <typename T, typename U>
concept lazy_tag = has_binary_plus<T> && requires(U& u, const T& t, int z) {
  { u.Apply(t, z) } -> std::same_as<void>;
};

}  // namespace internal

template <internal::monoid S, internal::lazy_tag<S> T>
class LazySegmentTree {
 public:
  explicit LazySegmentTree(int n)
      : tree_(std::bit_ceil(static_cast<uint32_t>(n)) * 2),
        tag_(std::bit_ceil(static_cast<uint32_t>(n))),
        n_(n),
        size_(static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n)))),
        log_(std::bit_width(static_cast<uint32_t>(size_) - 1)) {}

  template <typename U>
    requires std::is_assignable_v<S&, U>
  explicit LazySegmentTree(int n, const U& value) : LazySegmentTree(n) {
    std::ranges::fill(tree_, value);
    for (int i = size_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
    requires std::is_assignable_v<S&, U>
  explicit LazySegmentTree(const std::vector<U>& values)
      : LazySegmentTree(static_cast<int>(values.size())) {
    std::ranges::copy(values, tree_.begin() + size_);
    for (int i = size_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
    requires std::is_assignable_v<S&, U>
  void Set(int pos, const U& value) {
    pos += size_;
    for (int i = log_; i > 0; --i) {
      Push(pos >> i);
    }
    tree_[pos] = value;
    for (int i = 1; i <= log_; ++i) {
      Pull(pos >> i);
    }
  }

  void Apply(int pos, const T& tag) {
    pos += size_;
    for (int i = log_; i > 0; --i) {
      Push(pos >> i);
    }
    tree_[pos].Apply(tag, 1);
    for (int i = 1; i <= log_; ++i) {
      Pull(pos >> i);
    }
  }

  [[nodiscard]] S Get() const { return tree_[1]; }

  S Get(int pos) {
    pos += size_;
    for (int i = log_; i > 0; --i) {
      Push(pos >> i);
    }
    return tree_[pos];
  }

  void Apply(int first, int last, const T& tag) {
    first += size_;
    last += size_;
    DoPush(first, last);
    {
      int first_copy = first;
      int last_copy = last;
      while (first < last) {
        if ((first & 1) == 1) {
          ApplyAll(first++, tag);
        }
        if ((last & 1) == 1) {
          ApplyAll(--last, tag);
        }
        first >>= 1;
        last >>= 1;
      }
      first = first_copy;
      last = last_copy;
    }
    DoPull(first, last);
  }

  S Get(int first, int last) {
    first += size_;
    last += size_;
    DoPush(first, last);
    S res_left{};
    S res_right{};
    while (first < last) {
      if ((first & 1) == 1) {
        res_left = res_left + tree_[first++];
      }
      if ((last & 1) == 1) {
        res_right = tree_[--last] + res_right;
      }
      first >>= 1;
      last >>= 1;
    }
    return res_left + res_right;
  }

  using Predicate = std::function<bool(const S&)>;

  [[nodiscard]] int FindFirst(int first, int last, const Predicate& pred) {
    first += size_;
    last += size_;
    DoPush(first, last);
    int first_copy = first;
    int last_copy = last;
    while (first_copy < last_copy) {
      if ((first_copy & 1) == 1) {
        if (pred(tree_[first_copy])) {
          return Descent(first_copy, pred, DescentDirection::kLeft);
        }
        first_copy += 1;
      }

      first_copy >>= 1;
      last_copy >>= 1;
    }

    int height = std::bit_width(static_cast<uint32_t>(last_copy));
    for (int i = log_ - height; i >= 0; --i) {
      last_copy <<= 1;
      if (((last >> i) & 1) == 1) {
        if (pred(tree_[last_copy])) {
          return Descent(last_copy, pred, DescentDirection::kLeft);
        }
        last_copy += 1;
      }
    }

    return -1;
  }

  [[nodiscard]] int FindLast(int first, int last, const Predicate& pred) {
    first += size_;
    last += size_;
    DoPush(first, last);

    uint32_t mask = 1;
    int first_copy = first;
    int last_copy = last;
    while (first_copy < last_copy) {
      mask <<= 1;
      if ((first_copy & 1) == 1) {
        mask += 1;
        first_copy += 1;
      }

      if ((last_copy & 1) == 1) {
        last_copy -= 1;
        if (pred(tree_[last_copy])) {
          return Descent(last_copy, pred, DescentDirection::kRight);
        }
      }

      first_copy >>= 1;
      last_copy >>= 1;
    }

    while (!std::has_single_bit(mask)) {
      int z = std::countr_zero(mask) + 1;
      mask >>= z;
      first_copy <<= z;
      first_copy -= 1;
      if (pred(tree_[first_copy])) {
        return Descent(first_copy, pred, DescentDirection::kRight);
      }
    }

    return -1;
  }

 private:
  enum class DescentDirection {
    kLeft = 0,
    kRight = 1,
  };

  [[nodiscard]] int Descent(int k, const Predicate& pred,
                            DescentDirection direction) {
    while (k < size_) {
      Push(k);
      k <<= 1;
      k ^= static_cast<int>(direction);
      k ^= !pred(tree_[k]);
    }
    return k - size_;
  }

  void ApplyAll(int k, const T& tag) {
    int length = size_ >> (std::bit_width(static_cast<uint32_t>(k)) - 1);
    tree_[k].Apply(tag, length);
    if (length != 1) [[likely]] {
      tag_[k] = tag_[k] + tag;
    }
  }

  void DoPull(int first, int last) {
    const int last_common_bit = std::max(
        1, static_cast<int>(
               std::bit_width(static_cast<uint32_t>(first ^ (last - 1)))));
    int first_z = std::max(1, std::countr_zero(static_cast<uint32_t>(first)));
    for (int i = first_z; i < last_common_bit; ++i) {
      Pull(first >> i);
    }
    int last_z = std::max(1, std::countr_zero(static_cast<uint32_t>(last)));
    for (int i = last_z; i < last_common_bit; ++i) {
      Pull((last - 1) >> i);
    }
    for (int i = last_common_bit; i <= log_; ++i) {
      Pull(first >> i);
    }
  }

  void DoPush(int first, int last) {
    const int last_common_bit = std::max(
        1, static_cast<int>(
               std::bit_width(static_cast<uint32_t>(first ^ (last - 1)))));
    for (int i = log_; i >= last_common_bit; --i) {
      Push(first >> i);
    }
    int first_z = std::max(1, std::countr_zero(static_cast<uint32_t>(first)));
    for (int i = last_common_bit - 1; i >= first_z; --i) {
      Push(first >> i);
    }
    int last_z = std::max(1, std::countr_zero(static_cast<uint32_t>(last)));
    for (int i = last_common_bit - 1; i >= last_z; --i) {
      Push((last - 1) >> i);
    }
  }

  void Push(int k) {
    ApplyAll(k << 1, tag_[k]);
    ApplyAll(k << 1 | 1, tag_[k]);
    tag_[k] = T();
  }

  inline void Pull(int k) {
    tree_[k] = tree_[k << 1] + tree_[k << 1 | 1];
    int length = size_ >> (std::bit_width(static_cast<uint32_t>(k)) - 1);
    tree_[k].Apply(tag_[k], length);
  }

  std::vector<S> tree_;
  std::vector<T> tag_;
  const int n_;
  const int size_;
  const int log_;
};

}  // namespace lazy_segment_tree

}  // namespace ds

struct S {
  S() : value(0) {}
  S(int value) : value(value) {}

  void Apply(int tag, int) { value += tag; }

  friend S operator+(const S& lhs, const S& rhs) {
    return lhs.value < rhs.value ? rhs : lhs;
  }

  int value;
};

struct Searcher {
  explicit Searcher(int k) : k(k) {}

  bool operator()(const S& s) const { return s.value >= k; }
  int k;
};

void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  ds::lazy_segment_tree::LazySegmentTree<S, int> lazy_segment_tree(n);
  while (q--) {
    int t;
    std::cin >> t;
    if (t == 1) {
      int first;
      int last;
      int value;
      std::cin >> first >> last >> value;
      first = n - first;
      last = n - last;
      std::swap(first, last);
      lazy_segment_tree.Apply(first, last, value);
    } else {
      int k;
      int last;
      std::cin >> k >> last;
      int res = lazy_segment_tree.FindLast(0, n - last, Searcher(k));
      if (res != -1) {
        res = n - 1 - res;
      }
      std::cout << res << "\n";
    }
  }
}

void Main() {
  int testcases = 1;
  // std::cin >> testcases;
  for (int tt = 1; tt <= testcases; ++tt) {
    RunCase(tt);
  }
}

}  // namespace

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  Main();
  return 0;
}
