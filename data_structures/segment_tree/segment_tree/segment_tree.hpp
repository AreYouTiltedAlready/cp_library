#include <algorithm>
#include <bit>
#include <cstdint>
#include <functional>
#include <vector>

namespace ds {
namespace segment_tree {

namespace internal {

template <typename T>
constexpr bool has_binary_plus = requires(const T& lhs, const T& rhs) {
  { lhs + rhs } -> std::same_as<T>;
};

template <typename T>
concept monoid = std::is_default_constructible_v<T> && has_binary_plus<T>;

}  // namespace internal

template <internal::monoid S>
class SegmentTree {
 public:
  explicit SegmentTree(int n)
      : tree_(std::bit_ceil(static_cast<uint32_t>(n)) * 2),
        n_(n),
        size_(static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n)))),
        log_(std::bit_width(static_cast<uint32_t>(size_) - 1)) {}

  template <typename U>
    requires std::is_assignable_v<S&, U>
  explicit SegmentTree(int n, const U& value) : SegmentTree(n) {
    std::ranges::fill(tree_, value);
    for (int i = size_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
    requires std::is_assignable_v<S&, U>
  explicit SegmentTree(const std::vector<U>& values)
      : SegmentTree(static_cast<int>(values.size())) {
    std::ranges::copy(values, tree_.begin() + size_);
    for (int i = size_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
    requires std::is_assignable_v<S&, U>
  void Set(int pos, const U& value) {
    pos += size_;
    tree_[pos] = value;
    DoPull(pos);
  }

  void Apply(int pos, const S& value) {
    pos += size_;
    tree_[pos] = tree_[pos] + value;
    DoPull(pos);
  }

  [[nodiscard]] S Get() const { return tree_[1]; }

  [[nodiscard]] S Get(int pos) const {
    pos += size_;
    return tree_[pos];
  }

  [[nodiscard]] S Get(int first, int last) const {
    first += size_;
    last += size_;
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

  [[nodiscard]] int FindFirst(int first, int last,
                              const Predicate& pred) const {
    first += size_;
    last += size_;

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

  [[nodiscard]] int FindLast(int first, int last, const Predicate& pred) const {
    first += size_;
    last += size_;

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
                            DescentDirection direction) const {
    while (k < size_) {
      k <<= 1;
      k ^= static_cast<int>(direction);
      k ^= !pred(tree_[k]);
    }
    return k - size_;
  }

  void DoPull(int v) {
    for (int i = 1; i <= log_; ++i) {
      Pull(v >> i);
    }
  }

  inline void Pull(int k) { tree_[k] = tree_[k << 1] + tree_[k << 1 | 1]; }

  std::vector<S> tree_;
  int n_;
  int size_;
  int log_;
};

}  // namespace segment_tree

}  // namespace ds
