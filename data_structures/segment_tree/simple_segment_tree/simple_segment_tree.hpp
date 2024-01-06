#include <algorithm>
#include <type_traits>
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
  explicit SegmentTree(int n) : tree_(n * 2), n_(n) {}

  template <typename U>
  requires std::is_assignable_v<S&, U>
  explicit SegmentTree(int n, const U& value) : SegmentTree(n) {
    std::ranges::fill(tree_, value);
    for (int i = n_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
  requires std::is_assignable_v<S&, U>
  explicit SegmentTree(const std::vector<U>& values)
      : SegmentTree(static_cast<int>(values.size())) {
    std::ranges::copy(values, tree_.begin() + n_);
    for (int i = n_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
  requires std::is_assignable_v<S&, U>
  void Set(int pos, const U& value) {
    pos += n_;
    tree_[pos] = value;
    DoPull(pos);
  }

  void Apply(int pos, const S& value) {
    pos += n_;
    tree_[pos] = tree_[pos] + value;
    DoPull(pos);
  }

  [[nodiscard]] S Get() const { return tree_[1]; }

  [[nodiscard]] S Get(int pos) const {
    pos += n_;
    return tree_[pos];
  }

  [[nodiscard]] S Get(int first, int last) const {
    first += n_;
    last += n_;
    S res{};
    while (first < last) {
      if ((first & 1) == 1) {
        res = res + tree_[first++];
      }
      if ((last & 1) == 1) {
        res = res + tree_[--last];
      }
      first >>= 1;
      last >>= 1;
    }

    return res;
  }

 private:
  void DoPull(int v) {
    v /= 2;
    while (v > 0) {
      Pull(v);
      v /= 2;
    }
  }

  inline void Pull(int k) { tree_[k] = tree_[k << 1] + tree_[k << 1 | 1]; }

  std::vector<S> tree_;
  int n_;
};

}  // namespace segment_tree

}  // namespace ds
