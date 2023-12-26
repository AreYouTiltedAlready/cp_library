#include <type_traits>
#include <vector>

// Simple bottom-up segment tree
// Useful for easy invariants only
// Simple bottom-up segment tree
// Useful for easy invariants only
template <typename S, auto Op, auto E>
class SegmentTree {
  static_assert(std::is_convertible_v<decltype(Op), S (*)(S, S)>,
                "Op must work as S(S, S)");
  static_assert(std::is_convertible_v<decltype(E), S (*)()>,
                "E must work as S()");

 public:
  explicit SegmentTree(int n) : n_(n), tree_(n * 2, E()) {}

  template <typename U, typename std::enable_if_t<std::is_assignable_v<S&, U>,
                                                  void>* = nullptr>
  explicit SegmentTree(const std::vector<U>& values)
      : SegmentTree(static_cast<int>(values.size())) {
    std::copy(values.cbegin(), values.cend(), tree_.begin() + n_);
    for (int i = n_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  void Apply(int pos, const S& value) {
    tree_[pos += n_] = value;
    pos /= 2;
    while (pos > 0) {
      Pull(pos);
      pos /= 2;
    }
  }

  S Get(int left, int right) const {
    left += n_;
    right += n_;
    S res = E();
    while (left < right) {
      if (left % 2 == 1) {
        res = Op(res, tree_[left++]);
      }
      if (right % 2 == 1) {
        res = Op(res, tree_[--right]);
      }
      left /= 2;
      right /= 2;
    }
    return res;
  }

 private:
  inline void Pull(int k) { tree_[k] = Op(tree_[k * 2], tree_[k * 2 + 1]); }

  int n_;
  std::vector<S> tree_;
};
