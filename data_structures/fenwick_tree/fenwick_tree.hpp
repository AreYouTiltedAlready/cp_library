#include <bit>
#include <cstdint>
#include <vector>

namespace ds {
namespace fenwick_tree {

template <typename S>
class FenwickTree {
 public:
  explicit FenwickTree(int n) : tree_(n + 1), n_(n) {}

  void Apply(int pos, const S& value) {
    pos += 1;
    while (pos <= n_) {
      tree_[pos] += value;
      pos += pos & -pos;
    }
  }

  [[nodiscard]] S Get(int last) const {
    S res{};
    while (last > 0) {
      res += tree_[last];
      last -= last & -last;
    }
    return res;
  }

  [[nodiscard]] S Get(int first, int last) const {
    return Get(last) - Get(first);
  }

  [[nodiscard]] int LowerBound(S value) const {
    int res = 0;
    for (int i = std::bit_width(static_cast<uint32_t>(n_)); i >= 0; --i) {
      if (int next = res + (1 << i); next <= n_ && tree_[next] < value) {
        res = next;
        value -= tree_[res];
      }
    }
    return res;
  }

 private:
  std::vector<S> tree_;
  int n_;
};

}  // namespace fenwick_tree
}  // namespace ds
