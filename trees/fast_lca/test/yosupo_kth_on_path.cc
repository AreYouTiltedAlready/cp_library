// Problem: https://judge.yosupo.jp/problem/jump_on_tree
// Submission: https://judge.yosupo.jp/submission/181701

#include <algorithm>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

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
namespace rmq {
namespace internal {

template <typename T, typename Comp>
class RMQSolver {
 public:
  static constexpr int kBlockLength = 64;

  template <typename U,
            std::enable_if_t<std::is_same_v<std::decay_t<U>, std::vector<T>>,
                             void>* = nullptr>
  explicit RMQSolver(U&& values)
      : values_(std::forward<U>(values)),
        small_blocks_(values_.size()),
        comp_(),
        n_(static_cast<int>(values_.size())),
        blocks_count_((values_.size() + kBlockLength - 1) / kBlockLength) {
    // building last 64 elements of monotonic stack for each index
    {
      std::vector<int> stack;
      stack.reserve(n_);
      stack.push_back(0);
      small_blocks_[0] = 1;
      for (int i = 1; i < n_; ++i) {
        small_blocks_[i] = small_blocks_[i - 1] << 1;
        while (!stack.empty() && values_[stack.back()] >= values_[i]) {
          if (int d = i - stack.back(); d < kBlockLength) {
            small_blocks_[i] ^= static_cast<uint64_t>(1) << d;
          }
          stack.pop_back();
        }
        stack.push_back(i);
        small_blocks_[i] += 1;
      }
    }

    // Building the sparse table for blocks of size n / kBlockLength
    {
      const int blocks_log =
          std::bit_width(static_cast<uint32_t>(blocks_count_));
      blocks_table_.resize(blocks_log);
      for (int i = 0; i < blocks_log; ++i) {
        blocks_table_[i].resize(blocks_count_ - (1 << i) + 1);
      }
      for (int i = 0; i < blocks_count_; ++i) {
        blocks_table_[0][i] = i * kBlockLength;
      }
      for (int i = 0; i < n_; ++i) {
        blocks_table_[0][i / kBlockLength] =
            Merge(blocks_table_[0][i / kBlockLength], i);
      }
      for (int i = 1; i < blocks_log; ++i) {
        for (int j = 0; j < blocks_count_ - (1 << i) + 1; ++j) {
          blocks_table_[i][j] = Merge(blocks_table_[i - 1][j],
                                      blocks_table_[i - 1][j + (1 << (i - 1))]);
        }
      }
    }
  }

  [[nodiscard]] int GetIndex(int first, int last) const {
    auto [first_q, first_r] = std::div(first, kBlockLength);
    auto [last_q, last_r] = std::div(last, kBlockLength);

    if (first_q == last_q) {
      return GetSmallBlock(last - 1, last - first);
    }

    int result = first;
    if (first_r != 0) {
      result = Merge(result, GetSmallBlock(first + kBlockLength - 1 - first_r,
                                           kBlockLength - first_r));
      first_q += 1;
    }
    if (last_r != 0) {
      result = Merge(result, GetSmallBlock(last - 1, last_r));
    }
    if (first_q < last_q) {
      result = Merge(result, GetOnBlocks(first_q, last_q));
    }

    return result;
  }

  [[nodiscard]] T GetValue(int first, int last) const {
    return values_[GetIndex(first, last)];
  }

 private:
  [[nodiscard]] inline int GetSmallBlock(int right, int length) const {
    return right + 1 -
           std::bit_width(small_blocks_[right] &
                          ((static_cast<uint64_t>(1) << length) - 1));
  }

  [[nodiscard]] inline int GetOnBlocks(int first, int last) const {
    int level = std::bit_width(static_cast<uint32_t>(last - first)) - 1;
    return Merge(blocks_table_[level][first],
                 blocks_table_[level][last - (1 << level)]);
  }

  [[nodiscard]] inline int Merge(int lhs, int rhs) const {
    return comp_(values_[lhs], values_[rhs]) ? lhs : rhs;
  }

  std::vector<T> values_;
  std::vector<uint64_t> small_blocks_;
  std::vector<std::vector<int>> blocks_table_;

  const Comp comp_;

  int n_;
  int blocks_count_;
};

}  // namespace internal

template <typename T>
using RMQMinSolver = internal::RMQSolver<T, std::less<>>;

template <typename T>
using RMQMaxSolver = internal::RMQSolver<T, std::greater<>>;

}  // namespace rmq
}  // namespace ds

namespace trees {
namespace fast_lca {
// Okay, this is just a generalization of basic hld (on top of hld, we maintain
// euler tour for lca in $O(1))
// Note: similar to hld, one must call Build() before queries
// In case of construction from adjacency list, the Build() is called
// immediately All queries except of LA are O(1) now
class LcaForest {
 public:
  explicit LcaForest(int n)
      : in_(n),
        out_(n),
        euler_(n * 2),
        entry_(n),
        head_(n),
        tour_(n),
        size_(n),
        depth_(n),
        parent_(n),
        g_(n),
        rmq_solver_(nullptr),
        n_(n),
        tour_id_(0),
        euler_id_(0) {}

  explicit LcaForest(const std::vector<std::vector<int>>& g)
      : LcaForest(static_cast<int>(g.size())) {
    Build();
  }

  void AddEdge(int u, int v) {
    g_[u].push_back(v);
    g_[v].push_back(u);
  }

  void Build(const std::vector<int>& roots = std::vector<int>(1, 0)) {
    if (!roots.empty()) {
      for (int root : roots) {
        head_[root] = parent_[root] = root;
        SizeDfs(root);
      }
    } else {
      for (int i = 0; i < n_; ++i)
        if (size_[i] == 0) {
          head_[i] = parent_[i] = i;
          SizeDfs(i);
        }
    }
    for (int i = 0; i < n_; ++i) {
      std::erase(g_[i], parent_[i]);
      if (g_[i].empty()) {
        continue;
      }
      std::swap(
          *(g_[i].begin()),
          *std::max_element(g_[i].begin(), g_[i].end(),
                            [&](int u, int v) { return size_[u] < size_[v]; }));
    }
    if (!roots.empty()) {
      for (int root : roots) {
        TourDfs(root);
        euler_[euler_id_++] = -1;
      }
    } else {
      for (int i = 0; i < n_; ++i) {
        if (out_[i] == 0) {
          TourDfs(i);
          euler_[euler_id_++] = -1;
        }
      }
    }
    std::vector<int> euler_depths(2 * n_);
    for (int i = 0; i < 2 * n_; ++i) {
      if (euler_[i] < 0) {
        euler_depths[i] = -1;
      } else {
        euler_depths[i] = depth_[euler_[i]];
      }
    }
    rmq_solver_ = new ds::rmq::RMQMinSolver<int>(std::move(euler_depths));
  }

  // IsAncestor(u, u) is true for all u
  [[nodiscard]] bool IsAncestor(int u, int v) const {
    return in_[u] <= in_[v] && in_[v] < out_[u];
  }

  // Is z on the path from u to v
  [[nodiscard]] bool LiesOnPath(int z, int u, int v) const {
    return IsAncestor(Lca(u, v), z) && (IsAncestor(z, u) || IsAncestor(z, v));
  }

  // -1 if not connected
  [[nodiscard]] int Lca(int u, int v) const {
    if (entry_[u] > entry_[v]) {
      std::swap(u, v);
    }
    return euler_[rmq_solver_->GetIndex(entry_[u], entry_[v] + 1)];
  }

  // obviously, 0-indexed
  [[nodiscard]] int LevelAncestor(int v, int h) const {
    if (!(0 <= h && h <= depth_[v])) {
      return -1;
    }
    while (depth_[head_[v]] > h) {
      v = parent_[head_[v]];
    }
    return tour_[in_[head_[v]] + h - depth_[head_[v]]];
  }

  [[nodiscard]] int KthAncestor(int v, int k) const {
    return LevelAncestor(v, depth_[v] - k);
  }

  // 0-indexed
  // yields -1 if distance(u, v) > k
  [[nodiscard]] int KthNodeOnPath(int u, int v, int k) const {
    int z = Lca(u, v);
    int du = depth_[u] - depth_[z];
    int dv = depth_[v] - depth_[z];
    if (!(0 <= k && k <= du + dv)) {
      return -1;
    }
    if (k <= du) {
      return KthAncestor(u, k);
    }
    return KthAncestor(v, du + dv - k);
  }

  ~LcaForest() { delete rmq_solver_; }

 private:
  void SizeDfs(int v) {
    size_[v] = 1;
    for (int to : g_[v]) {
      if (to == parent_[v]) {
        continue;
      }
      parent_[to] = v;
      depth_[to] = depth_[v] + 1;
      SizeDfs(to);
      size_[v] += size_[to];
    }
  }

  void TourDfs(int v) {
    tour_[tour_id_] = v;
    in_[v] = tour_id_++;
    entry_[v] = euler_id_;
    euler_[euler_id_++] = v;
    bool heavy = true;
    for (int to : g_[v]) {
      head_[to] = heavy ? head_[v] : to;
      heavy = false;
      TourDfs(to);
      euler_[euler_id_++] = v;
    }
    out_[v] = tour_id_;
  }

  std::vector<int> in_;
  std::vector<int> out_;
  std::vector<int> euler_;
  std::vector<int> entry_;
  std::vector<int> head_;
  std::vector<int> tour_;
  std::vector<int> size_;
  std::vector<int> depth_;
  std::vector<int> parent_;
  std::vector<std::vector<int>> g_;
  ds::rmq::RMQMinSolver<int>* rmq_solver_;

  int n_;
  int tour_id_;
  int euler_id_;
};

}  // namespace fast_lca
}  // namespace trees

void RunCase([[maybe_unused]] int testcase) {
  int n, q;
  std::cin >> n >> q;

  trees::fast_lca::LcaForest lca_forest(n);
  for (int i = 1; i < n; ++i) {
    int u, v;
    std::cin >> u >> v;
    lca_forest.AddEdge(u, v);
  }

  lca_forest.Build();
  while (q--) {
    int u, v, k;
    std::cin >> u >> v >> k;
    std::cout << lca_forest.KthNodeOnPath(u, v, k) << "\n";
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
