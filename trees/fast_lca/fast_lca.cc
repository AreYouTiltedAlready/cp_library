#include <algorithm>
#include <cassert>
#include <memory>
#include <numeric>
#include <vector>

enum class RmqMode {
  kMax,
  kMin,
};

template <typename T, RmqMode mode>
class RmqSolver {
 public:
  template <typename U, typename std::enable_if_t<
                            std::is_same_v<std::decay_t<U>, std::vector<T>>,
                            void>* = nullptr>
  explicit RmqSolver(U values) : values_(values) {
    const int n = static_cast<int>(values.size());  // n = 0 is not allowed
    const int log = std::__lg(n * 2);
    matrix_.resize(log);
    for (int i = 0; i < log; ++i) { matrix_[i].resize(n - (1 << i) + 1); }
    std::iota(matrix_.front().begin(), matrix_.front().end(), 0);
    for (int i = 1; i < log; ++i) {
      for (int j = 0; j < n - (1 << i) + 1; ++j) {
        matrix_[i][j] =
            Merger(matrix_[i - 1][j], matrix_[i - 1][j + (1 << (i - 1))]);
      }
    }
  }

  // [first, last) : first == last is not allowed
  [[nodiscard]] T GetIndex(int first, int last) const {
    const auto level = std::__lg(last - first);
    return Merger(matrix_[level][first], matrix_[level][last - (1 << level)]);
  }

  [[nodiscard]] T GetValue(int first, int last) const {
    return values_[GetIndex(first, last)];
  }

 private:
  inline int Merger(int lhs, int rhs) const {
    if constexpr (mode == RmqMode::kMin) {
      return values_[lhs] < values_[rhs] ? lhs : rhs;
    }
    return values_[lhs] > values_[rhs] ? lhs : rhs;
  }

  std::vector<T> values_;
  std::vector<std::vector<int>> matrix_;
};

template <typename T>
using RmqSolverMax = RmqSolver<T, RmqMode::kMax>;

template <typename T>
using RmqSolverMin = RmqSolver<T, RmqMode::kMin>;

class LcaForest {
 public:
  explicit LcaForest(int n)
      : n_(n),
        tour_id_(0),
        euler_id_(0),
        in_(n),
        out_(n),
        tour_(n),
        size_(n),
        head_(n),
        depth_(n),
        entry_(n),
        parent_(n),
        euler_(n * 2),
        g_(n),
        rmq_ptr_(nullptr) {}

  explicit LcaForest(const std::vector<std::vector<int>>& g)
      : LcaForest(static_cast<int>(g.size())) {
    Build();
  }

  void AddEdge(int u, int v) {
    g_[u].push_back(v);
    g_[v].push_back(u);
  }

  void Build(std::vector<int> roots = std::vector<int>(1, 0)) {
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
      g_[i].erase(std::remove(g_[i].begin(), g_[i].end(), parent_[i]),
                  g_[i].end());
      if (g_[i].empty()) { continue; }
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
    rmq_ptr_ = std::make_unique<RmqSolverMin<int>>(std::move(euler_depths));
  }

  // IsAncestor(u, u) is true for all u
  bool IsAncestor(int u, int v) const {
    return in_[u] <= in_[v] && in_[v] < out_[u];
  }

  // Is z on the path from u to v
  bool LiesOnPath(int z, int u, int v) const {
    return IsAncestor(Lca(u, v), z) && (IsAncestor(z, u) || IsAncestor(z, v));
  }

  int Lca(int u, int v) const {
    assert(rmq_ptr_ != nullptr);
    if (entry_[u] > entry_[v]) { std::swap(u, v); }
    return euler_[rmq_ptr_->GetIndex(entry_[u], entry_[v] + 1)];
  }

  // obviously, 0-indexed
  int LevelAncestor(int v, int h) const {
    if (!(0 <= h && h <= depth_[v])) { return -1; }
    while (depth_[head_[v]] > h) { v = parent_[head_[v]]; }
    return tour_[in_[head_[v]] + h - depth_[head_[v]]];
  }

  int KthAncestor(int v, int k) const {
    return LevelAncestor(v, depth_[v] - k);
  }

  // 0-indexed
  // yields -1 if distance(u, v) > k
  int KthNodeOnPath(int u, int v, int k) const {
    int z = Lca(u, v);
    int du = depth_[u] - depth_[z];
    int dv = depth_[v] - depth_[z];
    if (!(0 <= k && k <= du + dv)) { return -1; }
    if (k <= du) { return KthAncestor(u, k); }
    return KthAncestor(v, du + dv - k);
  }

 private:
  void SizeDfs(int v) {
    size_[v] = 1;
    for (int to : g_[v]) {
      if (to == parent_[v]) { continue; }
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

  int n_;
  int tour_id_;
  int euler_id_;
  std::vector<int> in_;
  std::vector<int> out_;
  std::vector<int> euler_;
  std::vector<int> entry_;
  std::vector<int> head_;
  std::vector<int> tour_;
  std::vector<int> size_;
  std::vector<int> depth_;
  std::vector<int> parent_;
  std::unique_ptr<RmqSolverMin<int>> rmq_ptr_;
  std::vector<std::vector<int>> g_;
};