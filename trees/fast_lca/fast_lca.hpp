// #include "data_structures/rmq/farach-colton-bender.hpp"

#include <algorithm>
#include <vector>

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
