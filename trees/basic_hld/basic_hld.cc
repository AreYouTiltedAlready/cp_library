#include <algorithm>
#include <vector>

// Desription: Basic heavy-light decomposition stuff
// Build time is $O(n)$
// Lca, KthAncestor, LevelAncestor, KthNodeOnPath - all these functions work in
// $O(\log n)$
// IsAncestor is the only cheap check which works in $O(1)$
class BasicHLD {
 public:
  explicit BasicHLD(int n)
      : n_(n),
        in_(n),
        out_(n),
        head_(n),
        tour_(n),
        size_(n),
        depth_(n),
        parent_(n),
        g_(n) {}

  explicit BasicHLD(const std::vector<std::vector<int>>& g)
      : BasicHLD(static_cast<int>(g.size())) {
    Build();
  }

  void AddEdge(int u, int v) {
    g_[u].push_back(v);
    g_[v].push_back(u);
  }

  void Build(int root = 0) {
    parent_[root] = root;
    auto Dfs = [&](auto& self, int v) -> void {
      size_[v] = 1;
      for (int to : g_[v]) {
        if (to == parent_[v]) {
          continue;
        }
        parent_[to] = v;
        depth_[to] = depth_[v] + 1;
        self(self, to);
        size_[v] += size_[to];
      }
    };
    Dfs(Dfs, root);
    for (int i = 0; i < n_; ++i) {
      g_[i].erase(std::remove(g_[i].begin(), g_[i].end(), parent_[i]),
                  g_[i].end());
      if (g_[i].empty()) {
        continue;
      }
      std::swap(
          *(g_[i].begin()),
          *std::max_element(g_[i].begin(), g_[i].end(),
                            [&](int u, int v) { return size_[u] < size_[v]; }));
    }
    int time = 0;
    head_[root] = root;
    auto TourDfs = [&](auto& self, int v) -> void {
      tour_[time] = v;
      in_[v] = time++;
      bool heavy = true;
      for (int to : g_[v]) {
        head_[to] = heavy ? head_[v] : to;
        heavy = false;
        self(self, to);
      }
      out_[v] = time;
    };
    TourDfs(TourDfs, root);
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
    while (head_[u] != head_[v]) {
      if (depth_[head_[u]] > depth_[head_[v]]) {
        std::swap(u, v);
      }
      v = parent_[head_[v]];
    }
    if (depth_[u] > depth_[v]) {
      std::swap(u, v);
    }
    return u;
  }

  // obviously, 0-indexed
  int LevelAncestor(int v, int h) const {
    if (!(0 <= h && h <= depth_[v])) {
      return -1;
    }
    while (depth_[head_[v]] > h) {
      v = parent_[head_[v]];
    }
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
    if (!(0 <= k && k <= du + dv)) {
      return -1;
    }
    if (k <= du) {
      return KthAncestor(u, k);
    }
    return KthAncestor(v, du + dv - k);
  }

 private:
  std::vector<int> in_;
  std::vector<int> out_;
  std::vector<int> head_;
  std::vector<int> tour_;
  std::vector<int> size_;
  std::vector<int> depth_;
  std::vector<int> parent_;
  std::vector<std::vector<int>> g_;

  int n_;
};
