#include <functional>
#include <vector>

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun) : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

class HLD {
 public:
  explicit HLD(int n)
      : n_(n),
        in_(n),
        out_(n),
        tour_(n),
        size_(n),
        head_(n),
        depth_(n),
        parent_(n),
        g_(n) {}

  void AddEdge(int u, int v) {
    g_[u].push_back(v);
    g_[v].push_back(u);
  }

  void Build(int root = 0) {
    parent_[root] = root;
    y_combinator([&](auto&& dfs, int v) -> void {
      size_[v] = 1;
      for (int to : g_[v]) {
        if (to == parent_[v]) {
          continue;
        }
        parent_[to] = v;
        depth_[to] = depth_[v] + 1;
        dfs(to);
        size_[v] += size_[to];
      }
    })(root);

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
    y_combinator([&](auto&& dfs, int v) -> void {
      tour_[time] = v;
      in_[v] = time++;
      bool heavy_child = true;
      for (int to : g_[v]) {
        head_[to] = heavy_child ? head_[v] : to;
        heavy_child = false;
        dfs(to);
      }
      out_[v] = time;
    })(root);
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
  // returns v if k > distance(u, v)
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
  int n_;
  std::vector<int> in_;
  std::vector<int> out_;
  std::vector<int> head_;
  std::vector<int> tour_;
  std::vector<int> size_;
  std::vector<int> depth_;
  std::vector<int> parent_;
  std::vector<std::vector<int>> g_;
};
