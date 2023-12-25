#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

// Casual non-recursive dsu implementation (path compression + weight heuristic)
// $O(\alpha(n))$ per query (amortized)

class Dsu {
 public:
  explicit Dsu(int n) : weight_(n, 1), parent_(n) {
    std::iota(parent_.begin(), parent_.end(), 0);
  }

  int Find(int v) const {
    int u = v;
    while (v != parent_[v]) { v = parent_[v]; }
    while (u != v) { u = std::exchange(parent_[u], v); }
    return v;
  }

  int Weight(int v) const { return weight_[Find(v)]; }

  bool Same(int u, int v) const { return Find(u) == Find(v); }

  void Join(int u, int v) {
    u = Find(u);
    v = Find(v);
    if (u == v) { return; }

    if (weight_[u] > weight_[v]) { std::swap(u, v); }
    parent_[u] = v;
    weight_[v] += weight_[u];
  }

 private:
  std::vector<int> weight_;
  mutable std::vector<int> parent_;
};

// https://judge.yosupo.jp/problem/unionfind
// https://judge.yosupo.jp/submission/178571
int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n;
  int q;
  std::cin >> n >> q;

  Dsu dsu(n);
  for (int i = 0; i < q; ++i) {
    int t;
    int u;
    int v;
    std::cin >> t >> u >> v;
    if (t == 0) {
      dsu.Join(u, v);
    } else {
      std::cout << dsu.Same(u, v) << "\n";
    }
  }

  return 0;
}
