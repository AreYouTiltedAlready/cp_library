#include <algorithm>
#include <vector>

std::vector<int> FindSCC(const std::vector<std::vector<int>>& g) {
  const int n = static_cast<int>(g.size());
  std::vector<int> order;
  order.reserve(n);
  {
    std::vector<char> used(n);
    auto Dfs = [&](auto& self, int v) -> void {
      used[v] = true;
      for (int to : g[v]) {
        if (!used[to]) {
          self(self, to);
        }
      }
      order.push_back(v);
    };

    for (int i = 0; i < n; ++i) {
      if (!used[i]) {
        Dfs(Dfs, i);
      }
    }
  }
  const auto g_rev = [](const std::vector<std::vector<int>>& g)
      -> std::vector<std::vector<int>> {
    const int n = static_cast<int>(g.size());
    std::vector<std::vector<int>> g_rev(n);
    for (int i = 0; i < n; ++i) {
      for (int to : g[i]) {
        g_rev[to].push_back(i);
      }
    }
    return g_rev;
  }(g);
  int current_id = 0;
  std::vector<int> scc_id(n, -1);
  auto Dfs = [&](auto& self, int v) -> void {
    scc_id[v] = current_id;
    for (int to : g_rev[v]) {
      if (scc_id[to] == -1) {
        self(self, to);
      }
    }
  };
  std::reverse(order.begin(), order.end());
  for (int id : order) {
    if (scc_id[id] == -1) {
      Dfs(Dfs, id);
      current_id += 1;
    }
  }
  return scc_id;
}
