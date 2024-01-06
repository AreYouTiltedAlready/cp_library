// #include "graphs/graph_traits.hpp"

namespace graphs {
namespace dfs_forest {

template <graph_traits::undirected_graph Graph>
class DfsForest {
 public:
  template <typename U>
  requires std::is_same_v<std::decay_t<U>, Graph>
  explicit DfsForest(U&& graph)
      : graph_(std::forward<U>(graph)),
        low_(graph_.n()),
        depth_(graph_.n(), -1),
        par_edge_(graph_.n(), -1) {
    Build();
  }

  [[nodiscard]] bool IsBridge(int e_id) const {
    const auto& edge = graph_.edge(e_id);
    auto [u, v] = std::make_pair(edge.from, edge.to);
    if (depth_[u] > depth_[v]) {
      std::swap(u, v);
    }
    return par_edge_[v] == e_id && low_[v] > depth_[u];
  }

  [[nodiscard]] bool IsArticulationPoint(int v) const {
    int children = 0;
    bool dangling_child = false;
    for (int e_id : graph_.g()[v]) {
      const auto& edge = graph_.edge(e_id);
      int to = v ^ edge.from ^ edge.to;
      if (par_edge_[to] == e_id) {
        children += 1;
        dangling_child |= low_[to] >= depth_[v];
      }
    }
    if (par_edge_[v] == -1) {
      return children > 1;
    }
    return dangling_child;
  }

  std::vector<int> EdgeBiconnectivityComponents() {
    int last_id = 0;
    std::vector<int> id(graph_.n(), -1);
    auto Dfs = [&](auto& self, int v) -> void {
      for (int e_id : graph_.g()[v]) {
        const auto& edge = graph_.edge(e_id);
        int to = v ^ edge.from ^ edge.to;
        if (id[to] != -1) {
          continue;
        }
        if (low_[to] > depth_[v]) {
          id[to] = last_id++;
        } else {
          id[to] = id[v];
        }
        self(self, to);
      }
    };

    for (int i = 0; i < graph_.n(); ++i) {
      if (id[i] == -1) {
        id[i] = last_id++;
        Dfs(Dfs, i);
      }
    }

    return id;
  }

  [[nodiscard]] std::vector<int> VertexBiconnectivityComponents() const {
    int last_id = 0;
    std::vector<int> id(graph_.m(), -1);
    std::vector<bool> used(graph_.n());
    auto Dfs = [&](auto& self, int v, int p_id) -> void {
      used[v] = true;
      for (int e_id : graph_.g()[v]) {
        if (id[e_id] != -1) {
          continue;
        }
        const auto& edge = graph_.edge(e_id);
        int to = v ^ edge.from ^ edge.to;
        if (!used[to]) {
          if (low_[to] >= depth_[v]) {
            id[e_id] = last_id++;
          } else {
            id[e_id] = id[p_id];
          }
          self(self, to, e_id);
        } else {
          id[e_id] = id[p_id];
        }
      }
    };
    for (int i = 0; i < graph_.n(); ++i) {
      if (par_edge_[i] == -1) {
        Dfs(Dfs, i, -1);
      }
    }
    return id;
  }

 private:
  void Build() {
    for (int i = 0; i < graph_.n(); ++i) {
      if (depth_[i] == -1) {
        depth_[i] = 0;
        Dfs(i);
      }
    }
  }

  void Dfs(int v) {
    low_[v] = depth_[v];
    for (int e_id : graph_.g()[v]) {
      if (e_id == par_edge_[v]) {
        continue;
      }
      const auto& edge = graph_.edge(e_id);
      int to = v ^ edge.from ^ edge.to;
      if (depth_[to] == -1) {
        par_edge_[to] = e_id;
        depth_[to] = depth_[v] + 1;
        Dfs(to);
        low_[v] = std::min(low_[v], low_[to]);
      } else {
        low_[v] = std::min(low_[v], depth_[to]);
      }
    }
  }

  Graph graph_;
  std::vector<int> low_;
  std::vector<int> depth_;
  std::vector<int> par_edge_;
};

}  // namespace dfs_forest
} // namespace graphs