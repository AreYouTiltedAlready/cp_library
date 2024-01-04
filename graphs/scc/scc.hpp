// #include "graphs/graph_traits.hpp"

namespace graphs {

namespace scc {
namespace internal {

template <graph_traits::directed_graph Graph>
Graph Reversed(const Graph& graph) {
  Graph result(graph.n(), graph.m());
  for (auto edge : graph.edges()) {
    std::swap(edge.from, edge.to);
    result.AddEdge(std::move(edge));
  }
  return result;
}

}  // namespace internal

template <graph_traits::directed_graph Graph>
std::vector<int> FindSCC(const Graph& graph) {
  const int n = graph.n();
  std::vector<int> order;
  order.reserve(n);

  {
    const auto& g = graph.g();
    const auto& edges = graph.edges();
    std::vector<char> used(n);
    auto Dfs = [&](auto& self, int v) -> void {
      used[v] = true;
      for (int e_id : g[v]) {
        int to = edges[e_id].to;
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

  const Graph rev_graph = internal::Reversed(graph);
  const auto& rev_g = rev_graph.g();
  const auto& rev_edges = rev_graph.edges();

  int current_id = 0;
  std::vector<int> scc_id(n, -1);
  auto Dfs = [&](auto& self, int v) -> void {
    scc_id[v] = current_id;
    for (int e_id : rev_g[v]) {
      int to = rev_edges[e_id].to;
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

}  // namespace scc

}  // namespace graphs
