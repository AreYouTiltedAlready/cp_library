#include <optional>
// #include "graphs/graph_traits.hpp"
// #include "utils/simple_queue.hpp"

namespace graphs {
namespace topological_order {

template <graph_traits::directed_graph Graph>
std::optional<std::vector<int>> TopologicalOrder(const Graph& graph) {
  const int n = graph.n();
  std::vector<int> in_degree(n);
  for (const auto& edge : graph.edges()) {
    in_degree[edge.to] += 1;
  }

  ::utils::simple_queue<int> q;
  for (int i = 0; i < n; ++i) {
    if (in_degree[i] == 0) {
      q.push(i);
    }
  }

  std::vector<int> order;
  order.reserve(n);
  const auto& g = graph.g();
  const auto& edges = graph.edges();
  while (!q.empty()) {
    int v = q.poll();
    order.push_back(v);
    for (int e_id : g[v]) {
      const auto& edge = edges[e_id];
      if (--in_degree[edge.to] == 0) {
        q.push(edge.to);
      }
    }
  }

  if (order.size() < n) {
    return std::nullopt;
  }
  return order;
}

}  // namespace topological_order

}  // namespace graphs
