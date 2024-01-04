// #include <graphs/graph_traits.hpp>
// #include <utils/simple_queue.hpp>

#include <chrono>
#include <ext/pb_ds/priority_queue.hpp>
#include <random>

namespace graphs {

namespace shortest_paths {

template <typename Cost>
struct SSSPInfo {
  using cost_t = Cost;
  static constexpr cost_t kUnreachableSentinel =
      std::numeric_limits<cost_t>::max() / 2;
  explicit SSSPInfo(int n) : distance(n, kUnreachableSentinel), p_edge(n, -1) {}

  std::vector<Cost> distance;
  std::vector<int> p_edge;
};

template <typename Cost>
struct FordBellmanInfo {
  explicit FordBellmanInfo(int n) : sp_info(n) {}

  SSSPInfo<Cost> sp_info;
  std::vector<int> negative_cycle;
};

template <typename Distance, graph_traits::weighted_graph Graph>
  requires graph_traits::directed_graph<Graph>
SSSPInfo<Distance> Dijkstra(const Graph& graph,
                            const std::vector<int>& source) {
  using distance_t = Distance;
  static constexpr distance_t kUnreachableSentinel =
      std::numeric_limits<distance_t>::max() / 2;

  SSSPInfo<distance_t> result(graph.n());
  __gnu_pbds::priority_queue<std::pair<distance_t, int>, std::greater<>> heap;
  std::vector<typename decltype(heap)::point_iterator> its(graph.n(),
                                                           heap.end());

  for (int s : source) {
    result.distance[s] = 0;
    its[s] = heap.push(std::make_pair(0, s));
  }

  const auto& g = graph.g();
  const auto& edges = graph.edges();
  std::vector<distance_t>& dist = result.distance;
  std::vector<int>& p_edge = result.p_edge;

  while (!heap.empty()) {
    auto [d, v] = heap.top();
    its[v] = heap.end();
    heap.pop();
    for (int e_id : g[v]) {
      const auto& edge = edges[e_id];
      int to = v ^ edge.from ^ edge.to;
      if (distance_t new_cost = d + edge.cost; new_cost < dist[to]) {
        p_edge[to] = e_id;
        dist[to] = new_cost;
        auto new_state = std::make_pair(new_cost, to);
        if (its[to] == heap.end()) {
          its[to] = heap.push(new_state);
        } else {
          heap.modify(its[to], new_state);
        }
      }
    }
  }

  return result;
}

template <typename Distance, graph_traits::weighted_graph Graph>
  requires graph_traits::directed_graph<Graph>
SSSPInfo<Distance> Dijkstra(const Graph& graph, int source) {
  return Dijkstra<Distance, Graph>(graph, std::vector<int>(1, source));
}

template <typename Distance, graph_traits::weighted_graph Graph>
  requires graph_traits::directed_graph<Graph>
FordBellmanInfo<Distance> FordBellman(const Graph& graph,
                                      const std::vector<int>& source) {
  using distance_t = Distance;
  FordBellmanInfo<distance_t> result(graph.n());

  for (int s : source) {
    result.sp_info.distance[s] = 0;
  }

  auto& sp_info = result.sp_info;
  auto Relax = [&graph = std::as_const(graph), &sp_info](int e_id) -> bool {
    const auto& edge = graph.edge(e_id);
    if (distance_t new_cost = sp_info.distance[edge.from] + edge.cost;
        new_cost < sp_info.distance[edge.to]) {
      sp_info.p_edge[edge.to] = e_id;
      sp_info.distance[edge.to] = new_cost;
      return true;
    }
    return false;
  };

  static std::mt19937_64 rng(
      std::chrono::steady_clock::now().time_since_epoch().count());
  std::vector<int> order(graph.m());
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);

  bool any_relaxation = true;
  for (int i = 0; i < graph.n() - 1 && any_relaxation; ++i) {
    any_relaxation = false;
    for (int id : order) {
      any_relaxation |= Relax(id);
    }
  }

  if (!any_relaxation) {
    return result;
  }

  int last_edge = [&order = std::as_const(order), &Relax]() -> int {
    for (int id : order) {
      if (Relax(id)) {
        return id;
      }
    }
    return -1;
  }();

  if (last_edge == -1) {
    return result;
  }

  std::vector<char> seen(graph.m());
  seen[last_edge] = true;

  {
    std::vector<int> edges_stack;
    edges_stack.reserve(2 * graph.m());
    edges_stack.push_back(last_edge);
    int e_id = sp_info.p_edge[graph.edge(last_edge).from];
    while (!seen[e_id]) {
      seen[e_id] = true;
      edges_stack.push_back(e_id);
      e_id = sp_info.p_edge[graph.edge(e_id).from];
    }
    while (edges_stack.back() != e_id) {
      result.negative_cycle.push_back(edges_stack.back());
      edges_stack.pop_back();
    }
    result.negative_cycle.push_back(e_id);
  }

  return result;
}

template <typename Distance, graph_traits::weighted_graph Graph>
  requires graph_traits::directed_graph<Graph>
FordBellmanInfo<Distance> FordBellman(const Graph& graph, int source) {
  return FordBellman<Distance, Graph>(graph, source);
}

template <graph_traits::graph Graph>
SSSPInfo<int> Bfs(const Graph& graph, const std::vector<int>& source) {
  SSSPInfo<int> result(graph.n());
  ::utils::simple_queue<int> q;
  for (int s : source) {
    q.push(s);
    result.distance[s] = 0;
  }

  const auto& g = graph.g();
  const auto& edges = graph.edges();
  std::vector<int>& dist = result.distance;
  std::vector<int>& p_edge = result.p_edge;

  while (!q.empty()) {
    int v = q.poll();
    for (int e_id : g[v]) {
      const auto& edge = edges[e_id];
      int to = v ^ edge.from ^ edge.to;
      if (dist[to] == SSSPInfo<int>::kUnreachableSentinel) {
        q.push(to);
        p_edge[to] = e_id;
        dist[to] = dist[v] + 1;
      }
    }
  }

  return result;
}

template <graph_traits::graph Graph>
SSSPInfo<int> Bfs(const Graph& graph, int source) {
  return Bfs(graph, std::vector<int>(1, source));
}

}  // namespace shortest_paths

}  // namespace graphs
