// Problem: https://cses.fi/problemset/task/1132
// Submission: https://cses.fi/problemset/result/8100673/

#include <ext/pb_ds/priority_queue.hpp>
#include <iostream>
#include <optional>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun)  // NOLINT(*forwarding-reference*)
      : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

namespace utils {

template <typename T>
class simple_queue {
 public:
  simple_queue() : pos_(0), payload_({}) {}

  explicit simple_queue(int n) : simple_queue() { payload_.reserve(n); }

  void reserve(int n) { payload_.reserve(n); }

  void push(const T& value) { payload_.push_back(value); }

  void push(T&& value) { payload_.push_back(std::forward<T>(value)); }

  T poll() { return payload_[pos_++]; }

  [[nodiscard]] bool empty() const {
    return pos_ == static_cast<int>(payload_.size());
  }

  [[nodiscard]] int size() const {
    return static_cast<int>(payload_.size()) - pos_;
  }

 private:
  int pos_;
  std::vector<T> payload_;
};

}  // namespace utils

namespace graphs {
namespace graph_traits {
namespace internal {

struct BasicEdge {
  BasicEdge() : from(0), to(0) {}
  BasicEdge(int from, int to) : from(from), to(to) {}

  int from;
  int to;
};

template <typename T>
  requires std::is_arithmetic_v<T>
struct WeightedEdge : public BasicEdge {
  using cost_t = T;

  WeightedEdge() : BasicEdge(), cost(0) {}
  WeightedEdge(int from, int to, T cost) : BasicEdge(from, to), cost(cost) {}

  T cost;
};

template <typename T>
concept basic_edge = std::same_as<T, BasicEdge>;

template <typename T>
concept weighted_edge =
    std::is_base_of_v<BasicEdge, T> &&
    std::same_as<
        std::enable_if_t<std::is_arithmetic_v<typename T::cost_t>, void>, void>;

template <typename T>
concept edge = basic_edge<T> || weighted_edge<T>;

enum class GraphType { kDirected, kUndirected };

template <edge Edge, GraphType kGraphType>
class GraphBase {
 public:
  [[nodiscard]] int n() const { return n_; }

  [[nodiscard]] int m() const { return m_; }

  [[nodiscard]] const Edge& edge(int id) { return edges_[id]; }

  [[nodiscard]] const std::vector<Edge>& edges() const { return edges_; }

  [[nodiscard]] const std::vector<std::vector<int>>& g() const { return g_; }

 protected:
  explicit GraphBase(int n, int m) : g_(n), edges_(), n_(n), m_(0) {
    edges_.reserve(m);
  }

  std::vector<Edge> edges_;
  std::vector<std::vector<int>> g_;

  const int n_;
  int m_;
};

template <edge Edge>
class Digraph : public GraphBase<Edge, GraphType::kDirected> {
  using base = GraphBase<Edge, GraphType::kDirected>;
  using base::edges_;
  using base::g_;
  using base::m_;
  using base::n_;

 public:
  using edge_t = Edge;

  explicit Digraph(int n, int m) : base(n, m) {}

  template <typename... Args>
    requires std::is_constructible_v<edge_t, Args...>
  void AddEdge(Args&&... args) {
    int id = m_++;
    edges_.emplace_back(std::forward<Args>(args)...);
    g_[edges_.back().from].push_back(id);
  }
};

template <edge Edge>
class Graph : public GraphBase<Edge, GraphType::kUndirected> {
  using base = GraphBase<Edge, GraphType::kUndirected>;
  using base::edges_;
  using base::g_;
  using base::m_;
  using base::n_;

 public:
  using edge_t = Edge;
  explicit Graph(int n, int m) : base(n, m) {}

  template <typename... Args>
    requires std::is_constructible_v<edge_t, Args...>
  void AddEdge(Args&&... args) {
    int id = m_++;
    edges_.emplace_back(std::forward<Args>(args)...);
    g_[edges_.back().from].push_back(id);
    g_[edges_.back().to].push_back(id);
  }
};

template <typename T>
concept has_edge_t =
    std::same_as<std::enable_if_t<edge<typename T::edge_t>, void>, void>;

}  // namespace internal

template <typename T>
concept undirected_graph =
    internal::has_edge_t<T> &&
    std::is_base_of_v<internal::GraphBase<typename T::edge_t,
                                          internal::GraphType::kUndirected>,
                      T>;

template <typename T>
concept directed_graph =
    internal::has_edge_t<T> &&
    std::is_base_of_v<
        internal::GraphBase<typename T::edge_t, internal::GraphType::kDirected>,
        T>;

template <typename T>
concept weighted_graph =
    internal::has_edge_t<T> && internal::weighted_edge<typename T::edge_t>;

template <typename T>
concept graph = undirected_graph<T> || directed_graph<T>;

}  // namespace graph_traits

using Undigraph =
    graph_traits::internal::Graph<graph_traits::internal::BasicEdge>;
using Digraph =
    graph_traits::internal::Digraph<graph_traits::internal::BasicEdge>;

template <typename Cost>
using WeightedUndigraph =
    graph_traits::internal::Graph<graph_traits::internal::WeightedEdge<Cost>>;

template <typename Cost>
using WeightedDigraph =
    graph_traits::internal::Digraph<graph_traits::internal::WeightedEdge<Cost>>;

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

template <typename Distance, graph_traits::weighted_graph Graph>
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
SSSPInfo<Distance> Dijkstra(const Graph& graph, int source) {
  return Dijkstra<Distance, Graph>(graph, std::vector<int>(1, source));
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

void RunCase([[maybe_unused]] int testcase) {
  int n;
  std::cin >> n;

  graphs::Undigraph g(n, n - 1);
  for (int i = 1; i < n; ++i) {
    int u;
    int v;
    std::cin >> u >> v;
    --u;
    --v;
    g.AddEdge(u, v);
  }

  int first = [&](int v) -> int {
    auto [d, p] = graphs::shortest_paths::Bfs(g, v);
    return static_cast<int>(
        std::distance(d.cbegin(), std::max_element(d.cbegin(), d.cend())));
  }(0);

  auto [d, p] = graphs::shortest_paths::Bfs(g, first);
  int last = static_cast<int>(
      std::distance(d.cbegin(), std::max_element(d.cbegin(), d.cend())));

  std::vector<int> diameter = [&](int v) -> std::vector<int> {
    std::vector<int> res;
    res.reserve(d[v]);
    int id = p[v];
    while (id != -1) {
      res.push_back(v);
      const auto& edge = g.edge(id);
      v ^= edge.from ^ edge.to;
      id = p[v];
    }
    res.push_back(v);
    return res;
  }(last);

  const int d_size = static_cast<int>(diameter.size());
  std::vector<int> diam_id(n);
  for (int i = 0; i < d_size; ++i) {
    diam_id[diameter[i]] = i;
  }

  auto [d_all, p_all] = graphs::shortest_paths::Bfs(g, diameter);
  std::vector<int> ans(n);
  std::vector<int> diam_parent(n, -1);
  for (int v : diameter) {
    diam_parent[v] = v;
  }
  for (int i = 0; i < n; ++i) {
    int par{};

    {
      int u = i;
      while (diam_parent[u] == -1) {
        const auto& edge = g.edge(p_all[u]);
        u ^= edge.from ^ edge.to;
      }
      par = diam_parent[u];
    }

    {
      int u = i;
      while (diam_parent[u] == -1) {
        diam_parent[u] = par;
        const auto& edge = g.edge(p_all[u]);
        u ^= edge.from ^ edge.to;
      }
    }

    ans[i] = d_all[i] + std::max(diam_id[par], d_size - 1 - diam_id[par]);
  }

  for (int it : ans) {
    std::cout << it << " ";
  }
}

void Main() {
  int testcases = 1;
  // std::cin >> testcases;
  for (int tt = 1; tt <= testcases; ++tt) {
    RunCase(tt);
  }
}

}  // namespace

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  Main();
  return 0;
}
