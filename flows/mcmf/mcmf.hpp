#include <algorithm>
#include <vector>
#include <ext/pb_ds/priority_queue.hpp>

namespace flows {
namespace mcmf {
// Calculates MCMF in $O(VE + |F|E\log V)
// Uses Ford-Bellman for potentials initialization and Dijkstra on heap for
// finding $s \rightarrow t$ paths.
// Note: if initial network is acyclic, potentials can be calculated in $O(E)$
// with simple DAG dp
// Note: if graph is dense (e.g. Assignment problem),
// dijkstra with linear search is better
// Negative cycles are not allowed
// If you need some auxiliary stuff like decomposition, mincut etc., take a look
// on Dinic implementation - all these utils can be implemented in the same way
template <typename Flow, typename Cost>
class MinCostFlowGraph {
 public:
  using flow_t = Flow;
  using cost_t = Cost;

  struct Edge {
    Edge() : cap(0), flow(0), cost(0), from(0), to(0) {}
    Edge(int from, int to, flow_t cap, cost_t cost)
        : cap(cap), flow(0), cost(cost), from(from), to(to) {}

    flow_t cap;
    flow_t flow;
    cost_t cost;
    int from;
    int to;
  };

  explicit MinCostFlowGraph(int n, int m, int source, int sink)
      : edges_(m * 2),
        pot_(n),
        distance_(n),
        path_(n),
        sp_edge_(n),
        g_(n),
        heap_(),
        its_(n),
        n_(n),
        m_(m),
        source_(source),
        sink_(sink) {}

  int AddEdge(int from, int to, flow_t cap, cost_t cost) {
    int id = m_;
    edges_[m_++] = {from, to, cap, cost};
    edges_[m_++] = {to, from, 0, -cost};
    g_[from].push_back(id);
    g_[to].push_back(id + 1);
    return id;
  }

  [[nodiscard]] const Edge& GetEdge(int id) const { return edges_[id]; }

  [[nodiscard]] const std::vector<Edge>& Edges() const { return edges_; }

  void FindPath() {
    heap_.clear();
    std::fill(distance_.begin(), distance_.end(), kUnreachable);
    distance_[source_] = 0;

    std::fill(its_.begin(), its_.end(), heap_.end());
    its_[source_] = heap_.push(std::make_pair(0, source_));
    while (!heap_.empty()) {
      auto [d, v] = heap_.top();
      heap_.pop();
      for (int e_id : g_[v]) {
        const Edge& edge = edges_[e_id];
        if (edge.cap == edge.flow) {
          continue;
        }
        if (cost_t new_cost = -d + edge.cost + pot_[edge.from] - pot_[edge.to];
            new_cost < distance_[edge.to]) {
          sp_edge_[edge.to] = e_id;
          distance_[edge.to] = new_cost;
          auto new_pair = std::make_pair(-new_cost, edge.to);
          if (its_[edge.to] == heap_.end()) {
            its_[edge.to] = heap_.push(new_pair);
          } else {
            heap_.modify(its_[edge.to], new_pair);
          }
        }
      }
    }

    if (distance_[sink_] != kUnreachable) {
      std::replace(distance_.begin(), distance_.end(), kUnreachable,
                   static_cast<cost_t>(0));
      for (int i = 0; i < n_; ++i) {
        pot_[i] += distance_[i];
      }
    }
  }

  std::pair<flow_t, cost_t> MinCostFlow(
      flow_t flow_limit = std::numeric_limits<flow_t>::max()) {
    flow_t flow = 0;
    cost_t cost = 0;

    const bool negative_edges = [&]() -> bool {
      for (int i = 0; i < edges_.size() / 2; ++i) {
        if (edges_[i * 2].cost < 0) {
          return true;
        }
      }
      return false;
    }();

    if (negative_edges) {
      std::fill(pot_.begin(), pot_.end(), kUnreachable);
      pot_[source_] = 0;
      bool any = true;
      while (any) {
        any = false;
        for (int i = 0; i < static_cast<int>(edges_.size()); i += 2) {
          const Edge& edge = edges_[i];
          if (pot_[edge.from] == kUnreachable) {
            continue;
          }
          if (cost_t new_cost = pot_[edge.from] + edge.cost;
              new_cost < pot_[edge.to]) {
            any = true;
            pot_[edge.to] = new_cost;
          }
        }
      }
      std::replace(pot_.begin(), pot_.end(), kUnreachable,
                   static_cast<cost_t>(0));
    }

    FindPath();
    while (flow < flow_limit && distance_[sink_] != kUnreachable) {
      flow_t path_min = flow_limit - flow;
      int path_length = [&]() -> int {
        int id = 0;
        int v = sink_;
        while (v != source_) {
          int e_id = sp_edge_[v];
          path_[id++] = e_id;
          const Edge& edge = edges_[e_id];
          v = edge.from;
          path_min = std::min(path_min, edge.cap - edge.flow);
        }
        return id;
      }();

      cost_t additional = 0;
      flow += path_min;
      for (int i = 0; i < path_length; ++i) {
        Edge& edge = edges_[path_[i]];
        Edge& back_edge = edges_[path_[i] ^ 1];
        edge.flow += path_min;
        back_edge.flow -= path_min;
        additional += edge.cost;
      }

      cost += additional * path_min;
      FindPath();
    }

    return std::make_pair(flow, cost);
  }

 private:
  static constexpr cost_t kUnreachable = std::numeric_limits<cost_t>::max() / 2;

  std::vector<Edge> edges_;
  std::vector<cost_t> pot_;
  std::vector<cost_t> distance_;
  std::vector<int> path_;
  std::vector<int> sp_edge_;
  std::vector<std::vector<int>> g_;
  __gnu_pbds::priority_queue<std::pair<cost_t, int>> heap_;
  std::vector<typename decltype(heap_)::point_iterator> its_;

  int n_;
  int m_;
  const int source_;
  const int sink_;
};

}  // namespace mcmf
}  // namespace flows
