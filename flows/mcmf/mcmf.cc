#include <algorithm>
#include <ext/pb_ds/priority_queue.hpp>

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
template <typename F, typename C>
class MinCostFlowGraph {
 public:
  struct Edge {
    Edge() : cap(0), flow(0), cost(0), from(0), to(0) {}
    Edge(int from, int to, F cap, C cost)
        : cap(cap), flow(0), cost(cost), from(from), to(to) {}

    F cap;
    F flow;
    C cost;
    int from;
    int to;
  };

  explicit MinCostFlowGraph(int n, int m, int source, int sink)
      : n_(n),
        source_(source),
        sink_(sink),
        pot_(n),
        path_(n),
        distance_(n),
        edges_(),
        sp_edge_(n),
        g_(n),
        heap_(),
        its_(n) {
    edges_.reserve(m * 2);
  }

  int AddEdge(int from, int to, F cap, C cost) {
    int id = static_cast<int>(edges_.size());
    edges_.emplace_back(from, to, cap, cost);
    edges_.emplace_back(to, from, 0, -cost);
    g_[from].push_back(id);
    g_[to].push_back(id + 1);
    return id;
  }

  [[nodiscard]] const Edge& GetEdge(int id) const { return edges_[id]; }

  void FindPath() {
    heap_.clear();
    std::fill(distance_.begin(), distance_.end(), kUnreachable);
    distance_[source_] = 0;

    std::fill(its_.begin(), its_.end(), heap_.end());
    its_[source_] = heap_.push(std::make_pair(0, source_));
    while (!heap_.empty()) {
      auto [d, v] = heap_.top();
      heap_.pop();
      for (int eid : g_[v]) {
        const Edge& e = edges_[eid];
        if (e.cap == e.flow) {
          continue;
        }
        if (C new_cost = -d + e.cost + pot_[e.from] - pot_[e.to];
            new_cost < distance_[e.to]) {
          sp_edge_[e.to] = eid;
          distance_[e.to] = new_cost;
          auto new_pair = std::make_pair(-new_cost, e.to);
          if (its_[e.to] == heap_.end()) {
            its_[e.to] = heap_.push(new_pair);
          } else {
            heap_.modify(its_[e.to], new_pair);
          }
        }
      }
    }

    if (distance_[sink_] != kUnreachable) {
      std::replace(distance_.begin(), distance_.end(), kUnreachable,
                   static_cast<C>(0));
      for (int i = 0; i < n_; ++i) {
        pot_[i] += distance_[i];
      }
    }
  }

  std::pair<F, C> MinCostFlow(F flow_limit = std::numeric_limits<F>::max()) {
    F flow = 0;
    C cost = 0;

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
          const Edge& e = edges_[i];
          if (pot_[e.from] == kUnreachable) {
            continue;
          }
          if (C new_cost = pot_[e.from] + e.cost; new_cost < pot_[e.to]) {
            any = true;
            pot_[e.to] = new_cost;
          }
        }
      }
      std::replace(pot_.begin(), pot_.end(), kUnreachable, static_cast<C>(0));
    }

    FindPath();
    while (flow < flow_limit && distance_[sink_] != kUnreachable) {
      F path_min = flow_limit - flow;
      int path_length = [&]() -> int {
        int id = 0;
        int v = sink_;
        while (v != source_) {
          int eid = sp_edge_[v];
          path_[id++] = eid;
          const Edge& e = edges_[eid];
          v = e.from;
          path_min = std::min(path_min, e.cap - e.flow);
        }
        return id;
      }();

      C additional = 0;
      flow += path_min;
      for (int i = 0; i < path_length; ++i) {
        Edge& e = edges_[path_[i]];
        Edge& back = edges_[path_[i] ^ 1];
        e.flow += path_min;
        back.flow -= path_min;
        additional += e.cost;
      }

      cost += additional * path_min;
      FindPath();
    }

    return std::make_pair(flow, cost);
  }

 private:
  static constexpr C kUnreachable = std::numeric_limits<C>::max() / 2;

  const int n_;
  const int source_;
  const int sink_;
  std::vector<C> pot_;
  std::vector<int> path_;
  std::vector<C> distance_;
  std::vector<Edge> edges_;
  std::vector<int> sp_edge_;
  std::vector<std::vector<int>> g_;
  __gnu_pbds::priority_queue<std::pair<C, int>> heap_;
  std::vector<typename decltype(heap_)::point_iterator> its_;
};
