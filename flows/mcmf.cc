#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

template <typename F, typename C>
class MinCostFlowNetwork {
 public:
  struct Edge {
    Edge() : cap(0), flow(0), cost(0), from(0), to(0) {}
    Edge(int from, int to, F cap, C cost)
        : cap(cap), flow(0), cost(cost), from(from), to(to) {}

    F Residual() const { return cap - flow; }

    F cap;
    F flow;
    C cost;
    int from;
    int to;
  };

  explicit MinCostFlowNetwork(int n, int m, int source, int sink)
      : n_(n),
        source_(source),
        sink_(sink),
        pot_(n),
        path_(n),
        distance_(n),
        edges_(),
        sp_edge_(n),
        adj_list_(n) {
    edges_.reserve(m * 2);
  }

  int AddEdge(int from, int to, F cap, C cost) {
    int id = static_cast<int>(edges_.size());
    edges_.emplace_back(from, to, cap, cost);
    edges_.emplace_back(to, from, 0, -cost);
    adj_list_[from].push_back(id);
    adj_list_[to].push_back(id + 1);
    return id;
  }

  [[nodiscard]] const Edge& GetEdge(int id) const { return edges_[id]; }

  void PotentialsInit() {
    std::fill(pot_.begin(), pot_.end(), kUnreachable);
    pot_[source_] = 0;

    bool any = true;
    while (any) {
      any = false;
      for (int i = 0; i < static_cast<int>(edges_.size()); i += 2) {
        const Edge& e = edges_[i];
        if (C new_cost = pot_[e.from] + e.cost; new_cost < pot_[e.to]) {
          any = true;
          pot_[e.to] = new_cost;
        }
      }
    }

    std::replace(pot_.begin(), pot_.end(), kUnreachable, static_cast<C>(0));
  }

  bool Dijkstra() {
    std::fill(distance_.begin(), distance_.end(), kUnreachable);
    distance_[source_] = 0;

    std::priority_queue<std::pair<C, int>, std::vector<std::pair<C, int>>,
                        std::greater<>>
        heap;
    heap.emplace(0, source_);
    while (!heap.empty()) {
      auto [d, v] = heap.top();
      heap.pop();
      if (distance_[v] != d) {
        continue;
      }
      for (int edge_id : adj_list_[v]) {
        const Edge& e = edges_[edge_id];
        if (e.Residual() == 0) {
          continue;
        }
        if (C new_cost = d + e.cost + pot_[e.from] - pot_[e.to];
            new_cost < distance_[e.to]) {
          sp_edge_[e.to] = edge_id;
          distance_[e.to] = new_cost;
          heap.emplace(new_cost, e.to);
        }
      }
    }

    if (distance_[sink_] != kUnreachable) {
      std::replace(distance_.begin(), distance_.end(), kUnreachable,
                   static_cast<C>(0));
      for (int i = 0; i < n_; ++i) {
        pot_[i] += distance_[i];
      }
      return true;
    }

    return false;
  }

  std::pair<F, C> MinCostFlow(F flow_limit = std::numeric_limits<F>::max()) {
    F flow = 0;
    C cost = 0;

    while (flow < flow_limit && Dijkstra()) {
      F path_min = flow_limit - flow;
      int path_length = [&]() -> int {
        int id = 0;
        int v = sink_;
        while (v != source_) {
          int edge_id = sp_edge_[v];
          path_[id++] = edge_id;
          v = edges_[edge_id].from;
          path_min = std::min(path_min, edges_[edge_id].Residual());
        }
        return id;
      }();

      flow += path_min;
      for (int i = 0; i < path_length; ++i) {
        edges_[path_[i]].flow += path_min;
        edges_[path_[i] ^ 1].flow -= path_min;
        cost += edges_[path_[i]].cost * path_min;
      }
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
  std::vector<std::vector<int>> adj_list_;
};

