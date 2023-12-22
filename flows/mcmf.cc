#include <algorithm>
#include <limits>
#include <optional>
#include <queue>
#include <vector>

template <typename F, typename C>
class MCMFNetwork {
 public:
  struct Edge {
    Edge() : cap(0), flow(0), cost(0), from(0), to(0) {}
    Edge(int from, int to, F cap, C cost)
        : cap(cap), flow(0), cost(cost), from(from), to(to) {}
 
    F Residual() const { return cap - flow; };
 
    F cap;
    F flow;
    C cost;
    int from;
    int to;
  };
 
  static constexpr C kUnreachable = std::numeric_limits<C>::max() / 2;
  explicit MCMFNetwork(int n, int m, int source, int sink)
      : n_(n),
        source_(source),
        sink_(sink),
        pot_(n),
        distance_(n),
        parent_(n),
        adj_list_(n) {
    edges_.reserve(m * 2);
  }
 
  [[nodiscard]] const Edge& GetEdge(int id) const { return edges_[id]; }
 
  int AddEdge(int from, int to, F cap, C cost, F rev_cap = 0,
              std::optional<C> rev_cost = std::nullopt) {
    int id = static_cast<int>(edges_.size());
    edges_.emplace_back(from, to, cap, cost);
    edges_.emplace_back(to, from, rev_cap, rev_cost.value_or(-cost));
    adj_list_[from].push_back(id);
    adj_list_[to].push_back(id + 1);
    return id;
  }
 
  void PotentialsInit() {
    std::fill(pot_.begin(), pot_.end(), kUnreachable);
    pot_[source_] = 0;
 
    bool any = true;
    while (any) {
      any = false;
      for (int i = 0; i < static_cast<int>(edges_.size()); i += 2) {
        const Edge& e = edges_[i];
        if (pot_[e.from] + e.cost < pot_[e.to]) {
          any = true;
          pot_[e.to] = pot_[e.from] + e.cost;
        }
      }
    }
 
    std::replace(pot_.begin(), pot_.end(), kUnreachable, 0);
  }
 
  std::pair<F, C> MinCostFlow(F flow_limit = std::numeric_limits<F>::max()) {
    F flow = 0;
    C cost = 0;
    auto Walk = [&](int v, auto&& callback) -> void {
      while (v != source_) {
        int id = parent_[v];
        callback(id);
        v = edges_[id].from;
      }
    };
 
    while (flow < flow_limit && FindPath()) {
      C current_cost = 0;
      F path_min = flow_limit - flow;
      Walk(sink_, [&](int id) -> void {
        path_min = std::min(path_min, edges_[id].Residual());
      });
      Walk(sink_, [&](int id) -> void {
        Push(id, path_min);
        current_cost += edges_[id].cost;
      });
      flow += path_min;
      cost += current_cost * path_min;
      for (int i = 0; i < n_; ++i) {
        if (distance_[i] != kUnreachable) {
          pot_[i] += distance_[i];
        }
      }
    }
 
    return std::make_pair(flow, cost);
  }
 
 private:
  bool FindPath() {
    std::fill(distance_.begin(), distance_.end(), kUnreachable);
    distance_[source_] = 0;
 
    std::priority_queue<std::pair<C, int>, std::vector<std::pair<C, int>>,
                        std::greater<std::pair<C, int>>>
        heap;
    heap.emplace(0, source_);
    while (!heap.empty()) {
      const auto [distance, v] = heap.top();
      heap.pop();
      if (distance_[v] != distance) {
        continue;
      }
      for (int edge_id : adj_list_[v]) {
        const Edge& e = edges_[edge_id];
        if (e.Residual() == 0) {
          continue;
        }
        if (C new_cost = distance + e.cost + pot_[e.from] - pot_[e.to];
            new_cost < distance_[e.to]) {
          parent_[e.to] = edge_id;
          distance_[e.to] = new_cost;
          heap.emplace(new_cost, e.to);
        }
      }
    }
 
    return distance_[sink_] != kUnreachable;
  }
 
  void Push(int id, F flow) {
    edges_[id].flow += flow;
    edges_[id ^ 1].flow -= flow;
  }
 
  const int n_;
  const int source_;
  const int sink_;
  std::vector<C> pot_;
  std::vector<C> distance_;
  std::vector<int> parent_;
  std::vector<Edge> edges_;
  std::vector<std::vector<int>> adj_list_;
};
