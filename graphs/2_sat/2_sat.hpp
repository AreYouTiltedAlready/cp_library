// #include "graphs/graph_traits.hpp"
// #include "graphs/scc/scc.hpp"

#include <optional>

namespace graphs {
namespace two_sat {

class TwoSat {
 public:
  explicit TwoSat(int n, int m) : n_(n), g_(2 * n, 2 * m) {}

  void Implies(int u, int value_u, int v, int value_v) {
    g_.AddEdge(u * 2 + value_u, v * 2 + value_v);
  }

  void AddOrCondition(int u, int value_u, int v, int value_v) {
    Implies(u, value_u ^ 1, v, value_v);
    Implies(v, value_v ^ 1, u, value_u);
  }

  [[nodiscard]] std::optional<std::vector<char>> Solve() const {
    std::vector<char> result(n_);
    const std::vector<int> scc = scc::FindSCC(g_);
    for (int i = 0; i < n_; ++i) {
      if (scc[i * 2] == scc[i * 2 + 1]) {
        return std::nullopt;
      }
      result[i] = static_cast<char>(scc[i * 2] < scc[i * 2 + 1]);
    }
    return result;
  }

 private:
  int n_;
  Digraph g_;
};

}  // namespace two_sat

}  // namespace graphs
