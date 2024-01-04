// #include "modular/barrett.hpp"

namespace modular {
namespace modint {
namespace internal {

struct ModIntBase {};

}  // namespace internal

template <typename T>
concept modint = std::is_base_of_v<internal::ModIntBase, T>;

}  // namespace modint
} // namespace modular
