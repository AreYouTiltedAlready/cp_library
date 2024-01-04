#include <vector>

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
