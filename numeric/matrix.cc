#include <cstdint>
#include <type_traits>
#include <vector>

template <typename T>
class Matrix {
 public:
  explicit Matrix(int n, int m) : n_(n), m_(m), mat_(n, std::vector<T>(m)) {}

  std::vector<T>& operator[](int id) { return mat_[id]; }

  const std::vector<T>& operator[](int id) const { return mat_[id]; }

  Matrix<T>& operator+=(const Matrix<T>& other) {
    assert(n_ == other.n_ && m_ == other.m_);
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < m_; ++j) {
        (*this)[i][j] += other[i][j];
      }
    }
    return *this;
  }

  Matrix<T>& operator-=(const Matrix<T>& other) {
    assert(n_ == other.n_ && m_ == other.m_);
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < m_; ++j) {
        (*this)[i][j] -= other[i][j];
      }
    }
    return *this;
  }

  Matrix<T>& operator*=(T value) {
    for (std::vector<T>& row : mat_) {
      for (T& element : row) {
        element *= value;
      }
    }
    return *this;
  }

  Matrix<T>& operator/=(T value) {
    for (std::vector<T>& row : mat_) {
      for (T& element : row) {
        element /= value;
      }
    }
    return *this;
  }

  static Matrix<T> Identity(int n) {
    Matrix<T> result(n, n);
    for (int i = 0; i < n; ++i) {
      result[i][i] = T(1);
    }
    return result;
  }

  [[nodiscard]] int n() const noexcept { return n_; }
  [[nodiscard]] int m() const noexcept { return m_; }
  [[nodiscard]] bool IsSquare() { return n_ == m_; }

 private:
  int n_;
  int m_;
  std::vector<std::vector<T>> mat_;
};

template <typename T, typename U>
Matrix<std::common_type_t<T, U>> operator+(const Matrix<T>& lhs,
                                           const Matrix<U>& rhs) {
  using V = std::common_type_t<T, U>;
  assert(lhs.n() == rhs.n() && lhs.m() == rhs.m());
  Matrix<V> result(lhs.n(), lhs.m());
  for (int i = 0; i < lhs.n(); ++i) {
    for (int j = 0; j < lhs.m(); ++j) {
      result[i][j] = lhs[i][j] + rhs[i][j];
    }
  }
  return result;
}

template <typename T, typename U>
Matrix<std::common_type_t<T, U>> operator-(const Matrix<T>& lhs,
                                           const Matrix<U>& rhs) {
  using V = std::common_type_t<T, U>;
  assert(lhs.n() == rhs.n() && lhs.m() == rhs.m());
  Matrix<V> result(lhs.n(), lhs.m());
  for (int i = 0; i < lhs.n(); ++i) {
    for (int j = 0; j < lhs.m(); ++j) {
      result[i][j] = lhs[i][j] - rhs[i][j];
    }
  }
  return result;
}

template <typename T, typename U>
Matrix<std::common_type_t<T, U>> operator*(const Matrix<T>& lhs,
                                           const Matrix<U>& rhs) {
  using V = std::common_type_t<T, U>;
  assert(lhs.m() && rhs.n());
  Matrix<V> result(lhs.n(), rhs.m());
  for (int i = 0; i < lhs.n(); ++i) {
    for (int j = 0; j < rhs.m(); ++j) {
      for (int k = 0; k < lhs.m(); ++k) {
        result[i][j] += lhs[i][k] * rhs[k][j];
      }
    }
  }
  return result;
}

template <typename T>
Matrix<T> Power(Matrix<T> mat, uint64_t n) {
  assert(mat.IsSquare());
  Matrix<T> result = Matrix<T>::Identity(mat.n());

  while (n > 0) {
    if (n % 2 == 1) {
      result = result * mat;
    }
    mat = mat * mat;
    n /= 2;
  }

  return result;
}