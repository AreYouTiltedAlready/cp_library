#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun) : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

template <typename T, int kLog>
class FastFourierTransform {
 public:
  using Complex = std::complex<T>;

  FastFourierTransform() : root_({}), inverse_({}) {
    inverse_[0] = 0;
    for (int i = 1; i < 1 << kLog; ++i) {
      inverse_[i] = (inverse_[i >> 1] >> 1) + ((i & 1) << (kLog - 1));
    }
    root_[0] = root_[1] = 1;
    for (int i = 1; i < kLog; ++i) {
      const T angle = kPi / (1 << i);
      const Complex z(cosl(angle), sinl(angle));
      for (int j = 1 << (i - 1); j < 1 << i; ++j) {
        root_[j << 1] = root_[j];
        root_[j << 1 | 1] = root_[j] * z;
      }
    }
  }

  std::vector<int> ConvolutionMod(
      const typename std::vector<int>::const_iterator& lhs_first,
      const typename std::vector<int>::const_iterator& lhs_last,
      const typename std::vector<int>::const_iterator& rhs_first,
      const typename std::vector<int>::const_iterator& rhs_last, int mod) {
    const auto n = static_cast<int>(std::distance(lhs_first, lhs_last));
    const auto m = static_cast<int>(std::distance(rhs_first, rhs_last));
    int dft_size = 1;
    while (dft_size < std::max(n, m) * 2) { dft_size *= 2; }
    std::vector<Complex> dft_lhs(dft_size);
    for (int i = 0; i < n; ++i) {
      dft_lhs[i].real(lhs_first[i] % (1 << 15));
      dft_lhs[i].imag(lhs_first[i] >> 15);
    }
    std::vector<Complex> dft_rhs(dft_size);
    for (int i = 0; i < m; ++i) {
      dft_rhs[i].real(rhs_first[i] % (1 << 15));
      dft_rhs[i].imag(rhs_first[i] >> 15);
    }
    (*this)(dft_lhs.begin(), dft_lhs.end());
    (*this)(dft_rhs.begin(), dft_rhs.end());
    std::vector<Complex> edges(dft_size);
    std::vector<Complex> middle(dft_size);
    const Complex ratio(0.0, -0.25 / static_cast<int>(dft_size));
    for (int i = 0; i < dft_size; ++i) {
      int j = (dft_size - i) & (dft_size - 1);
      const Complex lhs_conj = std::conj(dft_lhs[j]);
      const Complex rhs_conj = std::conj(dft_rhs[j]);
      edges[j] = ((dft_lhs[i] + lhs_conj) * (dft_rhs[i] + rhs_conj)) *
                     Complex(0.25 / static_cast<T>(dft_size), 0.0) +
                 (dft_lhs[i] - lhs_conj) * (dft_rhs[i] - rhs_conj) * ratio;
      middle[j] = ((dft_lhs[i] - lhs_conj) * (dft_rhs[i] + rhs_conj) +
                   (dft_lhs[i] + lhs_conj) * (dft_rhs[i] - rhs_conj)) *
                  ratio;
    }
    (*this)(edges.begin(), edges.end());
    (*this)(middle.begin(), middle.end());
    std::vector<int> result(n + m - 1);
    for (int i = 0; i < n + m - 1; ++i) {
      result[i] = static_cast<int>(std::llroundl(edges[i].real()) % mod);
      if (result[i] += static_cast<int>(
              ((std::llroundl(middle[i].real()) % mod) << 15) % mod);
          result[i] >= mod) {
        result[i] -= mod;
      }
      if (result[i] += static_cast<int>(
              ((std::llroundl(edges[i].imag()) % mod) << 30) % mod);
          result[i] >= mod) {
        result[i] -= mod;
      }
    }
    return result;
  }

  template <typename U>
  std::vector<U> Convolution(
      const typename std::vector<int>::const_iterator& lhs_first,
      const typename std::vector<int>::const_iterator& lhs_last,
      const typename std::vector<int>::const_iterator& rhs_first,
      const typename std::vector<int>::const_iterator& rhs_last) {
    const auto n = static_cast<int>(std::distance(lhs_first, lhs_last));
    const auto m = static_cast<int>(std::distance(rhs_first, rhs_last));
    int dft_size = 1;
    while (dft_size < std::max(n, m) * 2) { dft_size *= 2; }
    std::vector<Complex> dft(dft_size);
    for (int i = 0; i < n; ++i) { dft[i].real(lhs_first[i]); }
    for (int i = 0; i < m; ++i) { dft[i].imag(rhs_first[i]); }
    (*this)(dft.begin(), dft.end());
    const Complex ratio(0.0, -0.25 / static_cast<T>(dft_size));
    std::vector<Complex> inv_dft(dft_size);
    for (int i = 0; i < dft_size; ++i) {
      int j = (dft_size - i) & (dft_size - 1);
      Complex conj = std::conj(dft[j]);
      inv_dft[j] = (dft[i] * dft[i] - conj * conj) * ratio;
    }
    (*this)(inv_dft.begin(), inv_dft.end());
    std::vector<U> result(n + m - 1);
    for (int i = 0; i < n + m - 1; ++i) {
      result[i] = static_cast<U>(std::llroundl(inv_dft[i].real()));
    }
    return result;
  }

  void operator()(const typename std::vector<Complex>::iterator& first,
                  const typename std::vector<Complex>::iterator& last) {
    const auto n = static_cast<int>(std::distance(first, last));
    {
      const auto shift = kLog - static_cast<int>(std::__lg(n));
      for (int i = 0; i < n; ++i) {
        if (int j = inverse_[i] >> shift; j < i) {
          std::swap(first[i], first[j]);
        }
      }
    }
    for (int i = 1; i < n; i *= 2) {
      for (int j = 0; j < n; j += i * 2) {
        for (int k = 0; k < i; ++k) {
          Complex z = root_[i + k] * first[i + j + k];
          first[i + j + k] = first[j + k] - z;
          first[j + k] = first[j + k] + z;
        }
      }
    }
  }

 private:
  static const inline T kPi = acosl(-1.0);
  std::array<Complex, 1 << kLog> root_;
  std::array<int, 1 << kLog> inverse_;
};

// https://judge.yosupo.jp/problem/convolution_mod
// https://judge.yosupo.jp/submission/179432
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int m;
  std::cin >> n >> m;

  std::vector<int> p(n);
  for (int& i : p) { std::cin >> i; }

  std::vector<int> q(m);
  for (int& i : q) { std::cin >> i; }

  FastFourierTransform<long double, 20> fft{};
  std::vector<int> C =
      fft.ConvolutionMod(p.cbegin(), p.cend(), q.cbegin(), q.cend(), 998244353);
  for (int c : C) { std::cout << c << " "; }
}

void Main() {
  int testcases = 1;
  // std::cin >> testcases;
  for (int tt = 1; tt <= testcases; ++tt) { RunCase(tt); }
}

}  // namespace

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  Main();
  return 0;
}
