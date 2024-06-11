#include <vector>

template<typename T>
constexpr T dividePolynomials(const std::vector<T>& a, const std::vector<T>& b, T x,
                       std::vector<T>& quotient, std::vector<T>& remainder) {
    T result = 0;
    for (auto it = a.rbegin(); it != a.rend(); ++it) {
        result = result * x + *it;
    }
    T result2 = 0;
    for (auto it = b.rbegin(); it != b.rend(); ++it) {
        result2 = result2 * x + *it;
    }
    return result/result2;
}
