#include <cmath>

class AAD22 {
private:
    double m_val; // value in some point
    double m_d1[2]; //first derivative in that point -  x, y
    double m_d2[3]; // differentiation for x^2, y^2, xy

public:
    AAD22() = delete;

    //from a const
    AAD22(double c) :
            m_val(0),
            m_d1{0, 0},
            m_d2{0, 0, 0} {}

private:
    constexpr AAD22(int i, double v) // i = 0 for x, i=1 for y
            : m_val(v),
              m_d1{(i == 0) ? 1.0 : 0.0, (i == 1) ? 1.0 : 0.0},
              m_d2{0, 0, 0} {}

public:
    constexpr static AAD22 X(double v) { return AAD22(0, v); }

    constexpr static AAD22 Y(double v) { return AAD22(1, v); }

    // ACCESSORS
    double getM_val() const {
        return m_val;
    }
    void setM_val(double value) {
        m_val = value;
    }
    double getM_d1(int index) const {
        if (index >= 0 && index < 2)
            return m_d1[index];
        else
            static_assert(1);
    }
    void setM_d1(int index, double value) {
        if (index >= 0 && index < 2)
            m_d1[index] = value;
        else
            static_assert(1);
    }
    double getM_d2(int index) const {
        if (index >= 0 && index < 3)
            return m_d2[index];
        else
            static_assert(1);
    }
    void setM_d2(int index, double value) {
        if (index >= 0 && index < 3)
            m_d2[index] = value;
        static_assert(1);
    }

    //op overloading
    //+, +=
    AAD22 &operator+(AAD22 const &right) {
        AAD22 res = *this;
        res.m_val += right.m_val;
        res.m_d1[0] += right.m_d1[0];
        res.m_d1[1] += right.m_d1[1];
        res.m_d2[0] += right.m_d2[0];
        res.m_d2[1] += right.m_d2[1];
        res.m_d2[2] += right.m_d2[2];
        return res;
    }

    AAD22 &operator+=(AAD22 const &right) {
        this->m_val += right.m_val;
        this->m_d1[0] += right.m_d1[0];
        this->m_d1[1] += right.m_d1[1];
        this->m_d2[0] += right.m_d2[0];
        this->m_d2[1] += right.m_d2[1];
        this->m_d2[2] += right.m_d2[2];
        return *this;
    }

    // -, -=
    AAD22 &operator-(AAD22 const &right) {
        AAD22 res = *this;
        res.m_val -= right.m_val;
        res.m_d1[0] -= right.m_d1[0];
        res.m_d1[1] -= right.m_d1[1];
        res.m_d2[0] -= right.m_d2[0];
        res.m_d2[1] -= right.m_d2[1];
        res.m_d2[2] -= right.m_d2[2];
        return res;
    }

    AAD22 &operator-=(AAD22 const &right) {
        this->m_val -= right.m_val;
        this->m_d1[0] -= right.m_d1[0];
        this->m_d1[1] -= right.m_d1[1];
        this->m_d2[0] -= right.m_d2[0];
        this->m_d2[1] -= right.m_d2[1];
        this->m_d2[2] -= right.m_d2[2];
        return *this;
    }

    //*, *=
    AAD22 &operator*(AAD22 const &right) {
        AAD22 res = *this;
        res.m_val *= right.m_val;
        res.m_d1[0] *= right.m_d1[0];
        res.m_d1[1] *= right.m_d1[1];
        res.m_d2[0] *= right.m_d2[0];
        res.m_d2[1] *= right.m_d2[1];
        res.m_d2[2] *= right.m_d2[2];
        return res;
    }

    AAD22 &operator*=(AAD22 const &right) {
        this->m_val *= right.m_val;
        this->m_d1[0] *= right.m_d1[0];
        this->m_d1[1] *= right.m_d1[1];
        this->m_d2[0] *= right.m_d2[0];
        this->m_d2[1] *= right.m_d2[1];
        this->m_d2[2] *= right.m_d2[2];
        return *this;
    }

    // /, /=
    AAD22 &operator/(AAD22 const &right) {
        AAD22 res = *this;
        res.m_val /= right.m_val;
        res.m_d1[0] /= right.m_d1[0];
        res.m_d1[1] /= right.m_d1[1];
        res.m_d2[0] /= right.m_d2[0];
        res.m_d2[1] /= right.m_d2[1];
        res.m_d2[2] /= right.m_d2[2];
        return res;
    }

    AAD22 &operator/=(AAD22 const &right) {
        this->m_val /= right.m_val;
        this->m_d1[0] /= right.m_d1[0];
        this->m_d1[1] /= right.m_d1[1];
        this->m_d2[0] /= right.m_d2[0];
        this->m_d2[1] /= right.m_d2[1];
        this->m_d2[2] /= right.m_d2[2];
        return *this;
    }

    //unary +, -
    AAD22 operator+() const {
        return *this;
    }

    AAD22 operator-() const {
        AAD22 res = *this;
        res.m_val = -res.m_val;
        res.m_d1[0] = -res.m_d1[0];
        res.m_d1[1] = -res.m_d1[1];
        res.m_d2[0] = -res.m_d2[0];
        res.m_d2[1] = -res.m_d2[1];
        res.m_d2[2] = -res.m_d2[2];
        return res;
    }

    // sin
    friend AAD22 sin(const AAD22 &x) {
        AAD22 res = x;
        res.m_val = std::sin(x.m_val);
        res.m_d1[0] = x.m_d1[0] * std::cos(x.m_val);
        res.m_d1[1] = x.m_d1[1] * std::cos(x.m_val);
        res.m_d2[0] = -std::sin(x.m_val) * x.m_d1[0] * x.m_d1[0] + std::cos(x.m_val) * x.m_d2[0];
        res.m_d2[1] = -std::sin(x.m_val) * x.m_d1[1] * x.m_d1[1] + std::cos(x.m_val) * x.m_d2[1];
        res.m_d2[2] = -std::sin(x.m_val) * x.m_d1[1] * x.m_d1[0] + std::cos(x.m_val) * x.m_d2[2];
        return res;
    }

    // cos
    friend AAD22 cos(const AAD22 &x) {
        AAD22 res = x;
        res.m_val = std::sin(x.m_val);
        res.m_d1[0] = -std::sin(x.m_val) * x.m_d1[0];
        res.m_d1[1] = -std::sin(x.m_val) * x.m_d1[1];
        res.m_d2[0] = -std::cos(x.m_val) * x.m_d1[0] * x.m_d1[0] - std::sin(x.m_val) * x.m_d2[0];
        res.m_d2[1] = -std::cos(x.m_val) * x.m_d1[1] * x.m_d1[1] - std::sin(x.m_val) * x.m_d2[1];
        res.m_d2[2] = -std::cos(x.m_val) * x.m_d1[1] * x.m_d1[0] - std::sin(x.m_val) * x.m_d2[2];
        return res;
    }

    // exp
    friend AAD22 exp(const AAD22 &x) {
        AAD22 res = x;
        res.m_val = std::exp(x.m_val);
        res.m_d1[0] = std::exp(x.m_val) * x.m_d1[0];
        res.m_d1[1] = std::exp(x.m_val) * x.m_d1[1];
        res.m_d2[0] = std::exp(x.m_val) * x.m_d1[0] * x.m_d1[0] + std::exp(x.m_val) * x.m_d2[0];
        res.m_d2[1] = std::exp(x.m_val) * x.m_d1[1] * x.m_d1[1] + std::exp(x.m_val) * x.m_d2[1];
        res.m_d2[2] = std::exp(x.m_val) * x.m_d1[1] * x.m_d1[0] + std::exp(x.m_val) * x.m_d2[2];
        return res;
    }
};
