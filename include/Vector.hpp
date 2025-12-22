#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <array>
#include <cmath>
#include <iostream>

template <size_t Dim>
struct Vector {
    std::array<double, Dim> data;

    Vector() {
        data.fill(0.0);
    }

    // Variadic constructor to allow Vector<3> v(1.0, 2.0, 3.0);
    template<typename... Args>
    Vector(Args... args) : data{static_cast<double>(args)...} {}

    double& operator[](size_t i) { return data[i]; }
    const double& operator[](size_t i) const { return data[i]; }

    // Vector addition
    Vector operator+(const Vector& other) const {
        Vector result;
        for (size_t i = 0; i < Dim; ++i) result[i] = data[i] + other[i];
        return result;
    }

    // Vector subtraction
    Vector operator-(const Vector& other) const {
        Vector result;
        for (size_t i = 0; i < Dim; ++i) result[i] = data[i] - other[i];
        return result;
    }

    // Scalar multiplication
    Vector operator*(double scalar) const {
        Vector result;
        for (size_t i = 0; i < Dim; ++i) result[i] = data[i] * scalar;
        return result;
    }

    // Dot product
    double dot(const Vector& other) const {
        double sum = 0.0;
        for (size_t i = 0; i < Dim; ++i) sum += data[i] * other[i];
        return sum;
    }

    // Squared norm
    double norm2() const {
        return dot(*this);
    }

    // Norm
    double norm() const {
        return std::sqrt(norm2());
    }
};

// Global scalar multiplication (scalar * vector)
template <size_t Dim>
Vector<Dim> operator*(double scalar, const Vector<Dim>& v) {
    return v * scalar;
}

#endif // VECTOR_HPP
