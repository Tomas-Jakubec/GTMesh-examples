#pragma once

#include <lib/NumericHelpers.hpp>

template<typename F>
class Gradient {

    F f;
    double delta;
    double invDelta;

    template <unsigned int Index, unsigned int Dim>
    double calculateGradientComponent(const Vertex<Dim>& x) {
        Vector<Dim> deltaVec = e<Dim, Index>() * delta;
        return (f(x + deltaVec) - f(x - deltaVec)) * invDelta;
    }

    template <unsigned int Dim>
    constexpr void calculateGradient(Vector<Dim>& result, const Vertex<Dim>& x, std::integral_constant<unsigned int, 0>) {
        result[0] = calculateGradientComponent<0>(x);
    }

    template <unsigned int Dim, unsigned int Index = 0>
    constexpr void calculateGradient(Vector<Dim>& result, const Vertex<Dim>& x, std::integral_constant<unsigned int, Index>) {
        result[Index] = calculateGradientComponent<Index>(x);
        calculateGradient(result, x, std::integral_constant<unsigned int, Index - 1>());
    }
public:

    Gradient(F f, double delta): f(f), delta(delta), invDelta(0.5 / delta) {}

    template <unsigned int Dim>
    Vector<Dim> operator ()(const Vertex<Dim>& x) {
        Vector<Dim> gradient;

        calculateGradient(gradient, x, std::integral_constant<unsigned int, Dim - 1>());

        return gradient;
    }
};

template <typename F>
Gradient<F> gradientOf(F&& f, double delta) {
    return Gradient<F>(f, delta);
}
