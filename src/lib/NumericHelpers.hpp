#pragma once

#include <GTMesh/NumericStaticArray/Vector.h>

template <unsigned int BasisIndex, unsigned int Index>
constexpr double basisComponent() {
    return (Index == BasisIndex) ? 1.0 : 0.0;
}

template <unsigned int Dim, unsigned int BasisIndex>
constexpr void fillBasis(Vertex<Dim>& basis, std::integral_constant<unsigned int, 0>) {
    basis[0] = basisComponent<BasisIndex, 0>();
}

template <unsigned int Dim, unsigned int BasisIndex, unsigned int Index>
constexpr void fillBasis(Vertex<Dim>& basis, std::integral_constant<unsigned int, Index>) {
    fillBasis<Dim, BasisIndex>(basis, std::integral_constant<unsigned int, Index - 1>());
    basis[Index] = basisComponent<BasisIndex, Index>();
}

template <unsigned int Dim, unsigned int Index>
constexpr Vertex<Dim> e() {
    Vertex<Dim> basis {};
    fillBasis<Dim, Index>(basis, std::integral_constant<unsigned int, Dim - 1>());
    return basis;
}


template <unsigned int Dim>
constexpr void unitMatrix(Vector<Dim,Vector<Dim>>& matrix, std::integral_constant<unsigned int, 0>) {
    matrix[0] = e<Dim, 0>();
}

template <unsigned int Dim, unsigned int Index>
constexpr void unitMatrix(Vector<Dim,Vector<Dim>>& matrix, std::integral_constant<unsigned int, Index>) {
    unitMatrix(matrix, std::integral_constant<unsigned int, Index -1>());
    matrix[Index] = e<Dim, Index>();
}


template <unsigned int Dim>
constexpr Vector<Dim,Vector<Dim>> unitMatrix() {
    Vector<Dim,Vector<Dim>> matrix;
    unitMatrix(matrix, std::integral_constant<unsigned int, Dim - 1>());
    return matrix;
}



Vector<2> multiply(const Vector<2,Vector<2>>& A, const Vector<2>& x){
    Vector<2> result{};
    for(unsigned int i = 0; i < 2; ++i) {
        for(unsigned int j = 0; j < 2; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

double multiply(const Vector<2>& x, const Vector<2>& y){
    double result{};
    for(unsigned int i = 0; i < 2; ++i) {
        result += x[i] * y[i];
    }
    return result;
}

Vector<1> multiply(const Vector<1,Vector<1>>& A, const Vector<1>& x){
    return Vector<1>{A[0][0] * x[0]};
}

double multiply(const Vector<1>& x, const Vector<1>& y){
    return x[0] * y[0];
}

Vector<2,Vector<2>> tenzorMultiply(const Vector<2>& x, const Vector<2>& y_t) {
    Vector<2,Vector<2>> result{{},{}};
    for(unsigned int i = 0; i < 2; ++i) {
        for(unsigned int j = 0; j < 2; ++j) {
            result[i][j] += x[i] * y_t[j];
        }
    }
    return result;
}

Vector<1,Vector<1>> tenzorMultiply(const Vector<1>& x, const Vector<1>& y_t) {
    return {{x[0] * y_t[0]}};
}

Vector<2,Vector<2>> operator*(double alpha, const Vector<2,Vector<2>>& A) {
    auto result = A;
    for(unsigned int i = 0; i < 2; ++i) {
        for(unsigned int j = 0; j < 2; ++j) {
            result[i][j] *= alpha;
        }
    }
    return result;
}

Vector<1,Vector<1>> operator*(double alpha, const Vector<1,Vector<1>>& A) {
    return Vector<1, Vector<1>>{{A[0][0] * alpha}};
}
