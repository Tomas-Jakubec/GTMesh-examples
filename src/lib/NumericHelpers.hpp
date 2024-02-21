#pragma once

#include <GTMesh/NumericStaticArray/Vector.h>

template <unsigned int Dim, unsigned int Index>
constexpr Vertex<Dim> e() {
    Vertex<Dim> basis{};
    std::get<Index>(basis) = 1.0;
    return basis;
}


template <unsigned int Dim>
constexpr void unitMatrix(Vector<Dim,Vector<Dim>>& matrix, std::integral_constant<unsigned int, 0>) {
    std::get<0>(matrix) = e<Dim, 0>();
}

template <unsigned int Dim, unsigned int Index>
constexpr void unitMatrix(Vector<Dim,Vector<Dim>>& matrix, std::integral_constant<unsigned int, Index>) {
    unitMatrix(matrix, std::integral_constant<unsigned int, Index -1>());
    std::get<Index>(matrix) = e<Dim, Index>();
}


template <unsigned int Dim>
constexpr Vector<Dim,Vector<Dim>> unitMatrix() {
    Vector<Dim,Vector<Dim>> matrix;
    unitMatrix(matrix, std::integral_constant<unsigned int, Dim - 1>());
    return matrix;
}


template<unsigned Dimension>
Vector<Dimension> multiply(const Vector<Dimension,Vector<Dimension>>& A, const Vector<Dimension>& x){
    Vector<Dimension> result{};
    for(unsigned int i = 0; i < Dimension; ++i) {
        for(unsigned int j = 0; j < Dimension; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

template<unsigned Dimension>
double multiply(const Vector<Dimension>& x, const Vector<Dimension>& y){
    double result{};
    for(unsigned int i = 0; i < Dimension; ++i) {
        result += x[i] * y[i];
    }
    return result;
}

template<unsigned Dimension>
Vector<Dimension,Vector<Dimension>> tenzorMultiply(const Vector<Dimension>& x, const Vector<Dimension>& y_t) {
    Vector<Dimension,Vector<Dimension>> result{{},{}};
    for(unsigned int i = 0; i < Dimension; ++i) {
        for(unsigned int j = 0; j < Dimension; ++j) {
            result[i][j] += x[i] * y_t[j];
        }
    }
    return result;
}

template<unsigned Dimension>
Vector<Dimension,Vector<Dimension>> operator*(double alpha, const Vector<Dimension,Vector<Dimension>>& A) {
    auto result = A;
    for(unsigned int i = 0; i < Dimension; ++i) {
        for(unsigned int j = 0; j < Dimension; ++j) {
            result[i][j] *= alpha;
        }
    }
    return result;
}

template <unsigned int Dim>
Vector<Dim,Vector<Dim>> calculateInvMatrix(const Vector<Dim,Vector<Dim>>& inA) {
    // Create augmented matrix
    auto A = inA;
    auto B = unitMatrix<Dim>();

    // Perform Gaussian elimination to transform augmented matrix into upper triangular form
    for (unsigned int i = 0; i < Dim; ++i) {
        // Find pivot row
        unsigned int pivot = i;
        for (unsigned int j = i + 1; j < Dim; ++j) {
            if (abs(A[j][i]) > abs(A[pivot][i])) {
                pivot = j;
            }
        }
        // Swap rows if necessary
        if (pivot != i) {
            swap(A[i], A[pivot]);
            swap(B[i], B[pivot]);
        }
        // Eliminate lower triangular elements
        for (unsigned int j = i + 1; j < Dim; ++j) {
            double factor = A[j][i] / A[i][i];
            for (unsigned int k = i; k < Dim; ++k) {
                A[j][k] -= factor * A[i][k];
                B[j][k] -= factor * B[i][k];
            }
        }
    }

    // Perform back substitution to transform upper triangular matrix into diagonal matrix
    for (unsigned int i = Dim - 1; i >= 0; --i) {
        for (unsigned int j = i - 1; j >= 0; --j) {
            double factor = B[j][i] / B[i][i];
            for (unsigned int k = i; k < Dim; ++k) {
                A[j][k] -= factor * A[i][k];
                B[j][k] -= factor * B[i][k];
            }
        }
    }

    // Perform Gaussian elimination to transform diagonal matrix into identity matrix
    for (unsigned int i = 0; i < Dim; ++i) {
        double factor = 1 / A[i][i];
        for (unsigned int j = i; j < Dim; ++j) {
            A[i][j] *= factor;
            B[i][j] *= factor;
        }
    }

    return B;
}