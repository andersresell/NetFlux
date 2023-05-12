#include "../include/Container.hpp"
namespace container {
    template<typename T, Index M, Index N, Index K>
    void matmult(const Container2D<T, M, N>& A,
                        const Container2D<T, N, K>& B,
                        Container2D<T, M, K>& C) {
        static_assert(A.rows() == C.rows() && A.columns() == B.rows(),
                      "Matrix dimensions are not compatible for multiplication.");

        // Perform matrix multiplication
        for (Index i = 0; i < A.rows(); ++i) {
            for (Index j = 0; j < B.columns(); ++j) {
                T sum = 0;
                for (Index k = 0; k < A.columns(); ++k) {
                    sum += A(i, k) * B(k, j);
                }
                C(i, j) = sum;
            }
        }
    }
}

int main() {
    using CMat = Container2D<double, 3, 5>;
    using EMat = Eigen::Matrix<double, 3, 5>;

    const Index M = 3, N = 5, K = 2;

    Container2D<double, M, N> A;
    Container2D<double, N, K> B;
    Container2D<double, M, K> C;

    container::matmult(A, B, C);

    return 0;
}