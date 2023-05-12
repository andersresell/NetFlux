#pragma once

#include "../includes.hpp"



/*Generic matrix type used directly or for building other data types*/
template<typename T, Index M, Index N>
class StaticContainer2D{
    using C2D = StaticContainer2D<T,M,N>;
    static constexpr size_t SIZE = M * N;
protected:
    T data[M * N];
public:
    StaticContainer2D() : data{} {}

    template<typename... Args>
    StaticContainer2D(Args... args) : data{static_cast<T>(args)...} {
        static_assert(sizeof...(args) == SIZE);
    }
    
    /*Stored in row major order*/
    T& operator()(Index i, Index j) {
        assert(i<M && j<N);
        assert(M > 1); //Use the [] for 1D arrays
        return data[i*N+j];
    }
    /*Stored in row major order*/
    const T& operator()(Index i, Index j) const{
        assert(i<M && j<N);
        return data[i*N+j];
    }

    T& operator[](Index i) {
        assert(i<SIZE);
        return data[i * N];
    }

    const T& operator[](Index i) const {
        assert(i<SIZE);
        return data[i * N];
    }
    
    C2D operator+(C2D rhs) const {
        for (Index i{0}; i<SIZE; i++) rhs.data[i] += data[i];
        return rhs;
    };

    C2D operator-(const C2D& rhs) const {
        C2D lhs{*this};
        for (Index i{0}; i<SIZE; i++) lhs.data[i] -= rhs.data[i];
        return lhs;
    };

    void operator-=(const C2D& rhs) {
        for (Index i{0}; i<SIZE; i++) data[i] -= rhs.data[i];
    };

    void operator+=(const C2D& rhs) {
        for (Index i{0}; i<SIZE; i++) data[i] += rhs.data[i];
    };

    void operator*=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] *= rhs;
    }

    void operator=(T rhs){
        for (Index i{0}; i<SIZE; i++) data[i] = rhs;
    }

    static constexpr Index size() {return SIZE;}
    static constexpr Index rows() {return M;}
    static constexpr Index cols() {return N;}
    
    friend std::ostream& operator<<(std::ostream& os, const C2D& rhs){
        for (Index i{0}; i<M; i++){ 
            for (Index j{0}; j<N; j++) 
                os << rhs(i,j) << " ";
            os << "\n";
        }
        return os << "\n";
    }
};

template<typename T, Index M>
class StaticContainer1D : public StaticContainer2D<T, M, 1>{
    
public:
    template<typename... Args>
    StaticContainer1D(Args&&... args) : StaticContainer2D<T,M,1>(std::forward<Args>(args)...) {}
};




namespace container {
    template<typename T, Index M, Index N, Index K>
    inline void matmult(const StaticContainer2D<T, M, N>& A,
                        const StaticContainer2D<T, N, K>& B,
                        StaticContainer2D<T, M, K>& C) {
        static_assert(A.rows() == C.rows() && A.cols() == B.rows(),
                      "Matrix dimensions are not compatible for multiplication.");

        for (Index i = 0; i < M; ++i) {
            for (Index j = 0; j < K; ++j) {
                T sum = 0;
                for (Index k = 0; k < N; ++k) {
                    sum += A(i, k) * B(k, j);
                }
                C(i, j) = sum;
            }
        }
    }

    template<typename T, Index M, Index N, Index K>
    inline void Vec3_mult(const StaticContainer2D<T, M, N>& A,
                        const Vec3& x,
                        StaticContainer2D<T, N_DIM, K>& b) {
        static_assert(A.rows() == b.rows(), "Different dimensions between A and b\n");
        static_assert(N==1 || K==1 && !(K==1 && N==1)); //Either N or K should be 1, but not both

        if constexpr (K == 1){
            for (Index i = 0; i < M; ++i) {
                b[i] = 0;
                for (Index j = 0; j < N_DIM; ++j) 
                    b[i] += A(i, j) * x[j];
            }
        }
        else{ //N == 1
            for (Index i = 0; i < M; ++i) {    
                for (Index j = 0; j < N_DIM; ++j) 
                    b(i,j) = A[i] * x[j];
            }
        }
    }



    
}


