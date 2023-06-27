#pragma once
#include "../Includes.hpp"
#include "StaticContainer.hpp"

/*Dynamic container used to hold solution variables
The type below is used to store vectors of matrices (such as gradients).
They are stored in this fashion [A0, A1 ... Ai ... AN], where Ai is a MxN matrix.
M will typically represent the number of solution variables, while
It is assumed that the number of colums of the matrices N is known at compile time*/

template <typename T, ShortIndex cols_>
class DynamicContainer3D
{

protected:
    using DC3D = DynamicContainer3D;

    Index size_{0}; // Length of outer vector
    Index rows_{0}; // number of rows of each variable,

    T *data_{nullptr};

public:
    DynamicContainer3D(Index size, Index rows) : size_{size}, rows_{rows}
    {
        data_ = new T[size_ * rows_ * cols_]{0};
    }

    T &operator()(Index l, Index i, Index j)
    {
        assert(l < size_ && i < rows_ && j < cols_);
        return data_[l * rows_ * cols_ + i * cols_ + j];
    }

    const T &operator()(Index l, Index i, Index j) const
    {
        assert(l < size_ && i < rows_ && j < cols_);
        return data_[l * rows_ * cols_ + i * cols_ + j];
    }

    /*Returns pointer to matrix l*/
    T *operator[](Index l)
    {
        assert(l < size_);
        return data_ + l * rows_ * cols_;
    }
    const T *operator[](Index l) const
    {
        assert(l < size_);
        return data_ + l * rows_ * cols_;
    }

    void operator*=(T rhs)
    {
        for (Index i{0}; i < size_ * rows_ * cols_; i++)
            data_[i] *= rhs;
    }

    void operator=(T rhs)
    {
        for (Index i{0}; i < size_ * rows_ * cols_; i++)
            data_[i] = rhs;
    }

    void set_zero()
    {
        for (size_t i{0}; i < size_ * rows_ * cols_; i++)
            data_[i] = 0;
    }

    void operator=(const DC3D &other)
    {
        assert(size() == other.size());
        std::copy(other.data_, other.data_ + size_ * rows_ * cols_, data_);
    }

    template <typename StaticEigenType>
    void set_variable(Index l, const StaticEigenType &variable)
    {
        assert(StaticEigenType::RowsAtCompileTime == rows_);
        static_assert(StaticEigenType::ColsAtCompileTime == cols_);
        std::copy(variable.data(), variable.data() + StaticEigenType::SizeAtCompileTime, (*this)[l]);
    }

    template <typename StaticEigenType>
    void set_constant_field_segment(const StaticEigenType &variable, Index first_index, Index last_index)
    {
        assert(StaticEigenType::RowsAtCompileTime == rows_);
        static_assert(StaticEigenType::ColsAtCompileTime == cols_);
        assert(first_index < last_index && last_index <= size_);
        for (Index l{first_index}; l < last_index; l++)
            set_variable<StaticEigenType>(l, variable);
    }

    template <typename StaticEigenType>
    void set_constant_field(const StaticEigenType &variable)
    {
        set_constant_field_segment(variable, 0, size_);
    }

    template <typename StaticEigenType>
    Eigen::Map<StaticEigenType> get_variable(Index l)
    {
        assert(StaticEigenType::RowsAtCompileTime == rows_ && StaticEigenType::ColsAtCompileTime == cols_);
        return Eigen::Map<StaticEigenType>(data_ + l * rows_ * cols_);
    }

    template <typename StaticEigenType>
    const Eigen::Map<StaticEigenType> get_variable(Index l) const
    {
        assert(StaticEigenType::RowsAtCompileTime == rows_ && StaticEigenType::ColsAtCompileTime == cols_);
        return Eigen::Map<StaticEigenType>(data_ + l * rows_ * cols_);
    }

    /*Checks all values for nan or inf*/
    bool values_are_valid() const
    {
        for (Index i{0}; i < size_ * rows_ * cols_; i++)
            if (!num_is_valid(data_[i]))
                return false;
        return true;
    }

public:
    Index size() const
    {
        return size_;
    }
    Index rows() const { return rows_; }
    constexpr Index cols() const { return cols_; }

    ~DynamicContainer3D() { delete[] data_; }
};

template <typename T>
class DynamicContainer2D : public DynamicContainer3D<T, 1>
{

public:
    using DynamicContainer3D<T, 1>::operator=;

    DynamicContainer2D() {}

    DynamicContainer2D(Index size, Index rows) : DynamicContainer3D<T, 1>(size, rows) {}

    T &operator()(Index l, Index i)
    {
        assert(l < this->size_ && i < this->rows_);
        return this->data_[l * this->rows_ + i];
    }

    const T &operator()(Index l, Index i) const
    {
        assert(l < this->size_ && i < this->rows_);
        return this->data_[l * this->rows_ + i];
    }
};