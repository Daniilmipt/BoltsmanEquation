//
// Created by daniil on 28.11.23.
//

#include <complex>
#include "Vector_3d.h"


template <typename T>
Vector_3d<T>::Vector_3d(): x(), y(), z() {}

template <typename T>
Vector_3d<T>::Vector_3d(const T& a): x(a), y(a), z(a) {}

template <typename T>
Vector_3d<T>::Vector_3d(const T& x0, const T& y0, const T& z0): x(x0), y(y0), z(z0) {}

template <typename T>
Vector_3d<T>::Vector_3d(const Vector_3d& vector) : x(vector.x), y(vector.y), z(vector.z) {}

template <typename T>
Vector_3d<T>::Vector_3d(Vector_3d&& vector) : x(std::move(vector.x)), y(std::move(vector.y)), z(std::move(vector.z)) {}

template <typename T>
T Vector_3d<T>::norm_squared() const {
    return x * x + y * y + z * z;
}

template <typename T>
T Vector_3d<T>::norm() const {
    return std::sqrt(this->norm_squared());
}

template <typename T>
T Vector_3d<T>::x_mult_y() const {
    return std::sqrt(x * x + y * y);
}


template <typename T>
Vector_3d<T>& Vector_3d<T>::operator = (const Vector_3d& vector) {
    if (this != &vector) {
        x = vector.x;
        y = vector.y;
        z = vector.z;
    }
    return *this;
}


template <typename T>
Vector_3d<T>& Vector_3d<T>::operator = (Vector_3d<T>&& vector) {
    if (this != &vector) {
        x = std::move(vector.x);
        y = std::move(vector.y);
        z = std::move(vector.z);
    }
    return *this;
}


template <typename T>
bool Vector_3d<T>::operator == (const Vector_3d<T>& vector) const {
    return (x == vector.x) && (y == vector.y) && (z == vector.z);
}


template <typename T>
bool Vector_3d<T>::operator != (const Vector_3d<T>& vector) const {
    return !(*this == vector);
}


template <typename T>
Vector_3d<T> Vector_3d<T>::operator + (const Vector_3d<T>& vector) const {
    Vector_3d res(x + vector.x, y + vector.y, z + vector.z);
    return res;
}


template <typename T>
Vector_3d<T>& Vector_3d<T>::operator += (const Vector_3d<T>& vector) {
    x += vector.x;
    y += vector.y;
    z += vector.z;
    return *this;
}


template <typename T>
Vector_3d<T> Vector_3d<T>::operator - (const Vector_3d<T>& vector) const {
    Vector_3d res(x - vector.x, y - vector.y, z - vector.z);
    return res;
}


template <typename T>
Vector_3d<T>& Vector_3d<T>::operator -= (const Vector_3d<T>& vector) {
    x -= vector.x;
    y -= vector.y;
    z -= vector.z;
    return *this;
}


template <typename T>
Vector_3d<T> Vector_3d<T>::operator * (const T& a) const {
    Vector_3d res(x * a, y * a, z * a);
    return res;
}


template <typename T>
Vector_3d<T> Vector_3d<T>::operator / (const T& a) const {
    Vector_3d res(x / a, y / a, z / a);
    return res;
}


template <typename T>
T & Vector_3d<T>::operator[] (int i) {  // TODO: может и нужна эта функия???
    switch (i)
    {
        case 0:
            return this->x;
        case 1:
            return this->y;
        case 2:
            return this->z;
        default:
            throw "Неправильный индекс. Размер вектора равен 3";
    }
}


template <typename T>
T Vector_3d<T>::operator[] (int i) const {
    switch (i)
    {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw "Неправильный индекс. Размер вектора равен 3";
    }
}
