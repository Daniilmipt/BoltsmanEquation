//
// Created by daniil on 28.11.23.
//

#ifndef TASK_RELAX_VECTOR_3D_H
#define TASK_RELAX_VECTOR_3D_H


#include <iostream>

template <typename T>
struct Vector_3d {
    T x, y, z;

    Vector_3d();
    Vector_3d(const T& a);
    Vector_3d(const T& x0, const T& y0, const T& z0);

    Vector_3d(const Vector_3d& vector);
    Vector_3d(Vector_3d&& vector);

    T norm_squared() const;
    T norm() const;
    T x_mult_y() const;

    Vector_3d& operator = (const Vector_3d& vector);

    Vector_3d& operator = (Vector_3d&& vector);

    bool operator == (const Vector_3d& vector) const;

    bool operator != (const Vector_3d& vector) const;

    Vector_3d operator + (const Vector_3d& vector) const;

    Vector_3d& operator += (const Vector_3d& vector);

    Vector_3d operator - (const Vector_3d& vector) const;

    Vector_3d& operator -= (const Vector_3d& vector);

    Vector_3d operator * (const T& a) const;

    Vector_3d operator / (const T& a) const;

    T & operator[] (int i);
    T operator[] (int i) const;

    friend std::ostream& operator<< (std::ostream& out, const Vector_3d& vec) {
        out << "{" << vec.x << " " << vec.y << " " << vec.z << "}";
        return out;
    }
};
#endif //TASK_RELAX_VECTOR_3D_H
