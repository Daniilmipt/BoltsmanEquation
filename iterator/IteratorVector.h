//
// Created by daniil on 30.11.23.
//

#ifndef TASK_RELAX_ITERATORVECTOR_H
#define TASK_RELAX_ITERATORVECTOR_H

#include <cstddef>
#include "../vector/Vector_3d.h"

// Итератор для трехмерного массива. *it = {i_x, i_y, i_z}. Чтобы получить само значение используйте Сontainer(*it) или Container[*it] 
template<std::size_t N_x, std::size_t N_y, std::size_t N_z, typename ValueType>
class IteratorVector : public std::iterator<std::input_iterator_tag, ValueType> {
private:
    ValueType ind;

public:
    IteratorVector(const ValueType& other_ind) : ind(other_ind) {} // TODO: лучше это отпрваить в private,наваерное
    IteratorVector(const IteratorVector& it) : ind(it.ind) {}
    bool operator== (const IteratorVector& other_it) const {
        return this->ind == other_it.ind;
    }
    bool operator!= (const IteratorVector& other_it) const {
        return this->ind != other_it.ind;
    }
    ValueType& operator* () {
        return ind;
    }
    IteratorVector& operator++ () {  // end = {Nx, 0, 0}
        if (ind == ValueType(N_x, 0, 0))
            return *this;

        ind.z++;
        if (ind.z == N_z) {
            ind.z = 0;
            ind.y++;
        }
        if (ind.y == N_y) {
            ind.y = 0;
            ind.x++;
        }
        return *this;
    }
};

#endif //TASK_RELAX_ITERATORVECTOR_H
