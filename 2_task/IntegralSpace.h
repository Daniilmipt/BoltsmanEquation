//
// Created by daniil on 28.11.23.
//

#ifndef TASK_RELAX_INTEGRALSPACE_H
#define TASK_RELAX_INTEGRALSPACE_H

#include "../vector/Vector_3d.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Vec_v = Vector_3d<double>;
using Vec_node = Vector_3d<int>;

template<std::size_t N_vx, std::size_t N_vy, std::size_t N_vz>
struct IntegralSpace {
public:
    double d_max = 1.0;  // в площадях d^2, где d - диаметр молекулы
    double v_cut = 5.0;
    Vec_v vec_dv = calc_vec_dv();
public:

    //семинар 2
    double v_iter(int iteration, double d_v) const {
        return -v_cut + (iteration + 0.5) * d_v;
    }

    //семинар 2 (используется в 2 и 3 семинаре, поэтому здесь)
    Vector_3d<double> iter_vec_v(const Vec_node &node, double d) const {
        Vec_v v;
        v.x = v_iter(node.x, this->vec_dv.x);
        v.y = v_iter(node.y, this->vec_dv.y);
        v.z = v_iter(node.z, this->vec_dv.z);
        return v;
    }

    //семинар 1
    Vec_v calc_vec_allowed_node(double x, double y, double z) {
        Vec_v v_grid = Vec_v(x, y, z);
        Vec_v v = v_grid * 2 * v_cut - v_cut;
        return v;
    }

    //семинар 1
    double calc_dv(int n) const {
        return 2 * v_cut / n;
    }

    //семинар 1
    Vec_v calc_vec_dv() const {
        return Vec_v(calc_dv(N_vx), calc_dv(N_vy), calc_dv(N_vz));
    }

    //семинар 1
    double get_v_cut() const {
        return this->v_cut;
    }

    //семинар 1
    void set_v_cut(double v_max) {
        this->v_cut = v_max;
        this->vec_dv = calc_dv();
    }

    //семинар 1
    std::pair<Vec_v, Vec_v> calc_new_v(double theta, double e, const Vec_v& v1,
                                       const Vec_v& v2){
        Vec_v v_relative = v2 - v1;

        double v_relative_norm = v_relative.norm();
        double v_relative_xy = v_relative.x_mult_y();

        Vec_v new_v_relative;
        if (v_relative_xy == 0) {
            new_v_relative.x = v_relative_norm * std::sin(e) * std::sin(theta);
            new_v_relative.y = v_relative_norm * std::cos(e) * std::sin(theta);
            new_v_relative.z = v_relative_norm * std::cos(theta);
        } else {
            new_v_relative.x = v_relative.x * std::cos(theta)
                               - (v_relative.x*v_relative.z/v_relative_xy) * std::cos(e) * std::sin(theta)
                               + v_relative_norm * v_relative.y / v_relative_xy * std::sin(e) * std::sin(theta);

            new_v_relative.y = v_relative.y * std::cos(theta)
                               - (v_relative.y * v_relative.z / v_relative_xy) * std::cos(e) * std::sin(theta)
                               - v_relative_norm * v_relative.x / v_relative_xy * std::sin(e) * std::sin(theta);

            new_v_relative.z = v_relative.z * std::cos(theta)
                               + v_relative_xy * std::cos(e) * std::sin(theta);
        }
        //новые разлетные скорости
        Vec_v u1 = (v1 + v2) / 2 - new_v_relative / 2;
        Vec_v u2 = (v1 + v2) / 2 + new_v_relative / 2;
        return std::make_pair(u1, u2);
    }

    IntegralSpace() = default;
    IntegralSpace(double v_max) : v_cut(v_max) {}
    IntegralSpace(double v_max, double d_max) : v_cut(v_max), d_max(d_max) {}
};


#endif //TASK_RELAX_INTEGRALSPACE_H
