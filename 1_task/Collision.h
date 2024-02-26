//
// Created by daniil on 28.11.23.
//

#ifndef TASK_RELAX_COLLISION_H
#define TASK_RELAX_COLLISION_H

#include "../vector/Vector_3d.h"
#include "../2_task/IntegralSpace.h"
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

template<std::size_t N_vx = 10, std::size_t N_vy = 10, std::size_t N_vz = 10>
class Collision : IntegralSpace<N_vx, N_vy, N_vz> {
public:
    Vec_node node_1;
    Vec_node node_2;
    double s;
    double e;
    bool isGood;

private:

    //семинар 2


    //семинар 2
    int get_center_node(double v, double dv) {
        //N = 2* v_cut / dv
        int nd_center = std::lround((v + this->get_v_cut()) / dv - 0.5);
        return nd_center;
    }

    //семинар 2
    Vec_node get_vec_center_node(const Vec_v& v) {
        Vec_node node_center;
        node_center.x = get_center_node(v.x, this->vec_dv.x);
        node_center.y = get_center_node(v.y, this->vec_dv.y);
        node_center.z = get_center_node(v.z, this->vec_dv.z);
        return node_center;
    }

    //семинар 2
    double get_relative_energy(const Vec_v& v, const Vec_v& v_centre) {
        return (v - v_centre).norm_squared();
    }

    //семинар 2
    int find_id_min_node(double v, double dv) {
        //N = 2* v_cut / dv
        int node_id = std::floor((v + this->get_v_cut()) / dv - 0.5);  // TODO: устранить повтор в find_near_i
        return node_id;
    }

    //семинар 2
    Vec_node find_min_node(const Vec_v v) {
        Vec_node min_node;
        min_node.x = find_id_min_node(v.x, this->vec_dv.x);
        min_node.y = find_id_min_node(v.y, this->vec_dv.y);
        min_node.z = find_id_min_node(v.z, this->vec_dv.z);
        return min_node;
    }

    //семинар 2
    //конце второго семинара, находим минимальный по расстоянию узел
    Vec_node find_second_node(Vec_v v1, Vec_v v_mean, const Vec_node& sum_node,
                                                           double energy_0, bool is_energy0_bigger) {
        Vec_node center = find_min_node(v1);  // center index
        Vec_node best_ind;
        double min = 1'000'000'000;

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    double vec_cut = this->get_v_cut();

                    Vec_node cur_ind = center + Vec_node(i, j, k);
                    Vec_v cur_vec = this->iter_vec_v(cur_ind, 0);

                    // после разлета узел
                    // суммарный вектор после разлеты должен быть тем же,
                    // поэтому так считаем paired_ind
                    Vec_node paired_ind = sum_node - cur_ind;
                    Vec_v paired_vec = this->iter_vec_v(paired_ind, 0);

                    if ((cur_vec.norm_squared() < vec_cut*vec_cut) &&
                        (paired_vec.norm_squared() < vec_cut*vec_cut)) {
                            double cur_energy = get_relative_energy(cur_vec, v_mean);

                            if (((energy_0 > cur_energy) == is_energy0_bigger) ||
                                (energy_0 == cur_energy)) {
                                    double distance = (v1 - cur_vec).norm_squared();
                                    if (distance < min) { // хотим найти min
                                        min = distance;
                                        best_ind = cur_ind;
                                    }
                            }
                    }
                }
            }
        }
        if (min == 1'000'000'000) {
            this->isGood = false;
        }
        return best_ind;
    }

    //семинар 2
    double find_r(double energy_0,
                                               Vec_node lambda_mu,
                                               Vec_node lambda_mu_new,
                                               Vec_v v_mean) {

        double energy_1 = get_relative_energy(this->iter_vec_v(lambda_mu, 0), v_mean);
        double energy_2 = get_relative_energy(this->iter_vec_v(lambda_mu_new, 0), v_mean);
        return (energy_0 - energy_1) / (energy_2 - energy_1);
    }

public:
    Collision() {}

    Collision(std::vector<double>& node_vector,
              double v_max,
              double d_max)
            : IntegralSpace<N_vx, N_vy, N_vz>(v_max, d_max) {

        node_1 = get_vec_center_node(
                this->calc_vec_allowed_node(
                        node_vector[0],
                        node_vector[1],
                        node_vector[2]
                )
        );
        node_2 = get_vec_center_node(
                this->calc_vec_allowed_node(
                        node_vector[3],
                        node_vector[4],
                        node_vector[5]
                )
        );

        if (this->iter_vec_v(node_1, 0).norm() < this->v_cut &&
            this->iter_vec_v(node_2, 0).norm() < this->v_cut) {
            isGood = true;
        } else {
            isGood = false;
        }
        s = node_vector[6] * this->d_max * this->d_max;
        e = node_vector[7] * 2 * M_PI;
    }

    //семинар 1-2
    std::tuple<Vec_node, Vec_node, Vec_node, Vec_node, double, int> get_new_v() {
        Vec_node lambda, lambda_new, mu, mu_new;
        double r;
        int count_ignore = 0;
        if (isGood) {
            // Находим скорости после столкновения
            double theta = 2.0 * std::acos(std::sqrt(s));
            Vec_v v1 = this->iter_vec_v(node_1, 0);
            Vec_v v2 = this->iter_vec_v(node_2, 0);

            //новые разлетные скорости
            std::pair<Vec_v, Vec_v> v_pair = this->calc_new_v(theta, e, v1, v2);
            Vec_v u1, u2;
            std::tie(u1, u2) = v_pair;

            // Находим аппроксимипующие скорости (индексы) и r
            std::tie(lambda, mu, lambda_new, mu_new, r, count_ignore) =
                    find_second_v_nodes(u1, u2, v1, v2, get_relative_energy(v1, (v1+v2)/2));
//            // TODO: cтереть после проверки
//             if (isGood) {
//                 double energy_0 = v1.norm_squared() + v2.norm_squared();
//                 double energy_1 = this->iter_vec_v(lambda,0).norm_squared() + this->iter_vec_v(mu,0).norm_squared();
//                 double energy_2 = this->iter_vec_v(lambda_new,0).norm_squared() + this->iter_vec_v(mu_new,0).norm_squared();
//                 double r_add = (energy_0 - energy_1) / (energy_2 - energy_1);
//                 Vec_v v_lambda = this->iter_vec_v(lambda,0);
//                 Vec_v v_mu = this->iter_vec_v(mu,0);
//                 Vec_v v_lambda_new = this->iter_vec_v(lambda_new,0);
//                 Vec_v v_mu_new = this->iter_vec_v(mu_new,0);
//                 std::cout << energy_0 << " " << energy_1 << " " << energy_2 << '\n';
//                 std::cout << v_lambda + v_mu - v1 - v2 << v_lambda_new + v_mu_new - v1 - v2<< '\n';
//                 std::cout << r << " " << r_add << "\n\n";
//             }
        }
        return std::make_tuple(lambda, mu, lambda_new, mu_new, r, count_ignore);
    }

    //семинар 2, находим вторую пару скорости и узлов
    std::tuple<Vec_node, Vec_node, Vec_node, Vec_node, double, int> find_second_v_nodes(Vec_v u1, Vec_v u2, Vec_v v1, Vec_v v2, double energy_0){
        Vec_node lambda, lambda_new, mu, mu_new;
        Vec_v v_mean_old = (v1 + v2) / 2;
        double r;
        // Находим аппроксимипующие скорости, узлы и r
        Vec_node node_center_1 = get_vec_center_node(u1);
        Vec_node node_center_2 = get_vec_center_node(u2);
        Vec_v u_center_1 = this->iter_vec_v(node_center_1, 0);
        Vec_v u_center_2 = this->iter_vec_v(node_center_2, 0);

        int count_ignore = 0;

        if ((u_center_1.norm() >= this->v_cut) || (u_center_2.norm() >= this->v_cut)) {
            isGood = false;
            ++count_ignore;
        } else {
            double energy_center = get_relative_energy(u_center_1, v_mean_old);
            if (energy_0 == energy_center) {
                lambda = node_center_1;
                lambda_new = node_center_1;
                mu = node_center_2;
                mu_new = node_center_2;
                r = 1;
            } else {
                if (energy_0 < energy_center) {
                    lambda_new = node_center_1;
                    mu_new = node_center_2;
                    lambda = this->find_second_node(
                            u1,
                            v_mean_old,
                            node_center_1 + node_center_2,
                            energy_0,
                            true
                    );
                    mu = lambda_new + mu_new - lambda;
                } else {
                    lambda = node_center_1;
                    mu = node_center_2;
                    lambda_new = this->find_second_node(
                            u1,
                            v_mean_old,
                            node_center_1 + node_center_2,
                            energy_0,
                            false);
                    mu_new = lambda + mu - lambda_new;
                }
                if (isGood) {
                    r = find_r(energy_0, lambda, lambda_new, v_mean_old);
                }
            }
        }
        return std::make_tuple(lambda, mu, lambda_new, mu_new, r, count_ignore);
    }

    bool operator == (const Collision& collision) const {
        if (this->node_2 == collision.node_2)
            return true;
        return false;
    }

    friend std::ostream& operator << (std::ostream& out, const Collision<N_vx, N_vy, N_vz>& collision) {
        out << collision->calc_v(collision.node_1) << "\t"
            << collision->calc_v(collision.node_2) << "\t"
            << collision.s << "\t"
            << collision.e;
        return out;
    }
};

#endif //TASK_RELAX_COLLISION_H
