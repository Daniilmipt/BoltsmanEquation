//
// Created by daniil on 21.01.24.
//

#ifndef TASK_RELAX_CONDUCTIVITY_TASK_H
#define TASK_RELAX_CONDUCTIVITY_TASK_H

#include <cstddef>
#include <complex>
#include <vector>
#include "vector/Vector_3d.h"

using Vec_v = Vector_3d<double>;
using Vec_node = Vector_3d<int>;

struct Conduct{
    Vec_node N_v = Vec_node(5, 5, 5);
    Vec_node N_c = Vec_node(3, 3, 3);
    std::size_t L = 50;

    double h = (double) L/N_c.x;
    Vec_v v_cut = {5, 5, 5};
    Vec_v dv= Vec_v(2*N_v.x/v_cut.x, N_v.y/v_cut.y, N_v.z/v_cut.z);

    uint N0 = 0;
    double tau = 2*h/v_cut.x;

    double T1 = 1;
    double T2 = 2;

    double T0 = 1;
    double n0 = 1/T0;

    //N_c.x + 2 на N0
    double** f;

    ~Conduct(){
        for (int i = 0; i < N_c.x + 2; ++i) {
            delete[] f[i];
        }
        delete[] f;
    }

    //task1
    Vec_v speed_greed_item(Vec_node index){
        Vec_v v;
        v.x = -v_cut.x + (index.x-0.5)*dv.x;
        v.y = (index.y-0.5)*dv.y;
        v.z = (index.z-0.5)*dv.z;
        return v;
    }

    //task1
    double coord_speed_grid(int index){
        return h*(index-0.5);
    }

    //task1
    void find_N0(){
        for (int i = 0; i < N_v.x; ++i) {
            for (int j = 0; j < N_v.y; ++j) {
                for (int k = 0; k < N_v.z; ++k) {
                    Vec_v speed = speed_greed_item(Vec_node(i, j, k));
                    if (speed.norm_squared() <= v_cut.x*v_cut.x){
                        ++N0;
                    }
                }
            }
        }

        f = new double*[N_c.x + 2];
        for (int i = 0; i < N_c.x + 2; ++i) {
            f[i] = new double[N0];
        }
    }

    //task2
    void init_f(){
        double c = 0;
        for (int i = 1; i <= N_c.x; ++i) {
            int p = 0;
            for (int j = 0; j < N_v.x; ++j) {
                for (int k = 0; k < N_v.y; ++k) {
                    for (int l = 0; l < N_v.z; ++l) {
                        Vec_v speed = speed_greed_item(Vec_node(j, k, l));
                        if (speed.norm_squared() <= v_cut.x*v_cut.x){
                            ++p;
                            double t = T0 + coord_speed_grid(i) * (T2 - T1) / L;
                            double n = 1 / t;

                            f[i][p] = n*std::exp(-0.5*speed.norm_squared()/t);
                            c += std::exp(-0.5*speed.norm_squared()/t);
                        }
                    }
                }
            }
        }
        c = 1/(4*c*dv.x*dv.y*dv.z);

        for (int i = 1; i <= N_c.x; ++i) {
            for (int j = 0; j < N0; ++j) {
                f[i][j] *= c;
            }
        }

        for (int j = 0; j < N0; ++j) {
            f[0][j] = std::max(0.0, 2 * f[1][j] - f[2][j]);
            f[N_c.x + 1][j] = std::max(0.0, 2 * f[N_c.x][j] - f[N_c.x - 1][j]);
        }

        for (int i = 0; i <= N_c.x + 1; ++i) {
            for (int p = 0; p < N0; ++p) {
                std::cout << f[i][p] << "\t";
            }
            std::cout << "\n";
        }
    }

    //task3
    double phy_func(double val){
        return std::max(0.0, std::min(2.0, std::min(2*val, (1+val)/2)));
    }

    //task3
    int signum(double val){
        return val > 0 ? 1 : -1;
    }

    //task3
    double get_new_f(double f_k, double f_i1, double f_i, Vec_v speed, double k_big, double theta){
        return f_k + 1/2 * signum(speed.x) * (1-k_big) * phy_func(theta) * (f_i1-f_i);
    }

    //task3
    //N_c.x на N0
    std::vector<std::vector<double>> find_theta(){
        std::vector<std::vector<double>> f_new(N_c.x);
        std::vector<std::vector<double>> f_evol(N_c.x);

        double tau_0 = tau / 5;

        for (int i = 1; i <= N_c.x; ++i) {
            int p = 0;
            std::vector<double> sub_vector(N0);
            std::vector<double> sub_vector_evol(N0);

            f_new[i-1] = sub_vector;
            f_evol[i-1] = sub_vector_evol;
            for (int j = 0; j < N_v.x; ++j) {
                for (int k = 0; k < N_v.y; ++k) {
                    for (int l = 0; l < N_v.z; ++l) {
                        Vec_v speed = speed_greed_item(Vec_node(j, k, l));
                        if (speed.norm() <= v_cut.x) {
                            double theta;
                            int k_value;
                            if (speed.x > 0) {
                                theta = (f[i][p] - f[i - 1][p]) / (f[i + 1][p] - f[i][p]);
                                k_value = i;
                            } else {
                                //theta = (f[i + 2][p] - f[i + 1][p]) / (f[i + 1][p] - f[i][p]);
                                theta = (f[i][p] - f[i - 1][p]) / (f[i + 1][p] - f[i][p]);
                                k_value = i + 1;
                            }
                            double k_big = speed.x * tau_0 / h;

                            double f_new_pos = get_new_f(f[k_value][p], f[i + 1][p], f[i][p], speed, k_big, theta);
                            f_new[i-1][p] = f_new_pos;

                            if (i > 1) {
                                double result = -speed.x * tau_0 / h * (f_new[i-1][p] - f_new[i - 2][p]) + f[i - 1][p];
                                f_evol[i-1][p] = result;
                            }
                            ++p;
                        }
                    }
                }
            }
        }
        return f_evol;
    }

    //task 7
    void print_temperature(std::vector<std::vector<double>> f_evol){
        std::cout << "\n\n";
        for (int i = 0; i < N_c.x; ++i) {
            int p = 0;

            double energy = 0;
            double n = 0;
            Vec_v v {0, 0, 0};
            for (int j = 0; j < N_v.x; ++j) {
                for (int k = 0; k < N_v.y; ++k) {
                    for (int l = 0; l < N_v.z; ++l) {
                        Vec_v speed = speed_greed_item(Vec_node(j, k, l));

                        if (speed.norm() <= v_cut.x) {
                            energy += speed.norm_squared() * f_evol[i][p] * dv.x_mult_y() * dv.z;
                            n += f_evol[i][p] * dv.x_mult_y() * dv.z;
                            v += speed.norm() * f_evol[i][p] * dv.x_mult_y() * dv.z;

                            ++p;
                        }
                    }
                }
            }

            v = v / n;
            auto val = energy/n - v.norm_squared();
            double temp = (energy/n - v.norm_squared()) / 3;
            std::cout << i << "\t" << temp << "\n";
        }
    }
};


#endif //TASK_RELAX_CONDUCTIVITY_TASK_H
