//
// Created by daniil on 28.11.23.
//

#ifndef TASK_RELAX_DISTRIBUTION_H
#define TASK_RELAX_DISTRIBUTION_H

#include "vector/Vector_3d.cpp"
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
const static double NORM_COEF = 1 / (2 * M_PI) / sqrt(2 * M_PI);

struct DistributionFunc {
    double temp = 1.0;
    DistributionFunc(double temperature) : temp(temperature) {}

    virtual double get (const Vec_v&) = 0;
};

struct Maxwell : DistributionFunc {
    Vec_v u;
    double temperature = 1;
    Maxwell(double temperature, Vec_v shift):
                        DistributionFunc(temperature),
                        u(std::move(shift)){}
    
    double get (const Vec_v& v) override {
        return std::exp((v - u).norm_squared() / (-2.0 * temperature));
    }
};

struct TwoGauss : DistributionFunc {
    Vec_v u;
    TwoGauss(double temperature = 1.0,
             Vec_v shift = {2.0, 0.0, 0.0}) :
                                                             DistributionFunc(temperature),
                                                             u(std::move(shift)) {}

    double get (const Vec_v& v) override {
        return NORM_COEF * 0.5 *
                        (
                                std::exp((v - u).norm_squared() / (-2.0 * temp))
                                + std::exp((v + u).norm_squared() / (-2.0 * temp))
                        ) / sqrt(std::pow(temp, 3));
    }
};

struct MaxwellBig : DistributionFunc {
    Vec_v u;
    double temperature = 1;
    MaxwellBig(double temperature, Vec_v shift):
            DistributionFunc(temperature),
            u(std::move(shift)){}

    double get (const Vec_v& v) override {
        return std::exp(-0.5 * v.norm_squared())
                + std::exp(v.norm_squared() / (-2.0 * temp)) / std::pow(temp, 3/2);
    }
};

#endif //TASK_RELAX_DISTRIBUTION_H
