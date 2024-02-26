//
// Created by daniil on 30.11.23.
//

#ifndef TASK_RELAX_DISTRIBUTION3_H
#define TASK_RELAX_DISTRIBUTION3_H

#include "../2_task/IntegralSpace.h"
#include "../iterator/IteratorVector.h"
#include "../1_task/Collision.h"
#include "../Grid/KorobovGrid.h"
#include "../Distribution.h"

using Vec_v = Vector_3d<double>;
using Vec_node = Vector_3d<int>;

const double NORM_COEF_3 = 1 / (2 * M_PI) / sqrt(2 * M_PI);

template<std::size_t N_vx = 10, std::size_t N_vy = 10, std::size_t N_vz = 10>
class Distribution3 : public IntegralSpace<N_vx, N_vy, N_vz>, public std::array<std::array<std::array<double, N_vz>, N_vy>, N_vx> {
private:
    // алгоритм - "решето Эратосфена" для нахождения след простого числа
    int next_prime(int n) {
        int size = 2 * n;
        bool is_prime[size];
        for (int i = 0; i < size; ++i)
            is_prime[i] = true;

        for (int i = 2; i * i < size; ++i) {
            if (is_prime[i]) {
                for (int j = i * i; j < size; j += i) {
                    is_prime[j] = false;
                }
            }
        }
        int next_prime = n;
        while (!is_prime[next_prime]) {  // от [n до 2n-1] точно есть простое, так что цикл не вечный
            next_prime++;
        }
        return next_prime;
    }

    //7 пункт 3 семинар, пункт Б)
    // првоеряем на неотрицательность полученных значений
    bool change_positive_value(double** speed_nodes, const double* mult_coeff, int n) {
        bool is_positive = false;

        for (int i = 0; i < n; ++i) {
            if (*(speed_nodes[i]) + mult_coeff[i] < 0) {
                is_positive = true;
                break;
            }
        }

        if (!is_positive) {
            for (int j = 0; j < n; ++j) {
                *(speed_nodes[j]) += mult_coeff[j];
            }
        }
        return is_positive;
    }

    int find_n0() {  // число узлов внутри сферы радиуса v_cut
        int n0 = 0;
        for (auto& node : *this) {
            if (this->iter_vec_v(node, 0).norm() < this->v_cut)
                n0++;
        }
        return n0;
    }

public:
    int n0 = find_n0();

    typedef IteratorVector<N_vx, N_vy, N_vz, Vec_node> iterator;
    typedef IteratorVector<N_vx, N_vy, N_vz, Vec_node> const_iterator;  // TODO: подумать над разницей между iterator и const_iterator

    iterator begin() {
        return iterator({0, 0, 0});
    }
    iterator end() {
        return iterator({N_vx, 0, 0});
    }
    const_iterator begin() const {
        return const_iterator({0, 0, 0});
    }
    const_iterator end() const {
        return const_iterator({N_vx, 0, 0});
    }

    //вычисляем значение функции в узле по оператору ()
    double& get (const Vec_node& node) {
        return (*this)[node.x][node.y][node.z];
    }
    double get (const Vec_node& node) const {
        return (*this)[node.x][node.y][node.z];
    }

    //dv_x*dv_y*dv_z
    double dv_xyz_multiply() const {
        return  this->vec_dv.x * this->vec_dv.y * this->vec_dv.z;
    }

    // нормировка, чтобы концентрация газа была равна 1
    double get_norm_coeff() const {
        double n = 0.0;
        for (auto& node : *this)
            n += (*this).get(node);
        double v = n * dv_xyz_multiply();
        return n * dv_xyz_multiply();
    }

    Vec_v momentum() const {
        Vec_v momentum(0.0, 0.0, 0.0);
        for (auto& node : *this)
            momentum += this->iter_vec_v(node, 0) * (*this).get(node);
        return momentum * dv_xyz_multiply();
    }

    double energy() const {
        double energy = 0.0;
        for (auto& node : *this)
            energy += this->iter_vec_v(node, 0).norm_squared() * (*this).get(node);
        return energy * dv_xyz_multiply();
    }

    double temperature() const {
        return (energy() / get_norm_coeff() - momentum().norm_squared()) / 3;
    }

    // M_{i,j} = integral(v_i * v_j * f * d^3v)
    double M_2(int alpha, int beta) {   // node = 1, 2 or 3
        double m = 0.0;
        for (auto& node : *this)
            m += this->v_iter(node[alpha], this->vec_dv[alpha])
                 * this->v_iter(node[beta], this->vec_dv[beta])
                 * (*this).get(node);
        return m * dv_xyz_multiply();
    }

    //нормализация
    void normalization() {
        double n = get_norm_coeff();
        for (auto& node : *this)
            (*this).get(node) /= n;
    }

    template<typename Functional>
    Distribution3(Functional obj_or_func, bool is_normalization = true) {  // obj_or_func должен быть вызываемым от v_t, т.е. obj_or_func(v) должно выдавать double
        for (Vec_node& node : *this) {
            Vec_v v = this->iter_vec_v(node, 0);
            if (v.norm() < this->get_v_cut())
                (*this).get(node) = obj_or_func(v);
            else
                (*this).get(node) = 0.0;
        }
        if (is_normalization)
            normalization();
    }

    Distribution3<N_vx, N_vy, N_vz> func_distr_maxwell(double t, Vec_v shift){
        Distribution3<N_vx, N_vy, N_vz> distribution = Distribution3<N_vx, N_vy, N_vz>();
        double n = 0;
        for (Vec_node& node : distribution) {
            distribution.get(node) = std::exp(-0.5 * (distribution.iter_vec_v(node,0)-shift).norm_squared()/t);
            n += distribution.get(node);
        }
        n = n * dv_xyz_multiply();

        for (Vec_node& node : distribution) {
            distribution.get(node) /= n;
        }
        return distribution;
    }

    Distribution3() = default;
    Distribution3<N_vx, N_vy, N_vz> get_distribution3(Distribution3 f1, Distribution3 f2){
        Distribution3<N_vx, N_vy, N_vz> distribution3 = Distribution3<N_vx, N_vy, N_vz>();
        for (Vec_node& node : distribution3) {
            distribution3.get(node) = f1.get(node) - f2.get(node);
        }
        return distribution3;
    }

    Distribution3<N_vx, N_vy, N_vz> get_distribution_time(Distribution3 f1, Distribution3 f2, double t){
        Distribution3<N_vx, N_vy, N_vz> distribution3 = Distribution3<N_vx, N_vy, N_vz>();
        for (Vec_node& node : distribution3) {
            distribution3.get(node) = (f1.get(node) - f2.get(node)) / t;
        }
        return distribution3;
    }

    Distribution3(double v_max) : IntegralSpace<N_vx, N_vy, N_vz>(v_max) {}
    Distribution3(double v_max, double d_max) : IntegralSpace<N_vx, N_vy, N_vz>(v_max, d_max) {}

    template<typename Functional>
    Distribution3(Functional VelocityDistFunc, double v_max = 5.0, double b_max_sq = 1.0, bool haveRationing = true) : IntegralSpace<N_vx, N_vy, N_vz>(v_max, b_max_sq) {
        init(VelocityDistFunc, haveRationing);
    }

    template<typename Functional>
    void init(Functional VelocityDistFunc, double temp = 1.0, bool is_normalized = true) {  // TODO: delete or make constructor
        for (Vec_node & node : *this) {
            Vec_v v = this->iter_vec_v(node, 0);
            if (v.norm() < this->get_v_cut())
                (*this).get(node) = VelocityDistFunc.get(v);
            else
                (*this).get(node) = 0.0;
        }
        if (is_normalized)
            normalization();
    }

    void add_header(std::ofstream& out) const {
        out << "# n_vx = " << N_vx << ", n_vy = " << N_vy << ", n_vz = " << N_vz << '\n';
    }
    void add_param(std::ofstream& out) const {
        out << "n = " << get_norm_coeff() << "\nT = " << temperature() << "\nE = " << energy() << "\nvec_p = " << momentum() << '\n';
    }

    void save_to_file_3dim(std::string& filename, bool haveParam = false) const {
        std::ofstream outfile(filename);
        outfile << "# f(v_x, v_y, v_z)\n";
        add_header(outfile);
        if (haveParam)
            add_param(outfile);
        for (Vec_node& node : *this) {
            Vec_v v = this->iter_vec_v(node, 0);
            outfile << v.x << ' ' << v.y << ' ' << v.z << ' ' << (*this).get(node) << '\n';
        }
    }
    void save_to_file_1dim(std::string filename, Vec_node print = Vec_node(-1, N_vy / 2, N_vz / 2), bool haveParam = false) const {
        std::ofstream outfile(filename);
        add_header(outfile);
        outfile << "# f(v_x, v_y, v_z)\n";
        int out_i = 0;
        for (int i : {1, 2, 3}) {
            if (print[i] >= 0) {
                out_i += i;
                outfile << ", where v_";
                switch (i)
                {
                    case 1:
                        outfile << "x";
                        break;
                    case 2:
                        outfile << "y";
                        break;
                    case 3:
                        outfile << "z";
                        break;

                    default:
                        break;
                }
                outfile << " = " << this->iter_vec_v(print[i], this->vec_dv[i]) << ' ';
            }
        }
        outfile << "\n\n";
        if (haveParam)
            add_param(outfile);

        out_i = 1 + 2 + 3 - out_i;
        int not_out1 = out_i % 3 + 1;  // 1, 2 or 3
        int not_out2 = (out_i + 1) % 3 + 1;
        for (auto& node : (*this)) {
            if (node[not_out1] == print[not_out1] && node[not_out2] == print[not_out2])
                outfile << this->iter_vec_v(node[out_i], this->dv[out_i]) << '\t' << (*this).get(node) << '\n';
        }
    }

    void make_iteration(double tau, int prime_number, bool haveShuffle = false, bool haveShifft = false) {
        KorobovGrid kor_grid(next_prime(prime_number), 8, haveShuffle, haveShifft);

        //6 пункт 3 семинар
        const double c = 1 / (4 * std::sqrt(2)) * this->d_max * dv_xyz_multiply() * tau * n0 * n0 / prime_number;

        int count = 0;
        int count_modified = 0;
        int count_ignore = 0;
        double maxDelta = 0;
        for (auto& row : kor_grid) {
            ++count;
            Collision<N_vx, N_vy, N_vz> collision(row, this->v_cut, this->d_max);
            Vec_node nd1 = collision.node_1;
            Vec_node nd2 = collision.node_2;
            if (collision.isGood) {
                double r;
                Vec_node lam, mu, lam_add, mu_add;

                //находим узлы и скорости
                int cnt = 0;
                std::tie(lam, mu, lam_add, mu_add, r, cnt) = collision.get_new_v();

                // TODO: cтереть после проверки
                // if (col.isGood && count == 56) {
                //     v_t v1 = iter_vec_v(col.node_1);
                //     v_t v2 = iter_vec_v(col.node_2);
                //     double energy_0 = v1.squared() + v2.squared();
                //     double energy_1 = iter_vec_v(lam).squared() + iter_vec_v(mu).squared();
                //     double energy_2 = iter_vec_v(lam_add).squared() + iter_vec_v(mu_add).squared();
                //     double r_add = (energy_0 - energy_1) / (energy_2 - energy_1);
                //     v_t v_lam = iter_vec_v(lam);
                //     v_t v_mu = iter_vec_v(mu);
                //     v_t v_lam_add = iter_vec_v(lam_add);
                //     v_t v_mu_add = iter_vec_v(mu_add);
                //     std::cout << energy_0 << " " << energy_1 << " " << energy_2 << '\n';
                //     std::cout << (v_lam + v_mu).x  << " " << (v_lam_add + v_mu_add).x << " " << (v1 + v2).x <<'\n';
                //     std::cout << r << " " << r_add << "\n\n";
                // }
                int v = 0;
                if (collision.isGood) {
                    auto val1 = (*this).get(lam) * (*this).get(mu);
                    auto val2 = std::pow((*this).get(lam_add) * (*this).get(mu_add) / (*this).get(lam) / (*this).get(mu), r);
                    auto val3 = (*this).get(collision.node_1) * (*this).get(collision.node_2);
                    auto val4 = (this->iter_vec_v(collision.node_1, 0) - this->iter_vec_v(collision.node_2, 0)).norm();
                    double omega = ((*this).get(lam) * (*this).get(mu)
                                    * std::pow(
                            (*this).get(lam_add) * (*this).get(mu_add) / (*this).get(lam) / (*this).get(mu), r)
                                    - (*this).get(collision.node_1) * (*this).get(collision.node_2))
                                   * (this->iter_vec_v(collision.node_1, 0) -
                                      this->iter_vec_v(collision.node_2, 0)).norm();
                    double delta = c * omega;

                    if (maxDelta < delta){
                        maxDelta = delta;
                        nd1 = collision.node_1;
                        nd2 = collision.node_2;
                    }
                    // if (omega != omega) {
                    //     std::cout << count  << ' ' << omega << '\t' << col << '\n';
                    //     std::cout << col.node_1 << ' ' << col.node_2 << ' ' << lam << ' ' << mu << lam_add << ' ' << mu_add << '\n';
                    //     std::cout << f(col.node_1) << ' ' << f(col.node_2) << ' ' << (*this)(lam) << ' ' << (*this)(mu) << (*this)(lam_add) << ' ' << (*this)(mu_add) << '\n';
                    // }

                    // 7 пункт 3 семинар, пункт А)
                    double* speed_nodes[6]{
                            &(*this).get(collision.node_1),
                            &(*this).get(collision.node_2),
                            &(*this).get(lam),
                            &(*this).get(mu),
                            &(*this).get(lam_add),
                            &(*this).get(mu_add)
                    };
                    double mult_coeff[]{delta, delta, (r - 1) * delta, (r - 1) * delta, -r * delta, -r * delta};
                    bool modified = change_positive_value(speed_nodes, mult_coeff, 6);
                    // if (change > max_delta)
                    //     max_delta = change;
                    if (modified) {
                        ++count_modified;
                    }
                    // count++;
                }
                else{
                    ++count_ignore;
                }
            }
            else{
                ++count_ignore;
            }
        }
        double perc = (count - count_ignore) * 100.0 / count;
        if (haveShifft)
            kor_grid.random_shift();

        // TODO: если не нужна таблица, то эти строчки не нужны
        // neg_f_table.push_back(count_neg);
        // neg_f_table.push_back(count);

        // return max_delta;
        // TODO: удалить после дебага
        // std::unordered_map<collision_t, int> table;
        // table.reserve(arr_col.size());
        // for (auto& col : arr_col) {
        //     // if (col.node_2.x == arr_col[5000].node_2.x && col.node_2.y == arr_col[5000].node_2.y && col.node_2.z == arr_col[5000].node_2.z)
        //     if (col.isGood)
        //         table[col]++;
        // }
        // int a = 0;
        // for (auto& el : table) {
        //     std::cout << el.first << "\t" << el.second << "\n";
        //     a += el.second;
        // }
        // std::cout << table.size() << " " << a << std::endl;

        // std::cout << count << "\n";
    }

    void make_iteration(double tau, bool haveShuffle = false, bool haveShift = false) {  // TODO: неоптимально, каждый раз надо считать n_col
        int n_col = next_prime(std::ceil(2 * tau * n0 * this->v_cut * this->d_max * n0 * dv_xyz_multiply() * NORM_COEF_3 / std::sqrt(2)));
        make_iteration(tau, n_col, haveShuffle, haveShift);
    }

    void simulate(int n_iter, double tau, double n_col) {
        for (int i = 0; i < n_iter; ++i) {
            make_iteration(tau, n_col, true, true);
        }
    }
    void simulate(double time, double tau, double n_col) {
        int n_iter = time / tau;
        simulate(n_iter, tau, n_col);
    }

    Vec_v do_norm_calc(Distribution3 f1, Distribution3 f2, Distribution3 f3){
        Vec_v v;
//        auto v1 = calc_norm(1, get_distribution3(f1, f2));
//        auto v2 = calc_norm(1, get_distribution3(f3, f2));
//        auto v3 = v1 / v2;
        v.x = -log(calc_norm(1, get_distribution3(f1, f2)) / calc_norm(1, get_distribution3(f3, f2)));
        v.y = -log(calc_norm(2, get_distribution3(f1, f2)) / calc_norm(2, get_distribution3(f3, f2)));
        v.z = -log(calc_norm(3, get_distribution3(f1, f2)) / calc_norm(3, get_distribution3(f3, f2)));
        return v;
    }

    double calc_norm(int k, Distribution3 y){
        double v = 0;
        double speed_squared = 0;
        int cnt = 0;
        for (auto& node : *this) {
            Vec_v vec = y.iter_vec_v(node, 0);
            speed_squared = (y.iter_vec_v(node, 0)).norm_squared();
            if (speed_squared <= k * this->v_cut/3 && speed_squared >= (k - 1) * this->v_cut/3) {
                v += std::pow(y.get(node), 2);
                ++cnt;
            }
        }
        return sqrt(v * this->dv_xyz_multiply());
    }


    double derivation_y_1_u(double temperature, Vec_v shift){
        double v2_x = 0.0;
        double v2 = 0;
        double v_x = 0;
        for (auto& node : *this) {
            v2_x += this->iter_vec_v(node, 0).norm_squared() * this->iter_vec_v(node, 0).x
                    * (*this).get(node);
            v2 += this->iter_vec_v(node, 0).norm_squared() * (*this).get(node);
            v_x += this->iter_vec_v(node, 0).x * (*this).get(node);
        }
        return (v2_x - v_x * v2) * dv_xyz_multiply() / (3 * temperature);
    }

    double derivation_y_2_u(double temperature, Vec_v shift){
        double v_x = 0;
        double v_x2 = 0;
        for (auto& node : *this) {
            v_x += this->iter_vec_v(node, 0).x * (*this).get(node);
            v_x2 += std::pow(this->iter_vec_v(node, 0).x, 2) * (*this).get(node);
        }
        return (v_x2 - v_x) * dv_xyz_multiply() / temperature;
    }

    double derivation_y_1_t(double temperature, Vec_v shift){
        double v2_x = 0.0;
        double v_x = 0.0;
        double v4 = 0;
        double v2 = 0;
        for (auto& node : *this) {
            v2_x += this->iter_vec_v(node, 0).norm_squared() * this->iter_vec_v(node, 0).x
                    * (*this).get(node);
            v4 += std::pow(this->iter_vec_v(node, 0).norm_squared(), 2) * (*this).get(node);
            v2 += this->iter_vec_v(node, 0).norm_squared() * (*this).get(node);
            v_x += this->iter_vec_v(node, 0).x * (*this).get(node);
        }
        return (v4 - 2*shift.x*v2_x - v2*v2 + 2*shift.x*v2*v_x) * dv_xyz_multiply() / (6 * temperature * temperature);
    }

    double derivation_y_2_t(double temperature, Vec_v shift){
        double v2_x = 0.0;
        double v_x = 0.0;
        double v_x2 = 0.0;
        double v4 = 0;
        double v2 = 0;
        for (auto& node : *this) {
            v2_x += this->iter_vec_v(node, 0).norm_squared() * this->iter_vec_v(node, 0).x
                    * (*this).get(node);
            v4 += std::pow(this->iter_vec_v(node, 0).norm_squared(), 2) * (*this).get(node);
            v2 += this->iter_vec_v(node, 0).norm_squared() * (*this).get(node);
            v_x += this->iter_vec_v(node, 0).x * (*this).get(node);
            v_x2 += std::pow(this->iter_vec_v(node, 0).x, 2) * (*this).get(node);
        }
        return (v2_x - 2*shift.x*v_x2 - v2*v_x + 2*shift.x*v_x*v_x) * dv_xyz_multiply() / (temperature * temperature);
    }

    Vec_node get_max_node(){
        Vec_node node;
        double max = -1;
        for (int i = 0; i < N_vx; ++i) {
            for (int j = 0; j < N_vy; ++j) {
                for (int k = 0; k < N_vz; ++k) {
                    if ((*this)[i][j][k] > max){
                        max = (*this)[i][j][k];
                        node.x = i;
                        node.y = j;
                        node.z = k;
                    }
                }
            }
        }
        return node;
    }

    double symmetric(Distribution3 f, Distribution3 f0, double time){
        Distribution3<N_vx, N_vy, N_vz> distr_I = get_distribution_time(f, f0, time);
        double value = 0;
        double squared = 0;
        for (int j = 1; j <= N_vx; ++j) {
            for (int k = 1; k <= N_vy; ++k) {
                for (int i = 1; i <= N_vz; ++i) {
                    value += std::pow(distr_I.get({i-1, j-1, k-1}) - distr_I.get({i-1, j-1, N_vz-k}), 2);
                    squared += std::pow(distr_I.get({i-1,j-1,k-1}), 2);
                }
            }
        }

        squared = std::sqrt(0.5 * squared);
        value = std::sqrt(value) / squared;
        return value;
    }
};


#endif //TASK_RELAX_DISTRIBUTION3_H
