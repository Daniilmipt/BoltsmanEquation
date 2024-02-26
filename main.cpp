#include <iostream>
#include "3_task/Distribution3.h"
#include "conductivity_task.h"

using Vec_v = Vector_3d<double>;
using Vec_node = Vector_3d<int>;

int main() {
    Conduct conduct;
    conduct.find_N0();
    conduct.init_f();

    std::vector<std::vector<double>> theta = conduct.find_theta();

    conduct.print_temperature(theta);


//    std::string data_dir = "data/";
//
//    double v_cut(5.0);
//    const std::size_t N1(20), N2(6), N3(8), N4(10), N5(20);
//    double temperature(0.975);
//    Vec_v u(2, 0, 0);
//    Distribution3<N1, N1, N1> f1(TwoGauss(temperature, u), v_cut, 1.0, true);
////    Distribution3<N2, N2, N2> f2(TwoGauss(temperature, u), v_cut, 1.0, true);
////    Distribution3<N3, N3, N3> f3(TwoGauss(temperature, u), v_cut, 1.0, true);
////    Distribution3<N4, N4, N4> f4(TwoGauss(temperature, u), v_cut, 1.0, true);
//    Distribution3<N5, N5, N5> f5(TwoGauss(temperature, u), v_cut, 1.0, true);
//
//    double w_min = f5.n0 * std::pow(f5.v_cut, 4) / (6 * std::sqrt(M_PI));
//
//    int n_col1(1000), n_col2(5000), n_col3(10000), n_col4(250000), n_col5(250000);
//    double time(2), tau(0.02);
//
//    int N_min = std::floor(w_min * tau);
//
////    std::cout << "1) n = " << f1.get_norm_coeff() << "\tT = " << f1.temperature() << "\tE = " << f1.energy() << "\tvec_p = " << f1.momentum() << '\n';
////    std::cout << "2) n = " << f2.get_norm_coeff() << "\tT = " << f2.temperature() << "\tE = " << f2.energy() << "\tvec_p = " << f2.momentum() << '\n';
////    std::cout << "3) n = " << f3.get_norm_coeff() << "\tT = " << f3.temperature() << "\tE = " << f3.energy() << "\tvec_p = " << f3.momentum() << '\n';
////    std::cout << "4) n = " << f4.get_norm_coeff() << "\tT = " << f4.temperature() << "\tE = " << f4.energy() << "\tvec_p = " << f4.momentum() << '\n';
////    std::cout << "5) n = " << f5.get_norm_coeff() << "\tT = " << f5.temperature() << "\tE = " << f5.energy() << "\tvec_p = " << f5.momentum() << '\n';
//
//
//    std::string filename = data_dir + "several_f.txt";
//    std::ofstream out(filename);
////    out << "# v_cut = " << v_cut << ", tau = " << tau << '\n';
////    out << "#1) N = " << N1 << ", n_col = " << n_col4 << '\n';
////    out << "#2) N = " << N2 << ", n_col = " << n_col4 << '\n';
////    out << "#3) N = " << N3 << ", n_col = " << n_col4 << '\n';
////    out << "#4) N = " << N4 << ", n_col = " << n_col4 << '\n';
////    out << "#5) N = " << N5 << ", n_col = " << n_col5 << '\n';
////
////    out << "0.0\t" << f1.M_2(1,1) << '\t'
////        << f2.M_2(1,1) << '\t'
////        << f3.M_2(1,1) << '\t'
////        << f4.M_2(1,1) << '\t'
////        << f5.M_2(1,1) << '\n';
//
//    std::cout << "Start\n";
//    Distribution3<N1, N1, N1> f1_0(MaxwellBig(temperature, u), v_cut, 1.0, true);
//    Distribution3<N2, N2, N2> f2_0(TwoGauss(temperature, u), v_cut, 1.0, true);
//    Distribution3<N3, N3, N3> f3_0(TwoGauss(temperature, u), v_cut, 1.0, true);
//    Distribution3<N4, N4, N4> f4_0(TwoGauss(temperature, u), v_cut, 1.0, true);
//    Distribution3<N5, N5, N5> f5_0(TwoGauss(temperature, u), v_cut, 1.0, true);
//    for (int i = 0; i < time / tau; ++i) {
//
//    }
//    out.close();
//
//    return 0;
}
