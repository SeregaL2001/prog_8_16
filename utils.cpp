#include "utils.hpp"

double u1_dxdx(double x, double y) {
    return y * y * (y - 1) * (6 * x - 2); 
}

double u1_dxdy(double x, double y) {
    return (3 * x * x - 2 * x) * (3 * y * y - 2 * y); 
}

double u1_dydy(double x, double y) {
    return x * x * (x - 1) * (6 * y - 2); 
}

double u2_dxdx(double x, double y) {
    return (y - 1) * std::sin(y) * (2 * std::cos(x) - (x - 1) * std::sin(x));
}

double u2_dxdy(double x, double y) {
    return (std::sin(y) + (y - 1) * std::cos(y)) * (std::sin(x) + (x - 1) * std::cos(x));
}

double u2_dydy(double x, double y) {
    return (x - 1) * std::sin(x) * (2 * std::cos(y) - (y - 1) * std::sin(y));
}

double f1(double x, double y) {
    return -2 * u1_dxdx(x, y) - u1_dydy(x, y) - u2_dxdy(x, y);
}

double f2(double x, double y) {
    return -2 * u2_dydy(x, y) - u2_dxdx(x, y) - u1_dxdy(x, y);
}

//------------------------------------------------------------
// double u1(double x, double y) {
//     return std::sin(M_PI*x)*std::sin(M_PI*y); 
// }

// double u2(double x, double y) {
//     return std::sin(M_PI*x)*std::sin(M_PI*y);
// }

// double f1(double x, double y) {
//     return 3*std::pow(M_PI,2)*std::sin(M_PI*x)*std::sin(M_PI*y) - std::pow(M_PI,2)*std::cos(M_PI*x)*std::cos(M_PI*y);
// }

// double f2(double x, double y) {
//     return 3*std::pow(M_PI,2)*std::sin(M_PI*x)*std::sin(M_PI*y) - std::pow(M_PI,2)*std::cos(M_PI*x)*std::cos(M_PI*y);
// }
//------------------------------------------------------------

double get_a_ij(int i, int j, int N, double h) {
    int i_f, j_f, i_u, j_u; 
    // Первая половина заполняет f_1, вторая f_2
    // Если j в первой половине, то умножается на u_1, иначе на u_2
    if (i <= (N - 1) * (N - 1)) {
        i_f = (i - 1) / (N - 1) + 1;
        j_f = (i - 1) % (N - 1) + 1;
        
        if (j <= (N - 1) * (N - 1)) {
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
            
            if ((i_f == i_u + 1 || i_f == i_u - 1) && j_f == j_u)
                return -2.0 / (h * h);
            if (i_f == i_u && j_f == j_u)
                return 6.0 / (h * h);
            if ((j_f == j_u + 1 || j_f == j_u - 1) && i_f == i_u)
                return -1.0 / (h * h);
        }
        else {
            j -= (N - 1) * (N - 1);
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
            
            if (i_f == i_u + 1 && j_f == j_u + 1)
                return -1.0 / (4 * h * h);
            if (i_f == i_u - 1 && j_f == j_u + 1)
                return 1.0 / (4 * h * h);
            if (i_f == i_u + 1 && j_f == j_u - 1)
                return 1.0 / (4 * h * h);
            if (i_f == i_u - 1 && j_f == j_u - 1)
                return -1.0 / (4 * h * h);
        }
    } 
    else {
        i -= (N - 1) * (N - 1);
        i_f = (i - 1) / (N - 1) + 1;
        j_f = (i - 1) % (N - 1) + 1;
        
        if (j <= (N - 1) * (N - 1)) {
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
            
            if (i_f == i_u + 1 && j_f == j_u + 1)
                return -1.0 / (4 * h * h);
            if (i_f == i_u - 1 && j_f == j_u + 1)
                return 1.0 / (4 * h * h);
            if (i_f == i_u + 1 && j_f == j_u - 1)
                return 1.0 / (4 * h * h);
            if (i_f == i_u - 1 && j_f == j_u - 1)
                return -1.0 / (4 * h * h);
        }
        else {
            j -= (N - 1) * (N - 1);
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
        
            if ((i_f == i_u + 1 || i_f == i_u - 1) && j_f == j_u)
                return -1.0 / (h * h);
            if (i_f == i_u && j_f == j_u)
                return 6.0 / (h * h);
            if ((j_f == j_u + 1 || j_f == j_u - 1) && i_f == i_u)
                return -2.0 / (h * h);
        }
    }
    
    return 0;
}

double get_f_i(int i, int N, std::vector<double> &x, std::vector<double> &y) {
    if (i <= (N - 1) * (N - 1)) {
        return f1(x[(i - 1) / (N - 1) + 1], y[(i - 1) % (N - 1) + 1]);  
    } 
    else {
        i -= (N - 1) * (N - 1);
        return f2(x[(i - 1) / (N - 1) + 1], y[(i - 1) % (N - 1) + 1]);    
    }
}