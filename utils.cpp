#include "utils.hpp"

double u1_dxdx(double x, double y) {
    return y * (y - 1) * 2; 
}

double u1_dxdy(double x, double y) {
    return (2 * std::pow(x,1) - 1) * (2 * std::pow(y,1) - 1); 
}

double u1_dydy(double x, double y) {
    return x * (x - 1) * 2; 
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

// извлечение нужного элемента матрицы A
double get_a_ij(int i, int j, int N, double h) {
    int i_f, j_f, i_u, j_u; 
    // Первая половина заполняет коэф соответствующие f1, вторая f2
    if (i <= (N - 1) * (N - 1)) {
        i_f = (i - 1) / (N - 1) + 1; // берем целую часть от деления 
        j_f = (i - 1) % (N - 1) + 1; // берем остаток от деления 
        
        // Если j в первой половине, то умножается на u1,
        if (j <= (N - 1) * (N - 1)) {
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
            
            if ((i_f == i_u + 1 || i_f == i_u - 1) && j_f == j_u) 
                return -2.0 / (h * h); // коэф при u1_i+1,j; u1_i-1,j
            if (i_f == i_u && j_f == j_u)
                return 6.0 / (h * h);  // коэф при u1_i,j
            if ((j_f == j_u + 1 || j_f == j_u - 1) && i_f == i_u)
                return -1.0 / (h * h);  // коэф при u1_i,j+1; u1_i,j-1
        }
        else { // иначе на u2
            j -= (N - 1) * (N - 1);
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
            
            if (i_f == i_u + 1 && j_f == j_u + 1)
                return -1.0 / (4 * h * h);  // коэф при u2_i+1,j+1 
            if (i_f == i_u - 1 && j_f == j_u + 1)
                return 1.0 / (4 * h * h); // коэф при u2_i-1,j+1
            if (i_f == i_u + 1 && j_f == j_u - 1)
                return 1.0 / (4 * h * h); // коэф при u2_i+1,j-1
            if (i_f == i_u - 1 && j_f == j_u - 1)
                return -1.0 / (4 * h * h); // коэф при u2_i-1,j-1
        }
    } 
    else { // -//- вторая - f2
        i -= (N - 1) * (N - 1);
        i_f = (i - 1) / (N - 1) + 1;
        j_f = (i - 1) % (N - 1) + 1;

        // -//- умножается на u1,
        if (j <= (N - 1) * (N - 1)) {
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
            
            if (i_f == i_u + 1 && j_f == j_u + 1)
                return -1.0 / (4 * h * h); // коэф при u1_i+1,j+1
            if (i_f == i_u - 1 && j_f == j_u + 1)
                return 1.0 / (4 * h * h); // коэф при u1_i-1,j+1
            if (i_f == i_u + 1 && j_f == j_u - 1)
                return 1.0 / (4 * h * h); // коэф при u1_i+1,j-1
            if (i_f == i_u - 1 && j_f == j_u - 1)
                return -1.0 / (4 * h * h); // коэф при u1_i-1,j-1
        }
        else { // иначе на u2
            j -= (N - 1) * (N - 1);
            i_u = (j - 1) / (N - 1) + 1;
            j_u = (j - 1) % (N - 1) + 1;
        
           if ((i_f == i_u + 1 || i_f == i_u - 1) && j_f == j_u)
                return -1.0 / (h * h); // коэф при u2_i,j; u2_i-1,j
            if (i_f == i_u && j_f == j_u)
                return 6.0 / (h * h); // коэф при u2_i,j
            if ((j_f == j_u + 1 || j_f == j_u - 1) && i_f == i_u)
                return -2.0 / (h * h); // коэф при u2_i,j+1; u2_i,j-1
        }
    }
    
    return 0;
}

// извлечение нужного элемента вектора fk, где k = 1,2
double get_f_k(int i, int N, std::vector<double> &x, std::vector<double> &y) {
    if (i <= (N - 1) * (N - 1)) {
        return f1(x[(i - 1) / (N - 1) + 1], y[(i - 1) % (N - 1) + 1]);  
    } 
    else {
        i -= (N - 1) * (N - 1);
        return f2(x[(i - 1) / (N - 1) + 1], y[(i - 1) % (N - 1) + 1]);    
    }
}