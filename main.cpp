#include <iostream>
#include <iomanip> 
#include <math.h>
#include <vector>
#include <unistd.h>
#include <string>
#include <fstream>
#include <armadillo>

#define EPS 1e-8
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace arma;
// g++ -std=c++17 -larmadillo main.cpp

double u1(double x, double y) {
    return x * x * y * y * (x - 1) * (y - 1); 
}

double u2(double x, double y) {
    return (x - 1) * (y - 1) * sin(x) * sin(y);
}

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
    return (y - 1) * sin(y) * (2 * cos(x) - (x - 1) * sin(x));
}

double u2_dxdy(double x, double y) {
    return (sin(y) + (y - 1) * cos(y)) * (sin(x) + (x - 1) * cos(x));
}

double u2_dydy(double x, double y) {
    return (x - 1) * sin(x) * (2 * cos(y) - (y - 1) * sin(y));
}

double f1(double x, double y) {
    return -2 * u1_dxdx(x, y) - u1_dydy(x, y) - u2_dxdy(x, y);
}

double f2(double x, double y) {
    return -2 * u2_dydy(x, y) - u2_dxdx(x, y) - u1_dxdy(x, y);
}


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


// Тут используется arma 
//Реализовывание метода Ричардсона, одношагового 
Col<double> richardson(mat &A, Col<double> &f) {
    Col<double> x(size(f), fill::zeros);    
    mat D  = diagmat(A);
    mat sumLR = A - D;
    mat invD = inv(diagmat(D));
    double tau = 0.7;
    Col<double> new_x(size(f), fill::zeros);
    while(true) {
        // Col<double> new_x(size(f), fill::zeros);
        new_x = (1.0 - tau)*x + tau * invD * (f - sumLR * x);
        auto nrm = norm(new_x - x); 
        if (nrm < EPS)
            break;
        x = new_x;
    }

    return x;
}

int main(int argc, char *argv[]) 
{
    std::ofstream output;
    int N;
    double h;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> u1_arr;
    std::vector<double> u2_arr;


    std::cout << "Введите число точек: " << std::endl;
    std::cin >> N;

    // // Парсинг параметров командной строки
    // int res = 0;
    // while ((res = getopt(argc, argv, "N:")) != -1){
    //     switch (res) {
    //         case 'N': 
    //             N = static_cast<unsigned int>(std::stoi(optarg)); 
    //             break;
    //         default:
    //             break;
    //     };
    // };
    
    std::cout << std::setprecision(12);
    
    // std::cout << "N = " << N << "\n";
    
    h = 1.0 / N;
        
    for (int i = 0; i <= N; i++) {
        x.push_back(static_cast<double>(i) / N);
        y.push_back(static_cast<double>(i) / N);
    }
    
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            u1_arr.push_back(u1(x[i], y[j]));
            u2_arr.push_back(u2(x[i], y[j]));
        }
    }
        
    
    // Заполняем матрицу A
    int m_size = 2 * (N - 1) * (N - 1);
    mat A(m_size, m_size, fill::zeros);
    for (int i = 0; i < m_size; i++) {
        for (int j = 0; j < m_size; j++) {
            // A начинает нумирацию с нуля, функция get_a_ij ожидает нумирацию с 1
            A(i, j) = get_a_ij(i + 1, j + 1, N, h);
        }
    }
                
    // Заполняем правую часть
    Col<double> f(m_size, fill::zeros);
    for (int i = 0; i < m_size; i++) {
        // f начинает нумирацию с нуля, функция get_f_i ожидает нумирацию с 1
        f(i) = get_f_i(i + 1, N, x, y);
    }
        
    Col<double> u_richardson = richardson(A, f);    
    Col<double> u_original(m_size, fill::zeros);
    for (int i = 0; i < 2 * (N - 1) * (N - 1); i++) {
        if (i < (N - 1) * (N - 1)) {
            u_original(i) = u1_arr[i];
        }
        else {
            u_original(i) = u2_arr[i - (N - 1) * (N - 1)];
        }
    }
        
    double n = norm(u_original - u_richardson);
    std::cout << "погрешность вычислений = " << n << std::endl;
        
    // output.open("./cppoutput/norms.txt", std::ios::app);
    // output << n << " ";
    // output.close();
    
    // output.open("./cppoutput/N.txt", std::ios::app);
    // output << N << " ";
    // output.close();
    
    return 0;
}