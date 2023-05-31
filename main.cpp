#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>

#define EPS 1e-8
#include "utils.hpp"
// using namespace arma;
// g++ -std=c++17 -larmadillo main.cpp utils.cpp

double u1(double x, double y) {
    return x * y * (x - 1) * (y - 1); 
}

double u2(double x, double y) {
    return std::sin(x) * std::sin(y) * (x - 1) * (y - 1);
}

// Реализовывание метода Ричардсона, одношагового
arma::Col<double> Richardson(arma::mat &A, arma::Col<double> &f, int &count) {
    using namespace arma;
    Col<double> x(size(f), fill::zeros);    
    mat D  = diagmat(A);
    mat LsumR = A - D;
    mat invD = inv(diagmat(D));
    double tau = 0.97;
    Col<double> new_x(size(f), fill::zeros);
    while(true) {
        new_x = (1.0 - tau)*x + tau * invD * (f - LsumR * x);
        double nrm = norm(new_x - x); 

        if (nrm < EPS)
            break;
        x = new_x;
        count+=1;
    }

    return x;
}

int main(int argc, char *argv[]) 
{
    using namespace arma;

    int N;
    std::cout << "Введите число точек: " << std::endl;
    std::cin >> N;
    double h = 1.0 / N;

    std::vector<double> x; 
    std::vector<double> y; 
    for (int i = 0; i <= N; i++) {
        x.push_back(static_cast<double> (i) / N);
        y.push_back(static_cast<double> (i) / N);
    }

    std::vector<double> u1_arr;
    std::vector<double> u2_arr;
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            u1_arr.push_back(u1(x[i], y[j]));
            u2_arr.push_back(u2(x[i], y[j]));
        }
    } 
    
    // Заполняем матрицу A элементами a_ij
    int f_size = 2 * (N - 1) * (N - 1);
    mat A(f_size, f_size, fill::zeros);
    for (int i = 0; i < f_size; i++) {
        for (int j = 0; j < f_size; j++) {
            // A начинает нумирацию с нуля, функция get_a_ij начинает нумирацию с 1 
            // => сдвиг по индексу в цикле  
            A(i, j) = get_a_ij(i + 1, j + 1, N, h);
        }
    }
    
    // A.print("A = ");
                
    // Заполняем f с помощью fk
    Col<double> f(f_size, fill::zeros);
    for (int i = 0; i < f_size; i++) {
        // f начинает нумирацию с нуля, функция get_f_i начинает нумирацию с 1 
        // => сдвиг по индексу в цикле
        f(i) = get_f_k(i + 1, N, x, y);
    }
    // f.print("f = ");
    // std:: cout << size(f) << std::endl;

    int count = 0;
    Col<double> u_richardson = Richardson(A, f, count);    
    Col<double> u_original(f_size, fill::zeros);
    for (int i = 0; i < f_size; i++) {
        if (i < (N - 1) * (N - 1)) {
            u_original(i) = u1_arr[i];
        }
        else {
            u_original(i) = u2_arr[i - (N - 1) * (N - 1)];
        }
    }
        
    double nrm = norm(u_original - u_richardson);
    std::cout << "погрешность вычислений = " << nrm << std::endl;
    
    return 0;
}