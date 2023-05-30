#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>

#define EPS 1e-8
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double u1_dxdx(double x, double y);
double u1_dxdy(double x, double y);
double u1_dydy(double x, double y);
double u2_dxdx(double x, double y);
double u2_dxdy(double x, double y);
double u2_dydy(double x, double y);
double f1(double x, double y);
double f2(double x, double y);

double get_a_ij(int i, int j, int N, double h);
double get_f_k(int i, int N, std::vector<double> &x, std::vector<double> &y);