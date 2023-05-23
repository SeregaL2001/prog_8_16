#include <cmath>
#include <fstream>
#include <iostream>

#include "linalg.hpp"

constexpr int N = 25;
// #define N 25

double f1(double x, double y) 
{ 
  return std::pow(M_PI, 2) * std::cos(M_PI * (x - y)); 
}

double f2(double x, double y) 
{
   return std::pow(M_PI, 2) * std::cos(M_PI * (x - y)); 
}
//~~~~~~~~~~~~~~~~~~~~~~~~

constexpr double h = 1.0 / (N - 1);
constexpr double h2 = h * h;

bool is_bound(int i, int j) { return i == 0 || i == N - 1 || j == 0 || j == N - 1; }

double f1(int i, int j) { return is_bound(i, j) ? 0.0 : h2 * f1(i * h, j * h); }
double f2(int i, int j) { return is_bound(i, j) ? 0.0 : h2 * f2(i * h, j * h); }

int pos(int i, int j) { return i + N * j; }

std::vector<double> f_k(int k) {
  std::vector<double> out(N * N);
  switch (k) {
    case 1: {
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) out[pos(i, j)] = f1(i, j);
      break;
    }
    case 2: {
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) out[pos(i, j)] = f2(i, j);
      break;
    }
    default: throw std::runtime_error("Unexpected index of f");
  }
  return out;
}

Matrix B() {
  Matrix out(N * N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      int p = pos(i, j);
      if (is_bound(i, j)) {
        out.at(p, p) = 1.0;
      } else {
        out.at(p, p) = 2.0;
        out.at(pos(i, j + 1), p) = -1.0;
        out.at(pos(i, j - 1), p) = -1.0;
      }
    }
  return out;
}

Matrix C() {
  Matrix out(N * N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      int p = pos(i, j);
      if (is_bound(i, j)) {
        out.at(p, p) = 1.0;
      } else {
        out.at(p, p) = 2.0;
        out.at(pos(i + 1, j), p) = -1.0;
        out.at(pos(i - 1, j), p) = -1.0;
      }
    }
  return out;
}

Matrix A() {
  Matrix out(N * N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      int p = pos(i, j);
      if (is_bound(i, j)) {
        out.at(p, p) = 1.0;
      } else {
        out.at(pos(i + 1, j + 1), p) = 0.25;
        out.at(pos(i - 1, j - 1), p) = 0.25;
        out.at(pos(i + 1, j - 1), p) = -0.25;
        out.at(pos(i - 1, j + 1), p) = -0.25;
      }
    }
  return out;
}
Matrix P() {
  Matrix out(N * N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) {
      int p = pos(i, j);
      if (is_bound(i, j)) {
        out.at(p, p) = 1.0;
      } else {
        out.at(p, p) = 4.0;
        out.at(pos(i, j + 1), p) = -1.0;
        out.at(pos(i, j - 1), p) = -1.0;
        out.at(pos(i + 1, j), p) = -1.0;
        out.at(pos(i - 1, j), p) = -1.0;
      }
    }
  return out.inv();
}

std::vector<double> operator+(std::vector<double>&& a, const std::vector<double>& b) {
  if (a.size() != b.size())
    throw std::runtime_error("incorrect vector size");
  for (int i = 0; i < a.size(); ++i) a[i] += b[i];
  return a;
}

std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size())
    throw std::runtime_error("incorrect vector size");
  for (int i = 0; i < a.size(); ++i) a[i] += b[i];
  return a;
}

std::vector<double> operator-(std::vector<double>&& a, const std::vector<double>& b) {
  if (a.size() != b.size())
    throw std::runtime_error("incorrect vector size");
  for (int i = 0; i < a.size(); ++i) a[i] -= b[i];
  return a;
}

std::vector<double>& operator-=(std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size())
    throw std::runtime_error("incorrect vector size");
  for (int i = 0; i < a.size(); ++i) a[i] -= b[i];
  return a;
}

std::vector<double> operator*(std::vector<double>&& a, double k) {
  for (int i = 0; i < a.size(); ++i) a[i] *= k;
  return a;
}

double norm(const std::vector<double>& a);

int main() {
  std::vector<double> u(N * N), v(N * N);
  Matrix _P(P()), PA(_P.dot(A())), PB(_P.dot(B())), PC(_P.dot(C()));
  std::vector<double> Pf1(_P.dot(f_k(1))), Pf2(_P.dot(f_k(2)));

  double tau = 1.25; // взяли подбором 
  // при N = 20 - 0.75
  // при N = 25 - 1.25 
  // при N = 50 - 
  double u_priv = 1, v_priv = 1;
  for (int counter = 1;; ++counter) {
    auto du = (PA.dot(v) + PB.dot(u) - Pf1) * tau;
    auto dv = (PA.dot(u) + PC.dot(v) - Pf2) * tau;
    double nu = norm(du);
    double nv = norm(dv);
    double error;
    if (counter % 100 == 0)
      // std::cout << "Iteration: " << counter << " |du|: " << nu << " |dv|: " << nv << " " << error << std::endl;
    if (nu + nv < 1e-10) {
      std::cout << "Iterations: " << counter << " |du|: " << nu << " |dv|: " << nv << " eror_i/error_i-1 = "<< error << std::endl;
      break;
    }
    u -= du;
    v -= dv;
    error = u_priv/nu;
    u_priv = nu;
    // std::cout << u_priv << std::endl;
    // std::cout << error << std::endl;
  }
  

  std::ofstream("u1.csv") << Matrix(std::move(u));
  std::ofstream("u2.csv") << Matrix(std::move(v));
}

double norm(const std::vector<double>& a) {
  double sum = 0;
  for (int i = 0; i < a.size(); ++i) sum += std::pow(a[i], 2);
  return std::sqrt(sum);
}