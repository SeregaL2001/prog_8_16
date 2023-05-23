#include "linalg.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

std::ostream& operator<<(std::ostream& os, const MatrixLike& m) {
  for (int j = 0; j < m.size(); ++j)
    for (int i = 0; i < m.size(); ++i) os << m.at(i, j) << (i == m.size() - 1 ? "\n" : ",");
  return os;
}

Matrix::Matrix(std::vector<double>&& data) : data_(data), size_(std::sqrt(data_.size())) {}
Matrix::Matrix(int size) : data_(size * size), size_(size) {}

int Matrix::size() const { return size_; }
double Matrix::at_impl(int i, int j) const { return data_[i + j * size_]; }
double& Matrix::at_impl(int i, int j) { return data_[i + j * size_]; }

Matrix Matrix::dot(const MatrixLike& a, const MatrixLike& b) {
  if (a.size() != b.size())
    throw std::runtime_error("incorrect matrix size");
  Matrix out(a.size());
  for (int i = 0; i < a.size(); ++i)
    for (int j = 0; j < a.size(); ++j) {
      auto& e = out.at(i, j);
      for (int k = 0; k < a.size(); ++k) e += a.at(k, j) * b.at(i, k);
    }
  return out;
}
Matrix Matrix::dot(const MatrixLike& other) const { return Matrix::dot(*this, other); }
std::vector<double> Matrix::dot(const MatrixLike& a, const std::vector<double>& vec) {
  if (a.size() != vec.size())
    throw std::runtime_error("incorrect matrix size");
  std::vector<double> out(a.size());
  for (int i = 0; i < a.size(); ++i) {
    auto& e = out[i];
    for (int k = 0; k < a.size(); ++k) e += a.at(k, i) * vec[k];
  }
  return out;
}

std::vector<double> Matrix::dot(const std::vector<double>& vec) const { return Matrix::dot(*this, vec); }
Matrix Matrix::inv(const MatrixLike& a) {
  LUDecomoposition d(a);
  Matrix iL(a.size()), iU(a.size());
  for (int i = 0; i < a.size(); ++i) {
    iL.at(i, i) = 1.0;
    for (int j = 0; j < i; ++j) {
      double s = 0;
      for (int k = j; k < i; ++k) s += iL.at(j, k) * d.L.at(k, i);
      iL.at(j, i) = -s * iL.at(i, i);
    }
    iU.at(i, i) = 1.0 / d.U.at(i, i);
    for (int j = 0; j < i; ++j) {
      double s = 0;
      for (int k = j; k < i; ++k) s += iU.at(k, j) * d.U.at(i, k);
      iU.at(i, j) = -s * iU.at(i, i);
    }
  }
  return iU.dot(iL);
}
Matrix Matrix::inv() const { return Matrix::inv(*this); }

LUDecomoposition::LMatrix::LMatrix(Matrix& data) : data_(data) {}

int LUDecomoposition::LMatrix::size() const { return data_.size(); }

double& LUDecomoposition::LMatrix::at_impl(int i, int j) {
  return (i - j < 0) ? data_.at(i, j) : throw std::runtime_error("accessing non writable part of L");
}
double LUDecomoposition::LMatrix::at_impl(int i, int j) const { return (i - j < 0) ? data_.at(i, j) : (i == j); }

LUDecomoposition::UMatrix::UMatrix(Matrix& data) : data_(data) {}

int LUDecomoposition::UMatrix::size() const { return data_.size(); }
double& LUDecomoposition::UMatrix::at_impl(int i, int j) {
  return (i - j < 0) ? throw std::runtime_error("accessing non writable part of U") : data_.at(i, j);
}
double LUDecomoposition::UMatrix::at_impl(int i, int j) const { return (i - j < 0) ? 0.0 : data_.at(i, j); }

LUDecomoposition::LUDecomoposition(const MatrixLike& matrix) : data_(matrix.size()), L(data_), U(data_) {
  int n = matrix.size();
  auto& lower = const_cast<LMatrix&>(L);
  auto& upper = const_cast<UMatrix&>(U);
  for (int i = 0; i < n; i++) {
    for (int k = i; k < n; k++) {
      double sum = 0;
      for (int j = 0; j < i; j++) sum += (L.at(j, i) * U.at(k, j));
      upper.at(k, i) = matrix.at(k, i) - sum;
    }
    for (int k = i + 1; k < n; k++) {
      double sum = 0;
      for (int j = 0; j < i; j++) sum += (L.at(j, k) * U.at(i, j));
      lower.at(i, k) = (matrix.at(i, k) - sum) / U.at(i, i);
    }
  }
}