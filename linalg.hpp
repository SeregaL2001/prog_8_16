#pragma once

#include <ostream>
#include <vector>

struct MatrixLike {
  double at(int i, int j) const { return at_impl(i, j); }
  double& at(int i, int j) { return at_impl(i, j); }
  virtual int size() const = 0;

 protected:
  virtual double at_impl(int i, int j) const = 0;
  virtual double& at_impl(int i, int j) = 0;
};

std::ostream& operator<<(std::ostream& os, const MatrixLike& m);

class Matrix : public MatrixLike {
 public:
  Matrix(std::vector<double>&& data);
  Matrix(int size = 0);

 public:
  int size() const override;
  static Matrix dot(const MatrixLike& a, const MatrixLike& b);
  Matrix dot(const MatrixLike& other) const;
  static std::vector<double> dot(const MatrixLike& a, const std::vector<double>& vec);
  std::vector<double> dot(const std::vector<double>& vec) const;
  static Matrix inv(const MatrixLike& a);
  Matrix inv() const;

 private:
  double at_impl(int i, int j) const override;
  double& at_impl(int i, int j) override;

 private:
  std::vector<double> data_;
  int size_;
};

class LUDecomoposition {
  class LMatrix : public MatrixLike {
   public:
    LMatrix(Matrix& data);

   public:
    int size() const override;

   private:
    double at_impl(int i, int j) const override;
    double& at_impl(int i, int j) override;

   private:
    Matrix& data_;
  };
  class UMatrix : public MatrixLike {
   public:
    UMatrix(Matrix& data);

   public:
    int size() const override;

   private:
    double at_impl(int i, int j) const override;
    double& at_impl(int i, int j) override;

   private:
    Matrix& data_;
  };

 public:
  LUDecomoposition(const MatrixLike& matrix);

 private:
  Matrix data_;

 public:
  const LMatrix L;
  const UMatrix U;
};