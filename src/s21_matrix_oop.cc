#include "s21_matrix_oop.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

using namespace std::literals;

S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ < 0 || cols_ < 0) {
    throw std::invalid_argument("rows/cols <= 0!"s);
  }

  try {
    matrix_ = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_]{};
    }
  }
  catch (const std::bad_alloc& e) {
    throw;
  }
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  rows_ = other.rows_;
  cols_ = other.cols_;

  try {
    matrix_ = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  catch (const std::bad_alloc& e) {
    throw;
  }  
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
  rows_ = std::move(other.rows_);
  cols_ = std::move(other.cols_);
  matrix_ = std::exchange(other.matrix_, nullptr);
}

S21Matrix::~S21Matrix() {
  DeleteMatrix();
  rows_ = 0;
  cols_ = 0;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("diffrent rows/cols counts"s);
  }
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result(i, j) = matrix_[i][j] + other(i, j);
    }
  }
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("diffrent rows/cols counts"s);
  }
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result(i, j) = matrix_[i][j] - other(i, j);
    }
  }
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  if (cols_ != other.rows_) {
    throw std::logic_error(
        "the number of columns of the first matrix is not equal to the number of rows of the second matrix"s);
  }
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < result.rows_; ++i) {
    for (int j = 0; j < result.cols_; ++j) {
      for (int k = 0; k < other.rows_; ++k)
        result(i, j) += matrix_[i][k] * other(k, j);
    }
  }
  return result;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result(i, j) = matrix_[i][j] * num;
    }
  }
  return result;
}

S21Matrix operator*(const double num, const S21Matrix& other) {
  S21Matrix result(other.rows_, other.cols_);
  for (int i = 0; i < other.rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      result(i, j) = other.matrix_[i][j] * num;
    }
  }
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) const noexcept {
  bool result = true;
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    result = false;
  } else {
    for (int i = 0; i < rows_ && result; ++i) {
      for (int j = 0; j < cols_ && result; ++j) {
        if (std::abs(matrix_[i][j] - other(i, j)) > 1e-7) {
          result = false;
        }
      }
    }
  }
  return result;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;
  }
  DeleteMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;

  try {
    matrix_ = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other(i, j);
      }
    }
  }
  catch (const std::bad_alloc& e) {
    throw;
  }  
  
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  DeleteMatrix();
  rows_ = std::move(other.rows_);
  cols_ = std::move(other.cols_);
  matrix_ = std::exchange(other.matrix_, nullptr);
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("diffrent rows/cols counts"s);
  }
  *this = *this + other;
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("diffrent rows/cols counts"s);
  }
  *this = *this - other;
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) noexcept {
  *this = *this * num;
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::logic_error(
        "the number of columns of the first matrix is not equal to the number of rows of the second matrix"s);
  }
  *this = *this * other;
  return *this;
}

int S21Matrix::GetRows() const noexcept { return rows_; }

int S21Matrix::GetCols() const noexcept { return cols_; }

void S21Matrix::SetRows(const int rows) {
  if (rows <= 0) {
    throw std::invalid_argument("rows <= 0!"s);
  }
  S21Matrix temp(rows, cols_);
  int min_rows = std::min(rows, rows_);
  for (int i = 0; i < min_rows; ++i) {
    for (int j = 0; j < cols_; ++j) {
      temp(i, j) = matrix_[i][j];
    }
  }
  *this = std::move(temp);
  rows_ = rows;
}

void S21Matrix::SetCols(const int cols) {
  if (cols <= 0) {
    throw std::invalid_argument("cols <= 0!"s);
  }
  S21Matrix temp(rows_, cols);
  int min_cols = std::min(cols, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < min_cols; ++j) {
      temp(i, j) = matrix_[i][j];
    }
  }
  *this = std::move(temp);
  cols_ = cols;
}

double& S21Matrix::operator()(const int i, const int j) const {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("index outside the matrix"s);
  }
  return matrix_[i][j];
}

void S21Matrix::SumMatrix(const S21Matrix& other) { *this += other; }

void S21Matrix::SubMatrix(const S21Matrix& other) { *this -= other; }

void S21Matrix::MulNumber(const double num) { *this *= num; }

void S21Matrix::MulMatrix(const S21Matrix& other) { *this *= other; }

bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  return *this == other;
}

S21Matrix S21Matrix::Transpose() const noexcept {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result(j, i) = matrix_[i][j];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::logic_error("the matrix is not square"s);
  } else if (rows_ == 0) {
    throw std::out_of_range("Matrix is 1x1"s);
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result(0, 0) = 1;
  } else {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        int k = (i + j) % 2 ? -1 : 1;
        result(i, j) = GetMinorMatrix(i, j).Determinant() * k;
      }
    }
  }
  return result;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::logic_error("the matrix is not square"s);
  }
  double result;
  if (rows_ == 1)
    result = matrix_[0][0];
  else if (rows_ == 2)
    result = matrix_[0][0] * matrix_[1][1] - matrix_[1][0] * matrix_[0][1];
  else {
    result = 0;
    for (int i = 0; i < rows_; ++i) {
      double k = ((i + cols_ - 1) % 2 ? -1. : 1.);
      double val = GetMinorMatrix(i, cols_ - 1).Determinant();
      result += k * matrix_[i][cols_ - 1] * val;
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  double det = Determinant();
  if (fabs(det) < 1e-6) {
    throw std::logic_error("the determinant of the matrix is 0"s);
  }
  return CalcComplements().Transpose() * (1. / det);
}

S21Matrix S21Matrix::GetMinorMatrix(const int row, const int col) const {
  if (row < 0 || col < 0) {
    throw std::out_of_range("row/col < 0!");
  } else if (rows_ != cols_) {
    throw std::logic_error("the matrix is not square"s);
  }
  S21Matrix result(rows_ - 1, cols_ - 1);
  for (int i = 0; i < result.rows_; ++i) {
    for (int j = 0; j < result.cols_; ++j) {
      int t_row = i >= row ? 1 : 0;
      int t_col = j >= col ? 1 : 0;
      result(i, j) = matrix_[i + t_row][j + t_col];
    }
  }
  return result;
}

void S21Matrix::DeleteMatrix() {
  if (matrix_) {
    for (int i = 0; i < rows_; ++i) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}