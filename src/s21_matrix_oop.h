#pragma once
class S21Matrix {
 public:
  S21Matrix() noexcept;
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;

  ~S21Matrix();

  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const double num) const;
  S21Matrix operator*(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other) const noexcept;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double num) noexcept;
  S21Matrix& operator*=(const S21Matrix& other);
  friend S21Matrix operator*(const double num, const S21Matrix& other);

  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(const int rows);
  void SetCols(const int cols);

  double& operator()(const int i, const int j) const;

  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  bool EqMatrix(const S21Matrix& other) const noexcept;
  S21Matrix Transpose() const noexcept;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;
  S21Matrix GetMinorMatrix(const int row, const int col) const;

 private:
  void DeleteMatrix();

  int rows_, cols_;
  double** matrix_;
};
