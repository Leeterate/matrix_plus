#include "matrix_oop.h"

Matrix::Matrix() {
  cols_ = 4;
  rows_ = 4;
  MakeMatrix();
}

Matrix::~Matrix() {
  if (this->matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
}
// конструктор с параметрами
Matrix::Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0)
    throw std::invalid_argument("Impossible operation: dimention of 0 or less");
  cols_ = cols;
  rows_ = rows;
  MakeMatrix();
}
// конструктор копирования
Matrix::Matrix(const Matrix& other)
    : Matrix::Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] = other.matrix_[i][j];
    }
  }
}
// конструктор переноса
Matrix::Matrix(Matrix&& other)
    : Matrix::Matrix(other.rows_, other.cols_) {
  std::swap(this->rows_, other.rows_);
  std::swap(this->cols_, other.cols_);
  std::swap(this->matrix_, other.matrix_);
}

void Matrix::SumMatrix(const Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Impossible operation: different size");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

bool Matrix::EqMatrix(const Matrix& other) {
  bool res = true;
  if (this->cols_ != other.cols_ || this->rows_ != other.rows_) {
    res = false;
  } else {
    for (int i = 0; i < this->rows_; i++) {
      if (res != true) {
        break;
      }
      for (int j = 0; j < this->cols_; j++) {
        if (fabs(this->matrix_[i][j] - other.matrix_[i][j]) > EPS) {
          res = false;
        }
      }
    }
  }
  return res;
}

void Matrix::SubMatrix(const Matrix& other) {
  if (this->cols_ != other.cols_ || this->rows_ != other.rows_)
    throw std::invalid_argument("Impossible operation: different size");
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_; j++) {
      this->matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= (double)num;
    }
  }
}

void Matrix::MulMatrix(const Matrix& other) {
  Matrix res(this->rows_, other.cols_);
  if (this->cols_ != other.rows_)
    throw std::invalid_argument(
        "Impossible operation: inappropriate sizes of arguments");
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      for (int k = 0; k < this->cols_; k++) {
        res.matrix_[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  SetCols(other.cols_);
  CopyMatrix(*this, res);
}

Matrix Matrix::Transpose() {
  Matrix res(this->cols_, this->rows_);
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      res.matrix_[i][j] = this->matrix_[j][i];
    }
  }
  return res;
}

Matrix Matrix::CalcComplements() {
  Matrix res(this->rows_, this->cols_);
  if (this->cols_ != this->rows_)
    throw std::invalid_argument("Impossible operation: unequal dimentions");
  for (int i = 0; i < res.rows_; i++) {
    for (int j = 0; j < res.cols_; j++) {
      Matrix tmp = Minor(i, j);
      double det = tmp.Determinant();
      res.matrix_[i][j] = det * (((i + j) % 2) ? -1.0 : 1.0);
    }
  }
  return res;
}

double Matrix::Determinant() {
  double res = 0.0;
  if (cols_ != rows_)
    throw std::invalid_argument("Impossible operation: unequal dimentions");
  if (rows_ < 3) {
    if (rows_ == 1) res = matrix_[0][0];
    if (rows_ == 2)
      res = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  }
  if (rows_ >= 3) {
    for (int i = 0; i < cols_; i++) {
      Matrix temp = Minor(0, i);
      double det = temp.Determinant();
      res += (pow(-1, i)) * matrix_[0][i] * det;
    }
  }
  return res;
}

Matrix Matrix::InverseMatrix() {
  if (this->rows_ != this->cols_ || this->Determinant() == 0)
    throw std::invalid_argument("Determinant = 0 or unequal dimentions");
  Matrix temp = this->CalcComplements();
  Matrix res = temp.Transpose();
  res.MulNumber(1.0 / this->Determinant());
  return res;
}

bool Matrix::operator==(const Matrix& other) {
  return this->EqMatrix(other);
}

Matrix& Matrix::operator=(const Matrix& other) {
  if (this != &other) CopyMatrix(*this, other);
  return *this;
}

Matrix Matrix::operator+(const Matrix& other) {
  Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

Matrix Matrix::operator-(const Matrix& other) {
  Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

Matrix Matrix::operator*(const Matrix& other) {
  Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

Matrix operator*(const double num, const Matrix& other) {
  Matrix res(other);
  res.MulNumber(num);
  return res;
}

Matrix operator*(const Matrix& other, const double num) {
  Matrix res(other);
  res.MulNumber(num);
  return res;
}

void Matrix::operator+=(const Matrix& other) { this->SumMatrix(other); }

void Matrix::operator-=(const Matrix& other) { this->SubMatrix(other); }

void Matrix::operator*=(const Matrix& other) { this->MulMatrix(other); }

void Matrix::operator*=(const double num) { this->MulNumber(num); }

double& Matrix::operator()(int i, int j) {
  if (i < 0 || i > rows_ || j < 0 || j > cols_)
    throw std::invalid_argument(
        "Input contains element with nonexistent index");
  return matrix_[i][j];
}

double Matrix::GetElement(int i, int j) const {
  if (i < 0 || i > rows_ || j < 0 || j > cols_)
    throw std::invalid_argument(
        "Input contains element with nonexistent index");
  return matrix_[i][j];
}

void Matrix::SetElement(int i, int j, double num) {
  if (i < 0 || i > rows_ || j < 0 || j > cols_)
    throw std::invalid_argument(
        "Input contains element with nonexistent index");
  matrix_[i][j] = num;
}

int Matrix::GetRows() const { return rows_; }

int Matrix::GetCols() const { return cols_; }

void Matrix::SetRows(int Rows) {
  if (this->rows_ != Rows) {
    Matrix tmp(Rows, this->cols_);
    if (Rows > this->rows_) {
      for (int i = 0; i < this->rows_; i++) {
        for (int j = 0; j < this->cols_; j++) {
          tmp.matrix_[i][j] = this->matrix_[i][j];
        }
      }
    } else if (Rows < this->rows_) {
      for (int i = 0; i < Rows; i++) {
        for (int j = 0; j < this->cols_; j++) {
          tmp.matrix_[i][j] = this->matrix_[i][j];
        }
      }
    }
    *this = std::move(tmp);
  }
  return;
}

void Matrix::SetCols(int cols) {
  if (this->cols_ != cols) {
    Matrix tmp(this->rows_, cols);
    if (cols > this->cols_) {
      for (int i = 0; i < this->rows_; i++) {
        for (int j = 0; j < this->cols_; j++) {
          tmp.matrix_[i][j] = this->matrix_[i][j];
        }
      }
    } else if (cols < this->cols_) {
      for (int i = 0; i < this->rows_; i++) {
        for (int j = 0; j < cols; j++) {
          tmp.matrix_[i][j] = this->matrix_[i][j];
        }
      }
    }
    *this = std::move(tmp);
  }
  return;
}

void Matrix::MakeMatrix() {
  matrix_ = new double*[rows_];
  if (matrix_) {
    for (int i = 0; i < rows_; i++) {
      matrix_[i] = new double[cols_];
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = 0;
      }
    }
  }
}

Matrix Matrix::Minor(int Irows, int Icols) {
  Matrix res((this->rows_ - 1), (this->cols_ - 1));
  int Rows2 = 0, Cols2 = 0;
  for (int i = 0; i < this->rows_; i++) {
    if (Irows == i) continue;
    for (int j = 0; j < this->cols_; j++) {
      if (Icols == j) continue;
      res.matrix_[Rows2][Cols2] = this->matrix_[i][j];
      Cols2++;
    }
    Rows2++;
    Cols2 = 0;
  }
  return res;
}

void Matrix::CopyMatrix(Matrix& Res, const Matrix& Src) {
  if (Res.matrix_) {
    for (int i = 0; i < Res.rows_; i++) delete[] Res.matrix_[i];
    delete[] Res.matrix_;
    Res.matrix_ = nullptr;
  }
  Res.rows_ = Src.rows_;
  Res.cols_ = Src.cols_;
  Res.MakeMatrix();
  for (int i = 0; i < Res.rows_; i += 1) {
    for (int j = 0; j < Res.cols_; j += 1) {
      Res.matrix_[i][j] = Src.matrix_[i][j];
    }
  }
}

void SetAllElements(Matrix* Src, double num) {
  for (int i = 0; i < Src->GetRows(); i++) {
    for (int j = 0; j < Src->GetCols(); j++) {
      Src->SetElement(i, j, num);
    }
  }
}
