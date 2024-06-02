#ifndef SRC_matrix_OOP_H_
#define SRC_matrix_OOP_H_

#include <cmath>
#include <iostream>

class Matrix {
 public:
  static constexpr double EPS = 1e-7;

  Matrix();   // Базовый конструктор
  ~Matrix();  // Деструктор
  Matrix(int rows, int cols);  // Параметризированный конструктор
  Matrix(const Matrix& other);  // Конструктор копирования
  Matrix(Matrix&& other);  // Конструктор переноса

  // Математические операции
  void SumMatrix(const Matrix& other);  // сумма матриц
  bool EqMatrix(const Matrix& other);   // равенство матриц
  void SubMatrix(const Matrix& other);  // вычитание матриц
  void MulNumber(const double num);  // умножение матрицы на число
  void MulMatrix(const Matrix& other);  // умножение на матрицу
  double Determinant();         // Возвращает детерминант
  Matrix Transpose();        // транспонирование матриц
  Matrix CalcComplements();  // Возвращает матрицу алгебраических дополнений
  Matrix InverseMatrix();    // обратная матрица

  // Операторы
  bool operator==(const Matrix& other);
  Matrix& operator=(const Matrix& other);
  Matrix operator+(const Matrix& other);
  Matrix operator-(const Matrix& other);
  Matrix operator*(const Matrix& other);
  void operator+=(const Matrix& other);  // увеличивает на аргумент
  void operator-=(const Matrix& other);
  void operator*=(const Matrix& other);
  void operator*=(const double num);
  double& operator()(int i, int j);  // Оператор индексирования
  friend Matrix operator*(const double num, const Matrix& other);
  friend Matrix operator*(const Matrix& other, const double num);

  // Вспомогательные функции
  double GetElement(int i, int j) const;
  void SetElement(int i, int j, double num);
  int GetRows() const;
  int GetCols() const;
  void SetCols(int cols);
  void SetRows(int rows);
  friend void SetAllElements(Matrix* Src, double num);

 private:
  int rows_, cols_;
  double** matrix_;
  void MakeMatrix();  // Запускается в конструкторах
  Matrix Minor(int row, int column);  // Часть детерминанта
  void CopyMatrix(Matrix& Res, const Matrix& Src);
};

#endif  // SRC_matrix_OOP_H_
