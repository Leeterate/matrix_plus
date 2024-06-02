#include <gtest/gtest.h>

#include "matrix_oop.h"

void CompareAllElems(Matrix test, double value) {
  for (int i = 0; i < test.GetRows(); i++) {
    for (int j = 0; j < test.GetCols(); j++) {
      EXPECT_NEAR(test(i, j), value, test.EPS);
    }
  }
}

TEST(testing_constructors, regular) {
  Matrix test;
  CompareAllElems(test, 0.00000001);
}

TEST(testing_constructors, pametrized) {
  Matrix test(9, 9);

  CompareAllElems(test, 0.00000001);
  EXPECT_EQ(test.GetCols(), 9);
  EXPECT_EQ(test.GetRows(), 9);
}

TEST(testing_constructors, copy) {
  Matrix test_init(5, 5);
  SetAllElements(&test_init, 2);

  Matrix test(test_init);

  CompareAllElems(test, 2);
  EXPECT_EQ(test.GetCols(), 5);
  EXPECT_EQ(test.GetRows(), 5);
}

TEST(testing_constructors, rvalue_mov_const) {
  Matrix test_init(3, 3);
  SetAllElements(&test_init, 2);

  Matrix test(std::move(test_init));

  CompareAllElems(test, 2);
  EXPECT_EQ(test.GetCols(), 3);
  EXPECT_EQ(test.GetRows(), 3);
}

TEST(testing_operators, Eq_Matrix_value_pos) {
  Matrix test, test2;

  SetAllElements(&test, 2.00000001);
  SetAllElements(&test2, 2);
  EXPECT_TRUE(test.EqMatrix(test2));
}

TEST(testing_operators, Eq_Matrix_size) {
  Matrix test(3, 3), test2(3, 5);

  SetAllElements(&test, 2.00000001);
  SetAllElements(&test2, 2);
  EXPECT_FALSE(test.EqMatrix(test2));
}

TEST(testing_operators, Eq_Matrix_value_neg) {
  Matrix test, test2;

  SetAllElements(&test, 2.00000001);
  SetAllElements(&test2, 2.000001);
  EXPECT_FALSE(test.EqMatrix(test2));
}

TEST(testing_operators, Sum_Matrix) {
  Matrix test, test2;

  SetAllElements(&test, 2.00000001);
  test.SumMatrix(test);
  SetAllElements(&test2, 4);
  EXPECT_TRUE(test.EqMatrix(test2));
}

TEST(testing_operators, Sub_Matrix) {
  Matrix test, test2;

  SetAllElements(&test, 6.00000001);
  SetAllElements(&test2, 3);
  test.SubMatrix(test2);
  EXPECT_TRUE(test.EqMatrix(test2));
}

TEST(testing_operators, mul_number) {
  Matrix test, test2;

  SetAllElements(&test, 4.00000001);
  SetAllElements(&test2, 20);
  test.MulNumber(5);
  EXPECT_TRUE(test.EqMatrix(test2));
}

TEST(testing_operators, Mul_Matrix) {
  Matrix test(3, 2), test2(2, 3), res(3, 3);

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 4);
  test.SetElement(1, 0, 2);
  test.SetElement(1, 1, 5);
  test.SetElement(2, 0, 3);
  test.SetElement(2, 1, 6);

  test2.SetElement(0, 0, 1);
  test2.SetElement(0, 1, -1);
  test2.SetElement(0, 2, 1);
  test2.SetElement(1, 0, 2);
  test2.SetElement(1, 1, 3);
  test2.SetElement(1, 2, 4);

  res.SetElement(0, 0, 9);
  res.SetElement(0, 1, 11);
  res.SetElement(0, 2, 17);
  res.SetElement(1, 0, 12);
  res.SetElement(1, 1, 13);
  res.SetElement(1, 2, 22);
  res.SetElement(2, 0, 15);
  res.SetElement(2, 1, 15);
  res.SetElement(2, 2, 27);

  test.MulMatrix(test2);
  EXPECT_TRUE(test.EqMatrix(res));
}

TEST(testing_operators, transpose) {
  Matrix test1(2, 3), test2(3, 2);

  test1.SetElement(0, 0, 1);
  test1.SetElement(0, 1, 2);
  test1.SetElement(0, 2, 3);
  test1.SetElement(1, 0, 4);
  test1.SetElement(1, 1, 5);
  test1.SetElement(1, 2, 6);

  test2.SetElement(0, 0, 1);
  test2.SetElement(0, 1, 4);
  test2.SetElement(1, 0, 2);
  test2.SetElement(1, 1, 5);
  test2.SetElement(2, 0, 3);
  test2.SetElement(2, 1, 6);

  Matrix r = test1.Transpose();
  EXPECT_TRUE(r.EqMatrix(test2));
}

TEST(testing_operators, complements) {
  Matrix test(3, 3), res(3, 3);

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 2);
  test.SetElement(0, 2, 3);
  test.SetElement(1, 0, 0);
  test.SetElement(1, 1, 4);
  test.SetElement(1, 2, 2);
  test.SetElement(2, 0, 5);
  test.SetElement(2, 1, 2);
  test.SetElement(2, 2, 1);

  res.SetElement(0, 0, 0);
  res.SetElement(0, 1, 10);
  res.SetElement(0, 2, -20);
  res.SetElement(1, 0, 4);
  res.SetElement(1, 1, -14);
  res.SetElement(1, 2, 8);
  res.SetElement(2, 0, -8);
  res.SetElement(2, 1, -2);
  res.SetElement(2, 2, 4);

  Matrix r = test.CalcComplements();
  EXPECT_TRUE(r.EqMatrix(res));
}

TEST(testing_operators, determinant) {
  Matrix test(3, 3);
  double det = 0.0;

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 2);
  test.SetElement(0, 2, 3);
  test.SetElement(1, 0, 0);
  test.SetElement(1, 1, 4);
  test.SetElement(1, 2, 2);
  test.SetElement(2, 0, 5);
  test.SetElement(2, 1, 2);
  test.SetElement(2, 2, 1);

  det = test.Determinant();
  EXPECT_TRUE(det - (-40) < test.EPS);
}

TEST(testing_operators, Inverse_Matrix) {
  Matrix test(3, 3), res(3, 3);

  test.SetElement(0, 0, 2);
  test.SetElement(0, 1, 6);
  test.SetElement(0, 2, 3);
  test.SetElement(1, 0, 4);
  test.SetElement(1, 1, -1);
  test.SetElement(1, 2, 3);
  test.SetElement(2, 0, 1);
  test.SetElement(2, 1, 3);
  test.SetElement(2, 2, 2);

  res.SetElement(0, 0, 0.8461538462);
  res.SetElement(0, 1, 0.2307692308);
  res.SetElement(0, 2, -1.6153846154);
  res.SetElement(1, 0, 0.3846153846);
  res.SetElement(1, 1, -0.0769230769);
  res.SetElement(1, 2, -0.4615384615);
  res.SetElement(2, 0, -1);
  res.SetElement(2, 1, 0);
  res.SetElement(2, 2, 2);

  Matrix r = test.InverseMatrix();

  EXPECT_TRUE(r.EqMatrix(res));
}

TEST(testing_operators, eqation) {
  Matrix test(3, 3);
  Matrix test2(3, 3);

  SetAllElements(&test, 9.2222);
  SetAllElements(&test2, 9.2222);
  EXPECT_TRUE(test == test2);
}

TEST(testing_operators, plus) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 3);
  SetAllElements(&res, 6);
  Matrix r = test + test;
  EXPECT_TRUE(r == res);
}
TEST(testing_operators, plus_small) {
  Matrix test(1, 1), res(1, 1);

  SetAllElements(&test, 11);
  SetAllElements(&res, 22);
  Matrix r = test + test;
  EXPECT_TRUE(r == res);
}

TEST(testing_operators, mul_num) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 7);
  SetAllElements(&res, 28);
  Matrix r = test * 4;
  EXPECT_TRUE(r == res);
}
TEST(testing_operators, mulnum_backw) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 7);
  SetAllElements(&res, 28);
  Matrix r = 4 * test;
  EXPECT_TRUE(r == res);
}

TEST(testing_operators, minus) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 14);
  SetAllElements(&res, 7);
  Matrix r = test - res;
  EXPECT_TRUE(r == res);
}

TEST(testing_operators, minus_new) {
  Matrix test(3, 3), test2(3, 3), res(3, 3);

  SetAllElements(&test, 6);
  SetAllElements(&test2, 7);
  SetAllElements(&res, -1);
  Matrix r = test - test2;
  EXPECT_TRUE(r == res);
}

TEST(testing_operators, multiply) {
  Matrix test(3, 2), test2(2, 3), res(3, 3);

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 4);
  test.SetElement(1, 0, 2);
  test.SetElement(1, 1, 5);
  test.SetElement(2, 0, 3);
  test.SetElement(2, 1, 6);

  test2.SetElement(0, 0, 1);
  test2.SetElement(0, 1, -1);
  test2.SetElement(0, 2, 1);
  test2.SetElement(1, 0, 2);
  test2.SetElement(1, 1, 3);
  test2.SetElement(1, 2, 4);

  res.SetElement(0, 0, 9);
  res.SetElement(0, 1, 11);
  res.SetElement(0, 2, 17);
  res.SetElement(1, 0, 12);
  res.SetElement(1, 1, 13);
  res.SetElement(1, 2, 22);
  res.SetElement(2, 0, 15);
  res.SetElement(2, 1, 15);
  res.SetElement(2, 2, 27);

  Matrix r = test * test2;
  EXPECT_TRUE(r.EqMatrix(res));
}

TEST(testing_operators, plus_eq) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 9);
  SetAllElements(&res, 18);
  test += test;
  EXPECT_TRUE(test == res);
}

TEST(testing_operators, minus_eq) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 14);
  SetAllElements(&res, 7);
  test -= res;
  EXPECT_TRUE(test == res);
}

TEST(testing_operators, multiply_eq_num) {
  Matrix test(3, 3), res(3, 3);

  SetAllElements(&test, 3);
  SetAllElements(&res, 12);
  test *= 4;
  EXPECT_TRUE(test == res);
}

TEST(testing_operators, copy) {
  Matrix test(3, 3), test2(3, 3), res(3, 3);

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 4);
  res = test;

  EXPECT_TRUE(test == res);
}

TEST(testing_operators, multiply_eq) {
  Matrix test(3, 2), test2(2, 3), res(3, 3);

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 4);
  test.SetElement(1, 0, 2);
  test.SetElement(1, 1, 5);
  test.SetElement(2, 0, 3);
  test.SetElement(2, 1, 6);

  test2.SetElement(0, 0, 1);
  test2.SetElement(0, 1, -1);
  test2.SetElement(0, 2, 1);
  test2.SetElement(1, 0, 2);
  test2.SetElement(1, 1, 3);
  test2.SetElement(1, 2, 4);

  res.SetElement(0, 0, 9);
  res.SetElement(0, 1, 11);
  res.SetElement(0, 2, 17);
  res.SetElement(1, 0, 12);
  res.SetElement(1, 1, 13);
  res.SetElement(1, 2, 22);
  res.SetElement(2, 0, 15);
  res.SetElement(2, 1, 15);
  res.SetElement(2, 2, 27);

  test *= test2;
  EXPECT_TRUE(test.EqMatrix(res));
}

TEST(testing_operators, indexation) {
  Matrix test(4, 4);

  SetAllElements(&test, 32.128);

  EXPECT_TRUE(test(0, 0) - 32.128 < test.EPS);
  EXPECT_TRUE(test(1, 0) - 32.128 < test.EPS);
  EXPECT_TRUE(test(0, 2) - 32.128 < test.EPS);
  EXPECT_TRUE(test(2, 0) - 32.128 < test.EPS);
  EXPECT_TRUE(test(1, 1) - 32.128 < test.EPS);
  EXPECT_TRUE(test(2, 2) - 32.128 < test.EPS);
  EXPECT_TRUE(test(3, 3) - 32.128 < test.EPS);
}

TEST(testing_operators, assign) {
  Matrix test(3, 2);

  test.SetElement(0, 0, 1);
  test.SetElement(0, 1, 70);
  test.SetElement(1, 0, 22);
  test.SetElement(1, 1, 5);
  test.SetElement(2, 0, 9);
  test.SetElement(2, 1, 16);

  Matrix test2 = test;
  EXPECT_TRUE(test.EqMatrix(test2));
}

TEST(testing_accessors, GetRows_GetCols) {
  Matrix test(3, 6);

  EXPECT_TRUE(test.GetCols() - 6 < test.EPS);
  EXPECT_TRUE(test.GetRows() - 3 < test.EPS);
}

TEST(testing_mutators, SetElement_GetElement) {
  Matrix test(3, 3);

  test.SetElement(2, 2, 16.5);
  EXPECT_TRUE(test.GetElement(2, 2) - 16.5 < test.EPS);
}

TEST(testing_mutators, SetRows_SetCols) {
  Matrix test(3, 3);

  test.SetRows(10);
  test.SetCols(5);
  SetAllElements(&test, 39.3);

  EXPECT_TRUE(test.GetCols() - 5 == 0);
  EXPECT_TRUE(test.GetRows() - 10 == 0);

  EXPECT_TRUE(test.GetElement(0, 4) - 39.3 < test.EPS);
  EXPECT_TRUE(test.GetElement(9, 4) - 39.3 < test.EPS);

  Matrix test2(4, 4);
  SetAllElements(&test2, 2);
  test2.SetRows(2);
  test2.SetCols(2);
  SetAllElements(&test2, 3);
  EXPECT_TRUE(test2.GetElement(0, 1) == 3);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
