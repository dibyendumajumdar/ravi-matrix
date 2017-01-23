#include <ravi_matrixlib.h>

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

int test_dgemm() {
  /* Tests here taken from GSL so GPL license applies */
  {
    bool transA = false;
    bool transB = true;
    double alpha = -0.3;
    double beta = 1;
    double A[] = {-0.795, 0.81, 0.388, 0.09};
    double B[] = {-0.847, 0.031, -0.938, 0.09, -0.286, -0.478, -0.981, 0.881};
    double C[] = {-0.242, -0.02};
    double C_expected[] = {-0.1562981, -0.0026243};
    ravi_matrix_t AA = {1, 4, A};
    ravi_matrix_t BB = {2, 4, B};
    ravi_matrix_t CC = {1, 2, C};
    ravi_matrix_multiply(&AA, &BB, &CC, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", CC.data[i], C_expected[i]);
      }
    }
  }

  {
    bool transA = true;
    bool transB = false;
    double alpha = -1;
    double beta = 0;
    double A[] = {-0.358, 0.224, -0.941, 0.513};
    double B[] = {-0.201, -0.159, -0.586, -0.016, -0.324, 0.411, 0.115, -0.229};
    double C[] = {0.558, 0.596};
    double C_expected[] = {-0.57956, 0.017636};
    ravi_matrix_t AA = {4, 1, A};
    ravi_matrix_t BB = {4, 2, B};
    ravi_matrix_t CC = {1, 2, C};
    ravi_matrix_multiply(&AA, &BB, &CC, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", CC.data[i], C_expected[i]);
      }
    }
  }

  {
    bool transA = true;
    bool transB = true;
    double alpha = -0.3;
    double beta = 1;
    double A[] = {-0.164, 0.522, 0.948, -0.624};
    double B[] = {-0.142, 0.778, 0.359, 0.622, -0.637, -0.757, -0.282, -0.805};
    double C[] = {-0.09, 0.183};
    double C_expected[] = {-0.0248334, 0.1884672};
    ravi_matrix_t AA = {4, 1, A};
    ravi_matrix_t BB = {2, 4, B};
    ravi_matrix_t CC = {1, 2, C};
    ravi_matrix_multiply(&AA, &BB, &CC, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", CC.data[i], C_expected[i]);
      }
    }
  }

  {
    bool transA = false;
    bool transB = false;
    double alpha = 0;
    double beta = 0;
    double A[] = {0.571, 0.081, 0.109, 0.988};
    double B[] = {-0.048, -0.753, -0.8, -0.89, -0.535, -0.017, -0.018, -0.544};
    double C[] = {-0.876, -0.792};
    double C_expected[] = {0.0, 0.0};
    ravi_matrix_t AA = {1, 4, A};
    ravi_matrix_t BB = {4, 2, B};
    ravi_matrix_t CC = {1, 2, C};
    ravi_matrix_multiply(&AA, &BB, &CC, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", CC.data[i], C_expected[i]);
      }
    }
  }

  return 0;
}

int test_outer_product() {
  double x[5] = {1, 2, 3, 4, 5};
  double y[5] = {6, 7, 8, 9, 10};
  double A[25] = {0};
  double expected[25] = {6.0,  12.0, 18.0, 24.0, 30.0, 7.0,  14.0, 21.0, 28.0,
                         35.0, 8.0,  16.0, 24.0, 32.0, 40.0, 9.0,  18.0, 27.0,
                         36.0, 45.0, 10.0, 20.0, 30.0, 40.0, 50.0};
  ravi_vector_outer_product(5, x, 5, y, A, 1.0);
  for (int i = 0; i < 25; i++) {
    if (expected[i] != A[i]) {
      fprintf(stderr, "mismatch at %d\n", i);
      return 1;
    }
  }

  double z1[2] = {71, 1};
  double z2[2] = {0, 45};
  double M[4] = {0};
  ravi_vector_outer_product(2, z1, 2, z2, M, 1.0);
  double expected2[4] = {0.0, 0.0, 3195.0, 45.0};
  for (int i = 0; i < 4; i++) {
    if (expected2[i] != M[i]) {
      fprintf(stderr, "mismatch at %d in matrix M\n", i);
      return 1;
    }
  }

  printf("outer_product OK\n");
  return 0;
}

int main(void) {
  int rc = test_dgemm();
  rc += test_outer_product();
  return rc != 0 ? 1 : 0;
}