#include <ravi_matrixlib.h>

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

int test_dgemm() {
  /* Tests here taken from GSL so GPL license applies */
  const ravi_matrix_ops_t *ops = ravi_matrix_get_implementation();
  {
    bool transA = false;
    bool transB = true;
    double alpha = -0.3;
    double beta = 1;
    double A[] = { -0.795, 0.81, 0.388, 0.09 };
    double B[] = { -0.847, 0.031, -0.938, 0.09, -0.286, -0.478, -0.981, 0.881 };
    double C[] = { -0.242, -0.02 };
    double C_expected[] = { -0.1562981, -0.0026243 };
    ops->multiply(1, 4, A, 2, 4, B, 1, 2, C, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", C[i], C_expected[i]);
      }
    }
  }

  {
    bool transA = true;
    bool transB = false;
    int order = 102;
    double alpha = -1;
    double beta = 0;
    double A[] = { -0.358, 0.224, -0.941, 0.513 };
    double B[] = { -0.201, -0.159, -0.586, -0.016, -0.324, 0.411, 0.115, -0.229 };
    double C[] = { 0.558, 0.596 };
    double C_expected[] = { -0.57956, 0.017636 };
    ops->multiply(4, 1, A, 4, 2, B, 1, 2, C, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", C[i], C_expected[i]);
      }
    }
  }

  {
    bool transA = true;
    bool transB = true;
    double alpha = -0.3;
    double beta = 1;
    double A[] = { -0.164, 0.522, 0.948, -0.624 };
    double B[] = { -0.142, 0.778, 0.359, 0.622, -0.637, -0.757, -0.282, -0.805 };
    double C[] = { -0.09, 0.183 };
    double C_expected[] = { -0.0248334, 0.1884672 };
    ops->multiply(4, 1, A, 2, 4, B, 1, 2, C, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", C[i], C_expected[i]);
      }
    }
  }

  {
    bool transA = false;
    bool transB = false;
    double alpha = 0;
    double beta = 0;
    double A[] = { 0.571, 0.081, 0.109, 0.988 };
    double B[] = { -0.048, -0.753, -0.8, -0.89, -0.535, -0.017, -0.018, -0.544 };
    double C[] = { -0.876, -0.792 };
    double C_expected[] = { 0.0, 0.0 };
    ops->multiply(1, 4, A, 4, 2, B, 1, 2, C, transA, transB, alpha, beta);
    {
      int i;
      for (i = 0; i < 2; i++) {
        printf("Got %f expected %f\n", C[i], C_expected[i]);
      }
    }
  }

  return 0;
}

int main(int argc, const char *argv[]) {
  int rc = test_dgemm();
	return rc != 0 ? 1 : 0;
}