#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const int size = 3;

int main(int argc, char *argv[]) {
    // 2*x1 + 1*x2 + 0*x3 = 4
    // 1*x1 + 2*x2 + 1*x3 = 7
    // 0*x1 + 1*x2 + 2*x3 = 1
    double A[3][3] = {
        {2, 1, 0},
        {1, 2, 1},
        {0, 1, 2}   
    };
    double F[3] = {4, 7, 1};
    double x[3] = {0, 0, 0};
    
    int i, j, k;
    int index1, index2, offset;
    double alpha, gamma;
    
    /* cycle reduction */
    for (i = 0; i < log2(size + 1) - 1; i++) {
        for (j = pow(2, i + 1) - 1; j < size; j = j + pow(2, i + 1)) {
            offset = pow(2, i);
            index1 = j - offset;
            index2 = j + offset;

            alpha = A[j][index1] / A[index1][index1];
            gamma = A[j][index2] / A[index2][index2];

            for (k = 0; k < size; k++) {
                A[j][k] -= (alpha * A[index1][k] + gamma * A[index2][k]);
            }
            
            F[j] -= (alpha * F[index1] + gamma * F[index2]);
        }
    }

    /* back substitution */
    int index = (size - 1) / 2;
    x[index] = F[index] / A[index][index];
    
    for (i = log2(size + 1) - 2; i >= 0; i--) {
        for (j = pow(2, i + 1) - 1; j < size; j = j + pow(2, i + 1)) {
            offset = pow(2, i);
            index1 = j - offset;
            index2 = j + offset;

            x[index1] = F[index1];
            x[index2] = F[index2];

            for (k = 0; k < size; k++) {
                if (k != index1) 
                    x[index1] -= A[index1][k] * x[k];
                if (k != index2)
                    x[index2] -= A[index2][k] * x[k];
            }

            x[index1] = x[index1] / A[index1][index1];
            x[index2] = x[index2] / A[index2][index2];
        }
    }

    for (i = 0; i < size; i++) {
        printf("x%d = %e\n", i, x[i]);
    }
    return 0;
}
