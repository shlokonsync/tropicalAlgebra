/**
 * @file test_week1.c
 * @brief Verification and test harness for Phase 1 & 2 of Tropical Algebra Library.
 * @author Embedded Systems & Applied Mathematics
 * @date 2026
 *
 * @details
 * Validates:
 * 1. Matrix multiplication and power in (max, +) semiring.
 * 2. Karp's algorithm host validation (step-by-step state vector D_k recursion via trop_matmul).
 * 3. Karp's cycle mean rational table evaluation yielding lambda = 3.
 * 4. Kleene Star metric closure and dependency elimination.
 */

#include <stdio.h>
#include <assert.h>
#include "tropical.h"

/**
 * @brief Nicely format and print a tropical matrix to stdout.
 *
 * @param name Identifier string.
 * @param M Pointer to matrix.
 */
static void print_trop_matrix(const char *name, const TropMatrix *M)
{
    printf("=== Matrix %s (%ux%u) ===\n", name, M->rows, M->cols);
    for (uint8_t i = 0; i < M->rows; ++i) {
        printf("  [ ");
        for (uint8_t j = 0; j < M->cols; ++j) {
            if (M->data[i][j] <= TROP_ZERO) {
                printf("%6s ", "eps"); /* eps / -inf */
            } else {
                printf("%6d ", M->data[i][j]);
            }
        }
        printf("]\n");
    }
    printf("\n");
}

/**
 * @brief Print a state column vector D_k formatted.
 *
 * @param k Step index.
 * @param D Pointer to column vector (n x 1).
 */
static void print_state_vector(uint8_t k, const TropMatrix *D)
{
    printf("  D_%u = [ ", k);
    for (uint8_t i = 0; i < D->rows; ++i) {
        if (D->data[i][0] <= TROP_ZERO) {
            printf("%4s ", "eps");
        } else {
            printf("%4d ", D->data[i][0]);
        }
    }
    printf("]^T\n");
}

int main(void)
{
    printf("===============================================================\n");
    printf("  STM32F407 Tropical (Max-Plus) Algebra Library - Test Suite   \n");
    printf("===============================================================\n\n");

    /*
     * 1. Define Test Matrix A (3x3):
     * Row 0: [ 2,  1,  -inf ]
     * Row 1: [ -inf,  3,  1 ]
     * Row 2: [ 1,  2,  2 ]
     */
    TropMatrix A;
    trop_matrix_init(&A, 3, 3);
    A.data[0][0] = 2;         A.data[0][1] = 1;         A.data[0][2] = TROP_ZERO;
    A.data[1][0] = TROP_ZERO; A.data[1][1] = 3;         A.data[1][2] = 1;
    A.data[2][0] = 1;         A.data[2][1] = 2;         A.data[2][2] = 2;

    print_trop_matrix("A (Original)", &A);

    /*
     * 2. Compute Matrix Powers: A^2 and A^3
     */
    TropMatrix A2, A3;
    trop_matpow(&A, 2, &A2);
    trop_matpow(&A, 3, &A3);

    print_trop_matrix("A^2 (A (x) A)", &A2);
    print_trop_matrix("A^3 (A^2 (x) A)", &A3);

    /* Verify A^2 analytical values */
    assert(A2.data[0][0] == 4 && A2.data[0][1] == 4 && A2.data[0][2] == 2);
    assert(A2.data[1][0] == 2 && A2.data[1][1] == 6 && A2.data[1][2] == 4);
    assert(A2.data[2][0] == 3 && A2.data[2][1] == 5 && A2.data[2][2] == 4);
    printf("[PASS] A^2 calculation verified.\n");

    /* Verify A^3 analytical values */
    assert(A3.data[0][0] == 6 && A3.data[0][1] == 7 && A3.data[0][2] == 5);
    assert(A3.data[1][0] == 5 && A3.data[1][1] == 9 && A3.data[1][2] == 7);
    assert(A3.data[2][0] == 5 && A3.data[2][1] == 8 && A3.data[2][2] == 6);
    printf("[PASS] A^3 calculation verified.\n\n");

    /*
     * 3. Karp's Algorithm Host Validation: Step-by-Step D_k Recursion
     *
     * The D_k recursion is repeated matrix-vector products:
     * D_k = A (x) D_{k-1} using trop_matmul
     */
    printf("===============================================================\n");
    printf("  Karp's Spectral Analysis: Step-by-Step Matrix-Vector Trace   \n");
    printf("===============================================================\n\n");

    TropMatrix D[4];
    /* D_0 = [0, eps, eps]^T */
    trop_matrix_init(&D[0], 3, 1);
    D[0].data[0][0] = 0;

    print_state_vector(0, &D[0]);

    /* Compute D_1, D_2, D_3 via trop_matmul */
    for (uint8_t k = 1; k <= 3; ++k) {
        trop_matmul(&A, &D[k - 1], &D[k]);
        print_state_vector(k, &D[k]);
    }
    printf("\n");

    /* Verify state vector trajectories */
    /* D_1 = [2, eps, 1]^T */
    assert(D[1].data[0][0] == 2 && D[1].data[1][0] == TROP_ZERO && D[1].data[2][0] == 1);
    /* D_2 = [4, 2, 3]^T */
    assert(D[2].data[0][0] == 4 && D[2].data[1][0] == 2 && D[2].data[2][0] == 3);
    /* D_3 = [6, 5, 5]^T */
    assert(D[3].data[0][0] == 6 && D[3].data[1][0] == 5 && D[3].data[2][0] == 5);
    printf("[PASS] D_k recurrence vectors D_0 -> D_3 verified.\n\n");

    /*
     * Rational cycle mean evaluation breakdown:
     * ratio(u, k) = (D_n(u) - D_k(u)) / (n - k)
     */
    printf("--- Karp's Rational Ratio Breakdown Table: (D_3(u) - D_k(u)) / (3 - k) ---\n");
    printf("  Node u |  k=0 (len 3)  |  k=1 (len 2)  |  k=2 (len 1)  |  min_k over u \n");
    printf("  -------+---------------+---------------+---------------+---------------\n");
    printf("    0    |  (6-0)/3 = 2  |  (6-2)/2 = 2  |  (6-4)/1 = 2  |      2        \n");
    printf("    1    |      eps      |      eps      |  (5-2)/1 = 3  |      3        \n");
    printf("    2    |      eps      |  (5-1)/2 = 2  |  (5-3)/1 = 2  |      2        \n");
    printf("  -------+---------------+---------------+---------------+---------------\n");
    printf("  Maximum Cycle Mean lambda = max_u { 2, 3, 2 } = 3\n\n");

    /* Run the library's trop_karp_lambda function */
    int16_t lambda = trop_karp_lambda(&A);
    printf("trop_karp_lambda(&A) returned: %d\n", lambda);
    assert(lambda == 3);
    printf("[PASS] Karp's eigenvalue assertion (lambda == 3) SUCCESSFUL!\n\n");

    /*
     * 4. Verify Kleene Star Closure and Dependency Elimination
     */
    printf("===============================================================\n");
    printf("  Kleene Star Closure & Discrete Event Dependency Elimination \n");
    printf("===============================================================\n\n");

    TropMatrix A0, A1, A_elim;
    trop_matrix_init(&A0, 2, 2);
    trop_matrix_init(&A1, 2, 2);

    /* A0: Direct internal dependency without positive cycles */
    A0.data[0][0] = TROP_ZERO; A0.data[0][1] = 2;
    A0.data[1][0] = TROP_ZERO; A0.data[1][1] = TROP_ZERO;

    /* A1: Event transitions */
    A1.data[0][0] = 5;         A1.data[0][1] = TROP_ZERO;
    A1.data[1][0] = 3;         A1.data[1][1] = 4;

    trop_eliminate_dependencies(&A0, &A1, &A_elim);
    print_trop_matrix("A_elim = A0* (x) A1", &A_elim);

    assert(A_elim.data[0][0] == 5 && A_elim.data[0][1] == 6);
    assert(A_elim.data[1][0] == 3 && A_elim.data[1][1] == 4);
    printf("[PASS] Dependency elimination (Kleene Star) verified.\n\n");

    printf("===============================================================\n");
    printf("  ALL TESTS PASSED! Production-ready for STM32F407 deployment. \n");
    printf("===============================================================\n");

    return 0;
}
