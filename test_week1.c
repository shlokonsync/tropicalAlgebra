/**
 * @file test_week1.c
 * @brief Verification and test harness for Phase 1 & 2 of Tropical Algebra Library.
 * @author Embedded Systems & Applied Mathematics
 * @date 2026
 *
 * @details
 * Validates:
 * 1. Matrix multiplication and power in (max, +) semiring.
 * 2. Kleene Star metric closure.
 * 3. Karp's algorithm for tropical eigenvalue (maximum cycle mean / asymptotic throughput).
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
    /*
     * Row 0: [ max(4, eps, eps)=4, max(3, 4, eps)=4, max(eps, 2, eps)=2 ]
     * Row 1: [ max(eps, eps, 2)=2, max(eps, 6, 3)=6, max(eps, 4, 3)=4 ]
     * Row 2: [ max(3, eps, 3)=3,   max(2, 5, 4)=5,   max(eps, 3, 4)=4 ]
     */
    assert(A2.data[0][0] == 4 && A2.data[0][1] == 4 && A2.data[0][2] == 2);
    assert(A2.data[1][0] == 2 && A2.data[1][1] == 6 && A2.data[1][2] == 4);
    assert(A2.data[2][0] == 3 && A2.data[2][1] == 5 && A2.data[2][2] == 4);
    printf("[PASS] A^2 calculation verified.\n");

    /* Verify A^3 analytical values */
    /*
     * Row 0: [ 6, 7, 5 ]
     * Row 1: [ 5, 9, 7 ]
     * Row 2: [ 5, 8, 6 ]
     */
    assert(A3.data[0][0] == 6 && A3.data[0][1] == 7 && A3.data[0][2] == 5);
    assert(A3.data[1][0] == 5 && A3.data[1][1] == 9 && A3.data[1][2] == 7);
    assert(A3.data[2][0] == 5 && A3.data[2][1] == 8 && A3.data[2][2] == 6);
    printf("[PASS] A^3 calculation verified.\n\n");

    /*
     * 3. Compute Tropical Eigenvalue (lambda) via Karp's Algorithm
     */
    printf("Running Karp's Maximum Cycle Mean Algorithm...\n");
    int16_t lambda = trop_karp_lambda(&A);
    printf("Computed Tropical Eigenvalue lambda = %d\n", lambda);

    /* Assert requirement: lambda must be 3 (governed by self-loop at node 1 with weight 3) */
    assert(lambda == 3);
    printf("[PASS] Assertion lambda == 3 SUCCESSFUL!\n\n");

    /*
     * 4. Verify Kleene Star Closure and Dependency Elimination
     */
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

    /*
     * A0* = I (+) A0 = [ [0, 2], [eps, 0] ]
     * A_elim = A0* (x) A1 = [ [max(5, 5), max(eps, 6)], [3, 4] ] = [ [5, 6], [3, 4] ]
     */
    assert(A_elim.data[0][0] == 5 && A_elim.data[0][1] == 6);
    assert(A_elim.data[1][0] == 3 && A_elim.data[1][1] == 4);
    printf("[PASS] Dependency elimination (Kleene Star) verified.\n\n");

    printf("===============================================================\n");
    printf("  ALL TESTS PASSED! Production-ready for STM32F407 deployment. \n");
    printf("===============================================================\n");

    return 0;
}
