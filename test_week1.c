/**
 * @file test_week1.c
 * @brief Verification and test harness for Phase 1 & 2 of Tropical Algebra Library.
 * @author Embedded Systems & Applied Mathematics
 * @date 2026
 *
 * @details
 * Validates:
 * 1. Matrix multiplication and power in (max, +) semiring.
 * 2. Karp's algorithm host validation (step-by-step state vector D_k recursion).
 * 3. Exact rational lambda evaluation (3/1 and fractional 7/2).
 * 4. Kleene Star metric closure and dependency elimination.
 * 5. Rejection of illegal positive cycles in Kleene Star (causality check).
 * 6. Overflow saturation in scalar tropical multiplication.
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
 * @brief Print a state column vector formatted.
 *
 * @param k Step index.
 * @param vec Array of values.
 * @param n Dimension.
 */
static void print_vector(uint8_t k, const int16_t *vec, uint8_t n)
{
    printf("  D_%u = [ ", k);
    for (uint8_t i = 0; i < n; ++i) {
        if (vec[i] <= TROP_ZERO) {
            printf("%4s ", "eps");
        } else {
            printf("%4d ", vec[i]);
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
    int8_t status_p2 = trop_matpow(&A, 2, &A2);
    int8_t status_p3 = trop_matpow(&A, 3, &A3);
    assert(status_p2 == TROP_SUCCESS);
    assert(status_p3 == TROP_SUCCESS);

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
     */
    printf("===============================================================\n");
    printf("  Karp's Spectral Analysis: Step-by-Step Matrix-Vector Trace   \n");
    printf("===============================================================\n\n");

    int16_t trace_D[4][3];
    trace_D[0][0] = 0;
    trace_D[0][1] = TROP_ZERO;
    trace_D[0][2] = TROP_ZERO;

    print_vector(0, trace_D[0], 3);

    for (uint8_t k = 1; k <= 3; ++k) {
        for (uint8_t i = 0; i < 3; ++i) {
            int16_t acc = TROP_ZERO;
            for (uint8_t j = 0; j < 3; ++j) {
                acc = trop_add(acc, trop_mul(A.data[i][j], trace_D[k - 1][j]));
            }
            trace_D[k][i] = acc;
        }
        print_vector(k, trace_D[k], 3);
    }
    printf("\n");

    /* Verify state vector trajectories */
    /* D_1 = [2, eps, 1]^T */
    assert(trace_D[1][0] == 2 && trace_D[1][1] == TROP_ZERO && trace_D[1][2] == 1);
    /* D_2 = [4, 2, 3]^T */
    assert(trace_D[2][0] == 4 && trace_D[2][1] == 2 && trace_D[2][2] == 3);
    /* D_3 = [6, 5, 5]^T */
    assert(trace_D[3][0] == 6 && trace_D[3][1] == 5 && trace_D[3][2] == 5);
    printf("[PASS] D_k recurrence vectors D_0 -> D_3 verified.\n\n");

    /* Rational table print */
    printf("--- Karp's Rational Ratio Breakdown Table: (D_3(u) - D_k(u)) / (3 - k) ---\n");
    printf("  Node u |  k=0 (len 3)  |  k=1 (len 2)  |  k=2 (len 1)  |  min_k over u \n");
    printf("  -------+---------------+---------------+---------------+---------------\n");
    printf("    0    |  (6-0)/3 = 2  |  (6-2)/2 = 2  |  (6-4)/1 = 2  |      2        \n");
    printf("    1    |      eps      |      eps      |  (5-2)/1 = 3  |      3        \n");
    printf("    2    |      eps      |  (5-1)/2 = 2  |  (5-3)/1 = 2  |      2        \n");
    printf("  -------+---------------+---------------+---------------+---------------\n");
    printf("  Maximum Cycle Mean lambda = max_u { 2, 3, 2 } = 3/1\n\n");

    /* Run the library's trop_karp_lambda function */
    int32_t lambda_num = 0;
    int32_t lambda_den = 1;
    int8_t karp_status = trop_karp_lambda(&A, &lambda_num, &lambda_den);
    printf("trop_karp_lambda(&A) returned: %d/%d (status: %d)\n", lambda_num, lambda_den, karp_status);
    assert(karp_status == TROP_SUCCESS);
    assert(lambda_num == 3 && lambda_den == 1);
    printf("[PASS] Karp's integer eigenvalue assertion (lambda == 3/1) SUCCESSFUL!\n\n");

    /*
     * 4. Validation of Fractional Lambda (Prevents Integer Truncation Loss)
     */
    printf("===============================================================\n");
    printf("  Fractional Eigenvalue Precision Verification                 \n");
    printf("===============================================================\n\n");

    TropMatrix A_frac;
    trop_matrix_init(&A_frac, 2, 2);
    /* 2-cycle with weights 5 and 2: total weight 7, cycle length 2 -> lambda = 7/2 = 3.5 */
    A_frac.data[0][0] = TROP_ZERO; A_frac.data[0][1] = 5;
    A_frac.data[1][0] = 2;         A_frac.data[1][1] = TROP_ZERO;

    int32_t f_num = 0, f_den = 1;
    int8_t f_status = trop_karp_lambda(&A_frac, &f_num, &f_den);
    printf("Fractional pipeline matrix lambda returned: %d/%d (float: %.2f, status: %d)\n",
           f_num, f_den, (float)f_num / (float)f_den, f_status);
    assert(f_status == TROP_SUCCESS);
    assert(f_num == 7 && f_den == 2);
    printf("[PASS] Fractional eigenvalue preserved exactly without integer truncation!\n\n");

    /*
     * 5. Verify Kleene Star Closure and Dependency Elimination
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

    int8_t elim_status = trop_eliminate_dependencies(&A0, &A1, &A_elim);
    assert(elim_status == TROP_SUCCESS);
    print_trop_matrix("A_elim = A0* (x) A1", &A_elim);

    assert(A_elim.data[0][0] == 5 && A_elim.data[0][1] == 6);
    assert(A_elim.data[1][0] == 3 && A_elim.data[1][1] == 4);
    printf("[PASS] Dependency elimination (Kleene Star) verified.\n\n");

    /*
     * 6. Positive Cycle Detection in Kleene Star (Causality Check)
     */
    printf("===============================================================\n");
    printf("  Kleene Star Positive Cycle Rejection (Causality Enforcement) \n");
    printf("===============================================================\n\n");

    TropMatrix A_bad;
    trop_matrix_init(&A_bad, 2, 2);
    /* Internal cycle with positive weight: 0 -> 1 (weight 2), 1 -> 0 (weight 1) => cycle weight +3 */
    A_bad.data[0][0] = TROP_ZERO; A_bad.data[0][1] = 2;
    A_bad.data[1][0] = 1;         A_bad.data[1][1] = TROP_ZERO;

    TropMatrix bad_closure;
    int8_t bad_kleene = trop_kleene_star(&A_bad, &bad_closure);
    printf("Kleene star on positive cycle returned status: %d (expected %d)\n",
           bad_kleene, TROP_ERR_POSITIVE_CYCLE);
    assert(bad_kleene == TROP_ERR_POSITIVE_CYCLE);

    int8_t bad_elim = trop_eliminate_dependencies(&A_bad, &A1, &A_elim);
    assert(bad_elim == TROP_ERR_POSITIVE_CYCLE);
    printf("[PASS] Positive cycle detected and successfully rejected by Kleene star!\n\n");

    /*
     * 7. Scalar Multiplication Overflow Saturation Verification
     */
    printf("===============================================================\n");
    printf("  Scalar Multiplication Saturation Protection                  \n");
    printf("===============================================================\n\n");

    int16_t saturated = trop_mul(20000, 20000);
    printf("trop_mul(20000, 20000) = %d (saturated at INT16_MAX = %d)\n", saturated, 32767);
    assert(saturated == 32767);
    printf("[PASS] Positive arithmetic overflow saturates safely without wraparound!\n\n");

    printf("===============================================================\n");
    printf("  ALL TESTS PASSED! Refactored & Hardened for Bare-Metal ARM. \n");
    printf("===============================================================\n");

    return 0;
}
