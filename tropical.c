/**
 * @file tropical.c
 * @brief Implementation of Max-Plus (Tropical) Semiring core operations.
 * @author Embedded Systems & Applied Mathematics
 * @date 2026
 *
 * @details
 * Bare-metal implementation compliant with MISRA C guidelines and tailored
 * for ARM Cortex-M4 (STM32F407) deployment.
 *
 * Memory Architecture:
 * - Scratch matrices are allocated as file-scope static variables in .bss to
 *   prevent stack overflow on resource-constrained embedded systems (such as
 *   the default 1 KB stack in Keil MDK startup files).
 * - Thread-safety / Reentrancy: Functions are non-reentrant. Safe for single-threaded
 *   bare-metal execution; must NOT be called concurrently or from ISRs.
 *
 * Mathematical Properties Verified in Implementation:
 * - a (+) b = max(a, b)
 * - a (x) b = a + b with absorbing zero epsilon = -32000 and positive saturation
 * - Matrix multiplication complexity: O(n^3)
 * - Matrix exponentiation complexity: O(n^3 log2 k) via binary exponentiation
 * - Kleene star closure: A* = \bigoplus_{k=0}^{n-1} A^k with positive cycle verification
 */

#include "tropical.h"

/* ========================================================================== */
/*                         STATIC SCRATCHPAD BUFFERS                          */
/* ========================================================================== */

/* File-scope static scratch buffers replace large stack temporaries */
static TropMatrix s_add_tmp;
static TropMatrix s_matmul_tmp;
static TropMatrix s_matpow_base;
static TropMatrix s_matpow_result;
static TropMatrix s_kleene_accum;
static TropMatrix s_kleene_power;
static TropMatrix s_elim_A0_star;

/* ========================================================================== */
/*                           SCALAR TROPICAL OPS                              */
/* ========================================================================== */

int16_t trop_add(int16_t a, int16_t b)
{
    /* In (max, +) semiring, addition is the maximum operator */
    return (a > b) ? a : b;
}

int16_t trop_mul(int16_t a, int16_t b)
{
    /*
     * In (max, +) semiring, multiplication is classical addition.
     * TROP_ZERO (-32000) acts as -\infty (absorbing element):
     * -\infty + x = -\infty for all x.
     */
    if ((a <= TROP_ZERO) || (b <= TROP_ZERO)) {
        return TROP_ZERO;
    }

    int32_t sum = (int32_t)a + (int32_t)b;

    /*
     * Positive saturation: Prevents signed 16-bit integer overflow (wraparound)
     * from producing false negative weights that corrupt max-plus calculations.
     */
    if (sum > INT16_MAX) {
        return INT16_MAX;
    }

    /* Headroom protection against values dropping to/below TROP_ZERO */
    if (sum <= TROP_ZERO) {
        return TROP_ZERO;
    }

    return (int16_t)sum;
}

/* ========================================================================== */
/*                           MATRIX TROPICAL OPS                              */
/* ========================================================================== */

int8_t trop_matrix_init(TropMatrix *M, uint8_t r, uint8_t c)
{
    if (M == 0) {
        return TROP_ERR_NULL_PTR;
    }

    if ((r == 0) || (r > MAX_DIM) || (c == 0) || (c > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    M->rows = r;
    M->cols = c;

    for (uint8_t i = 0; i < r; ++i) {
        for (uint8_t j = 0; j < c; ++j) {
            M->data[i][j] = TROP_ZERO;
        }
    }

    return TROP_SUCCESS;
}

int8_t trop_matrix_identity(TropMatrix *M, uint8_t n)
{
    if (M == 0) {
        return TROP_ERR_NULL_PTR;
    }

    if ((n == 0) || (n > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    int8_t status = trop_matrix_init(M, n, n);
    if (status != TROP_SUCCESS) {
        return status;
    }

    for (uint8_t i = 0; i < n; ++i) {
        M->data[i][i] = TROP_ONE; /* Multiplicative identity e = 0 */
    }

    return TROP_SUCCESS;
}

int8_t trop_matrix_add(const TropMatrix *A, const TropMatrix *B, TropMatrix *C)
{
    if ((A == 0) || (B == 0) || (C == 0)) {
        return TROP_ERR_NULL_PTR;
    }

    if ((A->rows != B->rows) || (A->cols != B->cols)) {
        return TROP_ERR_INVALID_DIM;
    }

    uint8_t r = A->rows;
    uint8_t c = A->cols;

    if ((r == 0) || (r > MAX_DIM) || (c == 0) || (c > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    s_add_tmp.rows = r;
    s_add_tmp.cols = c;

    for (uint8_t i = 0; i < r; ++i) {
        for (uint8_t j = 0; j < c; ++j) {
            s_add_tmp.data[i][j] = trop_add(A->data[i][j], B->data[i][j]);
        }
    }

    *C = s_add_tmp;
    return TROP_SUCCESS;
}

int8_t trop_matmul(const TropMatrix *A, const TropMatrix *B, TropMatrix *C)
{
    if ((A == 0) || (B == 0) || (C == 0)) {
        return TROP_ERR_NULL_PTR;
    }

    /* Dimension validation: cols of A must match rows of B */
    if (A->cols != B->rows) {
        return TROP_ERR_INVALID_DIM;
    }

    if ((A->rows == 0) || (A->rows > MAX_DIM) ||
        (A->cols == 0) || (A->cols > MAX_DIM) ||
        (B->cols == 0) || (B->cols > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    s_matmul_tmp.rows = A->rows;
    s_matmul_tmp.cols = B->cols;

    for (uint8_t i = 0; i < A->rows; ++i) {
        for (uint8_t j = 0; j < B->cols; ++j) {
            int16_t acc = TROP_ZERO;
            for (uint8_t k = 0; k < A->cols; ++k) {
                int16_t term = trop_mul(A->data[i][k], B->data[k][j]);
                acc = trop_add(acc, term);
            }
            s_matmul_tmp.data[i][j] = acc;
        }
    }

    *C = s_matmul_tmp;
    return TROP_SUCCESS;
}

int8_t trop_matpow(const TropMatrix *A, uint32_t power, TropMatrix *out)
{
    if ((A == 0) || (out == 0)) {
        return TROP_ERR_NULL_PTR;
    }

    if ((A->rows == 0) || (A->rows != A->cols) || (A->rows > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    uint8_t n = A->rows;

    /* Base case: A^0 = I */
    if (power == 0U) {
        return trop_matrix_identity(out, n);
    }

    /* Binary Exponentiation (O(n^3 * log2(power))) using static buffers */
    s_matpow_base = *A;
    trop_matrix_identity(&s_matpow_result, n);

    uint32_t p = power;
    while (p > 0U) {
        if ((p & 1U) != 0U) {
            trop_matmul(&s_matpow_result, &s_matpow_base, &s_matpow_result);
        }
        p >>= 1U;
        if (p > 0U) {
            trop_matmul(&s_matpow_base, &s_matpow_base, &s_matpow_base);
        }
    }

    *out = s_matpow_result;
    return TROP_SUCCESS;
}

int8_t trop_kleene_star(const TropMatrix *A, TropMatrix *A_star)
{
    if ((A == 0) || (A_star == 0)) {
        return TROP_ERR_NULL_PTR;
    }

    if ((A->rows == 0) || (A->rows != A->cols) || (A->rows > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    uint8_t n = A->rows;

    /* Base case: 1x1 matrix */
    if (n == 1) {
        if (A->data[0][0] > 0) {
            return TROP_ERR_POSITIVE_CYCLE;
        }
        trop_matrix_identity(A_star, 1);
        A_star->data[0][0] = trop_add(TROP_ONE, A->data[0][0]);
        return TROP_SUCCESS;
    }

    /*
     * Kleene Star definition for Max-Plus Semiring:
     * A* = I \oplus A \oplus A^2 \oplus ... \oplus A^(n-1)
     *
     * Convergence check:
     * A* exists if and only if there are no cycles of positive weight.
     * A positive cycle of length k manifests as a positive diagonal entry:
     * (A^k)_{i,i} > 0. If detected, abort and return error.
     */
    trop_matrix_identity(&s_kleene_accum, n); /* A^0 = I */
    s_kleene_power = *A;                      /* A^1 */

    /* Verify A^1 has no positive self-loops */
    for (uint8_t i = 0; i < n; ++i) {
        if (s_kleene_power.data[i][i] > 0) {
            return TROP_ERR_POSITIVE_CYCLE;
        }
    }

    /* Accumulate A^1: I \oplus A */
    trop_matrix_add(&s_kleene_accum, &s_kleene_power, &s_kleene_accum);

    /* Compute powers A^2 through A^n */
    for (uint8_t k = 2; k <= n; ++k) {
        trop_matmul(&s_kleene_power, A, &s_kleene_power); /* A^k = A^(k-1) (x) A */

        /* Check diagonal for positive cycles of length k */
        for (uint8_t i = 0; i < n; ++i) {
            if (s_kleene_power.data[i][i] > 0) {
                return TROP_ERR_POSITIVE_CYCLE;
            }
        }

        /* Kleene closure accumulates up to power n-1 */
        if (k < n) {
            trop_matrix_add(&s_kleene_accum, &s_kleene_power, &s_kleene_accum);
        }
    }

    *A_star = s_kleene_accum;
    return TROP_SUCCESS;
}

int8_t trop_eliminate_dependencies(const TropMatrix *A0, const TropMatrix *A1, TropMatrix *A)
{
    if ((A0 == 0) || (A1 == 0) || (A == 0)) {
        return TROP_ERR_NULL_PTR;
    }

    /*
     * In Discrete Event Systems (Timed Event Graphs):
     * Implicit equation: x(k) = A0 (x) x(k) (+) A1 (x) x(k-1)
     * Resolved form:     x(k) = (A0* (x) A1) (x) x(k-1)
     * Therefore:         A    = A0* (x) A1
     *
     * A0 must not contain positive cycles (physical causality constraint).
     */
    int8_t status = trop_kleene_star(A0, &s_elim_A0_star);
    if (status != TROP_SUCCESS) {
        return status;
    }

    return trop_matmul(&s_elim_A0_star, A1, A);
}
