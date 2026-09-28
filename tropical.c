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
 * Mathematical Properties Verified in Implementation:
 * - \f$ a \oplus b = \max(a, b) \f$
 * - \f$ a \otimes b = a + b \f$ with absorbing zero \f$ \varepsilon = -32000 \f$
 * - Matrix multiplication complexity: \f$ O(n^3) \f$ operations
 * - Matrix exponentiation complexity: \f$ O(n^3 \log_2 k) \f$ via binary exponentiation
 * - Kleene star closure: \f$ A^* = \bigoplus_{k=0}^{n-1} A^{\otimes k} \f$
 */

#include "tropical.h"

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
    return (int16_t)(a + b);
}

/* ========================================================================== */
/*                           MATRIX TROPICAL OPS                              */
/* ========================================================================== */

void trop_matrix_init(TropMatrix *M, uint8_t r, uint8_t c)
{
    if (M == 0) {
        return;
    }

    if (r > MAX_DIM) {
        r = MAX_DIM;
    }
    if (c > MAX_DIM) {
        c = MAX_DIM;
    }

    M->rows = r;
    M->cols = c;

    for (uint8_t i = 0; i < r; ++i) {
        for (uint8_t j = 0; j < c; ++j) {
            M->data[i][j] = TROP_ZERO;
        }
    }
}

void trop_matrix_identity(TropMatrix *M, uint8_t n)
{
    if (M == 0) {
        return;
    }

    trop_matrix_init(M, n, n);

    for (uint8_t i = 0; i < M->rows; ++i) {
        M->data[i][i] = TROP_ONE; /* Multiplicative identity e = 0 */
    }
}

void trop_matrix_add(const TropMatrix *A, const TropMatrix *B, TropMatrix *C)
{
    if ((A == 0) || (B == 0) || (C == 0)) {
        return;
    }

    /* Use local temporary to safely support in-place aliasing (C == A or C == B) */
    TropMatrix tmp;
    uint8_t r = (A->rows < B->rows) ? A->rows : B->rows;
    uint8_t c = (A->cols < B->cols) ? A->cols : B->cols;

    tmp.rows = r;
    tmp.cols = c;

    for (uint8_t i = 0; i < r; ++i) {
        for (uint8_t j = 0; j < c; ++j) {
            tmp.data[i][j] = trop_add(A->data[i][j], B->data[i][j]);
        }
    }

    *C = tmp;
}

void trop_matmul(const TropMatrix *A, const TropMatrix *B, TropMatrix *C)
{
    if ((A == 0) || (B == 0) || (C == 0)) {
        return;
    }

    /* Dimension validation: cols of A must match rows of B */
    if (A->cols != B->rows) {
        return;
    }

    /* Local stack buffer prevents corruption if C aliases A or B */
    TropMatrix tmp;
    tmp.rows = A->rows;
    tmp.cols = B->cols;

    for (uint8_t i = 0; i < A->rows; ++i) {
        for (uint8_t j = 0; j < B->cols; ++j) {
            int16_t acc = TROP_ZERO;
            for (uint8_t k = 0; k < A->cols; ++k) {
                int16_t term = trop_mul(A->data[i][k], B->data[k][j]);
                acc = trop_add(acc, term);
            }
            tmp.data[i][j] = acc;
        }
    }

    *C = tmp;
}

void trop_matpow(const TropMatrix *A, uint32_t power, TropMatrix *out)
{
    if ((A == 0) || (out == 0) || (A->rows != A->cols)) {
        return;
    }

    uint8_t n = A->rows;

    /* Base case: A^0 = I */
    if (power == 0U) {
        trop_matrix_identity(out, n);
        return;
    }

    /* Binary Exponentiation (O(n^3 * log2(power))) */
    TropMatrix base = *A;
    TropMatrix result;
    trop_matrix_identity(&result, n);

    uint32_t p = power;
    while (p > 0U) {
        if ((p & 1U) != 0U) {
            trop_matmul(&result, &base, &result);
        }
        p >>= 1U;
        if (p > 0U) {
            trop_matmul(&base, &base, &base);
        }
    }

    *out = result;
}

void trop_kleene_star(const TropMatrix *A, TropMatrix *A_star)
{
    if ((A == 0) || (A_star == 0) || (A->rows != A->cols)) {
        return;
    }

    uint8_t n = A->rows;

    /*
     * Kleene Star definition for Max-Plus Semiring:
     * A* = I \oplus A \oplus A^2 \oplus ... \oplus A^(n-1)
     */
    TropMatrix accum;
    trop_matrix_identity(&accum, n); /* A^0 = I */

    if (n <= 1) {
        *A_star = accum;
        return;
    }

    TropMatrix current_power = *A; /* A^1 */
    trop_matrix_add(&accum, &current_power, &accum); /* I \oplus A */

    for (uint8_t k = 2; k < n; ++k) {
        trop_matmul(&current_power, A, &current_power);  /* A^k = A^(k-1) \otimes A */
        trop_matrix_add(&accum, &current_power, &accum); /* accum = accum \oplus A^k */
    }

    *A_star = accum;
}

void trop_eliminate_dependencies(const TropMatrix *A0, const TropMatrix *A1, TropMatrix *A)
{
    if ((A0 == 0) || (A1 == 0) || (A == 0)) {
        return;
    }

    /*
     * In Discrete Event Systems (Timed Event Graphs):
     * Implicit equation: x(k) = A0 (x) x(k) (+) A1 (x) x(k-1)
     * Resolved form:     x(k) = (A0* (x) A1) (x) x(k-1)
     * Therefore:         A    = A0* (x) A1
     */
    TropMatrix A0_star;
    trop_kleene_star(A0, &A0_star);
    trop_matmul(&A0_star, A1, A);
}
