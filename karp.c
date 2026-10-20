/**
 * @file karp.c
 * @brief Karp's Minimum/Maximum Cycle Mean Algorithm for Tropical Semirings.
 * @author Embedded Systems & Applied Mathematics
 * @date 2026
 *
 * @details
 * Calculates the maximum cycle mean (\f$ \lambda \f$, the tropical eigenvalue)
 * of a Timed Event Graph (precedence matrix \f$ A \f$).
 *
 * Mathematical Formulation:
 * -------------------------
 * In classical graph theory, Karp's algorithm (1978) finds the minimum cycle mean.
 * In Max-Plus Algebra, with operations \f$ (\max, +) \f$, it computes the
 * MAXIMUM cycle mean:
 * \f[
 *   \lambda(A) = \max_{u \in V} \min_{0 \le k \le n-1} \frac{D_n(u) - D_k(u)}{n - k}
 * \f]
 * where:
 * - \f$ n \f$ is the number of states/nodes (\f$ n \le \text{MAX\_DIM} \f$).
 * - \f$ D_0 = [0, \varepsilon, \dots, \varepsilon]^T \f$ is the initial state vector.
 * - \f$ D_k = A \otimes D_{k-1} \f$ for \f$ 1 \le k \le n \f$.
 *
 * Preconditions:
 * --------------
 * Matrix @c A must be irreducible (strongly connected precedence graph).
 * When running from initial node 0, all cycles must be reachable from node 0.
 *
 * Embedded Safety:
 * ----------------
 * - Zero dynamic memory: state vectors D[k] are stored in a file-scope static 2D array
 *   (544 bytes in .bss). Stack footprint is under 48 bytes (safe for 1 KB Keil stack).
 * - Non-reentrant: Do not call concurrently or from an ISR.
 * - Exact rational arithmetic: All comparisons cross-multiply integer ratios to avoid
 *   floating-point or integer division rounding loss. Returns the exact reduced fraction.
 */

#include "tropical.h"

/* ========================================================================== */
/*                         STATIC VECTOR BUFFERS                              */
/* ========================================================================== */

/*
 * State trajectories: D[k][u] is the state value of node u after k steps.
 * Storing only 1D vectors consumes (17 * 16 * 2) = 544 bytes in .bss,
 * eliminating the 8.8 KB stack allocation of 17 full TropMatrix structs.
 */
static int16_t D[MAX_DIM + 1][MAX_DIM];

/* ========================================================================== */
/*                           HELPER FUNCTIONS                                 */
/* ========================================================================== */

/**
 * @brief Greatest Common Divisor (Euclidean algorithm) for fraction reduction.
 */
static int32_t calc_gcd(int32_t a, int32_t b)
{
    if (a < 0) {
        a = -a;
    }
    if (b < 0) {
        b = -b;
    }
    while (b != 0) {
        int32_t t = b;
        b = a % b;
        a = t;
    }
    return (a == 0) ? 1 : a;
}

/* ========================================================================== */
/*                         ALGORITHM IMPLEMENTATION                           */
/* ========================================================================== */

int8_t trop_karp_lambda(const TropMatrix *A, int32_t *lambda_num, int32_t *lambda_den)
{
    if ((A == 0) || (lambda_num == 0) || (lambda_den == 0)) {
        return TROP_ERR_NULL_PTR;
    }

    if ((A->rows == 0) || (A->rows != A->cols) || (A->rows > MAX_DIM)) {
        return TROP_ERR_INVALID_DIM;
    }

    const uint8_t n = A->rows;

    /*
     * Step 1: Initialize D_0 = [0, \varepsilon, \dots, \varepsilon]^T
     * Node 0 is the designated source vertex. (Requires A to be irreducible).
     */
    D[0][0] = TROP_ONE; /* 0 */
    for (uint8_t u = 1; u < n; ++u) {
        D[0][u] = TROP_ZERO;
    }

    /*
     * Step 2: Compute D_k = A (x) D_{k-1} for k = 1 ... n
     * Direct matrix-vector max-plus product avoids calling trop_matmul with
     * large matrix wrappers on the stack.
     */
    for (uint8_t k = 1; k <= n; ++k) {
        for (uint8_t i = 0; i < n; ++i) {
            int16_t acc = TROP_ZERO;
            for (uint8_t j = 0; j < n; ++j) {
                int16_t term = trop_mul(A->data[i][j], D[k - 1][j]);
                acc = trop_add(acc, term);
            }
            D[k][i] = acc;
        }
    }

    /*
     * Step 3: Evaluate lambda = max_u min_{0 <= k <= n-1} (D_n(u) - D_k(u)) / (n - k)
     *
     * Exact rational tracking: cross-multiplication avoids integer division rounding.
     */
    int32_t max_num = 0;
    int32_t max_den = 1;
    uint8_t has_valid_cycle = 0;

    for (uint8_t u = 0; u < n; ++u) {
        /* Node u must be reachable at step n */
        if (D[n][u] <= TROP_ZERO) {
            continue;
        }

        int32_t min_num_u = 0;
        int32_t min_den_u = 1;
        uint8_t has_valid_k = 0;

        for (uint8_t k = 0; k < n; ++k) {
            /* Node u must be reachable at step k */
            if (D[k][u] <= TROP_ZERO) {
                continue;
            }

            int32_t diff = (int32_t)D[n][u] - (int32_t)D[k][u];
            int32_t den  = (int32_t)n - (int32_t)k; /* Always > 0 since k < n */

            if (!has_valid_k) {
                min_num_u = diff;
                min_den_u = den;
                has_valid_k = 1;
            } else {
                /*
                 * Rational comparison: (diff / den) < (min_num_u / min_den_u)
                 * Equivalent to: diff * min_den_u < min_num_u * den
                 */
                if ((diff * min_den_u) < (min_num_u * den)) {
                    min_num_u = diff;
                    min_den_u = den;
                }
            }
        }

        if (has_valid_k) {
            if (!has_valid_cycle) {
                max_num = min_num_u;
                max_den = min_den_u;
                has_valid_cycle = 1;
            } else {
                /*
                 * Rational comparison: (min_num_u / min_den_u) > (max_num / max_den)
                 * Equivalent to: min_num_u * max_den > max_num * min_den_u
                 */
                if ((min_num_u * max_den) > (max_num * min_den_u)) {
                    max_num = min_num_u;
                    max_den = min_den_u;
                }
            }
        }
    }

    if (!has_valid_cycle) {
        return TROP_ERR_NO_CYCLE;
    }

    /* Standardize signs so denominator is strictly positive */
    if (max_den < 0) {
        max_num = -max_num;
        max_den = -max_den;
    }

    /* Reduce fraction to lowest canonical terms */
    int32_t g = calc_gcd(max_num, max_den);
    *lambda_num = max_num / g;
    *lambda_den = max_den / g;

    return TROP_SUCCESS;
}
