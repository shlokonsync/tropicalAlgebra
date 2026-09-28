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
 * Embedded Safety:
 * - Zero dynamic memory (all state vectors \f$ D_0 \dots D_n \f$ stored on stack).
 * - Exact rational comparison via cross-multiplication (\f$ \frac{a}{b} < \frac{c}{d} \iff a \cdot d < c \cdot b \f$)
 *   to avoid precision loss in embedded integer math.
 * - Strictly no standard library dependencies (<stdio.h>, <stdlib.h>, <assert.h>).
 */

#include "tropical.h"

int16_t trop_karp_lambda(const TropMatrix *A)
{
    if ((A == 0) || (A->rows == 0) || (A->rows != A->cols)) {
        return TROP_ZERO;
    }

    const uint8_t n = A->rows;

    /*
     * State trajectories: D[k][u] is the maximum weight path of length k to node u.
     * Dimensions: (MAX_DIM + 1) x MAX_DIM of int16_t
     * Memory consumption: 17 * 16 * 2 = 544 bytes on stack.
     */
    int16_t D[MAX_DIM + 1][MAX_DIM];

    /* Step 1: Initialize D_0 = [0, TROP_ZERO, ..., TROP_ZERO]^T */
    D[0][0] = TROP_ONE; /* 0 */
    for (uint8_t u = 1; u < n; ++u) {
        D[0][u] = TROP_ZERO;
    }

    /* Step 2: Iteratively compute D_k = A (x) D_{k-1} for k = 1 to n */
    for (uint8_t k = 1; k <= n; ++k) {
        for (uint8_t i = 0; i < n; ++i) {
            int16_t max_val = TROP_ZERO;
            for (uint8_t j = 0; j < n; ++j) {
                int16_t term = trop_mul(A->data[i][j], D[k - 1][j]);
                max_val = trop_add(max_val, term);
            }
            D[k][i] = max_val;
        }
    }

    /*
     * Step 3: Evaluate lambda = max_u min_{0 <= k <= n-1} (D_n(u) - D_k(u)) / (n - k)
     *
     * We track best rational numbers using integer numerators and denominators
     * to eliminate division rounding artifacts during comparison.
     */
    int32_t max_num = 0;
    int32_t max_den = 1;
    uint8_t has_valid_cycle = 0;

    for (uint8_t u = 0; u < n; ++u) {
        /* Node must be reachable at step n */
        if (D[n][u] <= TROP_ZERO) {
            continue;
        }

        int32_t min_num_u = 0;
        int32_t min_den_u = 1;
        uint8_t has_valid_k = 0;

        for (uint8_t k = 0; k < n; ++k) {
            /* Node must be reachable at step k */
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
        return TROP_ZERO;
    }

    /*
     * Integer evaluation of cycle time:
     * Standard integer division (or rounded fixed-point).
     */
    return (int16_t)(max_num / max_den);
}
