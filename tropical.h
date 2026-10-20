/**
 * @file tropical.h
 * @brief Bare-metal Max-Plus (Tropical) Semiring Algebra Library.
 * @author Embedded Systems & Applied Mathematics
 * @date 2026
 *
 * @details
 * This header defines the core data structures and operations for the
 * (max, +) tropical semiring, denoted as \f$ \mathbb{R}_{\max} = (\mathbb{R} \cup \{-\infty\}, \oplus, \otimes) \f$.
 *
 * Mathematical Foundations:
 * -------------------------
 * 1. Tropical Addition (\f$ \oplus \f$):
 *    \f[ a \oplus b \triangleq \max(a, b) \f]
 *    Identity element: \f$ \varepsilon = -\infty \f$ (represented by TROP_ZERO = -32000).
 *    Idempotency: \f$ a \oplus a = a \f$.
 *
 * 2. Tropical Multiplication (\f$ \otimes \f$):
 *    \f[ a \otimes b \triangleq a + b \f]
 *    Identity element: \f$ e = 0 \f$.
 *    Zero element: \f$ a \otimes \varepsilon = \varepsilon \otimes a = \varepsilon \f$.
 *    Positive saturation: Saturated at INT16_MAX (32767) to prevent signed overflow wraparound.
 *
 * 3. Matrix Multiplication (\f$ C = A \otimes B \f$):
 *    \f[ C_{i,j} = \bigoplus_{k=1}^{n} (A_{i,k} \otimes B_{k,j}) = \max_{1 \le k \le n} (A_{i,k} + B_{k,j}) \f]
 *
 * Target: STM32F407 (ARM Cortex-M4F, ARMv7E-M).
 * Memory Architecture:
 * - Zero dynamic memory allocation (no malloc/free).
 * - File-scope static scratch buffers (.bss) replace large local matrices, keeping stack
 *   consumption under 96 bytes per frame (fully safe for default 1 KB Keil MDK stacks).
 * - Non-reentrant: Designed for bare-metal main loop execution; DO NOT call from ISRs.
 * Core Dependencies: <stdint.h> only.
 */

#ifndef TROPICAL_H
#define TROPICAL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/*                              RETURN STATUS CODES                           */
/* ========================================================================== */

#define TROP_SUCCESS               ((int8_t)0)   /**< Operation completed successfully */
#define TROP_ERR_NULL_PTR          ((int8_t)-1)  /**< Null pointer passed to function */
#define TROP_ERR_INVALID_DIM       ((int8_t)-2)  /**< Dimension mismatch or exceeds MAX_DIM */
#define TROP_ERR_NO_CYCLE          ((int8_t)-3)  /**< No valid cycle found in precedence graph */
#define TROP_ERR_POSITIVE_CYCLE    ((int8_t)-4)  /**< Positive cycle detected (divergent Kleene star) */
#define TROP_ERR_OVERFLOW          ((int8_t)-5)  /**< Theoretical cumulative weight exceeds dynamic range */
#define TROP_ERR_NOT_REACHABLE     ((int8_t)-6)  /**< Graph not strongly connected / node unreachable from node 0 */

/* ========================================================================== */
/*                            CONSTANTS & CONFIG                              */
/* ========================================================================== */

/**
 * @brief Maximum supported matrix dimension (rows and columns).
 * @note 16x16 matrix of int16_t consumes 512 bytes of data payload.
 */
#define MAX_DIM 16

/**
 * @brief Tropical Zero (\f$ \varepsilon = -\infty \f$).
 *
 * @details Chosen as -32000 rather than INT16_MIN (-32768) to provide negative
 * headroom (~767 units). This prevents arithmetic underflow when negative edge
 * weights or offsets are added during intermediate semiring computations.
 */
#define TROP_ZERO ((int16_t)-32000)

/**
 * @brief Tropical Multiplicative Identity (\f$ e = 0 \f$).
 */
#define TROP_ONE  ((int16_t)0)

/**
 * @brief Safe upper bound on theoretical cumulative path weight (n * w_max).
 * If (n * w_max) exceeds this threshold, trop_karp_lambda returns TROP_ERR_OVERFLOW.
 */
#define TROP_WEIGHT_BOUND 30000

/**
 * @struct TropMatrix
 * @brief Matrix container for Max-Plus semiring calculations.
 *
 * Memory footprint: 2 + (16 * 16 * 2) = 514 bytes.
 */
typedef struct {
    uint8_t rows;                       /**< Number of active rows (1 <= rows <= MAX_DIM) */
    uint8_t cols;                       /**< Number of active columns (1 <= cols <= MAX_DIM) */
    int16_t data[MAX_DIM][MAX_DIM];     /**< 2D grid of tropical weights */
} TropMatrix;

/* ========================================================================== */
/*                           SCALAR TROPICAL OPS                              */
/* ========================================================================== */

/**
 * @brief Tropical Addition (\f$ a \oplus b = \max(a, b) \f$).
 *
 * @param a First operand in \f$ \mathbb{R}_{\max} \f$.
 * @param b Second operand in \f$ \mathbb{R}_{\max} \f$.
 * @return int16_t Maximum of @p a and @p b.
 */
int16_t trop_add(int16_t a, int16_t b);

/**
 * @brief Tropical Multiplication (\f$ a \otimes b = a + b \f$).
 *
 * @details Absorbing element rule: \f$ \varepsilon \otimes x = x \otimes \varepsilon = \varepsilon \f$.
 * Includes saturation check: sums exceeding INT16_MAX (32767) saturate to INT16_MAX to prevent
 * signed overflow wraparound from corrupting tropical maximum evaluations.
 *
 * @param a First operand in \f$ \mathbb{R}_{\max} \f$.
 * @param b Second operand in \f$ \mathbb{R}_{\max} \f$.
 * @return int16_t Classical addition @c (a + b) saturated, or @c TROP_ZERO if either operand is @c TROP_ZERO.
 */
int16_t trop_mul(int16_t a, int16_t b);

/* ========================================================================== */
/*                           MATRIX TROPICAL OPS                              */
/* ========================================================================== */

/**
 * @brief Initialize matrix filled with Tropical Zero (\f$ \varepsilon = -\infty \f$).
 *
 * @param[out] M Pointer to target matrix.
 * @param[in]  r Number of rows (1 <= r <= MAX_DIM).
 * @param[in]  c Number of columns (1 <= c <= MAX_DIM).
 * @return TROP_SUCCESS on success, or error code.
 */
int8_t trop_matrix_init(TropMatrix *M, uint8_t r, uint8_t c);

/**
 * @brief Construct the Tropical Identity Matrix (\f$ I_{\oplus} \f$).
 *
 * @details Diagonal elements are set to \f$ e = 0 \f$, and off-diagonal elements
 * are set to \f$ \varepsilon = \text{TROP\_ZERO} \f$.
 *
 * @param[out] M Pointer to target square matrix.
 * @param[in]  n Dimension (\f$ n \times n \f$, 1 <= n <= MAX_DIM).
 * @return TROP_SUCCESS on success, or error code.
 */
int8_t trop_matrix_identity(TropMatrix *M, uint8_t n);

/**
 * @brief Element-wise Tropical Matrix Addition (\f$ C = A \oplus B \f$).
 *
 * @details \f$ C_{i,j} = \max(A_{i,j}, B_{i,j}) \f$.
 * Safe against aliasing (@p C can be the same pointer as @p A or @p B).
 * Uses static scratch buffer (non-reentrant).
 *
 * @param[in]  A Pointer to first input matrix.
 * @param[in]  B Pointer to second input matrix.
 * @param[out] C Pointer to destination matrix.
 * @return TROP_SUCCESS on success, or error code.
 */
int8_t trop_matrix_add(const TropMatrix *A, const TropMatrix *B, TropMatrix *C);

/**
 * @brief Tropical Matrix Multiplication (\f$ C = A \otimes B \f$).
 *
 * @details Computes:
 * \f[ C_{i,j} = \max_{0 \le k < \text{cols}_A} (A_{i,k} + B_{k,j}) \f]
 * Safe against aliasing (@p C can be the same pointer as @p A or @p B).
 * Uses static scratch buffer (non-reentrant).
 *
 * @param[in]  A Pointer to left matrix (\f$ r_A \times c_A \f$).
 * @param[in]  B Pointer to right matrix (\f$ c_A \times c_B \f$).
 * @param[out] C Pointer to output matrix (\f$ r_A \times c_B \f$).
 * @return TROP_SUCCESS on success, or error code.
 */
int8_t trop_matmul(const TropMatrix *A, const TropMatrix *B, TropMatrix *C);

/**
 * @brief Tropical Matrix Power (\f$ \text{out} = A^{\otimes \text{power}} \f$).
 *
 * @details Implemented via binary exponentiation (repeated squaring) in \f$ O(n^3 \log k) \f$.
 * Safe against aliasing (@p out can be the same pointer as @p A).
 * Uses static scratch buffers (non-reentrant).
 *
 * @param[in]  A     Square input matrix.
 * @param[in]  power Exponent \f$ k \ge 0 \f$.
 * @param[out] out   Destination matrix.
 * @return TROP_SUCCESS on success, or error code.
 */
int8_t trop_matpow(const TropMatrix *A, uint32_t power, TropMatrix *out);

/**
 * @brief Tropical Kleene Star Closure (\f$ A^* \f$) with Positive Cycle Detection.
 *
 * @details Computes the subdiagonal/metric matrix closure:
 * \f[ A^* = \bigoplus_{k=0}^{n-1} A^{\otimes k} = I \oplus A \oplus A^{\otimes 2} \oplus \dots \oplus A^{\otimes (n-1)} \f]
 *
 * Precondition / Convergence Guarantee:
 * The Kleene star closure exists and is finite if and only if @p A contains no positive cycles
 * (i.e. maximum cycle mean \f$ \lambda(A) \le 0 \f$).
 * If any intermediate power \f$ A^k \f$ has a positive diagonal entry \f$ (A^k)_{i,i} > 0 \f$,
 * a positive cycle exists, causing the closure to diverge. This function halts and returns
 * @c TROP_ERR_POSITIVE_CYCLE.
 *
 * Uses static scratch buffers (non-reentrant).
 *
 * @param[in]  A      Square transition matrix (\f$ n \times n \f$).
 * @param[out] A_star Resulting Kleene star matrix.
 * @return TROP_SUCCESS on success, TROP_ERR_POSITIVE_CYCLE if positive cycle found, or error code.
 */
int8_t trop_kleene_star(const TropMatrix *A, TropMatrix *A_star);

/**
 * @brief Eliminate internal dependencies in Max-Plus linear systems.
 *
 * @details Solves the implicit discrete event dynamic equation:
 * \f[ x(k) = A_0 \otimes x(k) \oplus A_1 \otimes x(k-1) \f]
 * By Kleene's Star Lemma, this resolves to the explicit state evolution:
 * \f[ x(k) = (A_0^* \otimes A_1) \otimes x(k-1) \f]
 * This function computes \f$ A = A_0^* \otimes A_1 \f$.
 *
 * @param[in]  A0 Internal concurrent dependency matrix (must have no positive cycles).
 * @param[in]  A1 State transition matrix across consecutive stages/events.
 * @param[out] A  Equivalent reduced transition matrix.
 * @return TROP_SUCCESS on success, TROP_ERR_POSITIVE_CYCLE if A0 has positive cycle, or error code.
 */
int8_t trop_eliminate_dependencies(const TropMatrix *A0, const TropMatrix *A1, TropMatrix *A);

/* ========================================================================== */
/*                         SPECTRAL / CYCLE TIME OPS                          */
/* ========================================================================== */

/**
 * @brief Computes the Tropical Eigenvalue / Maximum Cycle Mean (\f$ \lambda \f$) using Karp's Algorithm.
 *
 * @details For a discrete event system modeled by \f$ x(k) = A \otimes x(k-1) \f$,
 * the maximum cycle mean \f$ \lambda(A) \f$ corresponds to the asymptotic cycle time
 * (throughput inverse):
 * \f[ \lambda = \max_{u \in V} \min_{0 \le k \le n-1} \frac{D_n(u) - D_k(u)}{n - k} \f]
 * where \f$ D_k = A \otimes D_{k-1} \f$.
 *
 * Exact Rational Representation:
 * The result is returned as an exact reduced fraction
 * \f$ \lambda = \frac{\text{lambda\_num}}{\text{lambda\_den}} \f$.
 *
 * Runtime Enforcements:
 * 1. Arithmetic Bound: Finds maximum finite weight w_max. If (n * w_max) > 30000,
 *    returns TROP_ERR_OVERFLOW to signal that weights risk overflow/saturation and
 *    the chosen time unit is too fine.
 * 2. Reachability Enforcement: Verifies every node u is reachable from node 0
 *    (has D[k][u] > TROP_ZERO for some k <= n). If any node is unreachable,
 *    returns TROP_ERR_NOT_REACHABLE.
 *
 * Embedded Resource Safety:
 * - State vectors \f$ D_k \f$ are stored in a file-scope static 2D array s_karp_D (544 bytes in .bss).
 * - Stack consumption is ~96 bytes (safe for 1 KB Keil stack).
 * - Non-reentrant: Do not call from ISR.
 *
 * @param[in]  A          Square adjacency/weight matrix of the Timed Event Graph (\f$ n \times n \f$).
 * @param[out] lambda_num Pointer to store the numerator of the exact cycle mean.
 * @param[out] lambda_den Pointer to store the denominator of the exact cycle mean.
 * @return TROP_SUCCESS on success,
 *         TROP_ERR_OVERFLOW if n * w_max > 30000,
 *         TROP_ERR_NOT_REACHABLE if any node is unreachable from node 0,
 *         TROP_ERR_NO_CYCLE if no valid cycle exists,
 *         or error code.
 */
int8_t trop_karp_lambda(const TropMatrix *A, int32_t *lambda_num, int32_t *lambda_den);

#ifdef __cplusplus
}
#endif

#endif /* TROPICAL_H */
