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
 *
 * 3. Matrix Multiplication (\f$ C = A \otimes B \f$):
 *    \f[ C_{i,j} = \bigoplus_{k=1}^{n} (A_{i,k} \otimes B_{k,j}) = \max_{1 \le k \le n} (A_{i,k} + B_{k,j}) \f]
 *
 * Target: STM32F407 (ARM Cortex-M4F, ARMv7E-M).
 * Memory: Static / Stack only (Zero Dynamic Memory).
 * Core Dependencies: <stdint.h> only.
 */

#ifndef TROPICAL_H
#define TROPICAL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
 * @struct TropMatrix
 * @brief Static matrix container for Max-Plus semiring calculations.
 *
 * Memory footprint: 2 + (16 * 16 * 2) = 514 bytes.
 * Safe for embedded stack and static BSS allocation on STM32F407.
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
 *
 * @param a First operand in \f$ \mathbb{R}_{\max} \f$.
 * @param b Second operand in \f$ \mathbb{R}_{\max} \f$.
 * @return int16_t Classical addition @c (a + b), or @c TROP_ZERO if either operand is @c TROP_ZERO.
 */
int16_t trop_mul(int16_t a, int16_t b);

/* ========================================================================== */
/*                           MATRIX TROPICAL OPS                              */
/* ========================================================================== */

/**
 * @brief Initialize matrix filled with Tropical Zero (\f$ \varepsilon = -\infty \f$).
 *
 * @param[out] M Pointer to target matrix.
 * @param[in]  r Number of rows.
 * @param[in]  c Number of columns.
 */
void trop_matrix_init(TropMatrix *M, uint8_t r, uint8_t c);

/**
 * @brief Construct the Tropical Identity Matrix (\f$ I_{\oplus} \f$).
 *
 * @details Diagonal elements are set to \f$ e = 0 \f$, and off-diagonal elements
 * are set to \f$ \varepsilon = \text{TROP\_ZERO} \f$.
 *
 * @param[out] M Pointer to target square matrix.
 * @param[in]  n Dimension (\f$ n \times n \f$).
 */
void trop_matrix_identity(TropMatrix *M, uint8_t n);

/**
 * @brief Element-wise Tropical Matrix Addition (\f$ C = A \oplus B \f$).
 *
 * @details \f$ C_{i,j} = \max(A_{i,j}, B_{i,j}) \f$.
 * Safe against aliasing (@p C can be the same pointer as @p A or @p B).
 *
 * @param[in]  A Pointer to first input matrix.
 * @param[in]  B Pointer to second input matrix.
 * @param[out] C Pointer to destination matrix.
 */
void trop_matrix_add(const TropMatrix *A, const TropMatrix *B, TropMatrix *C);

/**
 * @brief Tropical Matrix Multiplication (\f$ C = A \otimes B \f$).
 *
 * @details Computes:
 * \f[ C_{i,j} = \max_{0 \le k < \text{cols}_A} (A_{i,k} + B_{k,j}) \f]
 * Safe against aliasing (@p C can be the same pointer as @p A or @p B).
 *
 * @param[in]  A Pointer to left matrix (\f$ r_A \times c_A \f$).
 * @param[in]  B Pointer to right matrix (\f$ c_A \times c_B \f$).
 * @param[out] C Pointer to output matrix (\f$ r_A \times c_B \f$).
 */
void trop_matmul(const TropMatrix *A, const TropMatrix *B, TropMatrix *C);

/**
 * @brief Tropical Matrix Power (\f$ \text{out} = A^{\otimes \text{power}} \f$).
 *
 * @details Implemented via binary exponentiation (repeated squaring) in \f$ O(n^3 \log k) \f$.
 * Special cases:
 * - \f$ A^{\otimes 0} = I_{\oplus} \f$ (Tropical Identity)
 * - \f$ A^{\otimes 1} = A \f$
 *
 * @param[in]  A     Square input matrix.
 * @param[in]  power Exponent \f$ k \ge 0 \f$.
 * @param[out] out   Destination matrix.
 */
void trop_matpow(const TropMatrix *A, uint32_t power, TropMatrix *out);

/**
 * @brief Tropical Kleene Star Closure (\f$ A^* \f$).
 *
 * @details Computes the subdiagonal/metric matrix closure:
 * \f[ A^* = \bigoplus_{k=0}^{n-1} A^{\otimes k} = I \oplus A \oplus A^{\otimes 2} \oplus \dots \oplus A^{\otimes (n-1)} \f]
 *
 * For a system of discrete event dependencies without positive cycles,
 * \f$ A^*_{i,j} \f$ represents the maximum weight path between state @c j and state @c i.
 *
 * @param[in]  A      Square transition matrix (\f$ n \times n \f$).
 * @param[out] A_star Resulting Kleene star matrix.
 */
void trop_kleene_star(const TropMatrix *A, TropMatrix *A_star);

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
 */
void trop_eliminate_dependencies(const TropMatrix *A0, const TropMatrix *A1, TropMatrix *A);

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
 * @param[in] A Square adjacency/weight matrix of the Timed Event Graph (\f$ n \times n \f$).
 * @return int16_t Tropical eigenvalue \f$ \lambda \f$, or @c TROP_ZERO if no cycle exists.
 */
int16_t trop_karp_lambda(const TropMatrix *A);

#ifdef __cplusplus
}
#endif

#endif /* TROPICAL_H */
