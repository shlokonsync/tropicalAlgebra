# Tropical (Max-Plus) Algebra Library for Bare-Metal Embedded Systems

A hardened, zero-dynamic-memory C99 library implementing **Max-Plus (Tropical) Semiring Algebra** and **Spectral Graph Analysis (Karp's Cycle Mean Algorithm)** designed for resource-constrained ARM Cortex-M microcontrollers (STM32F407, ARMv7E-M).

---

## 📐 Mathematical Foundations

The $(\max, +)$ tropical semiring is denoted by $\mathbb{R}_{\max} = (\mathbb{R} \cup \{-\infty\}, \oplus, \otimes)$:

1. **Tropical Addition ($\oplus$):**
   $$a \oplus b \triangleq \max(a, b)$$
   - Identity element: $\varepsilon = -\infty$ (`TROP_ZERO = -32000`)
   - Idempotent: $a \oplus a = a$

2. **Tropical Multiplication ($\otimes$):**
   $$a \otimes b \triangleq a + b$$
   - Identity element: $e = 0$
   - Absorbing zero: $a \otimes \varepsilon = \varepsilon \otimes a = \varepsilon$
   - Positive saturation: Saturated at `INT16_MAX` ($32,767$) to eliminate signed overflow wraparound.

3. **Matrix Multiplication ($C = A \otimes B$):**
   $$C_{i,j} = \bigoplus_{k=1}^n (A_{i,k} \otimes B_{k,j}) = \max_{1 \le k \le n} (A_{i,k} + B_{k,j})$$

4. **Kleene Star Closure ($A^*$) with Positive Cycle Detection:**
   $$A^* = \bigoplus_{k=0}^{n-1} A^{\otimes k} = I \oplus A \oplus A^{\otimes 2} \oplus \dots \oplus A^{\otimes (n-1)}$$
   - **Precondition & Causality:** The closure exists if and only if $A$ contains no cycles with positive weight ($\lambda(A) \le 0$).
   - If any power $A^k$ ($1 \le k \le n$) has a positive diagonal entry $(A^k)_{i,i} > 0$, the closure diverges. The library halts and returns `TROP_ERR_POSITIVE_CYCLE`.
   - Resolves implicit discrete event system dependencies: $x(k) = A_0 x(k) \oplus A_1 x(k-1) \implies x(k) = (A_0^* A_1) x(k-1)$.

5. **Spectral Analysis (Karp's Maximum Cycle Mean Algorithm):**
   $$\lambda(A) = \max_{u \in V} \min_{0 \le k \le n-1} \frac{D_n(u) - D_k(u)}{n - k}$$
   where $D_k = A \otimes D_{k-1}$ and $D_0 = [0, \varepsilon, \dots, \varepsilon]^T$.
   - **Exact Rational Output:** Returns exact reduced fraction $\lambda = \frac{p}{q}$ (`lambda_num / lambda_den`) without integer truncation loss.
   - **Precondition:** Matrix $A$ must be irreducible (Timed Event Graph strongly connected), ensuring all cycles are reachable from the initial state vector.

---

## ⚡ Embedded Constraints & Memory Architecture

- **Zero Dynamic Memory Allocation:** No `malloc`, `calloc`, or `free`.
- **Minimal Stack Footprint (< 160 Bytes):** Large temporary matrices and vector trajectories are stored in file-scope static `.bss` buffers ($544\text{ B}$ for Karp vectors, $514\text{ B}$ per matrix scratchpad). This guarantees complete stack safety on default Keil MDK / CMSIS $1\text{ KB}$ startup stacks (worst-case call stack $\le 272\text{ bytes}$).
- **Concurrency & Non-Reentrancy:** Because file-scope static scratchpads are used, functions are **non-reentrant**. They are safe for single-threaded bare-metal main loop execution, but **must not be called from Interrupt Service Routines (ISRs)**.
- **Timing Units & Dynamic Range:**
  - With `int16_t`, $32,767\text{ cycles}$ at $168\text{ MHz}$ is only $\approx 195\ \mu\text{s}$.
  - For pipeline modeling with $n \le 16$, edge weights must satisfy $w \times n < 32,000$.
  - **Recommended Time Base:** Scale pipeline timings to $10\ \mu\text{s}$ ticks ($32,767\text{ ticks} \approx 327.67\text{ ms}$) or $1\ \mu\text{s}$ ticks, ensuring safe headroom for Cortex-M4 SIMD (`__QADD16`/`__SADD16`).
- **Strict ANSI C99:** Fully compatible with `-Wall -Wextra -Wpedantic -Werror -std=c99` and Keil MDK Arm Compiler 6 (`armclang`).
- **Zero Standard I/O in Core:** Core modules (`tropical.h`, `tropical.c`, `karp.c`) contain zero `printf` or standard library dependencies.

---

## 📁 Repository Structure

```
tropicalAlgebra/
├── tropical.h        # Core type definitions, error codes, and function prototypes
├── tropical.c        # Semiring scalar/matrix operations, powers, Kleene star
├── karp.c            # Karp's algorithm with static vector buffers and exact rational lambda
├── test_week1.c      # Verification suite (powers, rational lambda, positive cycle traps)
├── Makefile          # Host build & test automation (-Wpedantic enabled)
└── README.md         # Documentation & mathematical notes
```

---

## 🛠️ Building and Testing

### Host (GCC / Linux / Windows MinGW)
```bash
# Build the test harness with strict pedantic flags
gcc -Wall -Wextra -Wpedantic -Werror -O2 -std=c99 tropical.c karp.c test_week1.c -o test_week1

# Run the test suite
./test_week1
```

Or using `make`:
```bash
make
make test
```

### Keil MDK (Arm Compiler 6 / `armclang`)
1. Add `tropical.c` and `karp.c` to your project source group.
2. Add the repository directory to **Options for Target -> C/C++ (AC6) -> Include Paths**.
3. Include `"tropical.h"` in your STM32 firmware files.
4. Ensure routines are invoked from thread mode (main loop), not ISR context.

---

## 🧪 Test Results

```text
===============================================================
  STM32F407 Tropical (Max-Plus) Algebra Library - Test Suite   
===============================================================

=== Matrix A (Original) (3x3) ===
  [      2      1    eps ]
  [    eps      3      1 ]
  [      1      2      2 ]

=== Matrix A^2 (A (x) A) (3x3) ===
  [      4      4      2 ]
  [      2      6      4 ]
  [      3      5      4 ]

=== Matrix A^3 (A^2 (x) A) (3x3) ===
  [      6      7      5 ]
  [      5      9      7 ]
  [      5      8      6 ]

[PASS] A^2 calculation verified.
[PASS] A^3 calculation verified.

===============================================================
  Karp's Spectral Analysis: Step-by-Step Matrix-Vector Trace   
===============================================================

  D_0 = [    0  eps  eps ]^T
  D_1 = [    2  eps    1 ]^T
  D_2 = [    4    2    3 ]^T
  D_3 = [    6    5    5 ]^T

[PASS] D_k recurrence vectors D_0 -> D_3 verified.

--- Karp's Rational Ratio Breakdown Table: (D_3(u) - D_k(u)) / (3 - k) ---
  Node u |  k=0 (len 3)  |  k=1 (len 2)  |  k=2 (len 1)  |  min_k over u 
  -------+---------------+---------------+---------------+---------------
    0    |  (6-0)/3 = 2  |  (6-2)/2 = 2  |  (6-4)/1 = 2  |      2        
    1    |      eps      |      eps      |  (5-2)/1 = 3  |      3        
    2    |      eps      |  (5-1)/2 = 2  |  (5-3)/1 = 2  |      2        
  -------+---------------+---------------+---------------+---------------
  Maximum Cycle Mean lambda = max_u { 2, 3, 2 } = 3/1

trop_karp_lambda(&A) returned: 3/1 (status: 0)
[PASS] Karp's integer eigenvalue assertion (lambda == 3/1) SUCCESSFUL!

===============================================================
  Fractional Eigenvalue Precision Verification                 
===============================================================

Fractional pipeline matrix lambda returned: 7/2 (float: 3.50, status: 0)
[PASS] Fractional eigenvalue preserved exactly without integer truncation!

===============================================================
  Kleene Star Closure & Discrete Event Dependency Elimination 
===============================================================

=== Matrix A_elim = A0* (x) A1 (2x2) ===
  [      5      6 ]
  [      3      4 ]

[PASS] Dependency elimination (Kleene Star) verified.

===============================================================
  Kleene Star Positive Cycle Rejection (Causality Enforcement) 
===============================================================

Kleene star on positive cycle returned status: -4 (expected -4)
[PASS] Positive cycle detected and successfully rejected by Kleene star!

===============================================================
  Scalar Multiplication Saturation Protection                  
===============================================================

trop_mul(20000, 20000) = 32767 (saturated at INT16_MAX = 32767)
[PASS] Positive arithmetic overflow saturates safely without wraparound!

===============================================================
  ALL TESTS PASSED! Refactored & Hardened for Bare-Metal ARM. 
===============================================================
```
