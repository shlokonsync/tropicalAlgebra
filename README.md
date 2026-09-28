# Tropical (Max-Plus) Algebra Library for Bare-Metal Embedded Systems

A lightweight, zero-dynamic-memory C99 library implementing **Max-Plus (Tropical) Semiring Algebra** and **Spectral Graph Analysis (Karp's Cycle Mean Algorithm)** designed for resource-constrained ARM Cortex-M microcontrollers (STM32F407, ARMv7E-M).

---

##  Mathematical Foundations

The $(\max, +)$ tropical semiring is denoted by $\mathbb{R}_{\max} = (\mathbb{R} \cup \{-\infty\}, \oplus, \otimes)$:

1. **Tropical Addition ($\oplus$):**
   $$a \oplus b \triangleq \max(a, b)$$
   - Identity element: $\varepsilon = -\infty$ (`TROP_ZERO = -32000`)
   - Idempotent: $a \oplus a = a$

2. **Tropical Multiplication ($\otimes$):**
   $$a \otimes b \triangleq a + b$$
   - Identity element: $e = 0$
   - Absorbing zero: $a \otimes \varepsilon = \varepsilon \otimes a = \varepsilon$

3. **Matrix Multiplication ($C = A \otimes B$):**
   $$C_{i,j} = \bigoplus_{k=1}^n (A_{i,k} \otimes B_{k,j}) = \max_{1 \le k \le n} (A_{i,k} + B_{k,j})$$

4. **Kleene Star Closure ($A^*$):**
   $$A^* = \bigoplus_{k=0}^{n-1} A^{\otimes k} = I \oplus A \oplus A^{\otimes 2} \oplus \dots \oplus A^{\otimes (n-1)}$$
   Computes all-pairs longest paths in weighted event precedence graphs and resolves implicit state dependencies $x = A_0 x \oplus B \implies x = A_0^* B$.

5. **Spectral Analysis (Karp's Maximum Cycle Mean Algorithm):**
   $$\lambda(A) = \max_{u \in V} \min_{0 \le k \le n-1} \frac{D_n(u) - D_k(u)}{n - k}$$
   where $D_k = A \otimes D_{k-1}$ and $D_0 = [0, \varepsilon, \dots, \varepsilon]^T$. $\lambda(A)$ gives the maximum cycle mean, corresponding to the asymptotic cycle time (inverse throughput) of discrete event systems / Timed Event Graphs.

---

## ⚡ Embedded Constraints & Design

- **Zero Dynamic Memory Allocation:** No `malloc`, `calloc`, or `free`. Static structures and fixed stack buffers only.
- **Fixed-Width Types:** Pure `<stdint.h>` (`int16_t`, `uint8_t`, `uint32_t`) mapping directly to ARM Cortex-M4 DSP registers.
- **Arithmetic Safety:** `TROP_ZERO` is defined as `-32000` (leaving negative headroom above `INT16_MIN` = `-32768`) to prevent arithmetic underflow when negative offsets are added.
- **Strict ANSI C99:** Free of compiler-specific builtins, extensions, or VLAs. Fully portable between GCC and **Keil MDK Arm Compiler 6 (`armclang`)**.
- **No Standard I/O in Core:** Core modules (`tropical.h`, `tropical.c`, `karp.c`) contain zero `printf` or standard library calls.

---

## 📁 Repository Structure

```
tropicalAlgebra/
├── tropical.h        # Core type definitions and function prototypes
├── tropical.c        # Semiring scalar/matrix operations, powers, Kleene star
├── karp.c            # Karp's algorithm for tropical eigenvalue (cycle time)
├── test_week1.c      # Host test harness with formatted terminal output
├── Makefile          # Host build & test automation
└── README.md         # Documentation & mathematical notes
```

---

## 🛠️ Building and Testing

### Host (GCC / Linux / Windows MinGW)
```bash
# Build the test harness
make

# Run the test suite
make test
```

Or manually:
```bash
gcc -Wall -Wextra -Wpedantic -Werror -std=c99 tropical.c karp.c test_week1.c -o test_week1
./test_week1
```

### Keil MDK (Arm Compiler 6 / `armclang`)
1. Add `tropical.c` and `karp.c` to your project source group.
2. Add the repository directory to **Options for Target -> C/C++ (AC6) -> Include Paths**.
3. Include `"tropical.h"` in your STM32 firmware files.

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

Running Karp's Maximum Cycle Mean Algorithm...
Computed Tropical Eigenvalue lambda = 3
[PASS] Assertion lambda == 3 SUCCESSFUL!

=== Matrix A_elim = A0* (x) A1 (2x2) ===
  [      5      6 ]
  [      3      4 ]

[PASS] Dependency elimination (Kleene Star) verified.

===============================================================
  ALL TESTS PASSED! Production-ready for STM32F407 deployment. 
===============================================================
```
