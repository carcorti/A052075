# A052075 — Phase 1 Report: Validation, Heuristic Model, and Feasibility

**Project:** Computational extension of OEIS A052075  
**Sequence:** Primes *p* such that nextprime(*p*) is a substring of *p*³  
**Author:** Carlo Corti  
**Date:** March 2026  
**Status:** Phase 1 complete — ready for Phase 2 (C+GMP implementation)

---

## 1. Executive Summary

This report documents the results of Phase 1 of the A052075 extension project.
The key findings are:

1. **All 15 known terms verified** — the oracle set is confirmed correct.
2. **A probabilistic model** predicts ~0.78 terms per decimal decade, suggesting
   the sequence is infinite.
3. **The observed data** (15 terms through 12 decades) is consistent with the
   model, though slightly above expectation (ratio 1.43).
4. **A critical technical barrier** explains why the search stalled in 2018:
   `p³` exceeds `__uint128_t` at `p ≈ 6.98 × 10¹²`, exactly where Resta stopped.
5. **A novel algebraic identity** for trailing matches enables an extremely
   efficient modular filter (>99.9999% rejection rate).
6. **The project is feasible** on the target hardware (UM790 Pro) and
   scientifically worthwhile.

---

## 2. Sequence Definition and Known Data

**A052075** (OEIS, Patrick De Geest, 2000): Primes *p* such that the decimal
representation of nextprime(*p*) appears as a contiguous substring of *p*³.

### 2.1 Known terms (oracle set)

| n  | a(n)            | nextprime       | p³                                     | Match pos | Type     |
|----|-----------------|-----------------|----------------------------------------|-----------|----------|
| 1  | 11              | 13              | 1331                                   | 0         | LEADING  |
| 2  | 101             | 103             | 1030301                                | 0         | LEADING  |
| 3  | 2239            | 2243            | 11224377919                             | 2         | interior |
| 4  | 34297           | 34301           | 40343019516073                          | 2         | interior |
| 5  | 43789           | 43793           | 83964379378069                          | 4         | interior |
| 6  | 53549           | 53551           | 153551511228149                         | 1         | interior |
| 7  | 535487          | 535489          | 153548930496746303                      | 1         | interior |
| 8  | 59897017        | 59897039        | 214889691497505989703913                | 14        | interior |
| 9  | 430784719       | 430784737       | 79943078473759892945966959              | 3         | interior |
| 10 | 2549592677      | 2549592733      | 16573430415736921632549592733           | 19        | TRAILING |
| 11 | 2837138669      | 2837138677      | 22837138677705447754568672309           | 1         | interior |
| 12 | 97969345967     | 97969346063     | 940309072235302647342697969346063       | 22        | TRAILING |
| 13 | 100000000019    | 100000000057    | 1000000000570000000108300000006859      | 0         | LEADING  |
| 14 | 328096840219    | 328096840223    | 35318816603250425999328096840223459     | 20        | interior |
| 15 | 4110739763869   | 4110739763909   | 69464026243759115461642614110739763909  | 25        | TRAILING |

**Lower bound:** a(16) > 6.9 × 10¹² (Giovanni Resta, Jul 02 2018)

### 2.2 Oracle set verification

All 15 terms verified computationally (Python/sympy) on four properties:
primality of *p*, correctness of nextprime(*p*), correctness of *p*³, and
substring match at the expected position. Additionally, an exhaustive
completeness check up to 6 × 10⁵ confirmed that no terms are missing
from the oracle set in that range.

---

## 3. Probabilistic Model

### 3.1 Naive random model

For a *d*-digit prime *p*, we model the digits of *p*³ as independent uniform
random variables. Then:

- *p*³ has approximately *D* = 3*d* digits
- nextprime(*p*) has *d* digits
- Number of starting positions for a substring match: *D* − *d* + 1 = 2*d* + 1
- Probability of match at any given position: 10^(−*d*)
- Probability of at least one match: *P* ≈ (2*d* + 1) · 10^(−*d*)

The expected number of hits in the decade [10^(*d*−1), 10^*d*) is:

> *E_d* = π(10^*d*) − π(10^(*d*−1)) × *P*
>       ≈ (9 · 10^(*d*−1)) / (*d* · ln 10) × (2*d* + 1) / 10^*d*
>       = 9(2*d* + 1) / (10 · *d* · ln 10)

**Key result:** As *d* → ∞, *E_d* → 9/(5 ln 10) ≈ **0.7817**.

Since the expected hits per decade converge to a positive constant, the
sum Σ *E_d* diverges, suggesting A052075 is **infinite**.

### 3.2 Refined model

The refined model accounts for the fact that *p*³ can have 3*d* − 2, 3*d* − 1,
or 3*d* digits depending on the leading digits of *p*. The boundaries occur at
*p* = 10^(*d*−1) · ∛10 and *p* = 10^(*d*−1) · ∛100. This gives slightly higher
expected values (by ~5-20%) for small *d*, converging to the naive model for
large *d*.

### 3.3 Model vs. data

| Digits *d* | Range                  | Actual | Expected (refined) | Ratio |
|------------|------------------------|--------|-------------------|-------|
| 2          | [10, 100)              | 1      | 1.052             | 0.95  |
| 3          | [100, 1000)            | 1      | 0.950             | 1.05  |
| 4          | [1000, 10⁴)            | 1      | 0.904             | 1.11  |
| 5          | [10⁴, 10⁵)            | 3      | 0.878             | 3.42  |
| 6          | [10⁵, 10⁶)            | 1      | 0.861             | 1.16  |
| 7          | [10⁶, 10⁷)            | 0      | 0.849             | 0.00  |
| 8          | [10⁷, 10⁸)            | 1      | 0.840             | 1.19  |
| 9          | [10⁸, 10⁹)            | 1      | 0.834             | 1.20  |
| 10         | [10⁹, 10¹⁰)           | 2      | 0.828             | 2.41  |
| 11         | [10¹⁰, 10¹¹)          | 1      | 0.824             | 1.21  |
| 12         | [10¹¹, 10¹²)          | 2      | 0.820             | 2.44  |
| 13         | [10¹², 10¹³)          | 1      | 0.817             | 1.22  |
| **Total**  |                        | **15** | **10.458**        | **1.43** |

The observed count (15) exceeds the expected (10.5) by a factor of 1.43.
For a Poisson process with λ = 10.5, P(X ≥ 15) ≈ 0.10, so the excess is
within 2σ and not statistically significant. However, extending the data
will determine whether this ratio persists.

---

## 4. Match Position Analysis

### 4.1 Classification of matches

The 15 known matches fall into three structural categories:

- **LEADING** (3 cases: a(1), a(2), a(13)): nextprime(*p*) appears at the
  start of *p*³. This requires *p*³/nextprime(*p*) ≈ 10^*k* for some integer *k*.

- **TRAILING** (3 cases: a(10), a(12), a(15)): nextprime(*p*) appears at the
  end of *p*³. This requires *p*³ ≡ nextprime(*p*) (mod 10^*d*).

- **INTERIOR** (9 cases): nextprime(*p*) appears somewhere in the middle of *p*³.

### 4.2 Positional distribution

Under the random model, the relative position (0 = start, 1 = end) should
be uniformly distributed. The observed distribution is:

- First third [0.00, 0.33):  9/15 = **60%** (expected 33%)
- Middle third [0.33, 0.67): 1/15 = **7%** (expected 33%)
- Last third [0.67, 1.00]:   5/15 = **33%** (expected 33%)

The strong excess in the first third and deficit in the middle third is
noteworthy. With only 15 data points, this could be statistical noise,
but it suggests that matches near the beginning of *p*³ may be favored —
possibly because the leading digits of *p*³ are correlated with *p*
(and hence with nextprime(*p*)).

### 4.3 Bimodal structure

The distribution appears almost **bimodal**: matches cluster near the start
(positions 0–3) and near the end (positions corresponding to relative
positions > 0.85), with the middle almost empty. This is an interesting
pattern that Phase 2 data could confirm or refute.

---

## 5. Algebraic Structure of Trailing Matches

### 5.1 The identity

For a trailing match (nextprime(*p*) = *p* + *g* appears at the end of *p*³),
we have:

> *p*³ ≡ *p* + *g* (mod 10^*d*)

Since *p*³ − *p* = *p*(*p* − 1)(*p* + 1), this is equivalent to:

> *p*(*p* − 1)(*p* + 1) ≡ *g* (mod 10^*d*)

where *g* is the prime gap. This identity is verified for all three trailing
cases:

| a(n) | p              | gap *g* | *p*(*p*−1)(*p*+1) mod 10^*d* |
|------|----------------|---------|------------------------------|
| 10   | 2549592677     | 56      | 56                           |
| 12   | 97969345967    | 96      | 96                           |
| 15   | 4110739763869  | 40      | 40                           |

### 5.2 Modular obstruction filter

This identity enables an extremely efficient **pre-filter** for trailing matches:

1. For each prime *p* with *d* digits, compute *r* = *p*(*p*² − 1) mod 10^*d*
2. If *r* > (*ln p*)² (approximate maximal gap by Cramér's conjecture), reject
3. If *r* is even and *p* + *r* is not divisible by small primes, test whether
   *p* + *r* is prime

The rejection rate is astronomical: near *p* ~ 10¹⁴, the maximal likely gap
is ~1000, while 10^*d* = 10¹⁵, so >99.99999% of candidates are eliminated
without computing *p*³ at all.

**Important limitation:** This filter detects only trailing matches. Leading
and interior matches require the full *p*³ string computation.

---

## 6. Hardware and Feasibility Analysis

### 6.1 Target hardware

**MinisForum UM790 Pro**
- CPU: AMD Ryzen 9 7940HS (8 cores / 16 threads, Zen 4, boost to 5.2 GHz)
- RAM: 64 GB DDR5
- L3 cache: 16 MB
- Architecture: x86-64 with AVX-512

### 6.2 The __uint128_t barrier

The critical finding of this analysis: the limit of unsigned 128-bit integer
arithmetic for *p*³ is:

> *p* < 2^(128/3) ≈ 2^42.67 ≈ **6.98 × 10¹²**

Giovanni Resta's 2018 bound is a(16) > **6.9 × 10¹²** — essentially identical.
This strongly suggests Resta used `__uint128_t` for *p*³ and stopped at its
natural boundary.

### 6.3 GMP requirement

To extend the search beyond 7 × 10¹², the code must use GMP (`mpz_t`) or
equivalent multi-precision arithmetic. GMP's `mpz_pow_ui` and `mpz_get_str`
are the key operations for computing *p*³ and converting to a decimal string
for substring matching.

### 6.4 Time estimates

With a segmented sieve + GMP + OpenMP on the Ryzen 9 7940HS:

| Search limit | Primes to test | Estimated time  | Expected new terms |
|--------------|----------------|------------------|--------------------|
| 10¹³         | ~3.3 × 10¹¹   | ~1–2 hours       | ~0.8               |
| 10¹⁴         | ~3.1 × 10¹²   | ~14–20 hours     | ~0.8               |
| 10¹⁵         | ~2.9 × 10¹³   | ~6–8 days        | ~0.8               |
| 10¹⁶         | ~2.7 × 10¹⁴   | ~2–3 months      | ~0.8               |

A target of **10¹⁵** is realistic (about one week of computation) and would
add 2–3 decades of unexplored territory.

---

## 7. Conjectures

### Conjecture 1 (Infiniteness)
The sequence A052075 is infinite.

**Evidence:** The heuristic predicts E_d → C ≈ 0.782 hits per decade.
Since Σ E_d diverges, the expected total count grows without bound.
The 15 observed terms through 12 decades are consistent with this prediction.

### Conjecture 2 (Asymptotic density)
Let N(x) = #{p ≤ x : p ∈ A052075}. Then N(x) ~ C · log₁₀(x) where
C = 9/(5 ln 10) ≈ 0.7817.

**Status:** Heuristic. Observed average per decade: 1.25 (predicted: 0.78).
The slight excess may be statistically insignificant or may indicate that
the random model underestimates the true density.

### Conjecture 3 (Match position non-uniformity)
The relative position of the substring match in p³ is NOT uniformly
distributed; it is bimodal, with peaks near position 0 (leading) and
position 1 (trailing), and a deficit in the interior.

**Status:** Tentative. Based on only 15 data points (60% in first third,
7% in middle third, 33% in last third). Requires extended data to confirm.

### Open Question 1
Is there a number-theoretic explanation for the leading-match condition
*p*³/nextprime(*p*) ≈ 10^*k*?

### Open Question 2
For trailing matches, is the condition *p*(*p*² − 1) ≡ *g* (mod 10^*d*) — where
*g* is the prime gap — related to known results on the distribution of cubes
modulo powers of 10?

---

## 8. Phase 2 Plan

### 8.1 Deliverables
1. **C source code** with GMP for p³, segmented sieve with OpenMP, and
   optimized substring matching
2. **Trailing-match modular filter** as a fast pre-screening pass
3. **Checkpoint/restart** mechanism for multi-day runs
4. **Comprehensive logging** of all candidates, timings, and per-decade
   statistics

### 8.2 Algorithm design

```
FOR each segment [lo, hi) of the prime sieve:
  1. Generate all primes in [lo, hi) via segmented sieve
  2. FOR each consecutive pair (p, q) where q = nextprime(p):
     a. FAST PATH: compute p(p²-1) mod 10^d; if result is a
        plausible gap AND equals q-p, flag as trailing match candidate
     b. FULL PATH: compute p³ via GMP, convert to string,
        check if str(q) is substring of str(p³)
     c. If hit: log (p, q, p³, position, type, timestamp)
  3. Checkpoint: write segment progress to disk every N segments
```

### 8.3 Verification protocol
- Run the oracle set before and after every code change
- Cross-validate first 8 terms against the Python reference
- For any new term found: verify independently with Python/sympy

---

## 9. File Inventory

```
A052075/
├── oracle/
│   └── oracle.csv              # Frozen oracle set (15 terms)
├── python/
│   ├── a052075_reference.py    # Reference implementation + verification
│   └── a052075_heuristic.py    # Heuristic model + statistical analysis
├── analysis/                   # (Phase 2: computational results)
├── docs/
│   └── phase1_report.md        # This document
└── (Phase 2: C source, Makefile, logs, LaTeX)
```

---

## 10. References

1. OEIS Foundation, "A052075 — Primes p such that nextprime(p) is substring
   of p^3," https://oeis.org/A052075
2. P. De Geest, original author (2000)
3. G. Resta, extension to a(12)–a(15) and bound a(16) > 6.9×10¹² (2018)
4. C. W. Wu, definition clarification and Python program (2022)
5. H. Cramér, "On the order of magnitude of the difference between consecutive
   prime numbers," Acta Arithmetica, 1936
