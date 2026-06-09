# OEIS A052075 — Computational Extension

**Primes *p* such that nextprime(*p*) appears as a contiguous decimal
substring of *p*³**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19421914.svg)](https://doi.org/10.5281/zenodo.19421914)

---

## Main results

| Result          | Value                                                        |
| --------------- | ------------------------------------------------------------ |
| New term        | **a(16) = 287,902,832,031,253**                              |
| nextprime       | 287,902,832,031,277 (gap = 24)                               |
| p³              | 23863701656637949988435953863**287902832031277** (44 digits) |
| Match position  | 29 (trailing match)                                          |
| Lower bound     | **a(17) > 4 × 10¹⁴**                                         |
| Search interval | [6.9 × 10¹², 4 × 10¹⁴], 37 contiguous segments               |
| Primes tested   | ≈ 1.2 × 10¹³                                                 |
| Elapsed time    | ≈ 66 hours wall-clock (237,868 s total)                      |

The new term a(16) is the first addition to
[OEIS A052075](https://oeis.org/A052075) since Resta's 2018 lower bound.
Independent verification was performed with Python 3 / SymPy 1.14.0.

---

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── .gitignore
├── data/
│   ├── campaign_summary.csv
│   ├── campaign_summary.log
│   └── oracle.csv
├── docs/
│   └── phase1_report.md
├── logs/
│   ├── seg01.log … seg37.log
│   └── seg01.csv … seg37.csv
├── paper/
│   ├── A052075.tex
│   └── A052075.pdf
└── src/
    ├── a052075_v4.c
    └── Makefile
```

---

## Algorithm summary

The search program `a052075_v4.c` implements:

- **Segmented sieve of Eratosthenes** — candidate primes enumerated in
  macro-blocks of width 1.28 × 10⁸, with OpenMP parallelism across blocks
  (16 threads on the test hardware).
- **Exact arithmetic** — `__uint128_t` for p ≤ 6.98 × 10¹²; GMP (`mpz_t`)
  for larger primes, where p³ exceeds 128-bit capacity.
- **Modular filter for trailing matches** — an exact criterion based on
  the congruence p(p²−1) ≡ g (mod 10^d_q), proved as Lemma 1 in the paper.
  Candidates failing this filter are still fully scanned, so leading and
  interior matches are not missed.

Average throughput: ≈ 50 × 10⁶ prime evaluations per second
(≈ 3.1 × 10⁶ per thread), measured on an AMD Ryzen 9 7940HS
(8 cores / 16 threads), 64 GB DDR5 RAM, Linux Mint 22.3.

---

## Building and running

Requirements: GCC or Clang with GMP and OpenMP support.

```bash
# Install GMP (Debian/Ubuntu/Mint)
sudo apt install libgmp-dev

# Build
cd src
make

# Run a single segment (example: 280×10¹² to 292×10¹²)
./a052075_v4 280000000000000 292000000000000
```

The program writes results to stdout and appends hit records to
`hits.csv` in the working directory.

---

## Known terms

The complete oracle table for all 16 known terms is provided in
`data/oracle.csv`. Term 16 found in this work is:

| n   | p               | nextprime(p)    | gap | match type |
| --- | ---------------:| ---------------:| ---:| ---------- |
| 16  | 287902832031253 | 287902832031277 | 24  | trailing   |

Terms 1–15 are due to earlier OEIS contributors; term 16 was found in
this work.

---

## Paper

The accompanying manuscript is:

> Carlo Corti, *Computational Extension of OEIS Sequence A052075:
> A New Term and an Exhaustive Search to 4 × 10¹⁴*, 2026.

The LaTeX source (`paper/A052075.tex`) and compiled PDF
(`paper/A052075.pdf`) are included for transparency and reproducibility.
See the [LICENSE](LICENSE) for usage terms applicable to the paper files.

---

## Data and code availability

The complete dataset, including the oracle CSV file, per-segment logs,
campaign summary, and machine-readable results, is available in this
repository and in the corresponding Zenodo archive:

- **GitHub**: [github.com/carcorti/A052075](https://github.com/carcorti/A052075)
- **Zenodo**: [doi.org/10.5281/zenodo.19421914](https://doi.org/10.5281/zenodo.19421914)
- **OEIS**: [oeis.org/A052075](https://oeis.org/A052075)

The source file `src/a052075_v4.c`, the associated `src/Makefile`, and
the Python verification script `verify_a052075.py` are available in the
GitHub repository and Zenodo archive.

---

## License

This repository is released under the MIT License; see [LICENSE](LICENSE) for full terms.

Copyright (c) 2026 Carlo Corti.

---

## Citation

If you use data or code from this repository, please cite the Zenodo
archive:

```text
Carlo Corti (2026). OEIS A052075: Computational Extension.
Zenodo. https://doi.org/10.5281/zenodo.19421914
```

Citation metadata are provided in `CITATION.cff`.

---

## Author

**Carlo Corti** — independent researcher  
`carlo.corti@outlook.com`
