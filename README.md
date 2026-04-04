# OEIS A052075 — Computational Extension

**Primes *p* such that nextprime(*p*) appears as a contiguous decimal
substring of *p*³**

[![License: MIT (code)](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![License: CC0 (data)](https://img.shields.io/badge/data-CC0%201.0-lightgrey.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Main results

| Result | Value |
|--------|-------|
| New term | **a(16) = 287,902,832,031,253** |
| nextprime | 287,902,832,031,277 (gap = 24) |
| p³ | 23863701656637949988435953863**287902832031277** (44 digits) |
| Match position | 29 (trailing match) |
| Lower bound | **a(17) > 4 × 10¹⁴** |
| Search interval | [6.9 × 10¹², 4 × 10¹⁴], 37 contiguous segments |
| Primes tested | ≈ 1.2 × 10¹³ |
| Elapsed time | ≈ 66 hours wall-clock (237,868 s total) |

The new term a(16) is the first addition to
[OEIS A052075](https://oeis.org/A052075) since Resta's 2018 lower bound.
Independent verification was performed with Python 3 / SymPy 1.14.0.

---

## Repository structure

```
.
├── data/
│   ├── campaign_summary.csv   # per-segment statistics (machine-readable)
│   ├── campaign_summary.log   # per-segment statistics (human-readable)
│   └── oracle.csv             # all 16 known terms with full metadata
├── docs/
│   └── phase1_report.md       # internal campaign narrative
├── logs/
│   ├── seg01.log … seg37.log  # per-segment search logs
│   └── seg01.csv … seg37.csv  # per-segment hit records
├── paper/
│   ├── A052075_rev11.tex      # LaTeX source (article class)
│   ├── A052075_rev11.bbl      # bibliography (BibTeX output)
│   └── A052075_rev11.pdf      # compiled manuscript
└── src/
    ├── a052075_v4.c           # C search program (version 4)
    └── Makefile               # build instructions
```

---

## Algorithm summary

The search program `a052075_v4.c` implements:

- **Segmented sieve of Eratosthenes** — candidate primes enumerated in
  macro-blocks of width 1.28 × 10⁸, with OpenMP parallelism across blocks
  (16 threads on the test hardware).
- **Exact arithmetic** — `__uint128_t` for p ≤ 6.98 × 10¹²; GMP (`mpz_t`)
  for larger primes, where p³ exceeds 128-bit capacity.
- **Modular filter for trailing matches** — an exact early-exit criterion
  based on the congruence p(p²−1) ≡ g (mod 10^d_q), proved as
  Lemma 1 in the paper. This eliminates the vast majority of candidates
  before computing the full cube.

Average throughput: ≈ 50 × 10⁶ prime evaluations per second
(≈ 3.1 × 10⁶ per thread), measured on a MinisForum UM790 Pro
(AMD Ryzen 9 7940HS, 8 cores / 16 threads, 64 GB DDR5, Linux Mint 22.3).

---

## Building and running

Requirements: GCC (or Clang) with GMP and OpenMP support.

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

## Known terms (oracle)

| n  | p               | nextprime(p)    | gap | match type |
|----|-----------------|-----------------|-----|------------|
|  1 | 11              | 13              |   2 | leading    |
|  2 | 101             | 103             |   2 | leading    |
|  3 | 2239            | 2243            |   4 | interior   |
|  4 | 34297           | 34301           |   4 | interior   |
|  5 | 43789           | 43793           |   4 | interior   |
|  6 | 53549           | 53551           |   2 | interior   |
|  7 | 535487          | 535489          |   2 | interior   |
|  8 | 59897017        | 59897039        |  22 | interior   |
|  9 | 430784719       | 430784737       |  18 | interior   |
| 10 | 2549592677      | 2549592733      |  56 | trailing   |
| 11 | 2837138669      | 2837138677      |   8 | interior   |
| 12 | 97969345967     | 97969346063     |  96 | trailing   |
| 13 | 100000000019    | 100000000057    |  38 | leading    |
| 14 | 328096840219    | 328096840223    |   4 | interior   |
| 15 | 4110739763869   | 4110739763909   |  40 | trailing   |
| 16 | 287902832031253 | 287902832031277 |  24 | trailing   |

Terms 1–15 from De Geest (2000) and Handler (2006).
Term 16 found in this work (2026).

---

## Paper

The accompanying manuscript is:

> Carlo Corti, *Computational Extension of OEIS Sequence A052075:
> A New Term and an Exhaustive Search to 4 × 10¹⁴*,
> submitted to Journal of Integer Sequences, 2026.

The LaTeX source (`paper/A052075_rev11.tex`) and the compiled PDF
(`paper/A052075_rev11.pdf`) are included for transparency and
reproducibility. See the [LICENSE](LICENSE) for usage terms applicable
to the paper files.

---

## Data availability

- **GitHub**: this repository —
  [github.com/carcorti/OEIS-A052075](https://github.com/carcorti/OEIS-A052075)
- **Zenodo**: permanent archive —
  [doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)
  *(DOI will be finalised upon Zenodo release)*
- **OEIS**: [oeis.org/A052075](https://oeis.org/A052075)

---

## License

This repository uses a split license; see [LICENSE](LICENSE) for full terms.

- **Source code** (`src/`): MIT License
- **Data and logs** (`data/`, `logs/`): CC0 1.0 Universal (Public Domain)
- **Paper** (`paper/`): © 2026 Carlo Corti, all rights reserved

---

## Citation

If you use data or code from this repository, please cite the Zenodo
archive:

```
Carlo Corti (2026). OEIS A052075 – Computational Extension:
a(16) = 287,902,832,031,253 and exhaustive search to 4×10¹⁴.
Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
```

A BibTeX entry is provided below:

```bibtex
@misc{Corti2026_A052075,
  author       = {Carlo Corti},
  title        = {{OEIS A052075} – Computational Extension:
                  $a(16) = 287\,902\,832\,031\,253$ and exhaustive
                  search to $4 \times 10^{14}$},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

---

## Author

**Carlo Corti** — independent researcher, Kalinówka, Poland
`carlo.corti@outlook.com`
