# Exact moments of the Toeplitz and Hankel limiting spectral distributions

Code, data and paper for computing the moments of the limiting spectral
distributions of random symmetric Toeplitz and Hankel matrices (Bryc, Dembo and
Jiang 2006; Hammond and Miller 2005) in exact rational arithmetic, through
order 14.

| 2k | Toeplitz  | Hankel      |
|----|-----------|-------------|
| 2  | 1         | 1           |
| 4  | 8/3       | 2           |
| 6  | 11        | 11/2        |
| 8  | 908/15    | 281/15      |
| 10 | 415       | 2717/36     |
| 12 | 23840/7   | 1052/3      |
| 14 | 325719/10 | 1331087/720 |

The values at orders 10 to 14 are new. At order 8 our Toeplitz value differs
from the one in Hammond and Miller (964/15); the configuration count behind the
corrected value, and the checks, are in Section 5 of the paper.

Each moment is a sum of volumes of rational polytopes indexed by pair
partitions. The volumes are computed exactly (integer-only double description
plus triangulation, no external libraries) and cross-checked by an independent
vertex enumeration, Monte Carlo, eigenvalue simulation and an exact counting
route. The paper also proves that a pairing's polytope volume is 1 exactly when
the pairing is non-crossing, and at most 2/3 otherwise, in both ensembles.

## Layout

    paper/                  LaTeX source and PDF
    results/                per-orbit exact volumes, one JSON file per ensemble and order
    tk_core.py              pairings, constraint forms, dihedral orbits, Monte Carlo
    exact_volume.py         exact rational volume engine, with self-tests
    run_moments.py          production runner (parallel, checkpointed, resumable)
    verify_m8.py, mc_heavy_t8.py, brute_vertices_check.py,
    exact_count_check.py, sim_matrix.py, validate_more.py,
    tl_noncrossing.py       verification scripts
    analysis.py             deficit ratios, candidate exclusions, tail bounds
    fig_density.py          density figure
    references/             related code and documents

## Running

    pip install numpy scipy sympy
    python exact_volume.py                   # self-tests and published-value checks, about a minute
    python run_moments.py T 12 --mc 1000000  # exact m_12, Toeplitz, a couple of minutes on 14 cores
    python analysis.py

The exact engine itself is pure Python; numpy/scipy/sympy are only used by the
verification and analysis scripts. Toeplitz at order 14 (5283 orbits) takes
about 90 minutes on 14 cores; everything else is minutes.
