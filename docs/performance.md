# Performance

## How the kernels are built

- **One fused loop per physics operation.** `adm_rhs` computes the metric
  derivatives, the Christoffel symbols, the Ricci tensor and all 16
  right-hand sides point by point, in a single pass over memory.
- `@njit(parallel=True, fastmath=True, cache=True, error_model="numpy")`:
  threads over the outer `i` index, LLVM optimisations as in `-O3 -ffast-math`,
  and on-disk caching of the machine code.
- Direction-specialised stencils (`pynr.kernels.fd.dx, dxy, ...`). Numba
  then emits straight-line code with no branches on the direction index.
- The Ricci tensor needs only the two contractions of $\partial\Gamma$, not all
  162 components (see `pynr.kernels.adm`).
- Kreiss-Oliger dissipation loops with the contiguous index innermost, so it
  vectorises (SIMD).

## Measured on an Apple M3 Pro, 12 cores (6P+6E), 18 GB, macOS 27, Python 3.14, Numba 0.67, 12 threads

```{table} Measured kernel timings.
:name: tab-performance

| kernel | grid | time | per point |
|---|---|---|---|
| ADM RHS, 1 thread | $64^3$ | 131 ms | 670 ns |
| ADM RHS, 12 threads | $96^3$ | 78 ms | 107 ns (≈1.3 µs·thread) |
| KO dissipation, 16 vars, 12 threads | $96^3$ | 8.5 ms | 0.7 ns / var |
```

Reproduce with `python scripts/benchmark.py 96`.

A full RK4 step of problem 6 ($135^3$ points, ADM + dissipation + boundaries)
takes about 1 s on this machine: 576 steps in 593 s (the problem-set cost table, {numref}`tab-problem-cost`, lists all problems with their machine).

## Measured on an Intel i7-12700H laptop, and against the Einstein Toolkit

The same problems were run on a second machine, an Intel Core i7-12700H laptop (14 cores: 6 performance cores with
hyper-threading + 8 efficiency cores, 20 threads; 15 GB; Linux 6.8; Python 3.11, Numba 0.61, 20 Numba threads),
together with their Einstein Toolkit counterparts (ET_2026_05, Intel oneAPI 2023 build, 2 MPI ranks × 6 OpenMP
threads). Wall times include start-up, initial data, analysis and output.

```{table} Wall time per problem: pyNR on two machines and the Einstein Toolkit on the laptop. µs/(point·step): wall time divided by grid points (incl. ghosts) × time steps.
:name: tab-performance-machines

| problem | pyNR, M3 Pro | pyNR, i7-12700H | ET, i7-12700H | grid, steps (pyNR / ET) | µs/(point·step), pyNR / ET (i7) |
|---|---|---|---|---|---|
| 1 gauge wave | 3.2 s | 13 s ¹ | 58 s | $57{\times}9{\times}9$, 2000 / $56{\times}16{\times}16$, 2000 | — ² |
| 2 linear wave | 1.7 s | 2.0 s | 28 s | $57{\times}9{\times}9$, 1000 / $56{\times}16{\times}16$, 1000 | — ² |
| 3 geodesic slicing | 2.4 s (abort $t = 1$) | 4.4 s (abort $t = 1$) | 23 s (stop $t = 3.1$) | $68^3$, 20 / $68^3$, 62 | — ² |
| 4 1+log | 45 s (abort $t = 5$) | 82 s (abort $t = 5$) | 3.2 min ($t = 6$) | $108^3$, 100 / $108^3$, 120 | 0.65 (ADM) / 1.26 (BSSN) |
| 6 perturbed BH | 9.9 min | 18.7 min | 38 min | $136^3$, 576 / $135^3$, 640 | 0.77 (ADM) / 1.44 (BSSN) |
| 7 moving puncture | 5.6 min | 11.0 min | 13.0 min | $88^3$, 800 / $88^3$, 800 | **1.21 / 1.43 (both BSSN)** |
```

¹ First run of the session: includes the Numba compilation of the kernels. ² Too small for a meaningful per-point
cost: start-up and output dominate.

- **Problem 7 is the fair comparison:** both codes solve BSSN with the same gauge on the same grid for the same
  number of steps. pyNR's fused Numba kernels are about 15% faster per point and step than ML_BSSN (Kranc-generated
  C++) on this laptop. pyNR uses all 20 hardware threads (14 cores), the ET run 12 cores, so per core they are about equal.
- pyNR's ADM kernel costs half its BSSN kernel per point (fewer variables and no advection terms), which is why
  problems 4 and 6 look cheaper in pyNR: they use a different formulation there.
- For the 3D problems the laptop is 1.8–2.0× slower than the M3 Pro.
- The ET perturbed-black-hole cost includes the 2D elliptic solve of `IDAxiBrillBH` and a horizon search every 32
  steps.

## Where the remaining factor lies

A hand-tuned C++ ADM kernel would be roughly 3–5× faster per core. The gap
comes from the per-point scratch arrays (`g`, `gu`, `dg`, ...): LLVM cannot
prove they don't alias one another, so it keeps them in memory instead of
registers, and it cannot vectorise across grid points. The planned fix, the
first item on the roadmap, is to process a whole row of `k` at a
time, with every tensor operation as a short, SIMD-friendly loop over `k`.
This is the layout that code generators such as Kranc/McLachlan emit.
