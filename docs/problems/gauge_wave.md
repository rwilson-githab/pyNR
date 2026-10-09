# 1. Gauge wave

**File:** `par/gauge_wave.par`

## Physics

Minkowski spacetime in the coordinates

$$ ds^2 = -H\,dt^2 + H\,dx^2 + dy^2 + dz^2,\qquad H = 1 - A\sin\frac{2\pi(x-t)}{d}. $$ (eq-gauge-wave-metric)

Equation {eq}`eq-gauge-wave-metric` is flat space, but the lapse $\alpha = \sqrt{H}$, the metric
$\gamma_{xx} = H$ and the curvature $K_{xx} = -\frac{\pi A}{d}\cos\frac{2\pi(x-t)}{d}/\sqrt{H}$
all travel along $x$. The harmonic slicing $\partial_t\alpha = -\alpha^2 K$ is
consistent with this exact solution. It is test 1 of the "Apples with Apples"
suite (Alcubierre et al. 2004 {cite:p}`Alcubierre2004`;
[doi:10.1088/0264-9381/21/2/019](https://doi.org/10.1088/0264-9381/21/2/019),
[arXiv:gr-qc/0305023](https://arxiv.org/abs/gr-qc/0305023)).

(sec-gauge-wave-31)=
## The 3+1 form of the gauge wave

**Reading off the 3+1 variables.** Compare Eq. {eq}`eq-gauge-wave-metric` with the general 3+1 line element
$ds^2 = -\alpha^2dt^2 + \gamma_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)$. With $\varphi = 2\pi(x-t)/d$ and
$H = 1 - A\sin\varphi$,

$$
    \alpha = \sqrt H,\qquad \beta^i = 0,\qquad
    \gamma_{ij} = \operatorname{diag}(H, 1, 1),\qquad \gamma^{ij} = \operatorname{diag}(1/H, 1, 1).
$$ (eq-gw-31-vars)

Everything depends on $t$ and $x$ only through $x - t$, so $\partial_t H = -\partial_x H \equiv -H'$, with
$H' = -\frac{2\pi A}{d}\cos\varphi$.

**Extrinsic curvature.** With zero shift, the evolution equation for the metric,
$\partial_t\gamma_{ij} = -2\alpha K_{ij}$, *defines* $K_{ij}$:

$$
    K_{xx} = -\frac{\partial_t\gamma_{xx}}{2\alpha} = \frac{H'}{2\sqrt H}
           = -\frac{\pi A}{d}\,\frac{\cos\varphi}{\sqrt H},\qquad
    \text{all other } K_{ij} = 0,\qquad
    K = \gamma^{xx}K_{xx} = \frac{H'}{2H^{3/2}} .
$$ (eq-gw-31-K)

**The slice is flat.** The 3-metric $H(x)\,dx^2 + dy^2 + dz^2$ becomes $d\xi^2 + dy^2 + dz^2$ under
$d\xi = \sqrt H\,dx$, so $R_{ij} = 0$. All the dynamics is in $\alpha$ and $K_{ij}$: the slices are flat
hypersurfaces of Minkowski space that wave in time. This is why the test is a pure gauge test.

**What ADMEvolve integrates.** Only $\gamma_{xx}$, $K_{xx}$ and $\alpha$ change. The ADM system
(Eq. {eq}`eq-adm-evolution`) with harmonic slicing reduces to a 1+1 system:

$$
\begin{aligned}
    \partial_t\gamma_{xx} &= -2\alpha K_{xx},\\
    \partial_t K_{xx} &= -\left(\partial_x^2\alpha - \Gamma^x{}_{xx}\partial_x\alpha\right)
        + \alpha\left(\underbrace{R_{xx}}_{0} + K K_{xx} - 2K_{xx}K^x{}_x\right)
      = -\partial_x^2\alpha + \frac{\partial_x\gamma_{xx}}{2\gamma_{xx}}\partial_x\alpha
        - \alpha\frac{K_{xx}^2}{\gamma_{xx}},\\
    \partial_t\alpha &= -\alpha^2 K = -\alpha^2\frac{K_{xx}}{\gamma_{xx}},
\end{aligned}
$$ (eq-gw-31-system)

using $\Gamma^x{}_{xx} = \partial_x\gamma_{xx}/2\gamma_{xx}$ and $K = K^x{}_x = K_{xx}/\gamma_{xx}$.

**Checking that Eq. {eq}`eq-gw-31-vars` solves it.**

- *Metric:* $-2\alpha K_{xx} = -2\sqrt H\cdot H'/2\sqrt H = -H' = \partial_t H$. ✓
- *Lapse:* $\partial_t\sqrt H = -H'/2\sqrt H$, and $-\alpha^2K = -H\cdot H'/2H^{3/2} = -H'/2\sqrt H$. ✓
  More generally, with zero shift harmonic slicing and $\partial_t\ln\sqrt\gamma = -\alpha K$ give
  $\partial_t(\alpha/\sqrt\gamma) = 0$. The lapse is therefore locked to the volume element, $\alpha = \sqrt\gamma$,
  which is exactly the gauge-wave relation $\alpha^2 = \gamma_{xx} = H$.
- *Curvature:* with $\alpha = \sqrt H$, $\partial_x^2\alpha - \Gamma^x{}_{xx}\partial_x\alpha
  = \frac{H''}{2\sqrt H} - \frac{H'^2}{2H^{3/2}}$ and $\alpha K_{xx}^2/\gamma_{xx} = \frac{H'^2}{4H^{3/2}}$, so the
  right-hand side is $-\frac{H''}{2\sqrt H} + \frac{H'^2}{4H^{3/2}}$. The left-hand side is
  $\partial_t K_{xx} = -\partial_x\!\left(\frac{H'}{2\sqrt H}\right)$, which is the same expression. ✓

**Constraints.** $\mathcal H = R + K^2 - K_{ij}K^{ij} = 0 + K_{xx}^2/H^2 - K_{xx}^2/H^2 = 0$. For the momentum
constraint, $D_jK^j{}_x = \partial_x K^x{}_x + \Gamma^j{}_{jx}K^x{}_x - \Gamma^m{}_{xx}K^x{}_m = \partial_x K$, so
$\mathcal M_x = D_jK^j{}_x - \partial_x K = 0$.

**Why the constraints stay at round-off on the grid.** The grid is periodic with 3 points in $y$ and $z$, so all
$y$ and $z$ differences vanish exactly. In the discrete Ricci tensor Eq. {eq}`eq-disc-ricci` only $R^h_{xx}$
survives, and it equals $(D_x g^{-1})^{xx}\Gamma_{xxx} - \tfrac12(D_x g^{-1})^{xx}D_x g_{xx} = 0$ *identically*,
because $\Gamma_{xxx} = \tfrac12 D_x g_{xx}$ uses the same discrete $D_x$
(see {ref}`sec-discretisation`). $K^2 - K_{ij}K^{ij}$ cancels algebraically in the same way. What the test measures is
therefore the accuracy of Eq. {eq}`eq-gw-31-system` along $x$ (4th order), and whether ADM's weak hyperbolicity
lets the gauge mode grow (exercise 1).

## Setup

The grid is periodic in all directions and effectively 1D (3 points in $y$
and $z$). The amplitude is $A = 0.1$, the wavelength $d = 1$, and the run lasts
ten crossing times. {numref}`tab-settings-gauge-wave` lists the complete
settings of the reference run, read from its `parameters.par`.

## Expected results

```{include} ../figures/gauge_wave.md
```

- $\max\alpha$ stays at $\sqrt{1.1} = 1.0488$ and $\min\alpha$ at $\sqrt{0.9} = 0.9487$.
  The measured values drift by about $10^{-4}$ over 10 crossings.
- $H = 0$ to round-off, because the metric only depends on $x$.
- **Convergence:** the error in $\alpha$ against the exact solution falls by
  $\approx 16$ each time $\Delta x$ is halved (4th order). This is checked by
  `tests/test_evolution.py::test_gauge_wave_converges_4th_order` and shown in
  {numref}`fig-gauge-wave` (right) and `notebooks/01_gauge_wave.ipynb`.

## How to run

**1. Environment (once).** pyNR needs Python ≥ 3.10 with NumPy, Numba, SciPy and h5py (kuibit and
matplotlib for the plots). There is no compiler to configure and nothing to build: Numba compiles the kernels to
machine code through LLVM the first time they run and caches the result next to the sources, so only the first run
in a fresh checkout spends a few seconds compiling.

```bash
python3 -m venv .venv && source .venv/bin/activate     # or a conda env: see Installation
pip install -e ".[viz]"
```

| setting | how | default |
|---|---|---|
| threads | `export NUMBA_NUM_THREADS=8` | all cores |
| JIT cache location (clusters, read-only checkouts) | `export NUMBA_CACHE_DIR=$HOME/.numba_cache` | next to the sources |
| kernels | `pynr run … --backend numpy` (readable NumPy reference kernels, much slower) | Numba |
| output root | `export PYNR_OUTPUT_DIR=/scratch/runs` | `simulations/` of the checkout |

See [Installation](../installation.md) for conda, Jupyter and Codespaces.

**2. Run.**

```bash
pynr run par/gauge_wave.par                       # → simulations/gauge_wave/
pynr run par/gauge_wave.par --set Thorn::param=value   # change a parameter without editing the file
```

**3. Look at the output** (Carpet formats, readable by kuibit):

```python
from kuibit.simdir import SimDir
from pynr.paths import run_dir
sd = SimDir(run_dir("gauge_wave"))
```

**Cost.** 3.2 s on an Apple M3 Pro (12 threads); 13 s on an Intel i7-12700H laptop (20 threads; first run of the session, including Numba compilation).

**Einstein Toolkit counterpart:** {ref}`gauge wave <et-cmp-gauge-wave>` in the comparison with the Einstein Toolkit.

```{include} ../_generated/arch/workflow_gauge_wave.md
```

## Exercises

1. $A = 0.5$: run to $t = 100$. The ADM system is known to be unstable on
   large-amplitude gauge waves, so measure when the run fails.
2. Measure the convergence order with `MoL::ODE_Method = "RK2"`.
