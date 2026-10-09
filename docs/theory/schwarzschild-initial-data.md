(sec-schwarzschild-id)=
# Schwarzschild initial data in different coordinates

A single black hole of mass $M$ is one spacetime, but each coordinate system gives it different 3+1 data
$(\gamma_{ij}, K_{ij}, \alpha, \beta^i)$. Those data can be regular or singular at the horizon, and can cover or
miss the interior. This chapter derives the data for the four coordinate systems used in numerical relativity:
Schwarzschild (areal), isotropic, Kerr-Schild and Painlevé-Gullstrand. It shows how pyNR computes them in
{mod}`pynr.thorns.exact`, and how to check them. Background: chapter 1 and §6.1 of Alcubierre 2008
{cite:p}`Alcubierre2008`, and chapter 3 of Baumgarte & Shapiro 2010 {cite:p}`BaumgarteShapiro2010`.

Throughout, $G = c = 1$. $R$ is the areal radius (a sphere has area $4\pi R^2$), $r = \sqrt{x^2+y^2+z^2}$ is the
coordinate radius of the Cartesian grid, and $n_i = x_i/r$.

## 1. From a 4-metric to 3+1 data

Compare a given line element with the 3+1 form (Eq. {eq}`eq-adm-evolution` evolves its pieces):

$$
    ds^2 = -\alpha^2dt^2 + \gamma_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)
    \quad\Longleftrightarrow\quad
    g_{00} = -\alpha^2 + \beta_k\beta^k,\quad g_{0i} = \beta_i,\quad g_{ij} = \gamma_{ij}.
$$ (eq-sid-31-metric)

So the recipe has four steps:

1. $\gamma_{ij} = g_{ij}$ (the spatial block).
2. $\beta_i = g_{0i}$, then $\beta^i = \gamma^{ij}\beta_j$.
3. $\alpha = \sqrt{\beta_k\beta^k - g_{00}}$, or equivalently $\alpha = (-g^{00})^{-1/2}$.
4. $K_{ij}$ from the evolution equation of the metric, $\partial_t\gamma_{ij} = -2\alpha K_{ij} + D_i\beta_j + D_j\beta_i$.
   For a metric that does not depend on $t$ (every case below), this gives

   $$
       K_{ij} = \frac{1}{2\alpha}\left(D_i\beta_j + D_j\beta_i\right),
       \qquad D_i\beta_j = \partial_i\beta_j - \Gamma^k{}_{ij}\beta_k .
   $$ (eq-sid-K-stationary)

   Zero shift therefore means $K_{ij} = 0$: those slices are *time-symmetric*.

The data must satisfy the vacuum constraints (Eq. {eq}`eq-constraints-kernel`), $R + K^2 - K_{ij}K^{ij} = 0$ and
$D_j(K^{ij} - \gamma^{ij}K) = 0$. They hold automatically when the 4-metric is a vacuum solution, so evaluating them
on the grid tests the implementation and the resolution.

## 2. Schwarzschild (areal) coordinates

$$
    ds^2 = -\left(1 - \frac{2M}{R}\right)dt^2 + \frac{dR^2}{1 - 2M/R} + R^2 d\Omega^2 .
$$ (eq-sid-schwarzschild)

With $g_{0i} = 0$, the recipe gives

$$
    \alpha = \sqrt{1 - \frac{2M}{R}},\qquad \beta^i = 0,\qquad K_{ij} = 0,\qquad
    \gamma_{ij} = \delta_{ij} + \frac{2M}{R - 2M}\,n_i n_j \quad (r = R).
$$ (eq-sid-schwarzschild-31)

The Cartesian form follows from $dR = n_i\,dx^i$ and $\frac{dR^2}{1-2M/R} = dR^2 + \frac{2M}{R-2M}dR^2$.

These data are unusable on a grid that contains the horizon. $\gamma_{RR}\to\infty$ and $\alpha\to0$ at $R = 2M$,
and the slice does not continue inside: $R < 2M$ is not part of a $t = $ const slice. pyNR does not implement
them. They are the starting point for the other three systems, which are all coordinate changes of
Eq. {eq}`eq-sid-schwarzschild`.

## 3. Isotropic coordinates: the puncture

**Coordinate change.** Look for a radius $r$ in which the spatial metric is conformally flat,
$\gamma_{ij} = \psi^4\delta_{ij}$. Comparing the angular parts gives $R^2 = \psi^4 r^2$, so $R = r\psi^2$. The
radial parts then require $dR^2/(1 - 2M/R) = \psi^4 dr^2$. Try $\psi = 1 + M/2r$:

$$
    R = r\left(1 + \frac{M}{2r}\right)^2 = r + M + \frac{M^2}{4r},\qquad
    \frac{dR}{dr} = \psi\left(1 - \frac{M}{2r}\right),\qquad
    1 - \frac{2M}{R} = \frac{(1 - M/2r)^2}{\psi^2}.
$$ (eq-sid-iso-radius)

So $dR^2/(1-2M/R) = \psi^2(1 - M/2r)^2\,dr^2\cdot\psi^2/(1 - M/2r)^2 = \psi^4dr^2$, as required. The data are

$$
    \gamma_{ij} = \psi^4\delta_{ij},\quad \psi = 1 + \frac{M}{2r},\qquad K_{ij} = 0,\qquad
    \alpha = \frac{1 - M/2r}{1 + M/2r},\qquad \beta^i = 0.
$$ (eq-sid-isotropic)

This is Eq. {eq}`eq-schwarzschild-isotropic`, implemented by `schwarzschild_isotropic` in {mod}`pynr.thorns.exact`.

**Constraints.** With $K_{ij} = 0$, the momentum constraint is trivial. For a conformally flat metric,
$R = -8\psi^{-5}\nabla^2_{\rm flat}\psi$, so the Hamiltonian constraint becomes Laplace's equation,
$\nabla^2_{\rm flat}\psi = 0$. $\psi = 1 + M/2r$ is the monopole solution. A sum of such terms,
$\psi = 1 + \sum_A m_A/2|\mathbf x - \mathbf x_A|$, is again a solution: this is Brill-Lindquist data for several
black holes at rest (Brill & Lindquist 1963 {cite:p}`BrillLindquist1963`). Brandt & Brügmann 1997
{cite:p}`Brandt1997` added momentum and spin to it with the puncture method.

**Geometry: a wormhole.** $R(r)$ in Eq. {eq}`eq-sid-iso-radius` has a minimum at $r = M/2$, where $R = 2M$. This
minimal sphere is the horizon (the *throat*). The map $r \to M^2/4r$ leaves $R$ unchanged, so it is an isometry: the
region $r < M/2$ is a second copy of the exterior, and $r \to 0$ is a second spatial infinity, not a singularity.
The point $r = 0$ is called the **puncture**. Isotropic data never reach the physical singularity at $R = 0$.

**Lapse.** `ADMBase::initial_lapse` chooses between three options:

| keyword | $\alpha$ | at the throat $r = M/2$ | at $r\to0$ | used in |
|---|---|---|---|---|
| `exact` | $\frac{1 - M/2r}{1 + M/2r}$ (static) | 0 | $-1$ | stationary tests |
| `psi^-2` | $\psi^{-2}$ ("pre-collapsed") | $1/4$ | 0 | [problem 4](../problems/schwarzschild_1plog.md), the BSSN puncture run |
| `one` | 1 (geodesic slicing) | 1 | 1 | [problem 3](../problems/schwarzschild_geodesic.md) |

The static lapse is negative inside the throat, where time runs backwards on the other sheet. It is an exact
solution but useless for evolution. pyNR computes `psi^-2` as `gxx**-0.5`, which equals $\psi^{-2}$ only because
the data are conformally flat.

**On the grid.** $\psi$ diverges at $r = 0$. The parfiles place the puncture between grid points
(`xmin = -L - dx/2`), so the closest points sit at $r = \tfrac{\sqrt3}{2}\Delta x$, where
$\gamma_{xx} = \psi^4 \approx (M/2r)^4$ is large but finite. ADM cannot evolve these data for long
([problem 4](../problems/schwarzschild_1plog.md)). BSSN carries $W = \psi^{-2}\to0$ instead of $\psi$, so every
variable stays finite and the puncture can be evolved without excision. The slice then moves away from the second
infinity and settles on a *trumpet* (Hannam et al. 2007 {cite:p}`Hannam2007`).

## 4. Kerr-Schild coordinates (ingoing Eddington-Finkelstein)

**Coordinate change.** Replace the Schwarzschild time by $t = t_S + 2M\ln|R/2M - 1|$ and keep $R$. With
$f = 2M/R$, $dt_S = dt - \frac{f}{1-f}dR$, and

$$
    -(1-f)\left(dt - \tfrac{f}{1-f}dR\right)^2 + \frac{dR^2}{1-f}
    = -(1-f)\,dt^2 + 2f\,dt\,dR + (1+f)\,dR^2 ,
$$

so that

$$
    ds^2 = -\left(1 - \frac{2M}{R}\right)dt^2 + \frac{4M}{R}\,dt\,dR + \left(1 + \frac{2M}{R}\right)dR^2 + R^2d\Omega^2 .
$$ (eq-sid-ks-metric)

No component is singular at $R = 2M$. Ingoing light rays have $dR/dt = -1$, and the slices pass through the
horizon to the singularity at $R = 0$. In Cartesian coordinates ($r = R$) this is
$g_{\mu\nu} = \eta_{\mu\nu} + 2H\,l_\mu l_\nu$ with $H = M/r$ and $l_\mu = (1, n_i)$, the $a = 0$ case of
Eq. {eq}`eq-kerr-schild`.

**3+1 data.** From Eq. {eq}`eq-sid-31-metric` with $\gamma_{ij} = \delta_{ij} + 2Hn_in_j$ and $\beta_i = 2Hn_i$:

- The inverse metric follows from the Sherman-Morrison formula, because $n$ is a unit vector of the flat metric:
  $\gamma^{ij} = \delta^{ij} - \frac{2H}{1+2H}n^in^j$.
- Then $\beta^i = \frac{2H}{1+2H}n^i$, $\beta_k\beta^k = \frac{4H^2}{1+2H}$ and
  $\alpha^2 = 1 - 2H + \frac{4H^2}{1+2H} = \frac{1}{1+2H}$.

$$
    \gamma_{ij} = \delta_{ij} + \frac{2M}{r}n_in_j,\qquad
    \alpha = \frac{1}{\sqrt{1 + 2M/r}},\qquad
    \beta^i = \frac{2M/r}{1 + 2M/r}\,n^i .
$$ (eq-sid-ks-31)

Evaluating Eq. {eq}`eq-sid-K-stationary` gives the closed form

$$
    K_{ij} = \frac{2M\alpha}{r^2}\left[\delta_{ij} - \left(2 + \frac{M}{r}\right)n_in_j\right],\qquad
    K = \frac{2M\alpha^3}{r^2}\left(1 + \frac{3M}{r}\right).
$$ (eq-sid-ks-K)

Both were checked symbolically against Eq. {eq}`eq-sid-K-stationary`, together with the Hamiltonian constraint.
At the horizon $r = 2M$: $\alpha = 1/\sqrt2$, $\beta^r = 1/2$, and every quantity is finite. Only $r = 0$ is
singular, and there $H$, $\gamma_{ij}$ and $K_{ij}$ all diverge. pyNR evolves these data with excision
(`ADMEvolve::excision_radius` $\approx 0.9\,r_+$, [problem 5](../problems/kerr_schild.md)).

**How pyNR computes $K_{ij}$.** `kerr_schild` handles any spin $a$, for which the closed form is long. It
therefore evaluates Eq. {eq}`eq-sid-K-stationary` numerically from the analytic metric:
- It differentiates $\gamma_{ij}$ and $\beta_i$ with the 4th-order stencil of Eq. {eq}`eq-fd-stencils`, using a
  step $\epsilon = 10^{-4}$ in the *function arguments* (not on the grid).
- It builds $\Gamma_{kij}$ and $D_{(i}\beta_{j)}$ from those derivatives.

The truncation error is $\sim\epsilon^4 \approx 10^{-16}$ and the round-off error is $\sim 10^{-16}/\epsilon \approx 10^{-12}$
(relative, away from $r = 0$). Both are far below the grid's own $\mathcal O(\Delta x^4)$ error. For $a = 0$,
Eq. {eq}`eq-sid-ks-K` is an independent check (exercise 2).

## 5. Painlevé-Gullstrand coordinates

**Coordinate change.** Choose the time of observers who fall from rest at infinity,
$t = t_S + 2\sqrt{2MR} + 2M\ln\left|\frac{\sqrt{R/2M} - 1}{\sqrt{R/2M} + 1}\right|$. Then
$dt = dt_S + \frac{\sqrt{2M/R}}{1 - 2M/R}dR$, and

$$
    ds^2 = -dt^2 + \left(dR + \sqrt{\frac{2M}{R}}\,dt\right)^2 + R^2d\Omega^2 .
$$ (eq-sid-pg-metric)

**3+1 data.**

$$
    \gamma_{ij} = \delta_{ij},\qquad \alpha = 1,\qquad \beta^i = \sqrt{\frac{2M}{r}}\,n^i,\qquad
    K_{ij} = \sqrt{\frac{2M}{r^3}}\left(\delta_{ij} - \tfrac32 n_in_j\right),\qquad
    K = \tfrac32\sqrt{\frac{2M}{r^3}} .
$$ (eq-sid-pg-31)

The slices are **flat**. All the curvature of the black hole is in $K_{ij}$ and the shift: with $R = 0$, the
Hamiltonian constraint reads $K^2 = K_{ij}K^{ij}$. In units of $\sqrt{2M/r^3}$, $K^i{}_j$ has the eigenvalues
$-\tfrac12$ (radial) and $1, 1$ (tangential), so $K^2 = (\tfrac32)^2 = \tfrac14 + 1 + 1 = K_{ij}K^{ij}$
(exercise 3). Like Kerr-Schild, the slices are horizon-penetrating and end at $r = 0$. In both systems, outgoing
radial light moves at $dr/dt = \alpha\sqrt{\gamma^{rr}} - \beta^r$, which vanishes exactly at $r = 2M$. In
Painlevé-Gullstrand this happens because the shift reaches $\beta^r = 1 = \alpha$ there. pyNR does not implement
these data (exercise 4).

## 6. Comparison

| coordinates | $\gamma_{ij}$ | $\alpha$ | $\beta^i$ | $K_{ij}$ | horizon | covers interior | singular at | in pyNR |
|---|---|---|---|---|---|---|---|---|
| Schwarzschild | $\delta_{ij} + \frac{2M}{r-2M}n_in_j$ | $\sqrt{1-2M/r}$ | 0 | 0 | $r = 2M$, singular | no | $r = 2M$ | — |
| isotropic | $\psi^4\delta_{ij}$ | choice (§3) | 0 | 0 | $r = M/2$, minimal sphere | second exterior | $r = 0$ (puncture, an infinity) | `Schwarzschild/isotropic` |
| Kerr-Schild | $\delta_{ij} + \frac{2M}{r}n_in_j$ | $(1+2M/r)^{-1/2}$ | $\frac{2M/r}{1+2M/r}n^i$ | Eq. {eq}`eq-sid-ks-K` | $r = 2M$, regular | yes | $r = 0$ (physical) | `Kerr/Kerr-Schild` |
| Painlevé-Gullstrand | $\delta_{ij}$ | 1 | $\sqrt{2M/r}\,n^i$ | Eq. {eq}`eq-sid-pg-31` | $r = 2M$, regular | yes | $r = 0$ (physical) | — (exercise 4) |

The areal radius is $R = r$ in every row except isotropic, where $R = r(1 + M/2r)^2$.

(sec-sid-cartesian)=
## 7. Putting spherical data on the Cartesian grid

pyNR stores **only Cartesian components**. `ADMBase.U[0:6]` holds `gxx gxy gxz gyy gyz gzz`
(the table *Layout of the ADM state array* in [The 3+1 split](three-plus-one.md)) at the grid points $(x_i, y_j, z_k)$. There is no `grr`. The formulas in §2-§5 are written
with spherical components ($\gamma_{rr}$, $\gamma_{\theta\theta}$, $\beta^r$, $K_{rr}$, …), so they must be
**transformed as tensors** before they go into the grid. Copying $\gamma_{rr}$ into `gxx` would be wrong everywhere
except on the $x$ axis.

**Why Cartesian.** The evolution kernels ({ref}`sec-discretisation`) take every derivative with the same Cartesian
stencils. Spherical coordinates would add $1/r$ and $1/\sin\theta$ terms, which are singular at the origin and on the
axis. A Cartesian grid has neither problem, and it handles non-spherical data (spin, perturbations) without change.
The price is that spherical symmetry is no longer manifest in the components.

**The transformation.** Every spherically symmetric spatial metric can be written with two radial functions $A(r)$
and $B(r)$:

$$
    dl^2 = A(r)\,dr^2 + B(r)\,r^2 d\Omega^2 .
$$ (eq-sid-spherical-metric)

On the grid $x^i = r\,n^i$, so $dr = n_i\,dx^i$ and $r^2d\Omega^2 = \delta_{ij}dx^idx^j - dr^2$. Substituting gives

$$
    \gamma_{ij} = A\,n_in_j + B\,(\delta_{ij} - n_in_j)
    = B\,\delta_{ij} + (A - B)\,n_in_j .
$$ (eq-sid-cartesian-metric)

The two pieces are the radial projector $n_in_j$ and the tangential projector $\delta_{ij} - n_in_j$. Any symmetric
tensor with spherical symmetry, including $K_{ij}$, has the same form:

$$
    K_{ij} = K_{rr}\,n_in_j + K_T\,(\delta_{ij} - n_in_j),\qquad K_T = K_{\theta\theta}/r^2 .
$$ (eq-sid-cartesian-K)

Vectors follow from $\partial x^i/\partial r = n^i$: a radial shift $\beta^r$ becomes $\beta^i = \beta^r n^i$, and
$\beta_i = \beta_r n_i$. The lapse is a scalar and is copied unchanged.

Written out, the six stored metric components are

$$
    \texttt{gxx} = B + (A-B)\frac{x^2}{r^2},\quad
    \texttt{gxy} = (A-B)\frac{xy}{r^2},\quad
    \texttt{gyy} = B + (A-B)\frac{y^2}{r^2},\ \dots
$$ (eq-sid-gxx)

with the other components following the same pattern.

| coordinates | $A = \gamma_{rr}$ | $B = \gamma_{\theta\theta}/r^2$ | $\beta^r$ | $K_{rr}$ | $K_T = K_{\theta\theta}/r^2$ |
|---|---|---|---|---|---|
| Schwarzschild | $(1 - 2M/r)^{-1}$ | 1 | 0 | 0 | 0 |
| isotropic | $\psi^4$ | $\psi^4$ | 0 | 0 | 0 |
| Kerr-Schild | $1 + 2M/r$ | 1 | $\frac{2M/r}{1+2M/r}$ | $-\frac{2M\alpha}{r^2}\left(1 + \frac Mr\right)$ | $\frac{2M\alpha}{r^2}$ |
| Painlevé-Gullstrand | 1 | 1 | $\sqrt{2M/r}$ | $-\tfrac12\sqrt{2M/r^3}$ | $\sqrt{2M/r^3}$ |

Inserting these rows into Eqs. {eq}`eq-sid-cartesian-metric` and {eq}`eq-sid-cartesian-K` gives
Eqs. {eq}`eq-sid-schwarzschild-31`, {eq}`eq-sid-isotropic`, {eq}`eq-sid-ks-31`, {eq}`eq-sid-ks-K` and
{eq}`eq-sid-pg-31`.

Two special cases are worth knowing:

- **Isotropic data need no transformation** ($A = B$): $\gamma_{ij} = \psi^4\delta_{ij}$, so `gxx = gyy = gzz` $=\psi^4$
  and the off-diagonal components vanish. `schwarzschild_isotropic` sets exactly that. This is also why
  `initial_lapse = "psi^-2"` can be computed as `gxx**-0.5`.
- **Kerr-Schild is Cartesian by construction.** $g_{\mu\nu} = \eta_{\mu\nu} + 2H\,l_\mu l_\nu$ is already written
  with Cartesian $l_i$. For spin $a \ne 0$ the data are not spherically symmetric, $l_i$ is not $n_i$, and
  Eq. {eq}`eq-sid-cartesian-metric` does not apply. `_kerr_schild_metric` builds $\gamma_{ij} = \delta_{ij} + 2Hl_il_j$
  directly. Its $r$ is the Kerr radius, the root of $r^4 - (\rho^2 - a^2)r^2 - a^2z^2 = 0$ with
  $\rho^2 = x^2+y^2+z^2$; for $a = 0$ it equals $\rho$.

**Reading the output.** The same projections work backwards:

$$
    \gamma_{rr} = n^in^j\gamma_{ij},\qquad
    \gamma_{\theta\theta} = r^2\,\hat e_\theta^{\,i}\hat e_\theta^{\,j}\gamma_{ij},\qquad
    K_{rr} = n^in^jK_{ij},\qquad
    \beta^r = n_i\beta^i ,
$$ (eq-sid-project)

where $\hat e_\theta$ is the flat unit vector in the $\theta$ direction. On the positive $x$ axis, $n = (1,0,0)$ and
$\hat e_\theta = (0,0,-1)$, so

$$
    \texttt{gxx} = \gamma_{rr},\qquad \texttt{gyy} = \texttt{gzz} = \gamma_{\theta\theta}/r^2,\qquad
    \texttt{kxx} = K_{rr},\qquad \texttt{betax} = \beta^r .
$$ (eq-sid-xaxis)

1D output along $x$ (`gxx.x.dat`, `kxx.x.dat`, as plotted with kuibit) therefore shows the radial components directly.
Along any other line, use Eq. {eq}`eq-sid-project`. For example, on the diagonal $x = y = z$,
$\texttt{gxx} = \frac13\gamma_{rr} + \frac23\gamma_{\theta\theta}/r^2$.

**What `Exact.initial_data` does.** It evaluates the model on the full 3D coordinate arrays and writes the
components into the state array:

```python
X, Y, Z = grid.meshgrid()                    # Cartesian coordinates of every point, ghosts included
g, K, alp, beta = solution(0.0, X, Y, Z)     # g[6], K[6]: packed Cartesian components (PAIR order)
U[0:6] = g                                   # gxx gxy gxz gyy gyz gzz
U[6:12] = K                                  # kxx ... kzz
U[12] = alp                                  # or psi^-2 / one, by ADMBase::initial_lapse
U[13:16] = beta                              # contravariant beta^x, beta^y, beta^z, if initial_shift = "exact"
```

A new spherically symmetric model needs only $A$, $B$, $K_{rr}$, $K_T$, $\alpha$ and $\beta^r$ as functions of $r$:

```python
def spherical_to_cartesian(A, B, Krr, KT, alp, beta_r, x, y, z):
    r = np.sqrt(x * x + y * y + z * z)
    n = (x / r, y / r, z / r)
    g, K = np.empty((6, *r.shape)), np.empty((6, *r.shape))
    for c, (i, j) in enumerate(PAIR):                     # PAIR: (0,0) (0,1) (0,2) (1,1) (1,2) (2,2)
        d = 1.0 if i == j else 0.0
        g[c] = B * d + (A - B) * n[i] * n[j]              # Eq. (sid-cartesian-metric)
        K[c] = KT * d + (Krr - KT) * n[i] * n[j]          # Eq. (sid-cartesian-K)
    beta = np.array([beta_r * n[i] for i in range(3)])    # beta^i = beta^r n^i
    return g, K, alp, beta
```

With the Kerr-Schild row of the table and $M = 1$, this function reproduces `kerr_schild(spin=0)` at random points:
$\gamma_{ij}$ and $\beta^i$ to round-off, and $K_{ij}$ to $7\times10^{-12}$ (relative). The difference is the
finite-difference error of `kerr_schild`'s numerical $K_{ij}$ (§4).

Points with $r = 0$ give $n = 0/0$. Keep the centre between grid points, or inside the excision region, where
ADMEvolve replaces non-finite values with flat space.

## 8. Pseudocode: analytic metric to grid data

This is the generic version of what `kerr_schild` does. It works for any stationary metric given as functions
$\gamma_{ij}(\mathbf x)$, $\beta_i(\mathbf x)$ and $\alpha(\mathbf x)$:

```text
function ADM_DATA(metric, X, Y, Z, ε = 1e-4):
    γ, β_low, α ← metric(X, Y, Z)                         # γ_ij, β_i = g_0i, α
    β_up ← γ⁻¹ · β_low                                     # β^i
    for m in x, y, z:                                      # 4th-order derivative in the arguments
        ∂_m γ, ∂_m β_low ← Σ_s w_s · metric(X + s ε e_m) / (12 ε),  (s, w_s) ∈ {(−2,1), (−1,−8), (1,8), (2,−1)}
    for each pair (i ≤ j):
        Γβ  ← Σ_k ½ (∂_i γ_kj + ∂_j γ_ki − ∂_k γ_ij) β^k   # Γ_kij β^k = Γ^k_ij β_k
        K_ij ← (∂_i β_j + ∂_j β_i − 2 Γβ) / (2α)          # Eq. (sid-K-stationary)
    return γ, K, α, β_up
```

Then check that `ADMConstraints::H` and `ADMConstraints::M` converge to zero at 4th order as $\Delta x$ is halved,
away from excised or puncture points.

## Exercises

1. Show that $r\to M^2/4r$ maps $\gamma_{ij} = \psi^4\delta_{ij}$ to itself, and that the static lapse changes sign
   under it. Where is the minimal surface of the $t = $ const slice?
2. Run [problem 5](../problems/kerr_schild.md) with `Kerr_KerrSchild__spin = 0`. Compare `ADMBase::kxx` along the
   $x$ axis with Eq. {eq}`eq-sid-ks-K` ($n = (1,0,0)$ gives $K_{xx} = -\frac{2M\alpha}{r^2}(1 + \frac{M}{r})$).
3. Verify Eq. {eq}`eq-sid-pg-31`. Compute $D_{(i}\beta_{j)}$ for $\beta_i = \sqrt{2M/r}\,n_i$ on flat space, then
   check $K^2 = K_{ij}K^{ij}$ and $\partial_j K^{ij} = \partial^i K$.
4. Add a model `Schwarzschild/Painleve-Gullstrand` to {mod}`pynr.thorns.exact` using Eq. {eq}`eq-sid-pg-31`.
   Evolve it with excision and `lapse_evolution_method = "static"`, and compare its constraint violation and
   lifetime with Kerr-Schild at the same $\Delta x$.
5. Use the transformations of §4 and §5 to show that Kerr-Schild and Painlevé-Gullstrand data on the same slice
   $t = 0$ are *different slices* of the same spacetime. Hint: compare $t_{KS}$ and $t_{PG}$ as functions of $R$
   at fixed $t_S$.
