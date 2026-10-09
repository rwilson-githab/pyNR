(sec-discretisation)=
# Discretising the ADM update

This chapter follows one time step of an ADMEvolve run through the code and states, at each stage, the discrete
equation that is solved and its error. It covers the spatial operators ({mod}`pynr.kernels.fd`), the ADM right-hand
side ({mod}`pynr.kernels.adm`), Kreiss-Oliger dissipation ({mod}`pynr.kernels.dissipation`), the outer boundary
({mod}`pynr.kernels.boundary`) and the Runge-Kutta step ({mod}`pynr.utils.integrators`, {mod}`pynr.thorns.mol`).
The continuum equations are Eq. {eq}`eq-adm-evolution` in [The 3+1 split](three-plus-one.md). The textbook
treatment is chapter 9 of Alcubierre 2008 {cite:p}`Alcubierre2008`.

## 1. From PDE to ODE: the method of lines

Collect the 16 ADM variables (the table *Layout of the ADM state array* in [The 3+1 split](three-plus-one.md)) into $u = (\gamma_{ij}, K_{ij}, \alpha, \beta^i)$. The
continuum system is $\partial_t u = F[u]$, where $F$ contains first and second spatial derivatives of $u$. The
**method of lines** (MoL) discretises space first. Every variable becomes a grid function
$u_{\mathbf i}(t) = u(t, \mathbf x_{\mathbf i})$, every derivative becomes a finite-difference operator, and the PDE
becomes a large system of ordinary differential equations

$$
    \frac{d u_{\mathbf i}}{dt} = F_h[u]_{\mathbf i} + Q_h[u]_{\mathbf i},
    \qquad \mathbf i \in \text{interior},
$$ (eq-disc-semidiscrete)

closed by boundary conditions on the ghost points. $F_h$ is the ADM right-hand side with derivatives replaced by
stencils (§3), and $Q_h$ is the artificial dissipation (§4). Time is discretised separately, by a Runge-Kutta
method (§6). Space and time can therefore be analysed one at a time: the spatial error is $\mathcal O(h^4)$ and the
time error is $\mathcal O(\Delta t^4)$.

## 2. The grid

{mod}`pynr.cactus.grid` samples the physical domain $[x_{\min}, x_{\max}]$ with spacing $h_x$ (likewise for $y$
and $z$). It adds $n_g$ ghost points (`Driver::ghost_size`, default 3) on every side:

$$
    x_i = x_{\min} + (i - n_g)\,h_x,\qquad i = 0, \dots, n_x - 1,\qquad n_x = n_x^{\rm phys} + 2n_g .
$$ (eq-disc-grid)

All 16 variables are stored in one array `U[16, nx, ny, nz]`, with $z$ (index `k`) contiguous in memory. Kernels
receive the inverse spacings `idx` $= (1/h_x, 1/h_y, 1/h_z)$. Every RHS kernel loops over the **interior**
$n_g \le i < n_x - n_g$ (likewise for `j`, `k`). Stencils reach up to $n_g$ points beyond the point they update,
so they read ghost values but never write them.

| operator | half-width | needs |
|---|---|---|
| 4th-order first and second derivatives (§3) | 2 | $n_g \ge 2$ |
| Kreiss-Oliger dissipation (§4) | 3 | $n_g \ge 3$ (checked by `Dissipation.setup`) |

## 3. Finite-difference operators ({mod}`pynr.kernels.fd`)

### 3.1 First derivative

Let $D_x$ be the five-point centred operator of Eq. {eq}`eq-fd-stencils`,
$D_x f_i = (f_{i-2} - 8f_{i-1} + 8f_{i+1} - f_{i+2})/12h$. Taylor-expand
$f_{i\pm m} = \sum_n (\pm m h)^n f^{(n)}/n!$. Odd powers survive in $f_{i+m} - f_{i-m}$, so

$$
    8(f_{i+1} - f_{i-1}) - (f_{i+2} - f_{i-2})
    = 12 h f' + \underbrace{(16 - 16)}_{0}\frac{h^3}{6}f''' + (16 - 64)\frac{h^5}{120}f^{(5)} + \dots
$$

and

$$
    D_x f = f' - \frac{h^4}{30}f^{(5)} + \mathcal O(h^6).
$$ (eq-disc-d1-error)

The weights $(1,-8,8,-1)/12$ are the unique choice that cancels the $h^3$ term with half-width 2.

### 3.2 Second derivative

The compact five-point operator $D_{xx} f_i = (-f_{i-2} + 16f_{i-1} - 30f_i + 16f_{i+1} - f_{i+2})/12h^2$ uses the
even parts $f_{i+m} + f_{i-m}$. The same expansion gives

$$
    D_{xx} f = f'' - \frac{h^4}{90}f^{(6)} + \mathcal O(h^6).
$$ (eq-disc-d2-error)

Note that $D_{xx} \ne D_x D_x$. Applying $D_x$ twice would give a wider nine-point stencil that cannot see the
shortest grid mode (§3.4). The code always uses the compact $D_{xx}$, $D_{yy}$, $D_{zz}$ for pure second
derivatives.

### 3.3 Mixed derivatives

Mixed derivatives are tensor products of two first-derivative operators, $D_{xy} = D_x D_y$:

$$
    D_{xy} f_{ij} = \frac{1}{144\,h_x h_y}\sum_{p,q\in\{-2,-1,1,2\}} c_p c_q\, f_{i+p,\,j+q},
    \qquad c = (1, -8, 8, -1),
$$ (eq-disc-mixed)

a 16-point stencil in the $xy$ plane with error
$-\frac{h_x^4}{30}\partial_x^5\partial_y f - \frac{h_y^4}{30}\partial_x\partial_y^5 f$. The code builds it as
"$D_x$ applied to $D_y$" (`dxy` calls `dy` at four $x$ offsets).

The union of all stencils needed at one point (three axes, three mixed planes) covers
$1 + 3\cdot4 + 3\cdot16 = 61$ grid points per variable.

### 3.4 What the operators do to Fourier modes

On an infinite grid, insert $f_j = e^{ikx_j}$ and write $\theta = kh \in [0, \pi]$. Then $\theta = \pi$ is the
shortest mode the grid can carry, the "sawtooth" $(-1)^j$. The operators multiply the mode by their **symbols**:

$$
    D_x \to \frac{i}{h}\,\frac{8\sin\theta - \sin 2\theta}{6},\qquad
    D_{xx} \to -\frac{1}{h^2}\,\frac{30 - 32\cos\theta + 2\cos 2\theta}{12}.
$$ (eq-disc-symbols)

For small $\theta$ these are $ik(1 - \theta^4/30)$ and $-k^2(1 - \theta^4/90)$, which reproduces
Eqs. {eq}`eq-disc-d1-error` and {eq}`eq-disc-d2-error`. Two facts matter for stability:

- The symbol of $D_x$ is purely imaginary, so centred differences neither damp nor amplify a mode. It vanishes at
  $\theta = \pi$: **the first derivative does not see the sawtooth**. Error that collects there is not propagated
  away, and nonlinear terms such as $K_{ik}K^k{}_j$ keep feeding it. Section 4 adds a term that removes it.
- The largest magnitudes are $\max_\theta|D_x| = 1.372/h$ (at $\theta \approx 1.80$) and
  $\max_\theta|D_{xx}| = 16/3h^2$ (at $\theta = \pi$). They set the time-step limit (§6.2).

## 4. The ADM right-hand side ({func}`pynr.kernels.adm.adm_rhs`)

### 4.1 Which quantities are differenced

The kernel takes finite differences of **stored variables only**. Everything else, including the inverse metric,
Christoffel symbols, the Ricci tensor and the trace $K$, is built algebraically at the point from those
derivatives. At each interior point it evaluates:

| stored variable | operators | stencil evaluations |
|---|---|---|
| $\gamma_{ij}$ (6) | $D_m$, $D_{mn}$ | $6\times(3+6) = 54$ |
| $K_{ij}$ (6) | $D_m$ (advection) | $6\times3 = 18$ |
| $\alpha$ | $D_m$, $D_{mn}$ | $3 + 6 = 9$ |
| $\beta^i$ (3) | $D_m$ | $3\times 3 = 9$ |

That is 90 stencils per point. Symmetric tensors are computed for $a \le b$ only (the `PAIR` table) and mirrored.

### 4.2 Geometry at a point (`_geometry`)

Write $\partial$ for exact derivatives and $D$ for the operators of §3. With
$g_{ab} = \gamma_{ab}(\mathbf x_{\mathbf i})$, the kernel computes in order:

1. **Inverse metric** by cofactors, $g^{ab} = \operatorname{cof}(g)^{ab}/\det g$. This is exact; no
   differencing is involved.
2. **Christoffel symbols** from the differenced metric,

   $$
       \Gamma_{lab} = \tfrac12\left(D_a g_{lb} + D_b g_{la} - D_l g_{ab}\right),\qquad
       \Gamma^k{}_{ab} = g^{kl}\Gamma_{lab}.
   $$ (eq-disc-christoffel)

3. **Derivative of the inverse metric**, obtained algebraically and not by differencing $g^{kl}$:

   $$
       (D_m g^{-1})^{kl} := -g^{kp}g^{lq}\,D_m g_{pq}.
   $$ (eq-disc-dginv)

   In the continuum this is the identity $\partial_m\gamma^{kl} = -\gamma^{kp}\gamma^{lq}\partial_m\gamma_{pq}$.
   On the grid it is a *choice*: $D_m(g^{kl})$ would differ from it at $\mathcal O(h^4)$, and would need the inverse
   metric stored on the grid.
4. **Ricci tensor** from Eq. {eq}`eq-ricci`, with the two contractions of Eq. {eq}`eq-dgamma-contractions`
   discretised term by term:

   $$
       R^h_{ab} = (D_k g^{-1})^{kl}\,\Gamma_{lab}
         + \tfrac12 g^{kl}\left(D_{ka}g_{lb} + D_{kb}g_{la} - D_{kl}g_{ab} - D_{ab}g_{kl}\right)
         - \tfrac12 (D_b g^{-1})^{kl} D_a g_{kl}
         + \Gamma^l{}_{lk}\Gamma^k{}_{ab} - \Gamma^l{}_{bk}\Gamma^k{}_{la}.
   $$ (eq-disc-ricci)

   Second derivatives of the metric enter only through the $D_{mn}$ stencils, each applied once to a stored
   variable. The Christoffel symbols are never differenced, so this needs no extra storage and one pass over the
   grid.

Each term in Eq. {eq}`eq-disc-ricci` is a product of $\mathcal O(h^4)$-accurate factors, so
$R^h_{ab} = R_{ab} + \mathcal O(h^4)$. Because the formula is algebraic in the discrete derivatives, some
continuum identities still hold **exactly** on the grid. The gauge wave
({ref}`sec-gauge-wave-31`) is an example: there $R^h_{ab}$ is zero to
round-off at every resolution.

### 4.3 The evolution equations at a point

With $\beta^m$, $D_m\beta^a$, $\alpha$, $D_m\alpha$, $D_{mn}\alpha$ and $K_{ab}$ read at the point, the kernel
forms $K = g^{ab}K_{ab}$, $K^a{}_b = g^{am}K_{mb}$ and, for each of the 6 pairs $a \le b$,

$$
\begin{aligned}
    \frac{d\gamma_{ab}}{dt} &= -2\alpha K_{ab}
        + \beta^m D_m g_{ab} + g_{mb} D_a\beta^m + g_{am} D_b\beta^m, \\
    \frac{dK_{ab}}{dt} &= -\left(D_{ab}\alpha - \Gamma^m{}_{ab}D_m\alpha\right)
        + \alpha\left(R^h_{ab} + K K_{ab} - 2K_{am}K^m{}_b\right)
        + \beta^m D_m K_{ab} + K_{mb}D_a\beta^m + K_{am}D_b\beta^m .
\end{aligned}
$$ (eq-disc-adm-rhs)

This is Eq. {eq}`eq-adm-evolution` with every $\partial$ replaced by $D$. Two properties of this discretisation
matter later:

- **Advection is centred.** $\beta^m D_m$ uses the same centred operator as every other derivative. That is adequate
  here because the ADM runs in pyNR have zero or static, moderate shift.
- **No constraint is added.** The RHS is exactly the ADM system. This system is only weakly hyperbolic, and the
  growth of $\mathcal O(h^4)$ constraint errors is the instability that the ADM problems show.

The gauge equations are

$$
    \frac{d\alpha}{dt} =
    \begin{cases}
        0 & \texttt{static}\\
        -\alpha^2 K + [\beta^m D_m\alpha] & \texttt{harmonic}\\
        -2\alpha K + [\beta^m D_m\alpha] & \texttt{1+log}
    \end{cases},
    \qquad \frac{d\beta^i}{dt} = 0,
$$ (eq-disc-gauge)

where the bracket is present when `ADMEvolve::advect_lapse = yes`.

**Excision.** For $|\mathbf x_{\mathbf i}|^2 <$ `excision_r2` the kernel writes $du/dt = 0$ and skips the point.
The values there stay frozen at the (sanitised) initial data, and they are read by the stencils of nearby points
outside the excised region.

## 5. Kreiss-Oliger dissipation ({func}`pynr.kernels.dissipation.add_ko_dissipation`)

### 5.1 The operator

Let $D_+ f_i = (f_{i+1} - f_i)/h$ and $D_- f_i = (f_i - f_{i-1})/h$, so that $D_+D_- f_i = (f_{i+1} - 2f_i + f_{i-1})/h^2$
is the three-point Laplacian. Kreiss-Oliger dissipation of order $2r$ (Kreiss & Oliger 1973
{cite:p}`KreissOliger1973`) is

$$
    Q_h u = (-1)^{r+1}\,\frac{\epsilon}{2^{2r}}\sum_d h_d^{2r-1}\left(D_+^{(d)}D_-^{(d)}\right)^{r} u .
$$ (eq-disc-ko-general)

pyNR uses $r = 3$ (sixth differences), so $(-1)^{r+1} = +1$ and $2^{2r} = 64$. Expanding $(D_+D_-)^3$ along one
direction gives the binomial coefficients of $(1-1)^6$:

$$
    h^6 (D_+D_-)^3 u_i = u_{i-3} - 6u_{i-2} + 15u_{i-1} - 20u_i + 15u_{i+1} - 6u_{i+2} + u_{i+3},
$$ (eq-disc-delta6)

so that $Q_h u = \frac{\epsilon}{64}\sum_d \frac{1}{h_d}\,\delta_d^6 u$, which is Eq. {eq}`eq-kreiss-oliger`. In the
code these are the three sums `sx`, `sy`, `sz`, multiplied by `idx[d]` and by `c0 = eps/64`.

### 5.2 Fourier analysis: it only damps, and mostly the sawtooth

The symbol of $D_+D_-$ is $-\frac{4}{h^2}\sin^2\frac\theta2$, so

$$
    Q_h e^{ikx} = -\epsilon\sum_d \frac{1}{h_d}\sin^6\frac{\theta_d}{2}\; e^{ikx} .
$$ (eq-disc-ko-symbol)

Three properties follow:

- **Sign.** The symbol is real and $\le 0$. On its own, $\dot{\hat u} = Q_h\hat u$ decays and never grows.
- **Scale selectivity.** At $\theta = \pi$ the decay rate is $\epsilon/h$ per direction, which is the largest. For
  resolved modes, $\sin^6(\theta/2) \approx \theta^6/64$ and the term is $-\frac{\epsilon}{64}h^5k^6$.
- **Accuracy.** In real space $Q_h u = \frac{\epsilon}{64}\sum_d h_d^5\,\partial_d^6 u + \mathcal O(h^7)$. This is
  $\mathcal O(h^5)$, smaller than the $\mathcal O(h^4)$ truncation error of §3, so adding it does not lower the
  formal order of the scheme. (That holds only if the solution is smooth. Near a poorly resolved feature,
  dissipation is what dominates the error.)

Dissipation therefore removes exactly the modes that §3.4 showed the centred first derivative cannot see.

### 5.3 Where it acts

- It acts only on interior points $[n_g, n - n_g)$, and the 7-point stencil needs $n_g \ge 3$.
- It acts only on evolved variables: the `evolved_mask` that ADMEvolve registers with MoL excludes the shift and,
  under `static` slicing, the lapse.
- It is a MoL **RHS hook**. The Dissipation thorn registers `add` with {meth}`pynr.thorns.mol.MoL.add_rhs_hook`, so
  ADMEvolve does not know about it. The order of one RHS evaluation is: ADM RHS, then dissipation, then
  `rhs_final`. `rhs_final` (excision and outer boundary) comes last, so the excised region stays frozen and ghost
  RHS values are set by the boundary condition, not by dissipation.

## 6. Boundary conditions and the time step

### 6.1 Ghost points ({mod}`pynr.kernels.boundary`)

The ghosts must hold valid data at **every Runge-Kutta stage**, because interior stencils read them. pyNR fills
them in one of two ways:

- **RHS conditions**, set in `rhs_final` and integrated by the same RK method as the interior. With `static`,
  $du/dt = 0$ on ghost points. With `radiative`, the ghosts obey the Sommerfeld condition Eq. {eq}`eq-sommerfeld`,
  $du/dt = -[x^m\tilde D_m u + (u - u_\infty)]/r$. Here $\tilde D_m$ is a second-order derivative: centred
  $(f_{i+1}-f_{i-1})/2h$ where both neighbours exist, and one-sided inward $\pm(-3f_i + 4f_{i\pm1} - f_{i\pm2})/2h$
  in the ghost layer of that direction.
- **State conditions**, applied by `post(U)` after each stage. `periodic` copies the periodic image and `flat`
  copies the nearest interior value.

The boundary operators are second order (static is exact only for stationary data), so in a convergence test the
outer boundary can limit the measured order unless it is far away or periodic.

### 6.2 Runge-Kutta 4 and the CFL condition

$\Delta t =$ `Time::dtfac` $\times\min_d h_d$, with default `dtfac` $=\lambda = 0.25$. With
$\mathbf F(\mathbf u) = F_h + Q_h + \text{boundary}$ (all of Eq. {eq}`eq-disc-semidiscrete`), the classic RK4
step implemented in {func}`pynr.utils.integrators.rk4` is

$$
\begin{aligned}
    \mathbf k_1 &= \mathbf F(t, \mathbf u^n), &
    \mathbf u^{(1)} &= \mathbf u^n + \tfrac{\Delta t}{2}\mathbf k_1,\\
    \mathbf k_2 &= \mathbf F(t + \tfrac{\Delta t}{2}, \mathbf u^{(1)}), &
    \mathbf u^{(2)} &= \mathbf u^n + \tfrac{\Delta t}{2}\mathbf k_2,\\
    \mathbf k_3 &= \mathbf F(t + \tfrac{\Delta t}{2}, \mathbf u^{(2)}), &
    \mathbf u^{(3)} &= \mathbf u^n + \Delta t\,\mathbf k_3,\\
    \mathbf k_4 &= \mathbf F(t + \Delta t, \mathbf u^{(3)}), &
    \mathbf u^{n+1} &= \mathbf u^n + \tfrac{\Delta t}{6}(\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4),
\end{aligned}
$$ (eq-disc-rk4)

with `post` (periodic/flat ghosts) applied to each $\mathbf u^{(s)}$ and to $\mathbf u^{n+1}$. Its truncation
error is $\mathcal O(\Delta t^4) = \mathcal O(\lambda^4 h^4)$, the same order as the spatial error.

**Linear stability.** Linearise around flat space with $\alpha = 1$ and $\beta = 0$. One metric component then
obeys $\partial_t h = -2k$ and $\partial_t k = -\tfrac12\nabla^2 h$: a wave equation with speed 1. For a Fourier
mode, the semi-discrete eigenvalues times $\Delta t$ are

$$
    z(\boldsymbol\theta) = -\epsilon\lambda\sum_d\sin^6\frac{\theta_d}{2}
        \;\pm\; i\lambda\sqrt{\sum_d \frac{30 - 32\cos\theta_d + 2\cos2\theta_d}{12}} ,
$$ (eq-disc-eigen)

using Eqs. {eq}`eq-disc-symbols` and {eq}`eq-disc-ko-symbol` with equal spacings. The step is stable when every $z$
lies in the RK4 stability region $|1 + z + z^2/2 + z^3/6 + z^4/24| \le 1$. That region reaches $\pm 2.83i$ on the
imaginary axis and $-2.79$ on the real axis. Scanning all $\boldsymbol\theta \in [0,\pi]^3$ gives:

| | without dissipation | largest stable $\epsilon$ |
|---|---|---|
| $\lambda = 0.25$ (default) | stable ($\max\lvert z\rvert = 1$, amplification $0.994$) | $\epsilon \le 3.5$ |
| $\lambda = 0.5$ | stable | $\epsilon \le 1.15$ |
| limit | $4\lambda \le 2.83 \Rightarrow \lambda \le 0.71$ | — |

So `dtfac = 0.25` leaves a large margin, and the default `epsdis = 0.1` lies well inside it. These bounds are
necessary conditions for the linearised constant-coefficient problem. A black-hole spacetime has larger
characteristic speeds near the horizon and nonlinear terms, which is why production runs keep $\lambda \le 0.5$ and
$\epsilon \lesssim 0.5$.

## 7. One time step in pseudocode

The pseudocode follows the code: one MoL step at the top, then the RHS assembly, the pointwise ADM kernel and the
dissipation kernel. `U` is the state `U[16, nx, ny, nz]` and `dU` is a scratch array of the same shape.

```text
# ---------------------------------------------------------------- MoL::step  (EVOL bin)
procedure STEP(U, t, dt):                        # pynr.utils.integrators.rk4
    U0  ← copy(U)
    k   ← RHS(t,        U0);  acc ← k;           U ← U0 + dt/2 · k;  POST(U)
    k   ← RHS(t + dt/2, U );  acc ← acc + 2k;    U ← U0 + dt/2 · k;  POST(U)
    k   ← RHS(t + dt/2, U );  acc ← acc + 2k;    U ← U0 + dt   · k;  POST(U)
    k   ← RHS(t + dt,   U );  acc ← acc + k;     U ← U0 + dt/6 · acc; POST(U)

# ---------------------------------------------------------------- one RHS evaluation
function RHS(t, U) → dU:                         # MoL.step's closure
    ADM_RHS(U, dU)                               # ADMEvolve.rhs → kernels.adm.adm_rhs   (interior only)
    for each hook in MoL hooks:                  # Dissipation.add
        KO_DISSIPATION(U, dU, ε, evolved_mask)   # kernels.dissipation                   (interior only)
    RHS_FINAL(U, dU)                             # ADMEvolve.rhs_final
    return dU

procedure RHS_FINAL(U, dU):
    dU[:, excised points] ← 0                    # frozen region, overrides dissipation
    if grid not fully periodic and bound ≠ "flat":
        for each ghost point p, each variable c:
            if bound = "static":    dU[c,p] ← 0
            if bound = "radiative": dU[c,p] ← −( Σ_m x_m D̃_m U[c] + U[c,p] − u∞[c] ) / r_p   # 2nd order

procedure POST(U):                               # ADMEvolve.post
    if bound = "flat": copy nearest interior value into ghosts
    if any periodic direction: copy periodic images into ghosts

# ---------------------------------------------------------------- kernels.adm.adm_rhs
procedure ADM_RHS(U, dU):
    parallel for i in [ng, nx−ng):               # threads over i (prange)
      for j in [ng, ny−ng):
        for k in [ng, nz−ng):                    # k contiguous → streaming loads
          if r² < r_ex²: dU[:, i,j,k] ← 0; continue
          # --- geometry (_geometry) -------------------------------------------
          for each pair (a≤b): g_ab ← U[γ_ab]; dg[m,a,b] ← D_m U[γ_ab]; ddg[m,n,a,b] ← D_mn U[γ_ab]
          det ← det(g);   g^ab ← cofactor(g)^ab / det
          Γ_lab  ← ½ (dg[a,l,b] + dg[b,l,a] − dg[l,a,b]);    Γ^k_ab ← g^kl Γ_lab
          dgu[m,k,l] ← − g^kp g^lq dg[m,p,q]                 # algebraic, Eq. (disc-dginv)
          R_ab   ← Eq. (disc-ricci) from g^kl, dgu, dg, ddg, Γ
          # --- gauge and K -----------------------------------------------------
          α ← U[alp];  dα[m] ← D_m α;  ddα[m,n] ← D_mn α
          β^a ← U[β^a];  dβ[m,a] ← D_m β^a
          K_ab ← U[K_ab];  trK ← g^ab K_ab;  K^a_b ← g^am K_mb
          # --- evolution equations, Eq. (disc-adm-rhs) -------------------------
          for each pair (a≤b):
              dK[m] ← D_m U[K_ab]
              Lie_g ← Σ_m ( β^m dg[m,a,b] + g_mb dβ[a,m] + g_am dβ[b,m] )
              Lie_K ← Σ_m ( β^m dK[m]     + K_mb dβ[a,m] + K_am dβ[b,m] )
              DDα   ← ddα[a,b] − Σ_m Γ^m_ab dα[m]
              KK    ← Σ_m K_am K^m_b
              dU[γ_ab] ← −2α K_ab + Lie_g
              dU[K_ab] ← −DDα + α (R_ab + trK K_ab − 2 KK) + Lie_K
          adv ← (advect_lapse ? Σ_m β^m dα[m] : 0)
          dU[alp] ← { 0 | −α² trK + adv | −2α trK + adv }     # static | harmonic | 1+log
          dU[β^a] ← 0

# ---------------------------------------------------------------- kernels.dissipation
procedure KO_DISSIPATION(U, dU, ε, mask):
    c0 ← ε / 64
    parallel for i in [ng, nx−ng):
      for each variable c with mask[c]:
        for j in [ng, ny−ng):
          for k in [ng, nz−ng):
            δ6x ← U[c,i−3]−6U[c,i−2]+15U[c,i−1]−20U[c,i]+15U[c,i+1]−6U[c,i+2]+U[c,i+3]   # along x
            δ6y ← (same along j);   δ6z ← (same along k)
            dU[c,i,j,k] ← dU[c,i,j,k] + c0 · ( δ6x/h_x + δ6y/h_y + δ6z/h_z )
```

Per step, RK4 calls `ADM_RHS` four times. That kernel dominates the cost; see [Performance](../performance.md) for
the measured time per point on each machine.

## 8. Summary of errors

| part | discretisation | error |
|---|---|---|
| first/second/mixed derivatives | centred 5-point / compact 5-point / $4\times4$ product | $\mathcal O(h^4)$ |
| Ricci, Christoffels, $\partial g^{-1}$ | algebraic in the discrete derivatives | $\mathcal O(h^4)$ |
| shift advection | centred | $\mathcal O(h^4)$ |
| Kreiss-Oliger | $\frac{\epsilon}{64}h^5(D_+D_-)^3$ | $\mathcal O(\epsilon h^5)$ |
| outer boundary | static / Sommerfeld with 2nd-order one-sided differences | $\mathcal O(h^2)$ at the boundary |
| excision | frozen data | not convergent; it must lie inside the horizon |
| time | RK4, $\Delta t = \lambda h$ | $\mathcal O(\Delta t^4)$ |

On a periodic domain with smooth data the whole scheme is fourth order, which is what the gauge-wave test measures
(Eq. {eq}`eq-convergence-order`, {numref}`fig-gauge-wave`).

## Exercises

1. Derive Eq. {eq}`eq-disc-d2-error` and show that $D_xD_x$ has the symbol $-(8\sin\theta - \sin2\theta)^2/36h^2$.
   Which mode does it miss?
2. Show that $(-1)^{r+1}(D_+D_-)^r$ has a non-positive symbol for every $r$. Why would $r = 2$ (fourth differences)
   be a poor choice for a 4th-order scheme?
3. Run [problem 1](../problems/gauge_wave.md) with `Dissipation::epsdis = 4.0` and then `3.0`. Compare with the
   table in §6.2.
4. Replace Eq. {eq}`eq-disc-dginv` by differencing a stored $\gamma^{ij}$ (the inverse metric computed on the whole
   grid first). Does the Hamiltonian constraint of the gauge wave stay at round-off? (See {ref}`sec-gauge-wave-31`.)
