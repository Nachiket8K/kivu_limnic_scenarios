# Lake Kivu Limnic Eruption Scenario-Discovery Model — Mathematical Specification (paper-grounded)

This document extracts the *explicit mathematical formulas* and the *calculation workflow* needed to reproduce a **terrain-following dense CO₂ cloud** dispersion and its **human impact** (fatalities / hazard) as described in the core shallow-layer (TWODEE-family) literature, plus a standard regime-selection criterion used in recent probabilistic hazard workflows.

It is written so a separate implementation agent can translate it into code.

---

## 1) Core intent (as in the papers)

Your stated intent — “scenario discovery to see how the CO₂ flow will harm people and where” — maps directly to:  
1) **Dense-gas dispersion over topography** using a shallow-layer (depth-averaged) model (TWODEE family). citeturn28view1turn25view1  
2) **Vertical concentration and dosage at human-relevant heights** derived from the depth-averaged state variables. citeturn28view2  
3) **Impact / lethality mapping** using an empirical probability-of-death criterion depending on concentration and exposure duration. citeturn28view0  
4) Optional (but recommended): **source-regime selection** (passive vs dense) using a Richardson-number criterion. citeturn30view2  

---

## 2) Notation and state variables (paper-consistent)

All fields are functions of horizontal position \((x,y)\) and time \(t\).

**Terrain**
- \(e(x,y)\): terrain elevation. citeturn28view1turn25view1

**Ambient air**
- \(\rho_a\): air density (environment).
- \(\mathbf{u}_a=(v_a,w_a)\): ambient wind velocity vector (near-surface wind field). citeturn28view1

**Dense cloud (depth-averaged)**
- \(h(x,y,t)\): cloud depth (shallow-layer “height”). citeturn28view1turn25view0  
- \(\rho(x,y,t)\): depth-averaged cloud density. citeturn28view1  
- \(\mathbf{u}(x,y,t)=(u_x,u_y)\): depth-averaged cloud horizontal velocity. citeturn28view1  

**Gas properties**
- \(\rho_g\): density of pure gas (CO₂) at relevant conditions (used in mixing and concentration conversion). citeturn28view1turn28view2
- \(c_b\): background CO₂ concentration in air (ppm). citeturn28view2turn29view0

**Semi-empirical / closure parameters (appear explicitly)**
- \(S_1 \approx 0.5\): semi-empirical parameter in pressure/vertical-profile relations. citeturn28view1turn28view2  
- \(C_D\): drag coefficient for surface shear stress. citeturn28view1turn25view1  
- \(\kappa\): semi-empirical coefficient in leading-edge / interaction term. citeturn28view1turn28view2  
- \(\mathbf{F}\): turbulent shear stress force per unit area. citeturn28view1  
- \(u_e\): entrainment velocity modulus (ambient fluid entrainment into cloud). citeturn28view1turn25view1  
- \(u_s\): “gas inflow velocity modulus” used to represent the source term. citeturn28view1  

---

## 3) Shallow-layer governing equations (TWODEE-family)

### 3.1 Volume (cloud-depth) balance
\[
\frac{\partial h}{\partial t} + \nabla \cdot (h\,\mathbf{u}) = u_e + u_s.
\tag{1}
\]
This is the depth-integrated “volume” (cloud thickness) conservation. citeturn28view1

### 3.2 Mass (excess-density) balance
\[
\frac{\partial \big(h(\rho-\rho_a)\big)}{\partial t} + \nabla \cdot \big(h(\rho-\rho_a)\,\mathbf{u}\big)
= \rho_a u_e + \rho_g u_s.
\tag{2}
\]
This allows \(\rho\) (hence buoyancy) to vary in space/time through entrainment and source injection. citeturn28view1

### 3.3 Momentum balance (vector form)
\[
\frac{\partial(h\rho\mathbf{u})}{\partial t}
+ \nabla \cdot (h\rho\,\mathbf{u}\otimes\mathbf{u})
+ \frac{1}{2}S_1\,\nabla\big(g(\rho-\rho_a)h^2\big)
+ S_1 g(\rho-\rho_a)h\,\nabla e
+ \frac{1}{2}\rho C_D\lvert\mathbf{u}\rvert\mathbf{u}
+ \mathbf{F}
+ \kappa\rho_a\Big(\frac{\partial}{\partial t}
+ v_a\frac{\partial}{\partial x}
+ w_a\frac{\partial}{\partial y}\Big)\big(h(\mathbf{u}-\mathbf{u}_a)\big)
= \rho_a u_e\mathbf{u}_a.
\tag{3}
\]
Interpretation of terms (as given in the paper): local time derivative, convection, hydrostatic pressure-gradient, terrain-slope forcing, surface drag, turbulent shear stress force, and leading-edge / dense–ambient interaction terms. citeturn28view1turn28view2

> Practical note: TWODEE implements (1)–(3) numerically (in the original code base, via flux-corrected transport schemes), but the **equations above are the mathematical core** you’ll be matching. citeturn25view4  

---

## 4) Vertical profile + concentration conversion (from depth-averaged state)

The shallow-layer model predicts \(h,\rho,\mathbf{u}\). To compute **human-relevant CO₂ concentration at height** \(z\) (e.g., 0.5–2 m) you use the paper’s vertical profile.

### 4.1 Exponential density profile
Assuming an empirical exponential decay from the ground upward:
\[
\rho(z) = \rho_a + \frac{2}{S_1}(\rho-\rho_a)\exp\Big(-\frac{2}{S_1}\frac{z}{h}\Big),\quad 0\le z\le h.
\tag{4}
\]
citeturn28view2

### 4.2 Convert density profile to CO₂ concentration (ppm)
\[
c(z)= c_b + (10^6-c_b)\,\frac{\rho(z)-\rho_a}{\rho_g-\rho_a}.
\tag{5}
\]
This yields CO₂ concentration at height \(z\) in ppm. citeturn28view2turn25view0

### 4.3 Dosage over an exposure window
Given a toxicity exponent \(n\) (paper term), define dosage at height \(z\):
\[
D_o(t,z) = \int_0^{t} [c(z,t')]^{n}\,dt'.
\tag{6}
\]
citeturn28view2

> For your application (rapid lethal hazard), you can also work directly with \(c(z,t)\) and exposure duration \(d\) (minutes) for the fatality model (next section).

---

## 5) Impact criterion: probability of death vs concentration and exposure duration

Folch et al. (2017) introduce an “impact criterion” that outputs **percentage of fatalities** as a function of CO₂ concentration and exposure time (calibrated to HSE toxicology tables). citeturn28view0

### 5.1 Probability of death as a cumulative normal (error function)
Let:
- \(c\) = CO₂ concentration (in **% vol**, not ppm)
- \(d\) = exposure duration (minutes)
- \(P(c,d)\) = probability of death (0–1)

Then:
\[
P(c,d) = \frac{1}{2}\left[1+\operatorname{erf}\left(\frac{c-\mu}{\sqrt{2}\,\sigma}\right)\right].
\tag{7}
\]
with:
\[
\mu = a_0 + \frac{b_0}{1 + d^{c_0}},
\tag{8}
\]
\[
\sigma = a_1 + \frac{b_1}{1 + d^{c_1}}.
\tag{9}
\]
citeturn28view0

### 5.2 Parameter values (as calibrated in the paper)
After calibration (assuming SLOT at 3%), the paper reports:
- \(a_0=5.056\), \(b_0=17.885\), \(c_0=0.357\)  
- \(a_1=0.662\), \(b_1=2.421\), \(c_1=0.354\)  
citeturn28view0

### 5.3 Units and conversion (must be consistent)
Your dispersion model provides \(c(z,t)\) in **ppm** via Eq. (5). Convert to % vol for Eq. (7) as:
\[
c_{\%}(z,t)=\frac{c_{ppm}(z,t)}{10^4}.
\]
(Reason: 1% = 10,000 ppm.)

### 5.4 Scenario-discovery output you can compute from this
For each grid cell and height \(z\) (e.g., 1–2 m), compute:
- running exposure duration \(d\) (minutes) above chosen concentration thresholds (or simply total time in window)
- instantaneous \(P(c,d)\) per time step, and/or a “max over time” \(\max_t P\)
- expected fatalities = \(P\times\) local population (if you have population raster)

---

## 6) Regime selection: when a dense-gas model is needed (Ri criterion)

In probabilistic hazard workflows, the **flow Richardson number at the source** is used to decide whether dispersion is passive (wind-driven) or dense (gravity-current-like). citeturn30view2

### 6.1 Richardson number at source
\[
Ri = \frac{1}{v^2}\left(\frac{g'\,q}{r}\right)^{2/3},
\tag{10}
\]
where:
- \(v\): wind speed at source
- \(q\): **source volumetric flow rate**
- \(r\): effective source radius
- \(g'\): reduced gravity

citeturn30view2

### 6.2 Reduced gravity
\[
g' = g\,\frac{\rho_g-\rho_e}{\rho_e},
\tag{11}
\]
where \(\rho_e\) is environmental density and \(\rho_g\) is starting gas density. citeturn30view2

### 6.3 Regimes (as stated)
- If \(Ri<0.25\): passive dispersion dominated by wind advection/diffusion. citeturn30view2  
- If \(Ri>1\): dense-gas transport dominated by density contrast; gravity-current behaviour over topography. citeturn30view2  
- Intermediate (0.25–1): some workflows choose the dense model for caution (yields higher concentrations). citeturn30view2  

---

## 7) Recommended mathematical/spatial changes to your current POC (explicit mapping)

From your notebooks, your POC currently uses:
- slope-shaped “terrain velocity” + upwind advection,
- Laplacian diffusion,
- an explicit “loss” term (wind-driven decay),
- a circular initial disk source and metrics (centroid/front).  
(Notebook functions: `shape_slope_velocity`, `advect_upwind`, `diffuse_laplacian`, `apply_loss`, `H_to_ppm`, `disk_mask`, etc.)

Below are **paper-consistent upgrades** that keep your scenario-discovery focus but align the mathematics with the TWODEE-family framework.

### Change A — Replace heuristic advection with the shallow-layer PDE state: \(h,\rho,\mathbf{u}\)
**What to implement**
1. Promote your state from a single “CO₂ layer thickness / concentration” to the TWODEE set:
   - \(h\) (cloud depth)
   - \(\rho\) (depth-averaged density)
   - \(\mathbf{u}=(u_x,u_y)\) (depth-averaged velocity)

2. Evolve them with Eqs. (1)–(3). citeturn28view1  

**Why this matches your physics complaint**
- Terrain-following spread and pooling comes from the **pressure/slope forcing terms**
  \(\nabla(g(\rho-\rho_a)h^2)\) and \(h\nabla e\) in Eq. (3), not from wind alone. citeturn28view1  
- Wind still matters through \(\mathbf{u}_a\) and the interaction term in Eq. (3), but does not “own” the dynamics when density contrast is high. citeturn30view2  

**Discretisation guidance (process, not a specific scheme)**
At each time step \(t\to t+\Delta t\):
1. Compute closure terms (see Change C) at \(t\).
2. Update \(h\) via Eq. (1) (mass flux divergence + entrainment + source).
3. Update \(h(\rho-\rho_a)\) via Eq. (2), then recover \(\rho\).
4. Update momentum \(h\rho\mathbf{u}\) via Eq. (3), then recover \(\mathbf{u}\).
5. Enforce physically valid bounds: \(h\ge 0\); \(\rho\ge \rho_a\) in dense regime.

### Change B — Keep your “circular eruption footprint”, but express it as a source term \(u_s\) and \(\rho_g u_s\)
You asked for “erupts as a circle” with a radius and the volume stacked there.

Represent this as a **source region** \(\Omega_s\) (disk on the grid):
- define a **time-dependent** source injection rate (mass or volume)
- convert it to the terms required by Eqs. (1)–(2)

**Option B1 (volumetric source specified):**  
Choose a volumetric injection rate per unit area \(q_s(t)\) [m/s] over \(\Omega_s\).
Then set:
- \(u_s(x,y,t)=q_s(t)\) inside \(\Omega_s\), else 0.
- In Eq. (2), the source contributes \(\rho_g u_s\). citeturn28view1  

**Option B2 (total volume specified):**  
If you start from a total released volume \(V_{tot}\) over a duration \(T\),
- \(\dot V = V_{tot}/T\)
- if \(A_s\) is the disk area, \(q_s = \dot V/A_s\) so \(u_s=q_s\) on \(\Omega_s\).

This directly satisfies your “stack the whole volume in a circle” requirement while staying consistent with the governing equations. citeturn28view1  

### Change C — Explicitly compute near-ground hazard from \(h,\rho\) via Eqs. (4)–(9)
Your notebook already has `H_to_ppm`-style logic; replace/anchor it to the paper sequence:

Per grid cell and time step:
1. Use \(h\) and depth-averaged \(\rho\) to compute \(\rho(z)\) with Eq. (4). citeturn28view2  
2. Convert \(\rho(z)\to c(z)\) in ppm via Eq. (5). citeturn28view2  
3. Convert ppm to % (divide by 10,000).
4. Update exposure duration \(d\) (minutes) at that cell/height (e.g., cumulative time above some minimum CO₂ level).
5. Compute \(P(c,d)\) via Eqs. (7)–(9) and store either:
   - \(P\) time series, or
   - \(\max_t P\) and/or time-above-threshold maps.

This gives you paper-grounded “harm where” outputs without requiring full physiological modelling. citeturn28view0  

### Change D — Add a regime switch (optional, but aligns with scenario discovery)
For each scenario, compute source Ri using Eqs. (10)–(11). citeturn30view2  
- If \(Ri\gg 1\): use the dense shallow-layer equations (this document).
- If \(Ri\ll 0.25\): your simpler wind-advection + diffusion model is a reasonable approximation for a *passive* regime.
- If in-between: choose dense for conservatism (as done in hazard workflows). citeturn30view2  

### Change E — Time step requirement for lethality (your 1 h Δt is too coarse for impact)
The fatality criterion depends on exposure durations of a few minutes to tens of minutes. citeturn28view0turn29view2  
If you keep \(\Delta t=1\) hour for the physics (for speed), consider *sub-stepping* the impact calculation:
- interpolate \(c(t)\) within the hour (e.g., linear)
- compute \(d\) in minutes and update \(P\) at smaller increments (e.g., 1–5 min)

(Your user requirement for 1 h timestep is fine for POC dynamics, but not for exposure-based lethality outputs.)

---

## 8) Minimal end-to-end calculation workflow (for your other AI to implement)

**Inputs**
- DEM \(e(x,y)\)
- wind field \(\mathbf{u}_a(x,y,t)\) (at least near-ground)
- gas properties \(\rho_g\), \(\rho_a\), \(c_b\)
- closure parameters \(S_1, C_D, \kappa\) + closure laws for \(u_e\) and \(\mathbf{F}\) (from the TWODEE closure references)
- source geometry \(\Omega_s\) + emission schedule \(q_s(t)\) or mass flux schedule
- population raster (optional, for expected fatalities)

**Simulation loop**
For each time step:
1. Evaluate \(u_s(x,y,t)\) in \(\Omega_s\).
2. Compute \(u_e\), \(C_D\), \(\mathbf{F}\) (closure).
3. Solve Eqs. (1)–(3) to advance \(h,\rho,\mathbf{u}\).
4. For each chosen height \(z\) (e.g., 1 m and 2 m):
   - compute \(c(z,t)\) via (4)–(5)
   - update exposure duration \(d\) and compute \(P\) via (7)–(9)
5. Store outputs:
   - \(\max_t c(z,t)\), \(\max_t P\), time-above-threshold,
   - optionally expected fatalities = \(P\times\) population.

---

## 9) What this document intentionally does *not* specify
- The detailed *closure equations* for entrainment velocity \(u_e\), drag \(C_D\), and turbulent shear force \(\mathbf{F}\). Folch et al. (2017) explicitly points to the closure literature for these. citeturn28view1  
- Any PRIM / scenario discovery algorithm steps (you asked to keep PRIM aside for now).

If you want, I can add a second appendix that lists the most common closure options used with Hankin & Britter-style shallow-layer dense-gas models, but I’ll need the specific closure source(s) you want to treat as authoritative.

