# Heat and flow in 2D

A small C++ solver for transient incompressible flow and heat transfer,
with Python for visualization.

Built for teaching: explicit loops, simple rectangular domains,
and tests against analytical solutions.

## Build and run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/simulate heated_wall results
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 visualize.py results --field T
python3 visualize.py results --field T --save results/temperature.png
```

[Cases](src/cases.hpp) · [Parameters](src/solver.hpp) · [Heat tests](tests/heat_equation_test.cpp) · [Fluid tests](tests/fluid_test.cpp)

## Method

We solve incompressible flow and temperature transport in a 2D rectangle,
with constant material properties. In 2D, these are four equations for the
four unknowns $`u,v,P,T`$: continuity, two momentum components and heat.

**Continuity**

$$
\nabla\cdot\mathbf{u}=0
$$

**Momentum**

$$
\frac{\partial \mathbf{u}}{\partial t}
+(\mathbf{u}\cdot\nabla)\mathbf{u}
=-\frac{1}{\rho}\nabla P+\nu\nabla^2\mathbf{u}
+g\beta(T-T_{\mathrm{ref}})\mathbf{e}_y
$$

**Heat**

$$
\frac{\partial T}{\partial t}+\mathbf{u}\cdot\nabla T
=\alpha\nabla^2T
$$

Here $`\mathbf{u}=(u,v)`$ is velocity, $`t`$ time, $`P`$ pressure relative to the reference hydrostatic state, $`T`$ temperature, $`\rho`$ reference density, $`\nu`$ kinematic viscosity, and $`\alpha`$ thermal diffusivity. $`g`$ is the positive gravitational acceleration magnitude, $`\beta`$ the thermal expansion coefficient, and $`T_{\mathrm{ref}}`$ the reference temperature. $`\mathbf{e}_x=(1,0)`$ and $`\mathbf{e}_y=(0,1)`$ are unit vectors, with $`y`$ pointing upward; $`g\beta(T-T_{\mathrm{ref}})\mathbf{e}_y`$ is the Boussinesq buoyancy acceleration. $`\nabla`$ and $`\nabla^2`$ denote the gradient and Laplacian operators.

We use a uniform staggered grid: $`P,T`$ at cell centres, $`u`$ on vertical
faces and $`v`$ on horizontal faces. A 3 × 2 example (ghost cells omitted):

```text
            y, j
       e_y  ^
            |
 2dy        +--------v(1,2)---------+--------v(2,2)---------+--------v(3,2)---------+
            |                       |                       |                       |
 j=2     u(0,2)     P,T(1,2)     u(1,2)     P,T(2,2)     u(2,2)     P,T(3,2)     u(3,2)
            |                       |                       |                       |
 dy         +--------v(1,1)---------+--------v(2,1)---------+--------v(3,1)---------+
            |                       |                       |                       |
 j=1     u(0,1)     P,T(1,1)     u(1,1)     P,T(2,1)     u(2,1)     P,T(3,1)     u(3,1)
            |                       |                       |                       |
 (0,0)      O--------v(1,0)---------+--------v(2,0)---------+--------v(3,0)---------+ --> x, i
            0                      dx                     2dx                     3dx     e_x
```

The physical origin is the bottom-left corner; $`\mathbf{e}_x=(1,0)`$ points
right and $`\mathbf{e}_y=(0,1)`$ points upward. Indices $`i,j`$ increase in
those directions. The first cell centre, $`P,T(1,1)`$, is at
$`(\Delta x/2,\Delta y/2)`$; $`u(0,1)`$ and $`v(1,0)`$ lie on its left and bottom faces.

The code uses `pressure(i,j)` and `temperature(i,j)` for $`P,T`$.
`field(i,j)` accesses `field.values[j * field.nx + i]`: $`i`$ varies fastest,
and `field.nx` is the stored row width, including any ghost columns.
Row $`j=1`$ is near the bottom of the domain, even if a printed matrix lists it first.

When a velocity is needed where it is not stored, we interpolate by averaging.
For example, at the cell centre:

$$
u_c=\frac{u_{i-1,j}+u_{i,j}}{2},
\qquad
v_c=\frac{v_{i,j-1}+v_{i,j}}{2}.
$$

For momentum advection, we also need $`v`$ at a $`u`$ face and $`u`$ at a
$`v`$ face. We average the four surrounding values:

$$
v_{\text{at }u(i,j)}=
\frac{v_{i,j-1}+v_{i+1,j-1}+v_{i,j}+v_{i+1,j}}{4},
$$

$$
u_{\text{at }v(i,j)}=
\frac{u_{i-1,j}+u_{i,j}+u_{i-1,j+1}+u_{i,j+1}}{4}.
$$

### Boundary conditions and ghost cells

Boundary conditions prescribe velocity or heat transfer at the walls.
Ghost values outside the domain are calculated from the interior values to
enforce these conditions. At the bottom wall, $`u`$ is tangential and $`v`$ is normal:

```text
                    y
               e_y  ^
                    |
 dy/2, j=1       u(0,1)     T(1,1)     u(1,1)     T(2,1)     u(2,1)   fluid
                    |                       |                       |
 wall, y=0          O--------v(1,0)---------+--------v(2,0)---------+ --> x, e_x
                    :                       :                       :
-dy/2, j=0       u(0,0)     T(1,0)     u(1,0)     T(2,0)     u(2,0)   ghosts
```

The origin $`O=(0,0)`$ lies on the wall. The normal velocity is stored directly
there: impermeability gives $`v_{i,0}=0`$. The ghost values $`u_{i,0},T_{i,0}`$
are half a cell below the wall; their interior neighbours $`u_{i,1},T_{i,1}`$
are half a cell above it.

For a wall moving tangentially at speed $`U_{\mathrm{wall}}`$, no slip requires
the average at the wall to equal that speed:

$$
\frac{u_{i,0}+u_{i,1}}{2}=U_{\mathrm{wall}}
\quad\Longrightarrow\quad
u_{i,0}=2U_{\mathrm{wall}}-u_{i,1}.
$$

A fixed wall has $`U_{\mathrm{wall}}=0`$, so $`u_{i,0}=-u_{i,1}`$.

For an imposed wall temperature $`T_{\mathrm{wall}}`$, the same averaging gives:

$$
\frac{T_{i,0}+T_{i,1}}{2}=T_{\mathrm{wall}}
\quad\Longrightarrow\quad
T_{i,0}=2T_{\mathrm{wall}}-T_{i,1}.
$$

For an insulated wall, the normal temperature gradient is zero:

$$
\left.\frac{\partial T}{\partial y}\right|_{y=0}
\approx\frac{T_{i,1}-T_{i,0}}{\Delta y}=0
\quad\Longrightarrow\quad
T_{i,0}=T_{i,1}.
$$

The other walls use the same construction, with the corresponding normal
and tangential velocity components.

### Discretization and time stepping

Advection uses first-order upwind
differences; diffusion uses second-order central differences.
Time integration is first-order explicit, with a sequential thermal/fluid update.

Each call to [Solver::step](src/solver.cpp) follows this order, skipping disabled equations:

1. Check the time step using the current velocity and active diffusivities.
2. Advance temperature using the current velocity.
3. Predict velocity from advection, viscosity and forces, using the updated
   temperature for buoyancy.
4. Solve the pressure Poisson equation and subtract the pressure gradient
   to enforce incompressibility. Reapply the wall conditions at the new time.

The pressure equation uses Gauss–Seidel, optionally with relaxation.
Prescribed pressure boundaries set its reference; otherwise one cell fixes the arbitrary constant.
Channel cases use open ends with prescribed pressure.

The terminal shows Courant numbers, viscous/thermal Fourier numbers,
the combined explicit stability estimate and an estimated maximum time step.
The checks run every step; an excessive step or a nonfinite state stops the run.
These diagnostics use conservative velocity bounds including ghost cells.

## Transient heat conduction in water

Initially, the water is uniformly at 20°C. At t = 0 the left wall is set to
40°C; the right wall stays at 20°C. The top and bottom are insulated ($`\partial T/\partial y=0`$).
The 2 mm × 1 mm rectangle has no fluid motion or buoyancy: this isolates
thermal diffusion and gives a one-dimensional analytical solution.

![Water heated from one wall: transient numerical and analytical temperature](docs/media/heated_wall.gif)

12 s simulated in a 15 s video.
Colour shows temperature. Cyan points are the C++ solution, orange the analytical
transient, and dashed grey the steady linear profile.

```sh
./build/simulate heated_wall results_heated_wall
.venv/bin/python visualize.py results_heated_wall
.venv/bin/python visualize.py results_heated_wall --save results_heated_wall/heated_wall.gif
```

Water thermal diffusivity: 1.43e-7 m²/s. Mesh: 48 × 12, dt = 0.001 s.
[Transient test](tests/heated_wall_test.cpp): 0.1, 0.5, 1, 2, 5 and 12 s.
The maximum error must stay below 0.5% of the imposed 20 K temperature difference.

For the interior at t > 0, with initial temperature equal to the cold wall:

$$
T(x,t)=T_c+(T_h-T_c)\left(1-\frac{x}{L}\right)
-\frac{2(T_h-T_c)}{\pi}\sum_{n=1}^{\infty}
\frac{\sin(n\pi x/L)}{n}
\exp\left(-\frac{\alpha n^2\pi^2t}{L^2}\right)
$$

[Analytical method: heat equation with fixed end temperatures](https://www.math.purdue.edu/~pvankoug/math303/hw9solns.pdf)

## Poiseuille

The fluid starts at rest. A constant pressure gradient accelerates it, while viscosity
slows it down near the fixed walls. The velocity profile gradually approaches
the steady parabola.

![Transient Poiseuille flow: numerical and analytical velocity profiles](docs/media/poiseuille.gif)

30 s simulated, played at 2.5x speed.
White particles follow the numerical velocity; cyan points show the numerical
profile, orange the transient analytical solution, and dashed grey the steady limit.

```sh
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/simulate poiseuille results_poiseuille
.venv/bin/python visualize.py results_poiseuille
.venv/bin/python visualize.py results_poiseuille --save results_poiseuille/poiseuille.gif
```

Start from rest. Evolve for 30 s. Save every 0.1 s. Playback: 25 fps (2.5x physical speed).
The channel length is $`L=0.12\,\mathrm{m}`$, its height $`H=0.02\,\mathrm{m}`$,
and the density $`\rho=1000\,\mathrm{kg/m^3}`$.

We impose pressure directly at the left and right ends:

$$
P(0,y,t)=P_{\mathrm{in}}=1.2\,\mathrm{Pa},
\qquad
P(L,y,t)=P_{\mathrm{out}}=0\,\mathrm{Pa}.
$$

The pressure difference $`\Delta P=P_{\mathrm{in}}-P_{\mathrm{out}}`$ drives
the fluid from left to right. These are open boundaries: fluid crosses them,
and its velocity is calculated rather than prescribed. For this developed
channel flow, both ends use:

$$
\frac{\partial u}{\partial x}=0,\qquad v=0
\quad\text{at }x=0,L.
$$

The zero derivative extends the horizontal velocity profile through each end;
it does not set $`u`$ to zero. Fixed bottom and top walls impose no slip and
no penetration:

$$
\mathbf{u}(x,0,t)=\mathbf{u}(x,H,t)=\mathbf{0}.
$$

The pressure solve uses $`P_{\mathrm{in}}`$ and $`P_{\mathrm{out}}`$ as known
boundary values. The first and last cell centres are $`\Delta x/2`$ from these
boundaries. We use that same distance when correcting the velocities on the
inlet and outlet faces, so the adjacent cells also satisfy incompressibility.

For this case, the computed pressure is linear along the channel:

$$
P(x,y,t)=P_{\mathrm{in}}-\frac{\Delta P}{L}x,
\qquad
-\frac{1}{\rho}\frac{\partial P}{\partial x}=\frac{\Delta P}{\rho L}.
$$

The pressure gradient supplies the acceleration directly; no horizontal body
force is added. The output column `p` contains $`P`$ itself, in pascals relative
to the chosen outlet pressure.

[Poiseuille tests](tests/poiseuille_test.cpp): 0.1, 0.5, 2, 5, 8, 15, 30 s.
The final-profile test requires a maximum error below 0.2% of the steady
centerline velocity. With the default case:

| Physical time | Maximum difference from the steady parabola |
| --- | --- |
| 8 s | 14.30% |
| 30 s | 0.041% |

Differences are normalized by the theoretical steady centerline velocity.
These checks validate viscous diffusion, the channel boundaries, and the transient
profile. General two-dimensional nonlinear advection needs a separate benchmark.

$$
u(y,t)=\frac{\Delta P}{2\rho\nu L}y(H-y)
-\frac{4\Delta P H^2}{\rho\nu L\pi^3}
\sum_{n=1,3,5,\ldots}
\frac{\sin(n\pi y/H)}{n^3}
\exp\left(-\frac{\nu n^2\pi^2t}{H^2}\right)
$$

[Analytical reference: MIT](https://web.mit.edu/course/10/10.34/www/ch4_Poiseuille_startup.htm)

## Couette

The fluid starts at rest. The top wall suddenly moves at 0.05 m/s while the
bottom wall stays fixed ($`\mathbf{u}(x,0,t)=\mathbf{0}`$). Viscosity transfers momentum into the channel until
the velocity profile approaches a straight line.

![Transient Couette flow: numerical and analytical velocity profiles](docs/media/couette.gif)

30 s simulated, played at 2.5x speed.
White particles follow the numerical velocity; cyan points show the numerical
profile, orange the transient analytical solution, and dashed grey the steady limit.

```sh
./build/simulate couette results_couette
.venv/bin/python visualize.py results_couette
.venv/bin/python visualize.py results_couette --save results_couette/couette.gif
```

The left and right ends are open, with the same prescribed pressure:

$$
P(0,y,t)=P(L,y,t)=0.
$$

For this developed channel flow, both ends use:

$$
\frac{\partial u}{\partial x}=0,\qquad v=0
\quad\text{at }x=0,L.
$$

The horizontal velocity is calculated at both ends. The zero derivative means
the profile does not change along $`x`$ there; it still evolves with time.
Equal end pressures produce no driving pressure gradient: the moving top wall
sets the fluid in motion through viscosity.

The bottom wall is fixed and the top wall moves at $`U=0.05\,\mathrm{m/s}`$:

$$
\mathbf{u}(x,0,t)=(0,0),\qquad
\mathbf{u}(x,H,t)=(U,0),
$$

where $`L`$ is the channel length and $`H`$ its height.
Start from rest. Evolve for 30 s. Save every 0.1 s.

[Couette tests](tests/couette_test.cpp): 0.1, 0.5, 2, 5, 8, 15, 30 s.
Transient error must stay below 0.5% of wall speed; the final-profile test requires
an error below 0.2% of wall speed. With the default case:

| Physical time | Maximum difference from the steady line |
| --- | --- |
| 8 s | 8.85% |
| 30 s | 0.039% |

Differences are normalized by the moving wall speed. These checks validate
momentum diffusion from a moving wall and the transient approach to linear shear.

$$
u(y,t)=U\frac{y}{H}
+\frac{2U}{\pi}\sum_{n=1}^{\infty}
\frac{(-1)^n}{n}\sin\left(\frac{n\pi y}{H}\right)
\exp\left(-\frac{\nu n^2\pi^2t}{H^2}\right)
$$

[Analytical reference: transient Couette flow](https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node15.html)

## Stokes: oscillating wall in water

Water starts at rest. The bottom wall oscillates horizontally:
$`u_{\mathrm{wall}}(t)=U\sin(\omega t)`$, with $`U=0.01\,\mathrm{m/s}`$ and period 2 s.
The top wall stays fixed. The left and right ends are open, with equal pressure.
Temperature and buoyancy are disabled.

$$
P(0,y,t)=P(L,y,t)=0,\qquad
\frac{\partial u}{\partial x}=0,\qquad v=0
\quad\text{at }x=0,L.
$$

The horizontal velocity is calculated at both ends; the zero derivative is
spatial, so the velocity can still vary in time. Equal end pressures produce
no driving pressure gradient. The bottom wall drives the flow through viscosity,
and fluid can enter or leave either end as the wall reverses.
Here $`L`$ is the channel length.

![Oscillating wall in water: signed velocity, particles and analytical profiles](docs/media/stokes.gif)

Four cycles: 8 s of physical time shown in a 20 s animation.
The fluid is nearly black at rest. Warm tones show rightward velocity;
blue and cyan show leftward velocity. The moving particles reverse
direction; the cyan profile is numerical, orange is the analytical startup,
and dashed grey is the time-periodic solution.

```sh
./build/simulate stokes results_stokes
.venv/bin/python visualize.py results_stokes
.venv/bin/python visualize.py results_stokes --save results_stokes/stokes.gif
```

Water: density 998.2 kg/m³, kinematic viscosity 1.004e-6 m²/s.
Domain: 16 mm × 4 mm; mesh 24 × 64; dt = 0.0005 s.
This is the finite-channel version of Stokes' second problem; the analytical
comparison includes the fixed upper wall and the initial transient.

[Stokes tests](tests/stokes_test.cpp) check startup, reversal, incompressibility
and wall velocity at 0.1, 0.5, 1, 1.5, 2, 3.5, 5, 7.5 and 8 s.
Maximum startup error: below 0.5% of wall amplitude.
At 8 s the profile must also be within 0.3% of the time-periodic solution, which repeats every 2 s.

With $`k_n=n\pi/H`$ and $`\lambda_n=\nu k_n^2`$:

$$
u(y,t)=U\left(1-\frac{y}{H}\right)\sin(\omega t)
-\frac{2U\omega}{\pi}\sum_{n=1}^{\infty}
\frac{\sin(k_n y)}{n}
\frac{\lambda_n\cos(\omega t)+\omega\sin(\omega t)-\lambda_n e^{-\lambda_n t}}
{\lambda_n^2+\omega^2}
$$

[Background: Stokes' second problem and oscillatory Couette flow](https://arxiv.org/abs/2107.12894)

## Examples

Cases are defined explicitly in [src/cases.hpp](src/cases.hpp).
Choose a case with the first argument and a new output directory with the second.
Edit its dimensions, mesh, time step, duration and four walls to try another setup.

### Lid-driven cavity: water, no temperature

A 1 cm square cavity starts at rest. The top wall moves to the right at
0.01004 m/s; the other three walls stay fixed. Water properties:
density 998.2 kg/m³, kinematic viscosity 1.004e-6 m²/s. Reynolds number: $`\mathrm{Re}=UH/\nu=100`$, where $`U`$ is the lid speed and $`H`$ the cavity height.
Temperature and buoyancy are disabled.

![Water lid-driven cavity: transient velocity and circulation](docs/media/lid_cavity.gif)

Starts at rest; first 5 s at real speed, then 5x.
The colour shows speed; moving white particles reveal circulation.
Faint arrows show velocity direction and magnitude.
The cyan arrow indicates the motion of the lid.

```sh
cmake --build build -j
./build/simulate lid_cavity results_lid_cavity
.venv/bin/python visualize.py results_lid_cavity
.venv/bin/python visualize.py results_lid_cavity --save results_lid_cavity/lid_cavity.gif
```

32 × 32 cells, dt = 0.005 s, 8000 steps, one snapshot every 0.2 s.
This is a simulation example; no cavity benchmark validation is claimed.

### Buoyancy: air heated from below

A 4 cm × 2 cm rectangle starts at rest near 20°C. The bottom wall is held
at 40°C and the top wall at 20°C; the side walls are insulated ($`\partial T/\partial x=0`$, no heat crosses them). All four
walls are fixed and impermeable ($`\mathbf{u}=\mathbf{0}`$ at the walls: no slip or flow through them). Temperature and velocity evolve together
with the Boussinesq buoyancy force.

![Air heated from below: temperature, rising plumes and moving particles](docs/media/buoyancy.gif)

0–8 s of physical time, slowed down to a 15 s video.
Colours show temperature in °C; particles and arrows show the motion of the air.

```sh
cmake --build build -j
./build/simulate buoyancy results_buoyancy_demo
.venv/bin/python visualize.py results_buoyancy_demo
.venv/bin/python visualize.py results_buoyancy_demo --save results_buoyancy_demo/buoyancy.gif
.venv/bin/python visualize.py results_buoyancy_demo --field speed
```

32 × 16 cells, dt = 0.002 s, 4000 steps: 8 s simulated, one snapshot every 0.1 s.
Air: density 1.204 kg/m³, kinematic viscosity 1.516e-5 m²/s,
thermal diffusivity 2.14e-5 m²/s and expansion coefficient 3.43e-3 K⁻¹.
The initial temperature contains a deterministic 0.02 K sinusoidal perturbation
to break horizontal symmetry and seed convection; initial velocities are zero.
This is a qualitative example on a small grid.

[Background: natural convection in air](https://www.comsol.com/blogs/simulating-natural-convection-in-air)

### Rising thermal bubble: warm water in cold water

A circular patch of water near 35°C starts near the bottom of water at 20°C.
The entire domain starts at rest. Buoyancy makes the warm water rise and spread,
while the surrounding water recirculates. This is a temperature anomaly within
a single liquid, without a free surface, air bubble or surface tension.

![Warm water rising from rest, with temperature and tagged fluid particles](docs/media/thermal_bubble.gif)

4 s of physical time shown in a 15 s animation. Colour shows temperature in °C.
Dark particles mark water that was initially hot; white particles show the
surrounding circulation. All four walls are fixed and impermeable ($`\mathbf{u}=\mathbf{0}`$: no slip or flow through them), and insulated ($`\partial T/\partial n=0`$: no heat crosses them), where $`n`$ is the direction normal to each wall.

```sh
cmake --build build -j
./build/simulate thermal_bubble results_thermal_bubble_demo
.venv/bin/python visualize.py results_thermal_bubble_demo
.venv/bin/python visualize.py results_thermal_bubble_demo --save results_thermal_bubble_demo/thermal_bubble.gif
```

Domain: 20 mm × 30 mm. Initial patch: centre (10 mm, 6 mm), radius 3 mm,
with a smooth temperature transition at its edge.
Water: density 998.2 kg/m³, kinematic viscosity 1.004e-6 m²/s,
thermal diffusivity 1.43e-7 m²/s and expansion coefficient 2.07e-4 K⁻¹.
Mesh: 48 × 72; dt = 0.002 s; 2000 steps; snapshots every 0.04 s.
Pressure uses relaxation 1.7 to reduce iteration cost.
This is a qualitative convection example.

### Oscillating inlet jet: water in a channel

![Oscillating water jet: velocity and passive tracers over the first 45 seconds](docs/media/oscillating_jet.gif)

A parabolic jet enters a 60 mm × 30 mm channel through a 4.5 mm opening centred
at $`y_c=19.5\,\mathrm{mm}`$. Its peak speed is $`U=0.04462\,\mathrm{m/s}`$:

```math
U(y)=\begin{cases}
U\left[1-\left(\frac{y-y_c}{b}\right)^2\right], & |y-y_c|<b,\\
0, & |y-y_c|\ge b,
\end{cases}
\qquad b=2.25\,\mathrm{mm}.
```

The direction is $`\theta(t)=20^\circ\sin(2\pi t/\tau)`$, with $`\tau=4\,\mathrm{s}`$.
The inlet components are $`u(0,y,t)=U(y)\cos\theta(t)`$ and $`v(0,y,t)=U(y)\sin\theta(t)`$.
The top, bottom and the inlet outside the opening are fixed and impermeable.
The outlet has prescribed $`P=0`$; its horizontal velocity is calculated by the
momentum predictor and pressure correction. The initial velocity is projected
to satisfy incompressibility and the imposed inlet flux.
The inlet normal velocity for the new time is included before the pressure solve,
so each projected velocity field satisfies the changing volume flux.
Temperature and buoyancy are disabled. No additional perturbation is imposed.

Water: $`\nu=1.004\times10^{-6}\,\mathrm{m^2/s}`$;
$`\mathrm{Re}=U(2b)/\nu=200`$, based on peak inlet speed and opening width.
Mesh: 192 × 96 cells, $`\Delta t=0.0025\,\mathrm{s}`$, 50 s physical duration.
This is a qualitative laminar example.

```sh
cmake --build build -j
./build/simulate oscillating_jet results_oscillating_jet_20deg_re200_50s
.venv/bin/python visualize.py results_oscillating_jet_20deg_re200_50s --field speed --save results_oscillating_jet_20deg_re200_50s/oscillating_jet_particles.mp4
```

The video uses a fixed colour scale, velocity arrows and passive tracers;
The README GIF shows the first 45 s at real-time speed.
The full video shows 50 s of physical evolution in 50 s. Use `--field vorticity` for the vorticity field.
Velocity inlets currently support fluid-only cases with a pressure outlet on the right.
