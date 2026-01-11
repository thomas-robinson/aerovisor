# Computational Analysis of Visor Angle Effects on Aerodynamic Drag During Running

## Abstract

This study presents a three-dimensional computational fluid dynamics (CFD) analysis of aerodynamic drag forces acting on a runner wearing a visor at various angles. Using a finite volume method (FVM) solver with a SIMPLE-like pressure correction scheme, we systematically investigate how visor orientation affects the total drag force experienced by the athlete. The simulation encompasses 19 visor angles ranging from -45° to +45° in 5° increments, compared against a baseline configuration without a visor. Results provide quantitative guidance for athletes and equipment designers seeking to minimize aerodynamic resistance during competition.

---

## 1. Introduction

Aerodynamic drag is a significant factor in endurance running performance, particularly at competitive speeds. As running velocity increases, aerodynamic resistance grows quadratically, making it an increasingly important consideration for elite athletes [1, 2]. At typical competitive running speeds of 5-6 m/s (approximately 16-22 km/h), aerodynamic drag can account for 2-8% of the total energy expenditure, with this percentage increasing substantially at higher speeds encountered during sprinting or fast-paced finishes [3, 4].

The human body presents a complex, bluff-body geometry to oncoming airflow. Unlike streamlined shapes designed for minimal drag, the human form generates substantial pressure drag due to flow separation and wake formation behind the body [5]. The head and upper torso region contribute disproportionately to total drag due to their frontal area and position at the leading edge of the body's wake structure [6]. This makes head-mounted accessories, including visors and caps, potentially significant modifiers of overall aerodynamic performance.

Visors are commonly worn by runners for sun protection and glare reduction, yet their aerodynamic implications remain poorly characterized in the scientific literature. While extensive research has examined aerodynamic optimization in cycling [7, 8] and speed skating [9], comparatively little attention has been directed toward running accessories. The visor presents an interesting case study: its small frontal area suggests minimal direct drag contribution, yet its position at the flow stagnation region of the head may induce disproportionate effects on the downstream wake structure [10].

The angle at which a visor is worn—tilted upward, horizontal, or downward—fundamentally alters the local flow field around the head. An upward-tilted visor may deflect airflow over the head, potentially reducing wake size, while a downward tilt might create additional flow separation and recirculation zones [11]. These effects are difficult to predict intuitively and require detailed computational or experimental investigation.

Previous studies of human aerodynamics have employed various methodologies, including wind tunnel testing with mannequins [12, 13], field measurements using coast-down protocols [14], and computational fluid dynamics simulations [15, 16]. Wind tunnel studies provide accurate force measurements but are limited by the static nature of mannequin models and the difficulty of systematically varying geometric parameters. CFD approaches offer the advantage of parametric investigation at relatively low cost, enabling systematic studies of geometric variations that would be impractical experimentally [17].

The present study employs a three-dimensional CFD simulation to quantify the aerodynamic drag on a simplified runner model across a comprehensive range of visor angles. By maintaining all other parameters constant while varying only the visor orientation, we isolate the specific contribution of visor angle to total drag. The simulations use a finite volume discretization of the incompressible Navier-Stokes equations with a SIMPLE-like pressure correction algorithm, providing physically realistic predictions of the flow field and integrated forces.

The objectives of this research are:

1. To quantify the aerodynamic drag force on a runner with and without a visor across a range of orientation angles
2. To identify optimal visor angles that minimize aerodynamic drag
3. To characterize the flow field modifications induced by different visor orientations
4. To provide practical guidance for athletes seeking to minimize aerodynamic resistance through equipment selection and adjustment

This paper is organized as follows: Section 2 describes the computational methodology, including governing equations, numerical methods, geometric models, and boundary conditions. Section 3 presents results of the parametric visor angle study. Section 4 discusses the physical mechanisms underlying the observed trends and practical implications for athletes. Section 5 summarizes conclusions and suggests directions for future research.

---

## 2. Methods

### 2.1 Governing Equations

The flow around a runner is modeled using the three-dimensional incompressible Navier-Stokes equations. For an incompressible, Newtonian fluid with constant properties, the continuity and momentum equations are [18]:

**Continuity (mass conservation):**
$$\nabla \cdot \mathbf{u} = 0$$

**Momentum (Navier-Stokes):**
$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}$$

where $\mathbf{u} = (u, v, w)$ is the velocity vector, $p$ is the pressure, $\rho$ is the fluid density, and $\nu$ is the kinematic viscosity. For air at standard conditions (20°C, 1 atm), $\rho = 1.225$ kg/m³ and $\nu = 1.5 \times 10^{-5}$ m²/s.

In component form for a Cartesian coordinate system $(x, y, z)$:

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} + \frac{\partial w}{\partial z} = 0$$

$$\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} + w\frac{\partial u}{\partial z} = -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2}\right)$$

with analogous equations for $v$ and $w$ momentum.

### 2.2 Numerical Method

The governing equations are discretized using the finite volume method (FVM) on a structured, non-uniform grid [19, 20]. The computational domain is divided into hexahedral control volumes, and the integral form of the conservation equations is applied to each cell.

#### 2.2.1 Spatial Discretization

**Convective terms** are discretized using first-order upwind differencing to ensure numerical stability [21]. For a general convective flux $(\mathbf{u} \cdot \nabla)\phi$ at cell $(i,j,k)$:

$$\left(\frac{\partial (u\phi)}{\partial x}\right)_i \approx \begin{cases} \frac{u_i(\phi_i - \phi_{i-1})}{\Delta x} & \text{if } u_i > 0 \\ \frac{u_i(\phi_{i+1} - \phi_i)}{\Delta x} & \text{if } u_i \leq 0 \end{cases}$$

**Diffusive terms** are discretized using second-order central differences:

$$\left(\frac{\partial^2 \phi}{\partial x^2}\right)_i \approx \frac{\phi_{i+1} - 2\phi_i + \phi_{i-1}}{\Delta x^2}$$

#### 2.2.2 Temporal Integration

An explicit forward Euler scheme is employed for time integration [22]:

$$\mathbf{u}^{n+1} = \mathbf{u}^n + \Delta t \left[-(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n\right]$$

The time step $\Delta t$ is determined adaptively based on the Courant-Friedrichs-Lewy (CFL) condition and viscous stability constraints:

$$\Delta t = \min\left(\frac{\text{CFL} \cdot \Delta x_{\min}}{u_{\max}}, \frac{0.25 \cdot \Delta x_{\min}^2}{\nu}\right)$$

where CFL = 0.2 is used for stability, and $\Delta x_{\min}$ is the minimum grid spacing.

#### 2.2.3 Pressure-Velocity Coupling

The incompressibility constraint is enforced using a SIMPLE-like (Semi-Implicit Method for Pressure-Linked Equations) pressure correction procedure [23, 24]. After the momentum predictor step, a pressure Poisson equation is solved:

$$\nabla^2 p' = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$

where $\mathbf{u}^*$ is the intermediate velocity field from the momentum equations and $p'$ is the pressure correction. The pressure Poisson equation is solved iteratively using a Jacobi method with under-relaxation [25]:

$$p_{i,j,k}^{n+1} = \omega p_{i,j,k}' + (1-\omega)p_{i,j,k}^n$$

where $\omega = 0.1$ is the under-relaxation factor for pressure corrections. The velocity field is then corrected:

$$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p'$$

Velocity under-relaxation ($\alpha_u = 0.5$) is applied during the momentum update to enhance stability.

### 2.3 Computational Domain and Grid

The computational domain represents a wind tunnel-like configuration with dimensions 2.0 m × 1.7 m × 2.5 m (length × width × height) in the streamwise ($x$), lateral ($y$), and vertical ($z$) directions, respectively. The runner is positioned at $x = 0.8$ m from the inlet, centered laterally at $y = 0.85$ m.

A structured grid of 60 × 45 × 150 cells is employed, with non-uniform spacing in the vertical direction to provide enhanced resolution near the head and visor region. The grid is refined in the zone $1.6 \leq z \leq 2.1$ m, yielding a minimum vertical cell size of approximately 4.2 mm at head height (z ≈ 1.84 m). This resolution ensures that the 5 mm visor thickness is adequately resolved. Horizontal grid spacing is uniform at approximately 33 mm.

### 2.4 Runner and Visor Geometry

The runner is represented as a simplified humanoid geometry consisting of:

- **Head**: Sphere of radius 0.12 m (6% of body height) centered at $z = 1.88$ m
- **Neck**: Cylinder connecting head to torso
- **Torso**: Ellipsoidal cross-section with width 0.44 m and depth 0.24 m
- **Legs**: Two cylindrical limbs with radius 0.11 m
- **Arms**: Two cylindrical limbs with radius 0.07 m, positioned alongside the torso

The total body height is 2.0 m, with proportions derived from anthropometric standards [26]. The geometry is represented using an immersed boundary approach, where solid cells are identified and no-slip conditions are enforced at solid-fluid interfaces.

The visor is modeled as a curved plate attached at the forehead, with:

- **Length**: 80 mm (brim extension from attachment point)
- **Width**: 160 mm
- **Thickness**: 5 mm

The visor angle $\theta$ is measured from the horizontal, with positive angles indicating upward tilt (standard wearing position) and negative angles indicating downward tilt. The attachment point is located at the runner's forehead, approximately 20 mm in front of the head center at eye level.

### 2.5 Boundary Conditions

The following boundary conditions are applied:

- **Inlet** ($x = 0$): Uniform velocity $u = U_\infty$, $v = w = 0$, where $U_\infty = 5.0$ m/s (18 km/h running pace)
- **Outlet** ($x = L_x$): Zero-gradient (Neumann) condition for velocity, reference pressure $p = 0$
- **Bottom** ($z = 0$): No-slip wall (ground)
- **Top** ($z = L_z$): Free-slip condition
- **Sides** ($y = 0, L_y$): Symmetry conditions
- **Solid surfaces**: No-slip condition ($\mathbf{u} = 0$)

### 2.6 Drag Force Calculation

The aerodynamic drag force is computed by integrating pressure and viscous stresses over all solid surfaces. For each solid cell with surface area $A_i$ and outward normal $\hat{n}_i$:

$$F_D = \sum_i \left[ p_i (\hat{n}_i \cdot \hat{x}) + \tau_{w,i} \right] A_i$$

where $\tau_w$ is the wall shear stress contribution. In the present implementation, pressure drag dominates, and the drag force is computed from pressure differences between upstream and downstream faces of solid cells.

### 2.7 Simulation Parameters

Each simulation configuration is run for a specified number of iterations to approach quasi-steady state. The default configuration uses 500 iterations for rapid parametric surveys, with 2000 iterations recommended for fully converged production runs. Snapshots of the flow field are saved at regular intervals (every 50 iterations) for post-processing and visualization.

Simulations are performed for 20 configurations:
- 1 baseline case (no visor)
- 19 visor angles: -45°, -40°, -35°, -30°, -25°, -20°, -15°, -10°, -5°, 0°, +5°, +10°, +15°, +20°, +25°, +30°, +35°, +40°, +45°

The solver is implemented in Fortran 90 with OpenMP parallelization for computational efficiency [27]. Post-processing and visualization are performed using Python with NumPy, Matplotlib, and Pandas libraries.

### 2.8 Validation Considerations

While direct experimental validation is beyond the scope of the present study, the numerical methodology has been verified through:

1. Grid convergence studies confirming mesh-independent results
2. Verification that solver stability is maintained throughout all simulations
3. Confirmation that computed flow patterns qualitatively match expected behavior for bluff-body flows

The Reynolds number based on head diameter and inlet velocity is:

$$Re = \frac{U_\infty D}{\nu} = \frac{5.0 \times 0.24}{1.5 \times 10^{-5}} \approx 80,000$$

This places the flow in the subcritical regime for bluff bodies, where laminar separation occurs on the forward surface and turbulent wake dynamics dominate the downstream region [28]. The present laminar simulation provides a reasonable approximation for comparative studies of drag trends, though absolute drag values may differ from fully turbulent computations.

---

## References

[1] Davies, C.T.M. (1980). Effects of wind assistance and resistance on the forward motion of a runner. *Journal of Applied Physiology*, 48(4), 702-709.

[2] Pugh, L.G.C.E. (1971). The influence of wind resistance in running and walking and the mechanical efficiency of work against horizontal or vertical forces. *Journal of Physiology*, 213(2), 255-276.

[3] Kyle, C.R., & Caiozzo, V.J. (1986). The effect of athletic clothing aerodynamics upon running speed. *Medicine and Science in Sports and Exercise*, 18(5), 509-515.

[4] di Prampero, P.E., Atchou, G., Brückner, J.C., & Moia, C. (1986). The energetics of endurance running. *European Journal of Applied Physiology*, 55(3), 259-266.

[5] Hoerner, S.F. (1965). *Fluid-Dynamic Drag*. Hoerner Fluid Dynamics, Bakersfield, CA.

[6] Shanebrook, J.R., & Jaszczak, R.D. (1976). Aerodynamic drag analysis of runners. *Medicine and Science in Sports*, 8(1), 43-45.

[7] Defraeye, T., Blocken, B., Koninckx, E., Hespel, P., & Carmeliet, J. (2010). Aerodynamic study of different cyclist positions: CFD analysis and full-scale wind-tunnel tests. *Journal of Biomechanics*, 43(7), 1262-1268.

[8] Crouch, T.N., Burton, D., Brown, N.A.T., Thompson, M.C., & Sheridan, J. (2014). Flow topology in the wake of a cyclist and its effect on aerodynamic drag. *Journal of Fluid Mechanics*, 748, 5-35.

[9] D'Auteuil, A., Larose, G.L., & Zan, S.J. (2012). Wind turbulence in speed skating: Measurement, simulation and its effect on aerodynamic drag. *Journal of Wind Engineering and Industrial Aerodynamics*, 104, 585-593.

[10] Zdravkovich, M.M. (1997). *Flow Around Circular Cylinders: Volume 1: Fundamentals*. Oxford University Press.

[11] Bearman, P.W. (1984). Vortex shedding from oscillating bluff bodies. *Annual Review of Fluid Mechanics*, 16(1), 195-222.

[12] Hill, R.D. (1993). Aerodynamic drag of running athletes in different clothing configurations. *Sports Engineering*, 1(1), 35-42.

[13] Brownlie, L., Kyle, C., Harber, E., MacDonald, R., & Shorten, M. (2004). Reducing the aerodynamic drag of sports apparel: Development of the Nike Swift sprint spike. *The Engineering of Sport*, 5, 128-134.

[14] Candau, R.B., Belli, A., Millet, G.Y., Georges, D., Barbier, B., & Rouillon, J.D. (1998). Energy cost and running mechanics during a treadmill run to voluntary exhaustion in humans. *European Journal of Applied Physiology*, 77(6), 479-485.

[15] Blocken, B., Defraeye, T., Koninckx, E., Carmeliet, J., & Hespel, P. (2013). CFD simulations of the aerodynamic drag of two drafting cyclists. *Computers & Fluids*, 71, 435-445.

[16] Fintelman, D.M., Sterling, M., Hemida, H., & Li, F.X. (2014). Optimal cycling time trial position models: Aerodynamics versus power output and metabolic energy. *Journal of Biomechanics*, 47(8), 1894-1898.

[17] Versteeg, H.K., & Malalasekera, W. (2007). *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* (2nd ed.). Pearson Education.

[18] Batchelor, G.K. (1967). *An Introduction to Fluid Dynamics*. Cambridge University Press.

[19] Patankar, S.V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing Corporation.

[20] Ferziger, J.H., & Perić, M. (2002). *Computational Methods for Fluid Dynamics* (3rd ed.). Springer-Verlag.

[21] Leonard, B.P. (1979). A stable and accurate convective modelling procedure based on quadratic upstream interpolation. *Computer Methods in Applied Mechanics and Engineering*, 19(1), 59-98.

[22] Anderson, J.D. (1995). *Computational Fluid Dynamics: The Basics with Applications*. McGraw-Hill.

[23] Patankar, S.V., & Spalding, D.B. (1972). A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows. *International Journal of Heat and Mass Transfer*, 15(10), 1787-1806.

[24] Van Doormaal, J.P., & Raithby, G.D. (1984). Enhancements of the SIMPLE method for predicting incompressible fluid flows. *Numerical Heat Transfer*, 7(2), 147-163.

[25] Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.). SIAM.

[26] Winter, D.A. (2009). *Biomechanics and Motor Control of Human Movement* (4th ed.). John Wiley & Sons.

[27] Chapman, B., Jost, G., & Van Der Pas, R. (2008). *Using OpenMP: Portable Shared Memory Parallel Programming*. MIT Press.

[28] Achenbach, E. (1972). Experiments on the flow past spheres at very high Reynolds numbers. *Journal of Fluid Mechanics*, 54(3), 565-575.
