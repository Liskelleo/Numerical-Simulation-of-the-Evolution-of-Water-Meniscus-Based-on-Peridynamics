**Numerical-Simulation-of-the-Evolution-of-Water-Meniscus-Based-on-Peridynamics**

*Documents*

The `dissertation.pdf` and the accompanying presentation ppt have been uploaded, but the Particle Shifting Technology (PST) algorithm is not included. The recently uploaded `report.pdf` contains the details of the PST algorithm (see Chapter 2.7).

*Configuration*

- MATLAB 2015b (on Windows 10)

*Highlights*

- Efficient Operation

  - Create a single main program file containing all the calculations, instead of writing the main program in several sub-scripts and calling them again.

  - Perform a KD tree neighbourhood search for a given particle by calling the library function `rangesearch`.

  - Utilize library functions `alphaShape` (also `inpolygon`, etc.) to identify the free surface.

  - Reduce the number of nested `for` loops by transforming the data set to be computed into arrays and matrices.

- Realistic Physics

  - Consider interactions between ‘solid-liquid-gas’ systems.

    - Calculate the surface tension based on the meshless particle method.

    - Consider solid wall adhersion and apply anti-penetration methods.

*Issues*

- Not consider the evaporation effect.
  
- **Updated on January 9, 2026**: PST applied here (optional in `main - pst.m`) is refered to Xu et al. (2009) and is only applicable to the fluid flow without free surface. Accroding to Gao et al. (2020), if the free surface is involved, the improved version of PST by using Fick’s law (Lind et al., 2012) should be adopted. Due to this limitation, the current algorithm has limited effectiveness. If you don't want to apply the PST algorithm, please set the parameter `C_pst` to `0`. It should be borne in mind that, while Fick’s law acts to diffuse particles at the free surface, surface tension acts in an opposite sense, encouraging particle cohesion.

  - *References regarding PST*

    - Xu, R., Stansby, P., & Laurence, D. (2009). Accuracy and stability in incompressible SPH (ISPH) based on the projection method and a new approach. _Journal of computational Physics, 228_(18), 6703-6725.

    - Lind, S. J., Xu, R., Stansby, P. K., & Rogers, B. D. (2012). Incompressible smoothed particle hydrodynamics for free-surface flows: A generalised diffusion-based algorithm for stability and validations for impulsive flows and propagating waves. _Journal of Computational Physics, 231_(4), 1499-1523.

    - Gao, Y., & Oterkus, S. (2020). Multi-phase fluid flow simulation by using peridynamic differential operator. _Ocean Engineering, 216_, 108081.
