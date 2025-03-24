**Numerical-Simulation-of-the-Evolution-of-Water-Meniscus-Based-on-Peridynamics**

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

- Not apply an applicable particle shifting technology (PST).

- Not consider the evaporation effect.
