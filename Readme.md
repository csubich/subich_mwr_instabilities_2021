# Instabilities in the shallow-water system with a semi-Lagrangian, time-centered discretization

** Figure generating scripts **

Christopher Subich, July 2021

This repository contains the scripts necessary to regenerate the figures and data used in the above-titled paper, submitted to Monthly Weather Review in March 2021.

## Basic execution:

To generate all data and figures, run:

    # Generate data for growth-rate figures (about 72h on a desktop system)
    generate_data # produces saved data files

    # Generate figures
    mkdir paper_figs
    plot_paper

The ancillary files are:

* `cubic_instab.m`, `fourier_instab.m`, `linear_instab.m` -- calculates maximal growth rates of the shallow-water system over a wide range of Courant numbers and topographic wavelengths; parameters are modifiable in their respective preambles.
* `fd_ops.m` -- Helper function to generate finite-difference operators on the set of staggered grids
* `gen_ic.m` -- Helper function to generate a background steady state corresponding to a specific (and provided) $u_0$, $\Delta t$, $\phi_0$, and hill profile, via iteration of the linearized steady-state system
* `get_linearized_op.m` -- Helper function to provide the linearized shallow-water operators for a given background state
* `interesting_spectra.m` -- Collection of "interesting" points in Courant number / topographic wavenumber space, selected by the `CASE` variable
* `interp_cubic.m`, `interp_fourier.m`, `interp_linear.m` -- helper functions to compute the semi-Lagrangian operators for a given set of trajectories, using cubic, Fourier (trigonometric), or linear interpolation respectively.
* `label_eigvecs.m` -- Helper script to associate the eigenvectors of the advection or timestepping operator with a dominant physical mode
* `plot_growthrate.m` -- Helper script to give a contour plot of maximum growth rate over varying Courant number and hill wavenumbers, using the data produced by `*_instab.m`
* `plot_offcentering_graph.m` -- Separate script (called by `paper_figs.m`) to calculate and graph growth rates over varying hill amplitudes and off-centering parameters
* `plot_paper.m` -- Main plotting script, which generates all figures in the submitted paper
* `semilag_traj.m` -- Helper function to compute (iteratively) semi-Lagrangain departure points corresponding to a given velocity field
* `show_spectrum.m` -- Helper script to calculate and plot unstable modes in a particular realization of the shallow water operator
* `test_with_offcentering.m` -- Helper script to define an off-centered shallow water operator and compute its modes
*  time_convergence.m -- Script to generate a sample stable/unstable comparison, with reference solution (`ode45` and Fourier spectral method), from identical initial conditions
* `timestep.m` -- script to evaluate a shallow water case over time, to observe any instability growing
