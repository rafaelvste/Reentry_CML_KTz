# Reentry_CML_KTz
The codes in this repository reproduce all the simulations from the paper 'Cardiac reentry modeled by spatiotemporal chaos in a coupled map lattice', co-authored by Rafael V. Stenzinger and Marcelo H. R. Tragtenberg.

The paper was published in The European Physical Journal Special Topics - Complex Bio Rhythms (https://doi.org/10.1140/epjs/s11734-022-00473-1) and is also available on arXiv (https://arxiv.org/abs/2202.12406).

Most codes is standalone and has simple scripts to plot the results. The majority was written in Fortran 90, with a few in Matlab/Octave and Python. Fortran codes can be compiled using GFortran or ifort. They also require Gnuplot to open the plotting scripts (which are obviously optional).

It is important to note that numerical differences are expected in chaos simulations, given the sensitivity of the system to small numerical differences. Different compilers and hardware will lead to different numerical outcomes, even if the results are qualitatively the same.

Some of the parameters (like the number of iterations, number of realizations, etc) were simplified for quick running and testing. The values used in the paper are available in the caption of each figure.

Follows the description of which codes to use for each figure and some additional info:

KTz_cell.f90 is used to reproduce Fig.2a and the insets of Fig.3b. KTz_cell_phase_diagram_ISI.f90 is used to Fig.2b.

KTz_cell_bifurc_APD.f90 is used for Fig.3a and KTz_log_cell_Lyap_pacing.m for Fig.3b.

KTz_net_bifurc_APD.f90 is used for the bifurcation diagrams in Fig.4a and 4b.

KTz_log_net_Lyap_coupling.f90 is used for the Lyapunov exponent in Fig.4a and requires Lapack or Intel MKL libraries*. Adaptations of this code for Matlab/Octave and Python are also available and were written in attempts to calculate the Lyapunov spectrum of larger networks. Irrespective of the language, they should work on networks of arbitrary size, but require powerful hardware to run in a reasonable time. If using MKL, optimization flags when compiling are enormously helpful to speedup runtime. Flags can be obtained from the MKL Link Line Advisor page.

*The portion of this code that performs the LU factorization was adapted from an answer at https://stackoverflow.com/questions/40422172/is-there-a-command-or-subroutine-for-lu-factorization and credit goes to the author.

KTz_net_stand_dev_coupling.f90 is used for Fig.5a, KTz_net_stand_dev_size.f90 for Fig.5b and KTz_net_stand_dev_diagram.f90 for Fig.5c.

KTz_net_diff_coupl.f90 is used for Fig.6, Fig.8 and Fig.9b. The spatiotemporal portraits of Fig.7 and Fig.9a can be extracted from the last code or using KTz_net_diff_coupl_video.f90, used for the supplementary material videos (requires FFmpeg).

Videos.rar contains the videos that appear in the paper.
