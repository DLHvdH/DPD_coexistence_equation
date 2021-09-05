# DPD_coexistence_equation

A coexistence (or binodal) equation for mixing in Dissipative Particle Dynamics

DPD_binodal.M is a Wolfram Mathematica module that can be used
to determine binodal value for aij as function of φ, aii, ajj, Ni, and Nj.
This function is called aijVDH. Its definition appears at the end of the module.

The equation considers 
- Binary mixtures consisting of two species, i and j
- Combinations of polymer lengths Ni and Nj
- Self-self interactions aii and ajj within the range of 10 to 60
- The mass or number fraction φ of component i (the fraction of component j is given by φ-1)

The resulting output is then the interaction strength aij between particles i and j for which the system is in equilibrium (i.e. on the binodal).

The origin of the functional form and constants of the used equations is explained in:
van der Haven, Dingeman LH, et al. "Closed-Form Coexistence Equation for Phase Separation of Polymeric Mixtures in Dissipative Particle Dynamics."
The Journal of Physical Chemistry B 125.27 (2021): 7485-7498. https://doi.org/10.1021/acs.jpcb.0c11274
