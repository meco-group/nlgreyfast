# nlgreyfast
## Toolbox for identification of nonlinear state-space grey-box models

This toolbox supports estimating the parameters `theta` of nonlinear state-space models:

```
x[k+1]=f(x[k],u[k],theta)
y[k]=g(x[k],u[k],theta)
```

The idea is to derive such a discrete-time model using integration (RK4) from a continuous-time physics-based model of a system, and our goal is to know accurately the `theta` parameters of this system. 

This can be done by first making experiments on the system: applying known excitation and measuring the output, then afterwards applying an identification methods like those in this toolbox, in order to find out the parameters in the formulas accurately. We get a so-called grey-box model as a result. 

This toolbox supports multiple formulations of the parameter estimation problem (can also be regarded as a form of optimal control):
* single shooting (SS), 
* multiple shooting (MS), 
* partially constrained multiple shooting (PCMS),
* minimizing the N-step ahead prediction error (PEM).

As the objective of the optimization here is almost always non-convex, the initial guess of the parameters, that we give to the solver, plays an improtant role. The MS, PCMS and especially the PEM formulation is more robust against initial values of the decision variables for which the system becomes unstable, the simulation of the states becomes non-contractive.

Moreover, if the states can be measured or inferred, this toolbox can use this information in the estimation (required for PEM, can be used in PCMS/MS for initializing the decision variables).

For optimization this tool uses [CasADi](https://casadi.org), so that it can solve certain problems a magnitude faster than `nlgreyest` that ships with the System Identification Toolbox of MATLAB. This can be interesting for applications of online identificiation, e.g. identification combined with control.

A visual representation of the supported formulations:

![GitHubReadmeFigure](GitHubReadmeFigure.png)

## Quick start 

The basic steps of identification using this toolbox:

1. Design, build an perform experiment for system identification. Measure plant input and output. 
2. Define the grey-box model of your system in CasADi, as a continuous or discrete time nonlinear state-space model.
3. Use `nlgreyfast` to estimate the model parameters.

First, install the dependencies:
- MATLAB (tested with R2019b on Ubuntu),
- CasADi (tested with version 3.5.5), available [here](https://web.casadi.org/get/).
- Optional: MathWorks System Identification Toolbox (if not available, you can still use the `nlidcomb` interface).
- Optional: MATLAB Coder

Second, go through the script `nlgreyfast_emps_example.m` step by step.  
For an example of the `nlidcomb` interface, see: `nlidcomb_msd_example`.  

## Credits

Authors of the code in this package:
- András Retzler
- Joris Gillis -- `nlgreyfast` is based on Joris' CasADi implementation of the MS and SS method.

## Acknowledgement

This research is supported by Flanders Make through ICON project ID2CON: Integrated IDentification for CONtrol, and by the Research Foundation - Flanders (FWO - Flanders) through project G0A6917N.

Special thanks to Alexandre Janot, Mathieu Brunot, Maxime Gautier for their EMPS model and measurement data:

> A. Janot, M. Gautier and M. Brunot, Data Set and Reference Models of EMPS, 2019 Workshop on Nonlinear System Identification Benchmarks, Eindhoven, The Netherlands, April 10-12, 2019.

The source of data in `nlid_emps_sim_data_training.mat` and part of the code related to EMPS is the package available [here](https://www.nonlinearbenchmark.org/benchmarks/emps).

## Cite

If you use this toolbox in your research, please cite:

> A. Retzler, J. Swevers, J. Gillis and Zs. Kollár, "Shooting methods for identification of nonlinear state-space grey-box models". In _Proceedings of IEEE 17th International Conference on Advanced Motion Control_. Padova, Italy, 2022.
