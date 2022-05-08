# nlgreyfast
## Toolbox for identification of nonlinear state-space grey-box models

📜 You can find a **poster** about the toolbox below:

<a href="https://raw.githubusercontent.com/meco-group/nlgreyfast/master/flanders_make_conference_2022_retzler_nlgreyfast.pdf">![fmconf2022](flanders_make_conference_2022_retzler_nlgreyfast.png?raw=true)</a>

> A. Retzler, J. Swevers, J. Gillis, and Zs. Kollár, "Nlgreyfast: toolbox for nonlinear grey-box identification", *Flanders Make Conference*, Gent, Belgium, 2022.

📽 For a detailed introduction, watch the **talks** about this project here:

Talk at the [6th Workshop on Nonlinear System Identification Benchmarks](https://www.nonlinearbenchmark.org/history):

<a href="https://youtu.be/4na776RaUDs">![youtube](GitHubYouTubeBenchmarks.png?raw=true)</a>

A shorter talk at the [IEEE 17th International Conference on Advanced Motion Control](http://static.gest.unipd.it/AMC2022/):

<a href="https://youtu.be/J4RziJQEDDE">![youtube](GitHubYouTube.png?raw=true)</a>

🖼 We also have **visual documentation** that e.g. shows some relations between the formulas and the code, open this large poster, and 🔍 zoom into it with `Ctrl+mouse wheel`:

<a href="https://www.figma.com/file/bq76Z6uR22iQoa1iNRTzX6/nlgreyfast-overview-figure">![figma](GitHubFigmaZoomedOut.png?raw=true)</a>

## Compact summary

This toolbox supports estimating the parameters `theta` of nonlinear state-space models in the form:

```
x[k+1]=f(x[k],u[k],theta)
y[k]=g(x[k],u[k],theta)
```

Such models can often be derived based on physical knowledge of the system. If we have continuous-time ODEs describing the system, we can get a discrete-time model using integration (RK4) from the continuous-time model.

The parameter estimation task is ultimately about solving an optimization problem, of which this toolbox supports multiple formulations (these can also be regarded as a form of optimal control):
* single shooting (SS), 
* multiple shooting (MS), 
* partially constrained multiple shooting (PCMS),
* minimizing the N-step ahead prediction error (PEM).

The objective is always to minimize the least squares error between the noisy, measured system output and the simulated one, and is almost always non-convex, so the initial guess of the parameters play an improtant role. The MS, PCMS and especially the PEM formulations are more robust against initial values of the decision variables in a range where the system becomes unstable, in other words the simulation of the states becomes non-contractive.

Moreover, if the states can be measured or inferred, this toolbox can use this information in the estimation (it is required for PEM, but it can also be used in PCMS/MS for initializing the decision variables). In our experiments on a simple, simulated electromechanical system (see EMPS later), we have found that if we have the measured and inferred states, the combined PEM+SS method is the most robust against a non-contractive initialization point, and it is also the best choice with respect to running time.  

For optimization this tool uses [CasADi](https://casadi.org), so that it can solve certain problems a magnitude faster than `nlgreyest` that ships with the System Identification Toolbox of MATLAB. This can be interesting for applications of online identificiation, e.g. identification combined with control.

You can find more details about these [in our paper](https://ieeexplore.ieee.org/document/9729299).

## Quick start 

The basic steps of identification using this toolbox:

1. Design, build an perform experiment for system identification. Measure plant input and output. 
2. Define the grey-box model of your system in CasADi, as a continuous or discrete time nonlinear state-space model.
3. Choose formulation for the optimization problem, use `nlgreyfast` to estimate the model parameters.

All right, how to install it?  
First, install the dependencies:

- MATLAB (tested with R2019b on Ubuntu),
- CasADi (tested with version 3.5.5), available [here](https://web.casadi.org/get/).
- Optional: MathWorks System Identification Toolbox (if not available, you can still use the `nlidcomb` interface).
- Optional: MATLAB Coder

Second, go through the script `nlgreyfast_emps_example.m` step by step.  
For an example of the `nlidcomb` interface, see: `nlidcomb_msd_example.m`.  

## Credits

Authors of the code in this package:
- András Retzler
- Joris Gillis -- `nlgreyfast` is based on Joris' CasADi implementation of the MS and SS method.

## Acknowledgement

This research is supported by Flanders Make through ICON project ID2CON: Integrated IDentification for CONtrol, and by the Research Foundation - Flanders (FWO - Flanders) through project G0A6917N.

Special thanks to Alexandre Janot, Mathieu Brunot, Maxime Gautier for their EMPS model and measurement data:

> A. Janot, M. Gautier and M. Brunot, "Data Set and Reference Models of EMPS", *2019 Workshop on Nonlinear System Identification Benchmarks*, Eindhoven, The Netherlands, April 10-12, 2019.

The source of data in `nlid_emps_sim_data_training.mat` and part of the code related to EMPS is the package available [here](https://www.nonlinearbenchmark.org/benchmarks/emps).

## Cite

If you use this toolbox in your research, please [cite](https://zbib.org/8c5a12e1d89740609ec30d30e1e61a94):

> A. Retzler, J. Swevers, J. Gillis, and Zs. Kollár, "Shooting methods for identification of nonlinear state-space grey-box models", *IEEE 17th International Conference on Advanced Motion Control*, Padova, Italy, 2022.

