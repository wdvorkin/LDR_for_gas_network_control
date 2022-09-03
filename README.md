# GasLDR: Linear Decision Rules for Stochastic Control of Gas Networks

by Vladimir Dvorkin, Dharik Mallapragada, Audun Botterud, Jalal Kazempour and Pierre Pinson.
* * *

The repository contains data, code, and supplemental materials associated with this [preprint](https://arxiv.org/abs/2110.02824). 
<!-- If you find this preprint and code usefull for you research, please cite the preprint. -->

For the formulation of the stochastic OPF and network topology optimization problems, please refer to [Appendix.pdf](https://nbviewer.org/github/wdvorkin/LDR_for_gas_network_control/blob/main/Appendix.pdf).

The code is implemented in [Julia](https://julialang.org) (v1.6). To run the code, you will need the license for the [Mosek](https://www.mosek.com) solver, which is free for academic use.

To see the code in action, open a terminal and clone this repository by running
```
git clone https://github.com/wdvorkin/LDR_for_gas_network_control
```
Then, ```cd``` to project directory, and run the following command 
```
$ julia --project=@. main.jl 
```
where ```julia``` is an alias to Julia installation. This command optimizes the base stochastic control policy and stores the results in ```~/output```. The terminal output will look like this:
<img width="900" alt="Screen Shot 2021-10-06 at 10 12 32 AM" src="https://user-images.githubusercontent.com/31773955/136165304-d69c4f01-4714-49e5-b2de-0ea378606f41.png">

You can also specify different options, e.g., by typing
```
$ julia --project=@. main.jl -p 100
```
the program returns the variance-aware solution with the pressure variability penalty equal to 100. To see the list of all available options, type 
```
$ julia --project=@. main.jl --help
```

