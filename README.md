# Multi-Stage Linear Decision Rules for Stochastic Control of Natural Gas Networks with Linepack

by Vladimir Dvorkin, Dharik Mallapragada, Audun Botterud, Jalal Kazempour and Pierre Pinson.
* * *

The repository contains data, code, and supplemental materials associated with this [preprint](www.arxiv.org). 
<!-- If you find this preprint and code usefull for you research, please cite the preprint. -->

For the formulation of the stochastic OPF and network topology optimization problems, please refer to [Appendix.pdf](https://github.com/wdvorkin/LDR_for_gas_network_control/blob/main/Appendix.pdf).

The code is implemented in [Julia](https://julialang.org) (v1.6). To run the code, you will need the license for the [Mosek](https://www.mosek.com) solver, which is free for academic use

To see the code in action, open a terminal and clone this repository by running
```
git clone https://github.com/wdvorkin/LDR_for_gas_network_control
```
Then, ```cd``` to project directory, and run the following command 
```
$ julia --project=@. main.jl 
```
where ```julia``` is an alias to Julia installation. This command optimizes the base stochastic control policy and stores the results in ```~/output```. The terminal output will look like this:
<img width="677" alt="Screenshot 2020-10-08 at 10 57 57" src="https://user-images.githubusercontent.com/31773955/95437348-286c7a00-0955-11eb-9e77-8d7745f8c09f.png">


You can also specify different options, e.g., by typing
```
$ julia --project=@. main.jl -p 100
```
the program returns the variance-aware solution with the pressure variance penalty equal to 100. To see the list of all available options, type 
```
$ julia main.jl --help
```

