# Multi-Stage Linear Decision Rules for Stochastic Control of Natural Gas Networks with Linepack

by Vladimir Dvorkin, Dharik Mallapragada, Audun Botterud, Jalal Kazempour and Pierre Pinson.
* * *

Repository contains data, code and supplemental materials associated with this [preprint](www.arxiv.org). If you find this preprint and code usefull for you research, please cite the preprint.

For the formulation of the stochastic OPF and network topology optimization problems, please refer to [Appendix.pdf](https://github.com/wdvorkin/LDR_for_gas_network_control/blob/main/Appendix.pdf).

The code is implemented in [Julia](https://julialang.org) (v1.6). To run the code, you need to download, install and license (free for academic use) the [Mosek](https://www.mosek.com) solver to run the code. 

To see the code in action, first open a terminal and clone this repository by running
```
git clone https://github.com/wdvorkin/LDR_for_gas_network_control
```
Then, open Julia and activate the environment by running
```
$ julia 
julia> ]
pkg> activate .
pkg> instantiate
```
where ```julia``` is an alias to Julia installation.
