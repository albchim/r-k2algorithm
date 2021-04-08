# Learning the Topology of a Bayesian Network Using K2 Algorithm

### Prerequisites
Python versions supported:

[![](https://img.shields.io/badge/python-3.7+-blue.svg)](https://badge.fury.io/py/root_pandas)

### Authors:

- [Alberto Chimenti](https://github.com/albchim) (University of Padova)
- [Paolo Frazzetto](https://github.com/paolofraz) (University of Padova)
- [Vincenzo Schimmenti](https://github.com/Schimmenti) (University of Padova)


# Brief Introduction
A Byesian belief-network structure $B_s$ is a directed acyclic graph in which nodes represent domain variables and arcs between nodes represent probabilistic dependanceies. A variable in such a structure may be continuous or discrete. In the following approach we will consider the discrete case.
A BBN structure gives us informations about $\texttt{dependency relationships}$ among the nodes and can be used for $\texttt{probabilistic inference}$ over the system, providing a direct way to compute joint probabilities, given the conditional ones associated with the graph's arcs. Since the joint probability of a given instantiation of the system variables can be computed as:
$$P(X_1,X_2,...,X_n)=\prod_{i=1}^{n}P(X_{i}|\pi_{i})$$
where $\pi_{i}$ are the "parent" nodes of the variable $x_{i}$

Using a Bayesian approach such a network can be constructed starting from a database which presents several records of the values combination of the system variables.
One can find the most probable belief-network structure, given the dataset.

The `K2 Algorithm` gives a computationally optimized way to search for the most probable structure.

Instead of maximizing the probability over all structures, one assumes that a node has no parents and incrementally adds that parent whose addition most increases the probability of the resulting structure.

Therefore one computes the scoring function:
$$g(i,\pi_i)=\prod_{j=1}^{q_{i}}\frac{(r_{i}-1)!}{(N_{ij}+r_{i}-1)}\prod_{k=1}^{r_{i}}N_{ijk}!\;\;\;\;\;\;\;\;\; (4)$$
where $\pi_i$ is the set of parents for node i which, following the given ordering, is composed by all the variables candidate to be parents of the $i$th one. 
One incrementally adds one parent node to $\pi_i$ at each iteration and if $g(i,\pi_i')>g(i,\pi_i)$ the node is permanently added to the parents and the structure probability $g(i,\pi_i)$ is updated.



A complete explanation of the methods and results obtained can be found in the notebook Project3_CFS_final.ipynb


### Sources:
* [G.F.Cooper and E. Herskovits](https://link.springer.com/article/10.1007/BF00994110)
* [C. Ruiz](http://web.cs.wpi.edu/~cs539/s05/Projects/k2_algorithm.pdf)
* [A. Franzin et al.](https://cran.r-project.org/web/packages/bnstruct/vignettes/bnstruct.pdf)
* [`bnstruct`](https://github.com/sambofra/bnstruct)
* [`bnlearn`](http://www.bnlearn.com/)