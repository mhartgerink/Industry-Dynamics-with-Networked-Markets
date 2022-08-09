# Industry-Dynamics-with-Networked-Markets

This repository includes the code necessary for the computation of equilibrium values in the thesis "Industry Dynamics with Networked Markets". The most important files for the computation of equilibrium values are:
1. Primitive Matrices.jl
2. Capacity Constrained Networked Cournot Competition.jl
3. MPE_solver.jl

#
# Primitive Matrices
This file contains the functions used, in order to construct the matrices necessary for the equilibrium calculation. Most important are the LC-matrix, which contains all individual combinations of capacity and link choices, and the RSS-matrix, which is the reduced state space for two firms. Besides these two matrices, additional matrices are constructed of which most relate to identifying states with similar graph structures, or similar capacity structures. 

#
# Capactiy Constrained Networked Cournot Competition
This file contains a method that approximates the solution of the capacity constrained networked Cournot-competition game. The idea behind this method is to have the optimality conditions imposed, such that the solution to the unconstrained problem is simply the solution in the null-space. However, we only approximate the exact solution by using a non-negativite least-squares solver. After the unconstrained optimization problem is considered, additional constraints are imposed if the quantities supplied to the markets exceed one's capacity. For these constraints, we again use a non-negative least-squares solver and it is here that the error is being generated. The main reason for opting for this choice, is that we can already use an existing library, instead of putting additional constraints. Furthermore, the error remains relatively small. For these results, see the Appendix of the thesis.

#
# Markov-perfect Equilibrium solver
This file contains an algorithm that calculates the Markov-perfect Equilibrium (MPE) for our model. Using the optimal behaviour derived in the thesis, the expected value for each stage will be calculated for a given strategy of other players, and a given value function. In order to avoid cycling, we take the initial strategies of other players as the long-term average of all past iterations, furthermore it appears that the solution converges uniformly. At its core, the algorithm presented here is simply a reinforcement learning (RI)-algorithm, in which players gradually learn the optimal policies of other players. 
#

