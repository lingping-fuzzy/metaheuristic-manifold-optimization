# metaheuristic-manifold-optimization
 This is the implementation of paper: [From Constraints Fusion to Manifold Optimization: A New Directional Transport Manifold Metaheuristic Algorithm]()

<details>
  <summary><b>Why we study Manifold metaheuristic algorithms </b></summary>
Metaheuristic algorithms are suited for combinatorial optimization problems, given that, although they are not usually guaranteed to find the optimal global solution, they can often find a sufficiently good solution in a decent amount of time. 
Meta-heuristics can also be easily applied to many problems because they are not problem-specific and often incorporate some form of randomness to escape from local minima. The generalization and implementation of numerical optimization algorithms to the manifolds have been well studied and successfully applied to actual problems from science, engineering, and robotics.
</details>

In this paper, we experiment with six problem sets include the dominant invariant $p$-subspace, Procrustes problem, Semidefinite programs (SDP), Truncated singular value decomposition (SVD) problem , Thomson problem, and robot manipulation (stiffness simulation). 


<details>
 <summary><b>Contents of repository</b></summary>
- *preliminary introduction* To understand this paper, you should have a basic understanding the manifold. 
 
- *Source code*
 
 - *Understanding the repository* code review

 - * Directional Transport operator* (Understand the operator)
    
 - *Data*
 
- *How to reproduce the results*

- *How to cite*
- 
</details>

### preliminary introduction
### Source code
### Understanding the repository

# Understand the operator
Most implementations closely follow traditional metaheuristic algorithms by converting operations to tangent spaces without utilizing operations that specifically take advantage of manifold structures, as opposed to traditional Euclidean spaces. We propose an operation tailored for manifold learning. The figure below shows that after the directional transport of motion v, the decomposed $v_2$ moves out of the manifold (surface) into the complementary space. The retraction operation then pulls $v_2$ back to the original point $x$, regardless of how far the point travels in this direction.


### Data
### How to reproduce the results


### How to cite


![see the figures for ](https://github.com/lingping-fuzzy/metaheuristic-manifold-optimization/figs/DTMA.png)


