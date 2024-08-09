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

 - *Directional Transport operator* (Understand the operator)
    
 - *Data*
 
- *How to reproduce the results*

- *How to cite*
- 
</details>

### preliminary introduction
A metaheuristic manifold algorithm shares a similar structure to a traditional algorithm in many operators and differs in the implementation that guarantees the population's movement on the manifold. To propose a metaheuristic manifold algorithm, the population movement must be constrained in a specific `manner' so that each individual of the population remains on the manifold during the evolution process. Essentially, each individual is an element ( or point) of the manifold. For instance, the solution (point or element) is on Sphere manifolds. 
A very good and detailed introduction, implementation can be found [here](https://github.com/NicolasBoumal/manopt) and [here](https://www.manopt.org/)
### Source code
Our implementation is based on repository ManOpt [here](https://github.com/NicolasBoumal/manopt), thanks Manopt group. 

```matlab-usage
runexample.m - an example code.
```

### Understanding the repository

# Understand the operator
To understand the Directional Transport operator, you just need to consider the decomposition of signal/image/ or any other familar data. 
<details>
Most implementations closely follow traditional metaheuristic algorithms by converting operations to tangent spaces without utilizing operations that specifically take advantage of manifold structures, as opposed to traditional Euclidean spaces. We propose an operation tailored for manifold learning. The figure below shows that after the directional transport of motion v, the decomposed $v_2$ moves out of the manifold (surface) into the complementary space. The retraction operation then pulls $v_2$ back to the original point $x$, regardless of how far the point travels in this direction.
</details>

### Data
 We experiment with six problem sets including: the dominant invariant $p$-subspace, $p \in \mathbb{R}$, Procrustes problem, Semidefinite programs, Truncated singular value decomposition (SVD) problem, Thomson problem , and robot manipulation. We randomly generate five initial settings (if one problem involve an initial matrix, we randomly generate one as baseline) to construct five problem variants for each problem set, termed dataset1, dataset2,$\ldots$, to dataset5.

 Robotic data is from paper [2]
### How to Reproduce the Results

It is the same as running any one metaheuristic algorithm. Use the following command:

```matlab-usage
[X, xcost, info, ~] = mDTMA(problem, pn, itmax, 0.1, []);
```

Here is a breakdown of the parameters:

- `problem`: Defines the cost, manifold, and all basic settings.
- `pn`: Initial population size.
- `itmax`: Maximum iteration number.
- `0.1`: A parameter that can be tuned for better performance, referred to as `w0` in the paper.
- `mDTMA`: THE PROPOSED METHOD NAME

``` Robot arm experiment, 
    shared code by paper [1][2], and if you have problems access their, write to me. 
```


### citation
 [1] Saveriano, M., Abu-Dakka, F.J. and Kyrki, V., 2023. Learning stable robotic skills on Riemannian manifolds. Robotics and Autonomous Systems, 169, p.104510.
 [2] Abu-Dakka, F.J., Rozo, L. and Caldwell, D.G., 2018. Force-based variable impedance learning for robotic manipulation. Robotics and Autonomous Systems, 109, pp.156-167.
### How to cite
[link](https://www.sciencedirect.com/science/article/pii/S1566253524003749?dgcid=author)

```
@article{snavsel2024constraints,
  title={From constraints fusion to manifold optimization: A new directional transport manifold metaheuristic algorithm},
  author={Sn{\'a}{\v{s}}el, V{\'a}clav and Kong, Lingping and Das, Swagatam},
  journal={Information Fusion},
  pages={102596},
  year={2024},
  publisher={Elsevier}
}
```
If you have questions, please write to lingping_kong@yahoo.com



![see the figures for ](https://github.com/lingping-fuzzy/metaheuristic-manifold-optimization/figs/DTMA.png)

