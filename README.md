# metaheuristic-manifold-optimization
 This is the implementation of paper: [From Constraints Fusion to Manifold Optimization: A New Directional Transport Manifold Metaheuristic Algorithm]()


<details>
 <summary><b>contents </b></summary>
- *preliminary introduction*
- *source code*
- *How to reproduce the results*
- *Data*
- *Understanding the repository* code review
- * Directional Transport operator* (Understand the operator)
- *How to cite*
</details>

# Understand the operator
Most implementations closely follow traditional metaheuristic algorithms by converting operations to tangent spaces without utilizing operations that specifically take advantage of manifold structures, as opposed to traditional Euclidean spaces. We propose an operation tailored for manifold learning. The figure below shows that after the directional transport of motion v, the decomposed $v_2$ moves out of the manifold (surface) into the complementary space. The retraction operation then pulls $v_2$ back to the original point $x$, regardless of how far the point travels in this direction.



![see the figures for ](https://github.com/lingping-fuzzy/metaheuristic-manifold-optimization/figs/DTMA.png)


