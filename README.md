# metaheuristic-manifold-optimization
 This study pioneered the proposal and implementation of a metaheuristic manifold optimization, introducing a novel directional transport operator to rectify convergence and local minima issues


# Understand the operator
Most implementations closely follow traditional metaheuristic algorithms by converting operations to tangent spaces, without utilizing operations that specifically take advantage of manifold structures, as opposed to traditional Euclidean spaces. We propose an operation tailored for manifold learning. The figure below shows that after the directional transport of motion v, the decomposed $v_2$ moves out of the manifold (surface) into the complementary space. The retraction operation then pulls $v_2$ back to the original point $x$, regardless of how far the point travels in this direction.






