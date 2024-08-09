function basicexample()
   
% https://math.stackexchange.com/questions/1823891/what-is-the-meaning-of-subtracting-from-the-identity-matrix
% https://math.stackexchange.com/questions/1182679/basis-of-tangent-space-of-a-sphere
% https://math.stackexchange.com/questions/554156/the-boundary-of-an-n-manifold-is-an-n-1-manifold

    % Verify that Manopt was indeed added to the Matlab path.
    if isempty(which('spherefactory'))
        error(['You should first add Manopt to the Matlab path.\n' ...
		       'Please run importmanopt.']);
    end
    
    % Generate the problem data.
    n = 30;
    A = randn(n);
    A = .5*(A+A'); % create a symmetric matrix
    
    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;
    
    % Define the problem cost function and its gradient.
    problem.cost  = @(x) -x'*(A*x);
    problem.egrad = @(x) -2*A*x;
    problem.ehess = @(x, xdot) -2*A*xdot;
    
    % Numerically check gradient and Hessian consistency.
%     figure;
%     checkgradient(problem);
%     figure;
%     checkhessian(problem);
 
    % Solve.
%     [x, xcost, info] = trustregions(problem);          %#ok<ASGLU>
%     [x, xcost, info] = pso(problem);          %#ok<ASGLU>
      [x, xcost, info] = gao(problem, 10, 100);          %#ok<ASGLU>
    % Display some statistics.
%     figure;
%     semilogy([info.iter], [info.gradnorm], '.-');
%     xlabel('Iteration #');
%     ylabel('Gradient norm');
%     title('Convergence of the trust-regions algorithm on the sphere');
    
end
