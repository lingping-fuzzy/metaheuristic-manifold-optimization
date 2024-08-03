
function [xcost, bestall] = basicexamplerun(n, method, w0)
    % an example data is uploaded in /data/ folder. 
   filename = strcat('data\','eign-', int2str(n), '.mat');
   
   % method- 
    % Verify that Manopt was indeed added to the Matlab path.
    if isempty(which('spherefactory'))
        error(['You should first add Manopt to the Matlab path.\n' ...
		       'Please run importmanopt.']);
    end
    
    % Generate the problem data.
    A=load(filename).A;

    % Create the problem structure.
    manifold = spherefactory(n);
    problem.M = manifold;
    
    % Define the problem cost function and its gradient.
    problem.cost  = @(x) -x'*(A*x);
    problem.egrad = @(x) -2*A*x;
    problem.ehess = @(x, xdot) -2*A*xdot;
    
    % Numerically check gradient and Hessian consistency.

    % Solve.
     pn = 20; 
     itmax = 100; % parameters user defined
     [x, xcost, ~] = mDTMA(problem, pn, itmax, w0, []);          %#ok<ASGLU>

    
end
