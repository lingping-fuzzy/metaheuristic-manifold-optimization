%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%François Brancaleone                              %
%Master's thesis: Global optimization on manifolds %
%Advisor: P.-A. Absil                              %
%Readers: L. Jacques and P.-Y. Gousenbourger       %
%June 2018                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbest, fbest, info, options] = simulatedannealing(problem, pn, itmax, x, options)
% Simulated Annealing (SA) method for derivative-free minimization.
%
% function [x, cost, info, options] = simulatedannealing(problem)
% function [x, cost, info, options] = simulatedannealing(problem, x0)
% function [x, cost, info, options] = simulatedannealing(problem, x0, options)
% function [x, cost, info, options] = simulatedannealing(problem, [], options)
%
% To specify options whilst not specifying an initial guess, give x0 as []
% (the empty matrix).
%
% Apply the Simulated Annealing minimization algorithm to
% the problem defined in the problem structure, starting with the
% population x0 if it is provided (otherwise, a random population on the
% manifold is generated). A population is a cell containing points on the
% manifold. The number of elements in the cell must match the parameter
% options.populationsize.
%
% The outputs xbest and fbest are the best reached point on the manifold and its
% cost. The struct-array "info" contains information about the iterations:
%   iter : the iteration number
%   cost : cost value
%   time : elapsed time in seconds
%   gradnorm : Riemannian norm of the gradient
%   costevals : number of evaluations of the cost function
%   And possibly additional information logged by options.statsfun.
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%   maxcostevals (10000)
%       The algorithm terminates if maxcostevals is reached.
%   maxiter (500000)
%       The algorithm terminates if maxiter iterations have been executed.
%   maxtime (Inf)
%       The algorithm terminates if maxtime seconds elapsed.
%   populationsize (40)
%       The size of the population of points.
%   maxtemp (2000)
%       Maximal temperature (parameter of the SA method).
%   mintemp (10e-10)
%       Minimal/final temperature (parameter of the SA method).
%   iter_temp
%       Number of iterations for each temperature (parameter of the SA method).
%
% Original author: François Brancaleone, June 10, 2018.


% Verify that the problem description is sufficient for the solver.
if ~canGetCost(problem) %can the cost function be computed?
    warning('manopt:getCost', ...
        'No cost provided. The algorithm will likely abort.');
end

if ~canGetGradient(problem) % can the gradient function be computed?
    warning('manopt:getGradient', ...
        'No gradient provided.');
    indicator_grad=0; % no gradient provided
else
    indicator_grad=1; % the gradient is provided
end


% Dimension of the manifold
dim = problem.M.dim();


% Set local defaults here
localdefaults.storedepth = 0;      % no need for caching
localdefaults.maxcostevals = max(10000, 2*dim);

    localdefaults.maxiter = max(itmax, 10*dim);
    localdefaults.maxiter=min(2000, localdefaults.maxiter);
    localdefaults.populationsize = min(pn, 4*dim);
    
    
    localdefaults.maxcostevals = localdefaults.maxiter*localdefaults.populationsize;
% parameters of the SA algorithm
localdefaults.maxtemp=2000;
localdefaults.mintemp=10e-10;
localdefaults.iter_temp=40;


% Merge global and local defaults, then merge w/ user options, if any.
% getGlobalDefaults(): Returns a structure with default option values for Manopt.
% mergeOptions(): Merges two options structures with one having precedence over the other.
localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

if ~isfield(problem.M, 'log') % BM ==> error if the problem contains no 'log'
    error(['The manifold problem.M must provide a logarithmic map, ' ...
        'M.log(x, y). An approximate logarithm will do too.']);
end


% Start timing for initialization
timetic = tic();


% If no initial population x is given by the user,
% generate one at random.
if ~exist('x', 'var') || isempty(x)
    x = cell(options.populationsize, 1);
    for i = 1 : options.populationsize
        x{i} = problem.M.rand(); % random agent of the population
    end
else
    if ~iscell(x)
        error('The initial guess x0 must be a cell (a population).');
    end
    if length(x) ~= options.populationsize
        options.populationsize = length(x);
        warning('manopt:pso:size', ...
            ['The option populationsize was forced to the size' ...
            ' of the given initial population x0.']);
    end
end


% Create a store database and a key for each point x{i}
storedb = StoreDB(options.storedepth);%% key for each point????
xkey = cell(size(x));
for i = 1 : numel(x)
    xkey{i} = storedb.getNewKey();
end

key = storedb.getNewKey();
key_prop = storedb.getNewKey();


% Compute cost for each particle xi,
% initialize personal best costs,
% and setup a function evaluations counter.
costs = zeros(size(x));
for i = 1 : numel(x)
    costs(i) = getCost(problem, x{i}, storedb, xkey{i});
end
costevals = options.populationsize; %cost evaluation of the agents


% Identify the best agent and store its cost/position
[fbest, imin] = min(costs);
xbest = x{imin};
xbestkey = xkey{imin}; %#ok<NASGU>
if indicator_grad==1
    fgradx = getGradient(problem, xbest, storedb, key);
    norm_grad = problem.M.norm(xbest, fgradx);
end


% Iteration counter (at any point, iter is the number of fully executed
% iterations so far)
iter = 0;


% Save stats in a struct array info, and preallocate.
% savestats will be called twice for the initial iterate (number 0),
% which is unfortunate, but not problematic.
stats = savestats();
info(1) = stats;
info(min(10000, options.maxiter+1)).iter = [];



% Start SA algorithm
max_temp=options.maxtemp;
min_temp=options.mintemp;
iter_temp=options.iter_temp;

k_final=options.maxcostevals/iter_temp;
alpha=exp((log(min_temp/max_temp)/k_final));

temp_current=max_temp;
k=0;

% Start iterating until stopping criterion triggers
while true && temp_current > min_temp
    
    stats=savestats();
    info(iter+1)=stats;
    iter=iter+1;
    
    % Make sure we don't use too much memory for the store database
    storedb.purge();
    
    % Log / display iteration information here.
    if options.verbosity >= 2
        %fprintf('Cost evals: %7d\tBest cost: %+.8e\n', costevals, fbest);
    end
    
    % Start timing this iteration
    timetic = tic();
    
    % BM: Run standard stopping criterion checks.
    % BM: Stop if any particle triggers a stopping criterion.
    [stop, reason] = stoppingcriterion(problem, xbest, options, info, iter);
    if stop
        fprintf([reason '\n']);
        break;
    end
    
    temp_current=alpha*temp_current;
    
    for i=1:iter_temp
        x_new=problem.M.rand();
        cost_x_new=getCost(problem, x_new);
        costevals = costevals + 1; 
        if cost_x_new <= fbest
            xbest=x_new;
            fbest=cost_x_new;
        elseif exp(-(fbest-cost_x_new)/temp_current) > rand(1)
            x{1}=x_new;
        end
    end
    key=key_prop;
    if indicator_grad==1
        fgradx = getGradient(problem, xbest, storedb, key);
        norm_grad = problem.M.norm(xbest, fgradx);
    end
    k=k+1;
end

info = info(1:iter);

% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = fbest;
        stats.costevals = costevals;
        if indicator_grad==1
            stats.gradnorm=norm_grad;
        end
        stats.x = x;
        stats.xbest = xbest;
        if iter == 0
            stats.time = toc(timetic);
        else
            stats.time = info(iter).time + toc(timetic);
        end
        
        % BM: Begin storing user defined stats for the entire population
        num_old_fields = size(fieldnames(stats), 1);
        trialstats = applyStatsfun(problem, x{1}, storedb, xkey{1}, options, stats);% BM
        new_fields = fieldnames(trialstats);
        num_new_fields = size(fieldnames(trialstats), 1);
        num_additional_fields =  num_new_fields - num_old_fields; % User has defined new fields
        for jj = 1 : num_additional_fields % New fields added
            tempfield = new_fields(num_old_fields + jj);
            stats.(char(tempfield)) = cell(options.populationsize, 1);
        end
        for ii = 1 : options.populationsize % Adding information for each element of the population
            tempstats = applyStatsfun(problem, x{ii}, storedb, xkey{ii}, options, stats);
            for jj = 1 : num_additional_fields
                tempfield = new_fields(num_old_fields + jj);
                tempfield_value = tempstats.(char(tempfield));
                stats.(char(tempfield)){ii} = tempfield_value;
            end
        end
        % BM: End storing
        
    end
end
