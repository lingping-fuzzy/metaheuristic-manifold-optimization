%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%François Brancaleone                              %
%Master's thesis: Global optimization on manifolds %
%Advisor: P.-A. Absil                              %
%Readers: L. Jacques and P.-Y. Gousenbourger       %
%June 2018                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbest, fbest, info, options] = ChaoticDiffevolution(problem,x, options)
% Chaotic Differential Evolution (CDE) method.
%
% function [x, cost, info, options] = ChaoticDiffevolution(problem)
% function [x, cost, info, options] = ChaoticDiffevolution(problem, x0)
% function [x, cost, info, options] = ChaoticDiffevolution(problem, x0, options)
% function [x, cost, info, options] = ChaoticDiffevolution(problem, [], options)
%
% To specify options whilst not specifying an initial guess, give x0 as []
% (the empty matrix).
%
% Apply the Chaotic Differential Evolution minimization algorithm to
% the problem defined in the problem structure, starting with the
% population x0 if it is provided (otherwise, a random population on the
% manifold is generated). The population x0 has the particularity to be a
% "chaotic" population initialized by the function ChaosPopulation.
% A population is a cell containing points on the manifold.
% The number of elements in the cell must match the parameter
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
%   K (400)
%       maximal number of iterations to construct the chaso population
%   x_max
%       extreme point on the manifold and on the domain of the objective
%       function
%   x_min
%       extreme point on the manifold and on the domain of the objective
%       function
% Original author: François Brancaleone, June 10, 2018.


% Verify that the problem description is sufficient for the solver.
if ~canGetCost(problem) %can the cost function be computed?
    warning('manopt:getCost', ...
        'No cost provided. The algorithm will likely abort.');
end

if ~canGetGradient(problem) %can the gradient function be computed?
    warning('manopt:getGradient', ...
        'No gradient provided.');
    localdefaults.indicator_grad=0; % no gradient provided
else
    localdefaults.indicator_grad=1; % the gradient is provided
end


% Dimension of the manifold
dim = problem.M.dim();


% Set local defaults here
localdefaults.storedepth = 0;      % no need for caching
localdefaults.maxcostevals = max(30000, 2*dim);
localdefaults.maxiter = max(500000, 4*dim);

localdefaults.populationsize = min(40, 10*dim);


% parameters of the DE algorithm
localdefaults.Kconst=100;
localdefaults.x_max=problem.M.rand();
localdefaults.x_min=problem.M.rand();


% Merge global and local defaults, then merge w/ user options, if any.
% getGlobalDefaults(): Returns a structure with default option values for Manopt.
% mergeOptions(): Merges two options structures with one having precedence over the other.
localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);
indicator_grad=options.indicator_grad;
if ~isfield(problem.M, 'log') % BM
    error(['The manifold problem.M must provide a logarithmic map, ' ...
        'M.log(x, y). An approximate logarithm will do too.']);
end


% Start timing for initialization
timetic = tic();


% If no initial population x is given by the user,
% generate one at random.
KK=options.Kconst;
popsize=options.populationsize;
xxmin=options.x_min;
xxmax=options.x_max;
if ~exist('x', 'var') || isempty(x)
    x = ChaosPopulation(KK, popsize, problem,xxmin ,xxmax);
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
%key_prop = storedb.getNewKey();


% Initialize trial population
u = cell(size(x));
for i = 1 : numel(x)
    u{i} = problem.M.rand();
end


% Compute cost for each particle xi,
% initialize personal best costs,
% and setup a function evaluations counter.
costs = zeros(size(x));
for i = 1 : numel(x)
    costs(i) = getCost(problem, x{i}, storedb, xkey{i});
end
costevals = options.populationsize;


% Identify the best agent and store its cost/position
[fbest, imin] = min(costs);
xbest = x{imin};%xbbest;
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

K=options.Kconst;
x_min=options.x_min;
x_max=options.x_max;


% start the CDE algorithm

while true
    stats = savestats();
    info(iter+1) = stats;
    iter = iter + 1;
    
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
    for i = numel(x)
        [stop, reason] = stoppingcriterion(problem, x{i}, options, info, iter);
        if stop
            break;
        end
    end
    
    if stop
        if options.verbosity >= 1
            fprintf([reason '\n']);
        end
        break;
    end
    
    
    % call DE
    optionsDE.maxcostevals=1000;
    optionsDE.schema=7;
    [x_new, f_new, ~, ~] = diffevolution(problem, x, optionsDE);
    costevals=costevals+optionsDE.maxcostevals;
    if indicator_grad==1 %call the gradient method if gradient available
        optionsgrad.maxiter=400;
        optionsgrad.tolgradnorm=1e-9;
        [x_new, f_new, infograd, ~]=trustregions(problem, x_new, optionsgrad);
        costevals=costevals+length([infograd.iter]);
        
        fgradx = getGradient(problem, x_new, storedb, key);
        norm_grad = problem.M.norm(x_new, fgradx);
    end
    
    if f_new < fbest
        xbest=x_new;
        fbest=f_new;
    end
    
    if canGetGradient(problem)
    fgradx = getGradient(problem, x_new, storedb, key);
    norm_grad = problem.M.norm(x_new, fgradx);
    end
    
    
    %chaos population of size populationsize/2
    pophalf=options.populationsize/2;
    x_sub = ChaosPopulation(K, pophalf,problem, x_min,x_max);
    new_costs2=zeros(1,options.populationsize/2);
    for i=1:options.populationsize/2
        new_costs2(i) = getCost(problem, x_sub{i});
    end
    
    costevals=costevals+options.populationsize/2;
    
    % create a new population that is a mix of best agents of the initial population
    % and the new chaos population
    [new_costs1,I] = sort(costs);
    
    for i=1:options.populationsize/2
        x{i}=x{I(i)};
        costs(i)=new_costs1(i);
        
        x{options.populationsize/2+i}=x_sub{i};
        costs(options.populationsize/2+i)=new_costs2(i);
    end
    
    % replace the worst agent with xbest
    [~, index]=max(costs);
    x{index}=xbest;
    costs(index)=fbest;

end

info = info(1:iter);
fbest=getCost(problem, xbest);

% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = fbest;
        stats.costevals = costevals;
        if indicator_grad==1
            stats.gradnorm=norm_grad;
        end
        stats.x = x;
        %stats.v = v;
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










