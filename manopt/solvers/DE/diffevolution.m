%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%François Brancaleone                              %
%Master's thesis: Global optimization on manifolds %
%Advisor: P.-A. Absil                              %
%Readers: L. Jacques and P.-Y. Gousenbourger       %
%June 2018                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbest, fbest, info, options] = diffevolution(problem, pn, itmax, x, options)%,xbbest, fbbest)
% Differential Evolution (DE) method for derivative-free minimization.
%
% function [x, cost, info, options] = diffevolution(problem)
% function [x, cost, info, options] = diffevolution(problem, x0)
% function [x, cost, info, options] = diffevolution(problem, x0, options)
% function [x, cost, info, options] = diffevolution(problem, [], options)
%
% To specify options whilst not specifying an initial guess, give x0 as []
% (the empty matrix).
%
% Apply the Differential Evolution minimization algorithm to
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
%   schema (1,2,3,4,5,6 ==> by default 2)
%       Scheme for the mutation step (parameter of the DE method).
%   CR (0.6)
%       Crossover rate (parameter of the DE method).
%   lambda (0.95)
%       (parameter of the DE method).
%
% Original author: François Brancaleone, June 10, 2018.


% Verify that the problem description is sufficient for the solver.
if ~canGetCost(problem) %can the cost function be computed?
    warning('manopt:getCost', ...
        'No cost provided. The algorithm will likely abort.');
end

if ~canGetGradient(problem) %can the gradient function be computed?
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
localdefaults.maxcostevals = max(50000, 2*dim);
    localdefaults.maxiter = max(itmax, 10*dim);
    localdefaults.maxiter=min(800, localdefaults.maxiter);
    localdefaults.populationsize = min(pn, 4*dim);

% parameters of the DE algorithm
localdefaults.schema=2;
localdefaults.CR=0.6;
localdefaults.lambda=0.95;


% Merge global and local defaults, then merge w/ user options, if any.
% getGlobalDefaults(): Returns a structure with default option values for Manopt.
% mergeOptions(): Merges two options structures with one having precedence over the other.
localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

if ~isfield(problem.M, 'log') % BM
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
        x{i} = problem.M.rand();
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
x = load('specturmFigdraw.x.mat').x;
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



% Start DE algorithm
CR=options.CR;
lambda=options.lambda;
schema=options.schema;

% Start iterating until stopping criterion triggers
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
    
    for i=1:numel(x) %start looping on the agents of the population
        xx=x{i}; % agent i
        uu=u{i};
        r1=randi(numel(x));
        r2=randi(numel(x));
        r3=randi(numel(x));
        r4=randi(numel(x));
        r5=randi(numel(x));
        
        % mutation
        if options.schema==7
            schema=randi(6);
        end
        
        if schema==1 %DE/rand/1
            vv=problem.M.exp(x{r1}, problem.M.transp(x{r2},x{r1}, problem.M.lincomb(x{r2}, lambda, problem.M.log(x{r2},x{r3}))));
            
        elseif schema==2 %DE/best/1
            vv=problem.M.exp(xbest, problem.M.transp(x{r2},xbest, problem.M.lincomb(x{r2}, lambda, problem.M.log(x{r2},x{r3}))));
            
        elseif schema==3 %DE/rand to best/1
            vv1=problem.M.exp(xbest, problem.M.transp(x{r2},xbest, problem.M.lincomb(x{r2}, lambda, problem.M.log(x{r2},x{r3}))));
            vv=problem.M.exp(vv1, problem.M.transp(xbest,vv1, problem.M.lincomb(xbest, lambda, problem.M.log(xbest,x{r1}))));
            
        elseif schema==4 %DE/current to best/1
            vv1=problem.M.exp(x{i}, problem.M.transp(x{r2},x{i}, problem.M.lincomb(x{r2}, lambda, problem.M.log(x{r2},x{r3}))));
            vv=problem.M.exp(vv1, problem.M.transp(xbest,vv1, problem.M.lincomb(xbest, lambda, problem.M.log(xbest,x{i}))));
            
        elseif schema==5 %DE/ current to rand/1
            vv1=problem.M.exp(x{i}, problem.M.transp(x{r2},x{i}, problem.M.lincomb(x{r2}, lambda, problem.M.log(x{r2},x{r3}))));
            vv=problem.M.exp(vv1, problem.M.transp(x{1},vv1, problem.M.lincomb(x{r1}, rand(1), problem.M.log(x{r1},x{i}))));
            
        elseif schema==6 %DE/rand/2
            vv1=problem.M.exp(x{r1}, problem.M.transp(x{r2},x{r1}, problem.M.lincomb(x{r2}, lambda, problem.M.log(x{r2},x{r3}))));
            vv=problem.M.exp(vv1, problem.M.transp(x{r4},vv1, problem.M.lincomb(x{r4}, lambda, problem.M.log(x{r4},x{r5}))));
        end
        
        if isstruct(x{i})
            elems = fieldnames(x{i});
            nelems = numel(elems);
        else
            nn=length(x{i}); %dimension
            [nn, ~] =size(x{i}); %dimension
        end
        
       %%%%% recombination
        if isstruct(x{i})
            Irand=randi(nelems);
            for j=1:nelems
                randd=rand(1);
                if randd <= CR || j==Irand
                    uu.(elems{j})=vv.(elems{j});
                else
                    uu.(elems{j})=x{i}.(elems{j});
                end
            end
        elseif iscell(x{i})
            Irand=randi(nn);
            for j=1:nn
                randd=rand(1);
                if randd <= CR || j==Irand
                    uu{j}=vv{j};
                else
                    uu{j}=xx{j};
                end
            end
        elseif ismatrix(x{i}) % ca rentre ici aussi quand c'est un vecteur..
            Irand=randi(nn);
            for j=1:nn
                randd=rand(1);
                if randd <= CR || j==Irand
                    uu(j,:)=vv(j,:);
                else
                    uu(j,:)=xx(j,:);
                end
            end
            if nn==2
            uu=uu/norm(uu);
            end
         end
        
        u{i}=uu;
        
        % selection
        costuu=getCost(problem, u{i});
        costevals=costevals+1;
        
        if  costuu < costs(i)
            x{i}=u{i};
            costs(i)=costuu;
            if costuu < fbest
                xbest=x{i};
                fbest=costuu;
            end
        end
        if indicator_grad==1
            fgradx = getGradient(problem, xbest, storedb, key);
            norm_grad = problem.M.norm(xbest, fgradx);
        end
    end
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
        
        
        