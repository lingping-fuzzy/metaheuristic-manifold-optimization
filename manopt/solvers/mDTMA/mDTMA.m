function [xbest, fbest, info, options] = mDTMA(problem, pn, itmax, w0, x, options)
% Particle swarm optimization (PSO) for derivative-free minimization.
%
% function [x, cost, info, options] = mDTMA(problem)
% function [x, cost, info, options] = mDTMA(problem, x0)
% function [x, cost, info, options] = mDTMA(problem, x0, options)
% function [x, cost, info, options] = mDTMA(problem, [], options)
%
% Apply the genetic minimization algorithm to
% the problem defined in the problem structure, starting with the
% population x0 if it is provided (otherwise, a random population on the
% manifold is generated). A population is a cell containing points on the
% manifold. The number of elements in the cell must match the parameter
% options.populationsize.
%
% To specify options whilst not specifying an initial guess, give x0 as []
% (the empty matrix).
%
% None of the options are mandatory. See in code for details.
%
% Based on the original GA, PSO description in
% https://mse.redwoods.edu/darnold/math45/Activities/ProjectionMatrix/proj.pdf

% This file USE part of Manopt: www.manopt.org.
% Original author: Lingping KONG, April. 09, 2023.
% Contributors: Aug. 1, 2014.
% Change log:
%
   
    
    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem)
        warning('manopt:getCost', ...
            'No cost provided. The algorithm will likely abort.');
    end
    
    % Dimension of the manifold
    dim = problem.M.dim();
    
    % Set local defaults here
    localdefaults.storedepth = 0;                   % no need for caching

    localdefaults.maxiter = max(itmax, 10*dim);
    localdefaults.maxiter=min(800, localdefaults.maxiter);
    localdefaults.populationsize = min(pn, 4*dim);
    
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
   
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
    
    % Create a store database and a key for each point x{i}
    storedb = StoreDB(options.storedepth);
    xkey = cell(size(x));
    for i = 1 : numel(x)
        xkey{i} = storedb.getNewKey();
    end


    % Initialize personal best positions to the initial population
    y = cell(options.populationsize, 1);
    ykey = cell(size(y));
    
    
    % Initialize velocities for each particle 
    v = cell(size(x));
    
    % Compute cost for each particle xi,
    % initialize personal best costs,
    % and setup a function evaluations counter.
    costs = zeros(size(x));
    for i = 1 : numel(x)
        costs(i) = getCost(problem, x{i}, storedb, xkey{i});
    end
    costevals = options.populationsize;
    
    % Identify the best particle and store its cost/position
    [fbest, imin] = min(costs);
    xbest = x{imin};
    xbestkey = xkey{imin}; %#ok<NASGU>
    
    % Iteration counter (at any point, iter is the number of fully executed
    % iterations so far)
    iter = 0;
    
    % Save stats in a struct array info, and preallocate.
    % savestats will be called twice for the initial iterate (number 0),
    % which is unfortunate, but not problematic.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % Start iterating until stopping criterion triggers
    while true
        
        stats = savestats();
        info(iter+1) = stats; %#ok<AGROW>
        iter = iter + 1;
        
        % Make sure we don't use too much memory for the store database
        storedb.purge();
        
        % Log / display iteration information here.
        if options.verbosity >= 2
            fprintf('Cost evals: %7d\tBest cost: %+.8e\n', costevals, fbest);
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
        
        
        % Compute the inertia factor

        w = w0 + 0.1*(1-iter/options.maxiter); % KONG-
        % tournament selection
        MatingPool = TournamentSelection(2,numel(x),costs);
        [proC,disC,proM,disM] = deal(1,20,1,20);
        % storeb and its key-> attention
        velocityspring  = OperatorGA(problem, x(MatingPool), {proC,disC,proM,disM}); 
        
        if ~iscell(velocityspring)
            % just for double array, not struct
            for i = 1 : numel(x)
                xi1 = x{MatingPool(i)};
                vi1 = velocityspring(i,:)';
                vtang = (vi1.'*xi1/(xi1.'*xi1)).*xi1;
                vtemp = vtang - xi1;
                firstgo = problem.M.lincomb(xi1, w , vtemp);
                x{MatingPool(i)} = problem.M.retr(x{MatingPool(i)}, firstgo);
                
                vi1 = vi1-vtang;   %get the tangential (project vector).  
                v{MatingPool(i)} = problem.M.lincomb(x{MatingPool(i)}, w , vi1); % some point on tangent plane away from xi1. 

            end
        elseif size(velocityspring{1}, 1)~= 1
           for i = 1 : numel(x)
               xi1 = x{MatingPool(i)};
               vi1 = velocityspring{i};
               osize = size(vi1);
               sval = reshape(vi1, 1, []); dval = reshape(xi1, 1, []);
               % here we transfer a matrix to vector than do the
               % projection, --try use matrix project to matrix. may produce better result---###
                vtang = (sval*dval.'/(dval*dval.')).*dval;
                vtemp = reshape(vtang, osize) - xi1;
                firstgo = problem.M.lincomb(xi1, w , vtemp);
                x{MatingPool(i)} = problem.M.retr(xi1, firstgo);
                
                temp = sval-vtang;   %get the tangential (project vector).          
                vi1 = reshape(temp, osize);               
               v{MatingPool(i)} = problem.M.lincomb(x{MatingPool(i)}, w , vi1); % some point on tangent plane away from xi1.
            end
        else
             elems = fieldnames(x{1});
             nelems = numel(elems);
             
            for i = 1 : numel(x)
                xi1 = x{MatingPool(i)};
                vi1 = velocityspring{i};

                for idn = 1 : nelems
                    sval = vi1.(elems{idn}); dval = xi1.(elems{idn});
                    osize = size(sval);
                    sval = reshape(sval, 1, []); dval = reshape(dval, 1, []);
                    vtemp = (sval*dval.'/(dval*dval.'+1e-10)).*dval;
                    vtang.(elems{idn}) = reshape(vtemp, osize) - xi1.(elems{idn});
                end      
                firstgo = problem.M.lincomb(xi1, w , vtang); % some point on tangent plane away from xi1. 
                x{MatingPool(i)} = problem.M.retr(xi1, firstgo);
                
                for idn = 1 : nelems
                    sval = vi1.(elems{idn}); dval = xi1.(elems{idn});
                    osize = size(sval);
                    sval = reshape(sval, 1, []); dval = reshape(dval, 1, []);
                    temp = sval-(sval*dval.'/(dval*dval.'+1e-10)).*dval;
                    vi1.(elems{idn}) = reshape(temp, osize);
                end      
                v{MatingPool(i)} = problem.M.lincomb(x{MatingPool(i)}, w , vi1); % some point on tangent plane away from xi1. 
                
            end
        end
        
       ycosts = zeros(size(velocityspring, 1), 1);
       for i = 1 : numel(x)
            % compute new position of particle i
            id = MatingPool(i);
            y{i} = problem.M.retr(x{id}, v{id});
            ykey{i} = storedb.getNewKey();
            ycosts(i) = getCost(problem, y{i}, storedb, ykey{i});
            costevals = costevals + 1;    
       end
       

        tempx = [x; y];
        tempcost = [costs; ycosts];
        [~,idcost] = sort(tempcost);
        tempx = tempx(idcost);
        x = tempx(1: numel(x));
        tempcost =tempcost(idcost);
        costs = tempcost(1: numel(x));
        keys = [xkey;ykey];  keys = keys(idcost);
        xkey = keys(1: numel(x));

            % update self-best if necessary
        if costs(1) < fbest
             fbest = costs(1);
             xbest = x{1};
             xbestkey = xkey{1}; %#ok<NASGU>
        end    
        
    end
    
    
    info = info(1:iter);
     
    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = fbest;
        stats.costevals = costevals;
        stats.x = x;
%         stats.v = v;
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
