function Offspring = OperatorGA(Problem,Parent,Parameter)
%OperatorGA - Crossover and mutation operators of genetic algorithm.
%
%   Off = OperatorGA(Pro,P) uses genetic operators to generate offsprings
%   for problem Pro based on parents P. If P is an array of SOLUTION
%   objects, then Off is also an array of SOLUTION objects. While if P is a
%   matrix of decision variables, then Off is also a matrix of decision
%   variables, i.e., the offsprings are not evaluated. P is split into two
%   subsets P1 and P2 with the same size, where each object or row of P1
%   and P2 is used to generate two offsprings. Different operators are used
%   for different encoding schemes.
%
%   Off = OperatorGA(Pro,P,{proC,disC,proM,disM}) specifies the parameters
%   of operators, where proC is the probability of crossover, disC is the
%   distribution index of simulated binary crossover, proM is the
%   expectation of the number of mutated variables, and disM is the
%   distribution index of polynomial mutation.
%
%   Example:
%       Offspring = OperatorGA(Problem,Parent)
%       Offspring = OperatorGA(Problem,Parent.decs,{1,20,1,20})
%
%   See also OperatorGAhalf

%------------------------------- Reference --------------------------------
% [1] K. Deb, K. Sindhya, and T. Okabe, Self-adaptive simulated binary
% crossover for real-parameter optimization, Proceedings of the Annual
% Conference on Genetic and Evolutionary Computation, 2007, 1187-1194.
% [2] K. Deb and M. Goyal, A combined genetic adaptive search (GeneAS) for
% engineering design, Computer Science and informatics, 1996, 26: 30-45.
% [3] L. Davis, Applying adaptive algorithms to epistatic domains,
% Proceedings of the International Joint Conference on Artificial
% Intelligence, 1985, 162-164.
% [4] D. B. Fogel, An evolutionary approach to the traveling salesman
% problem, Biological Cybernetics, 1988, 60(2): 139-144.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if nargin > 2
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
    end
      
    Parent1   = Parent(1:floor(end/2),:);
    Parent2   = Parent(floor(end/2)+1:floor(end/2)*2,:);
    
    if ~isstruct(Parent{1}) % for there is named struct
        
        if size(Parent{1}, 1) ~= 1  % for solution is not a vector
            Offspring = cell(size(Parent, 1), 1);
            for pop = 1 : size(Parent1, 1)
                id = 2*pop-1;
                popOne = Parent1(pop);
                popTwo = Parent2(pop);
                input1 = cell2mat(popOne);
                input2 = cell2mat(popTwo);
                Offs = GArealstruct(input1, input2,...
                    proC,disC,proM,disM);
                n = size(Offs, 1);
                Offs = reshape(Offs, n,[]);
                Offspring{id,1} = reshape(Offs(1:n/2,:), size(input1));
                Offspring{id+1,1} = reshape(Offs(n/2+1:end,:), size(input1));
            end    
        else
            lower = zeros(1, Problem.M.dim()+1) - 1;
            upper = zeros(1, Problem.M.dim()+1) + 1;
            Offspring = GAreal(Parent1,Parent2,lower, upper, proC,disC,proM,disM);
        end
        
    else
        elems = fieldnames(Parent{1});
        nelems = numel(elems);
        
        assert(nelems >= 1, ...
            'elements must be a structure with at least one field.');
        
        Offspring = cell(size(Parent, 1), 1);
        for pop = 1 : size(Parent1, 1)
            id = 2*pop-1;
            for i = 1 : nelems
                popOne = Parent1(pop);
                popTwo = Parent2(pop);
                input1 = cell2mat(popOne).(elems{i});
                input2 = cell2mat(popTwo).(elems{i});
                Offs = GArealstruct(input1, input2,...
                    proC,disC,proM,disM);
                n = size(Offs, 1);
                Offs = reshape(Offs, n,[]);
                Offspring{id,1}.(elems{i}) = reshape(Offs(1:n/2,:), size(input1));
                Offspring{id+1,1}.(elems{i}) = reshape(Offs(n/2+1:end,:), size(input1));
            end
        end
    end
end


function data = cell2dataArray(celldata)
    data = zeros(numel(celldata), length(celldata{1}));
    for i = 1:numel(celldata)
       data(i,:) = cell2mat(celldata(i));
    end
end


function Offspring = GAreal(Parent1,Parent2,lower,upper,proC,disC,proM,disM)
% Genetic operators for real and integer variables

    %% Simulated binary crossover
    Parent1= cell2dataArray(Parent1);
    Parent2= cell2dataArray(Parent2);
    [N,D] = size(Parent1);
    beta  = zeros(N,D);
    mu    = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
%     Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
%                  (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    Offspring = [beta.*(Parent1-Parent2)/2
                 -beta.*(Parent1-Parent2)/2];             
    %% Polynomial mutation
    Lower = repmat(lower,2*N,1);
    Upper = repmat(upper,2*N,1);
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end


function Offspring = GArealstruct(Parent1,Parent2,proC,disC,proM,disM)
% Genetic operators for real and integer variables

    %% Simulated binary crossover
    sval = size(Parent1);
    if length(sval) == 2
      [N, D] =feval(@(x) x{:}, num2cell(sval));
       beta  = zeros(N,D);
       mu    = rand(N,D);
       wdim = 2;
    else
       [N, D, M] =feval(@(x) x{:}, num2cell(sval));
        beta  = zeros(N, D, M);
        mu    = rand(N,D, M);
        wdim = 3;
    end

    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    
    lower = zeros(1, D) - 1;
    upper = zeros(1, D) + 1;
    
    if  wdim == 2
        beta = beta.*(-1).^randi([0,1],N,D);
        beta(rand(N,D)<0.5) = 1;
        beta(repmat(rand(N,1)>proC,1,D)) = 1;
        Lower = repmat(lower,2*N,1);
        Upper = repmat(upper,2*N,1);
        Site  = rand(2*N,D) < proM/D;
        mu    = rand(2*N,D);
    elseif wdim == 3
        beta = beta.*(-1).^randi([0,1],N,D, M);
        beta(rand(N,D, M)<0.5) = 1;
        beta(repmat(rand(N,1, M)>proC,1,D)) = 1;
        Lower = repmat(lower,2*N,1, M);
        Upper = repmat(upper,2*N,1, M);
        Site  = rand(2*N,D, M) < proM/D;
        mu    = rand(2*N,D, M);
    end
  
%     Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
%                  (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    Offspring = [beta.*(Parent1-Parent2)/2
                 -beta.*(Parent1-Parent2)/2];             
    %% Polynomial mutation

    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end

