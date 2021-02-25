%% MVMC
% A multiobjective evolutionary algorithm for multiview clustering
% 
% Some general characteristics of this mvmc version:
%   + normalization of the objectives as fn=f-zref/znad-zref
%   + medoid-based representation
%   + uniform crossover (Cr=0.05) and 
%   + uniform mutation (Cm=0.03)
%
% Qingfu Zhang and Hui Li, 2007. MOEA/D: A Multiobjective Evolutionary 
% based on Decomposition optimization. IEEE Transactions on Evolutionary Computation.
% Vol 11, No. 6, 712-731.
%
function OUT = mvmc(Dataviews, K, MVMCparams, Display)

global Views_norm idealpoint nadirpoint;

MAXGEN  = MVMCparams.MAXGEN;    % Maximum number of generations.
NP      = MVMCparams.Popsize;   % Population size.
Nobj    = MVMCparams.NOBJ;      % Number of objectives.
Niche   = MVMCparams.Niche;     % The neighboursize  
Dmethod = MVMCparams.Method;    % The decomposition method, 'ws' or 'ts'.

T       = MVMCparams.Labels;    % TEMP

objective = 'wgs';
[cvifun,opt] = cviconfig(objective);

Views_norm = cell(1,Nobj);
for m = 1:Nobj
    Views_norm{m} = minmaxnorm(Dataviews{m});
end
N = size(Dataviews{1},1);

pC = 0.5;   % The probability of crossing over
pM = 0.03;  % The mutation probability
%% Compute the Nadir Point

MN = zeros(1,Nobj);
nadirpoint = zeros(Nobj,1);
for m = 1:Nobj
    nadirpoint(m) = min(sum(Dataviews{m},2),[],1);
end

%% Initial random population
idealpoint = ones(Nobj,1)*inf;

Subproplems = subproblems(NP, Niche, MVMCparams.W);
NP = length(Subproplems);

individual = struct('var',[], 'obj',[], 'nobj',[], 'w',[], 'clrs',[]);
individuals = repmat(individual,NP,1);

ARIstat = NaN(NP,2); ARItemp = NaN(NP,1);
%% Evalution of the initial population
for i=1:NP
    % Generation of medoids Pi in {1,...,N}
    Pi = datasample(1:N,K,'Replace',false); % NOTE: This function returns unique medoids

    % Assignation of each observation to one cluster
    Wi = (Subproplems(i).weight)';
    Clr = assignmentX(Pi,Wi,Nobj);
    
    % Evaluation of a sulution using different dissimilarity measures
    for m = 1:Nobj
        fit(m) = feval(cvifun,Clr,K,Dataviews{m},nadirpoint(m)); 
    end
    
    nfit = minmaxnorm(fit,[MN;nadirpoint']);
    individuals(i).obj = fit'; individuals(i).nobj = nfit';
    individuals(i).var = Pi';
    individuals(i).w = Wi';
    individuals(i).clrs = Clr;

    Subproplems(i).individual = individuals(i);
    ARItemp(i)= pairwiseindex(T,Clr);
end
Fit = [individuals.nobj];
if strcmpi(opt,'mn')
    idealpoint = min(idealpoint, min(Fit,[],2));
end

if Display
    [~,idx]=max(ARItemp);
    disp(['Iteration 0 | iARIavg=' num2str(mean(ARItemp)) ' iARImax=' num2str(max(ARItemp)) ' ID=' num2str(idx)]);
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
end
%% Evolution process
for g=1:MAXGEN
    
    masksc = rand(NP,K) < pC;
    masksm = rand(NP,K) < pM;
    for i=1:NP
        % Selection of parents from its neighborhood
        nindex = Subproplems(i).neighborhood;
        si = nindex(randperm(length(nindex)));
        Pj = [Subproplems(si(1:2)).individual];
        
        % Position-based crossover:
        % Gilbert Syswerda. Schedule Optimization using Genetic Algorithms. 
        % In Lawrence Davis, editor, Handbook of Genetic Algorithms, chapter 21, pages 332?349.
        medidx = NaN(1,K); parent1 = Pj(1).var';  parent2 = Pj(2).var'; % medidx is a new vector (child1)
        medidx(masksc(i,:)) = parent1(masksc(i,:));  % generate an offspring from P1
        
        % Remove from Parent2 the elements in medidx (Parent1)
        ids = parent1(masksc(i,:));
        for j=1:numel(ids) 
            parent2(parent2==ids(j))=[];
        end
        % Insert the elements of parent2 into medidx 
        medidx(~masksc(i,:)) = parent2(1:sum(~masksc(i,:)));
        if isnan(medidx); error('ERROR: vector having NaN values'); end

        % Uniform mutation
        if sum(masksm(i,:))
            %Pt = datasample(1:N,K,'Replace',false);
            Pt = setdiff(randperm(N),medidx,'stable');
            medidx(masksm(i,:)) = Pt(masksm(i,:));
        end
        % NOTE: the uniform crossover and mutation operator may produce duplicated
        % medoids. Thus, the "repair method" (validation of unique medoids)
        % is performed in the assigmentX() function (see pp. 718 of MOEA/D)
        if numel(unique(medidx)) ~= K
            medidx = repair_method(medidx,N,K);
        end

        % Assignation of each observation to one cluster
        Wi = Subproplems(i).weight';
        Clr = assignmentX(medidx,Wi,Nobj);

        % Evaluation of the sulution using the different data views
        for m = 1:Nobj
            fit(m) = feval(cvifun,Clr,K,Dataviews{m},nadirpoint(m)); 
        end
        
        nfit = minmaxnorm(fit,[MN;nadirpoint']);
        ind = struct('var', medidx','obj', fit','nobj', nfit', 'w', Wi', 'clrs', Clr);

        if strcmpi(opt,'mn')
            idealpoint = min(idealpoint, ind.nobj); 
        end
        % Update the neighbours        
        nindex = Subproplems(i).neighborhood;
        neighbours = Subproplems(nindex);
        newobj = subobjective([neighbours.weight], ind.nobj, idealpoint, Dmethod);
        oops = [neighbours.individual]; 
        oldobj = subobjective([neighbours.weight], [oops.nobj], idealpoint, Dmethod);

        C = newobj < oldobj;
        [neighbours(C).individual] = deal(ind);
        Subproplems(nindex) = neighbours;
        
        ARItemp(i)= pairwiseindex(T,Clr);
    end
    Pop = [Subproplems.individual];
    nPFront = [Pop.nobj];
    if Display
        ARIstat(g,1) = mean(ARItemp); ARIstat(g,2) = max(ARItemp);
        [~,idx]=max(ARItemp);
        
        disp(['Iteration ' num2str(g) '| iARIavg=' num2str(ARIstat(g,1)) ' iARImax=' num2str(ARIstat(g,2)) ' ID=' num2str(idx)]);
        disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        PrinterDisplay(nPFront',Nobj); % To print results on screen
    end
end

%% Report solutions
PSet    = [Pop.var];
PFront  = [Pop.obj];
nPFront = [Pop.nobj];
PClrs   = [Pop.clrs];
W       = [Pop.w];


Idx = DominanceFilter(nPFront',opt);
%Idx = ones(NP,1);
OUT.PSet    = PSet(:,Idx')';
OUT.nPFront = nPFront(:,Idx')'; % Normalized Pareto front approximation
OUT.PFront  = PFront(:,Idx')';   % Pareto front approximation
OUT.PClrs   = PClrs(:,Idx');
OUT.W       = W(:,Idx')';
OUT.ARIstat = ARIstat;

if Display
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    disp('Red  asterisks : Set Calculated.')
    disp('Black circles  : Pareto Front.')
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

    F=OUT.nPFront;
    for i=1:size(F,1)

        figure(123); hold on;
        if Nobj == 2
            plot(F(i,1),F(i,2),'ok','MarkerFaceColor','k');
        elseif Nobj == 3    
            plot3(F(i,1),F(i,2),F(i,3),'ok','MarkerFaceColor','k');
        end
        grid on; hold on;
    end
end
end
%% Print and Display information
function PrinterDisplay(PFront, Nobj)
figure(123);
%hold on;
if Nobj == 2
    plot(PFront(:,1),PFront(:,2),'*r'); grid on; 
    xlabel('View 1'); ylabel('View 2');
elseif Nobj == 3
    plot3(PFront(:,1),PFront(:,2), PFront(:,3),'*r'); grid on;
    xlabel('View 1'); ylabel('View 2'); zlabel('View 3');
end
end
%% Selection of unique medoids
% ---------------------------------------------------------------------
function Midx = repair_method(Midx,N,K)

[meds,p] = unique(Midx,'stable');
rmedoids = setdiff(randperm(N),meds,'stable');
rmedoids = rmedoids(1:K-numel(p));% true 

unimeds = zeros(1,numel(Midx)); unimeds(p) = 1; 
Midx(~unimeds) = rmedoids;

end
%% Decodification of a clustering solution
% ---------------------------------------------------------------------
function Clr = assignmentX(Midx,Wi,Nobj)
global Views_norm;

% Computation of the new weighted distance matrix
Dws = (Wi(1) .* Views_norm{1});
for m = 2:Nobj
    Dws = Dws + (Wi(m) .* Views_norm{m});
end

% Assignation of points to its clossest medoid
DXdist = Dws(:,Midx);
[~,Clr] = min(DXdist,[],2);  
end
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Dominance Filter
function Indx = DominanceFilter(F,opt)

N = size(F,1);
Indx = false(N,1);
for xpop = 1:N
    Dominado=0;
    for compara = 1:N
        if F(xpop,:)==F(compara,:)
            if xpop > compara
                Dominado=1;
                break;
            end
        else
            if strcmpi(opt,'mn')
                if F(xpop,:) >= F(compara,:)
                    Dominado=1;
                    break;
                end
            elseif strcmpi(opt,'mx')
                if F(xpop,:) <= F(compara,:)
                    Dominado=1;
                    break;
                end
            end
        end
    end
    
    if Dominado == 0
        Indx(xpop) = true;
    end
end
end

function Indx = DominanceFilter2(F,opt)

N = size(F,1);
Indx = false(N,1);
for xpop = 1:N
    Dominado=0;
    for compara = 1:N

            if strcmpi(opt,'mn')
                if F(xpop,:) > F(compara,:)
                    Dominado=1;
                    break;
                end
            elseif strcmpi(opt,'mx')
                if F(xpop,:) <= F(compara,:)
                    Dominado=1;
                    break;
                end
            end
        %end
    end
    
    if Dominado == 0
        Indx(xpop) = true;
    end
end
end
