%% Example file to run experiments using mvmv algorithm
% Look at the demo_mvmc.m file to run a single experiment
% ------------------------------------------------------------------------
% A multiobjective evolutionary algorithm for multiview clustering (MVCM).
% This algorithm using MOEA/D as the underlaying optimizer
% This toolbox was developed with MATLAB R2014a.
%
% Developed by
%   A José-García (adanjoga@gmail.com, adan.jose@cinvestav.mx)
%
% Please, cite the following paper where this algorithm was introduced:
%
%   A. José-García, J. Handl, W. Gómez-Flores, and M. Garza-Fabre.
%   "An Evolutionary Many-objective Approach to Multiview Clustering 
%   Using Feature and Relational Data". Applied Soft Computing (in press)
%
% NOTE: Modify the number number of experiments (N_EXPERIMENTS=21) and the
% number of generations (MVMCparams.MAXGEN=100) as needed. The results 
% obtained are saved in the results folder.
%------------------------------------------------------------------------
clear all; close all; clc;

addpath([pwd '/mvmc']); 
addpath([pwd '/datasets/synthetic']);
addpath([pwd '/datasets/real']);
addpath([pwd '/objectives']);
addpath([pwd '/model_selection']);
addpath([pwd '/utils']);

%% Preparation of input variables

% List of datasets provided
DSnames = {'Data_4_3','Sizes5','Spirals2','matlab_iris','matlab_wine'};
N_EXPERIMENTS = 21;
N_DATASETS = numel(DSnames);

% Number of objectives (or data views)
Nobjs  = 2;
% ------------------------------------------------------------------------
% Input precomputed weights
Weights   = load('weights.mat');
MVMCparams.W = Weights.W2D100;
if Nobjs == 2
    MOEADat.W = Weights.W2D100;
elseif Nobjs == 3
    MOEADat.W = Weights.W3D150;
elseif Nobjs == 4
    MOEADat.W = Weights.W4D165;    
elseif Nobjs == 5
    Dat.W = Weights.W5D210;
end
% ------------------------------------------------------------------------
% Variables regarding the multiobjective algorithm
MVMCparams.MAXGEN      = 100;                      % Number of generations
MVMCparams.Popsize     = size(MVMCparams.W,1);     % Population size
MVMCparams.Niche       = 10;                       % The neighboursize T = 20
MVMCparams.Method      = 'tch';                    % The decomposition method {tch}
MVMCparams.pC          = 0.5;                      % The probability of crossing over
MVMCparams.pM          = 0.03;                     % The mutation probability
MVMCparams.NOBJ        = Nobjs;

% ------------------------------------------------------------------------
for iDS = 1:N_DATASETS

% Variables for each data problem (dataset)    
ARIbest = NaN(1,N_EXPERIMENTS);  % best ARI values
CLRs    = cell(1,N_EXPERIMENTS); % Clustering solutions
PFAs   = cell(1,N_EXPERIMENTS);  % Pareto front approximations
    
disp(['Computing dataset (' DSnames{iDS} ')']);

% Load the input dataset
Data    = load(DSnames{iDS});
X       = Data.data(:,1:end-1);
Labels  = Data.data(:,end);
K       = numel(unique(Labels));

% Input multiple data views in the form of dissimilarity matrices
DXeuc  = pdist(X,'euclidean');  DXeuc = squareform(DXeuc);
DXcos  = pdist(X,'cosine');     DXcos = squareform(DXcos);
Dataviews = {DXeuc, DXcos};

for iExp = 1:N_EXPERIMENTS
    disp([' Computing dataset ' DSnames{iDS} '.' num2str(iExp)]);
    
    rng(iExp)
    % Run the mvmc algorithm
    OUT = mvmc(Dataviews, K, MVMCparams);
    % OUT.nPFront   returns the normalized Pareto fromt
    % OUT.PClrs     returns the consensus clustering solutions in the PF
    % OUT.W         returns the weights associated to each solution

    % Supervised Model Selection using ARI
    [Clr,ARIb,idx1,ARIvalues] = cmsARI(OUT.PClrs, Labels);
    
    % Save data for each experiment
    ARIbest(iExp) = ARIb;
    CLRs(iExp)    = {Clr};
    PFAs(iExp)    = {OUT.nPFront};
end

save(['results/' DSnames{iDS}],'ARIbest','CLRs','PFAs');
end