%% DEMO: Test the mvmv algorithm on different feature and relational data 
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
%------------------------------------------------------------------------
clear all; close all; clc;

addpath([pwd '/mvmc']); 
addpath([pwd '/datasets/synthetic']);
addpath([pwd '/datasets/real']);
addpath([pwd '/objectives']);
addpath([pwd '/model_selection']);
addpath([pwd '/utils']);

%% Preparation of input variables

% Variables regarding the clustering problem
% List of datasets provided
iDS     = 1;
DSnames = {'Data_4_3','Sizes5','Spirals','matlab_iris','matlab_wine'};

% Load the input dataset
Data    = load(DSnames{iDS});
X       = Data.data(:,1:end-1);
Labels  = Data.data(:,end);
K       = numel(unique(Labels));

% Input multiple data views in the form of dissimilarity matrices
DXeuc  = pdist(X,'euclidean');  DXeuc = squareform(DXeuc);
DXcos  = pdist(X,'cosine');     DXcos = squareform(DXcos);
Dataviews = {DXeuc, DXcos};
Nobjs  = 2;

% Input precomputed weights
Weights   = load('weights.mat');
MVMCparams.W = Weights.W2D100;

% ------------------------------------------------------------------------
% Variables regarding the multiobjective algorithm
MVMCparams.MAXGEN      = 100;                      % Number of generations
MVMCparams.Popsize     = size(MVMCparams.W,1);     % Population size
MVMCparams.Niche       = 10;                       % The neighboursize T = 20
MVMCparams.Method      = 'tch';                    % The decomposition method {tch}
MVMCparams.pC          = 0.5;                      % The probability of crossing over
MVMCparams.pM          = 0.03;                     % The mutation probability
MVMCparams.NOBJ        = Nobjs;

%% Run the mvmc algorithm
isDemo = true;

OUT = mvmc(Dataviews, K, MVMCparams, isDemo);

% OUT.nPFront   returns the normalized Pareto fromt
% OUT.PClrs     returns the consensus clustering solutions in the PF
% OUT.W         returns the weights associated to each solution

%% Supervised Model Selection using ARI

[Clr,ARIb,idx1,ARIvalues] = cmsARI(OUT.PClrs, Labels);
disp(['Best ARI = ' num2str(ARIb) ' | Best ID = ' num2str(idx1) ' | PFA size = ' num2str(size(OUT.PClrs,2))]);

% plotPFA(OUT.nPFront,Nobjs,idx1)

%% (OPTIONAL STEP) Unsupervised Model Selection using Silhouette index 

[Clr,SILb,idx2,SILvalues,Wi] = cmsSILw(OUT.PClrs, OUT.W, Dataviews, Nobjs);
disp(['Best Sil = ' num2str(SILb) ' | Best ID = ' num2str(idx2) ' | Best Sil ARI = ' num2str(ARIvalues(idx2))]);

plotPFA(OUT.nPFront,Nobjs,idx1, idx2)