%% DEMO: Test the mvmv algorithm on different feature and relational data 
% ------------------------------------------------------------------------
% A multiobjective evolutionary algorithm for multiview clustering (MVCM).
% This algorithm using MOEA/D as the underlaying optimizer
% This toolbox was developed with MATLAB R2014a.
%
% Developed by
%   Adan Jose-Garcia (adanjoga@gmail.com, adan.jose@cinvestav.mx)
%
% Please, cite the following paper where this algorithm was introduced:
%
%   [PENDING]...
%
%------------------------------------------------------------------------
clear all; close all; clc;

addpath([pwd '/mvmc']); 
addpath([pwd '/datasets/synthetic']);
addpath([pwd '/datasets/real']);
addpath([pwd '/objectives']);
addpath([pwd '/utils']);

%% Preparation of input variables

% Variables regarding the clustering problem
% List of datasets provided
iDS     = 4;
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
MVMCparams.NOBJ        = size(Dataviews,2);

MVMCparams.Labels      = Labels;
%% Run the mvmc algorithm
isDemo = true;

OUT = mvmc(Dataviews, K, MVMCparams, isDemo);

%% Plots and Stats

nPareto = size(OUT.nPFront,1);
ARIs = zeros(nPareto,1);
for i = 1:nPareto
    Yb = OUT.PClrs(:,i);
    ARIs(i) = pairwiseindex(MVMCparams.Labels,Yb);
end
[ARI,idx] = max(ARIs,[],1);
disp(['Best ARI = ' num2str(ARI) '| Best ID = ' num2str(idx) ' | PFA size = ' num2str(nPareto)]);

% Plot the best APF's solution and the best optimal Sil values
H = figure(123); hold on; set(H,'Color',[1, 1 ,1])
if MVMCparams.NOBJ == 2
    plot(OUT.nPFront(idx,1),OUT.nPFront(idx,2),'sr','MarkerFaceColor',[1 .6 .6],'MarkerSize',12);
else
    plot3(OUT.nPFront(idx,1),OUT.nPFront(idx,2),OUT.nPFront(idx,3),'sr','MarkerFaceColor',[1 .6 .6],'MarkerSize',12);
end

%% Selection of the best solution using Silhouette index
W = OUT.W; CLRs = OUT.PClrs; PFsize = size(OUT.nPFront,1);
nDX1 = minmaxnorm(Dataviews{1});
nDX2 = minmaxnorm(Dataviews{2});
SILs = NaN(PFsize,1);
for i =1:PFsize
    Wi = W(i,:); Ci = CLRs(:,i); Ki=numel(unique(Ci));
    Dws = (Wi(1).* nDX1) + (Wi(2).* nDX2);
    
    SILs(i) = ws_silindex(Ci, Ki, Dws);
end
[sARI, sidx] = max(SILs);
disp(['Sil ARI = ' num2str(ARIs(sidx)) '| Sil ID = ' num2str(sidx)]);
%figure; plot(SILs,'DisplayName','Sils') ; figure; area(W,'DisplayName','W')

% Plot the best APF's solution and the best optimal Sil values
H = figure(123); hold on; set(H,'Color',[1, 1 ,1])
if MVMCparams.NOBJ == 2
    plot(OUT.nPFront(sidx,1),OUT.nPFront(sidx,2),'sb','MarkerFaceColor',[1 .6 .6],'MarkerSize',12);
else
    plot3(OUT.nPFront(sidx,1),OUT.nPFront(sidx,2),OUT.nPFront(idx,3),'sb','MarkerFaceColor',[1 .6 .6],'MarkerSize',12);
end