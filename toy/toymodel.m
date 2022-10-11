initCobraToolbox

% Metabolic model
load('~/srv/mamba/MAMBA/toy/toy.mat')

% Gene expression
log_fc = dlmread('~/srv/mamba/MAMBA/toy/logfc.txt', ',');
pvals = dlmread('~/srv/mamba/MAMBA/toy/pval.txt', '\t');
fid = fopen('~/srv/mamba/MAMBA/toy/gpr.txt');
txt = textscan(fid,'%s','delimiter','\n');
txt;
gene_names = txt{1,1};

% Metabolomics
mets_ms = dlmread('~/srv/mamba/MAMBA/toy/mets_lfc.txt', '\t');
fid = fopen('~/srv/mamba/MAMBA/toy/mets.txt');
txt = textscan(fid,'%s','delimiter','\n');
txt;
mets_names = txt{1,1};

% ChIP-seq data
chip_ana = 'true';
chip_logfc = dlmread('~/srv/mamba/MAMBA/toy/chip_ms.txt', '\t');
fid = fopen('~/srv/mamba/MAMBA/toy/chip_names.txt');
txt = textscan(fid,'%s','delimiter','\n');
txt;
chip_names = txt{1,1};

% MAMBA
[model, gene_states,genes,sol,chip_out] = mamba(toymip,log_fc,pvals, ...
    mets_ms, mets_names, 1/3, ...
    chip_ana, ...
    'chip_logfc',chip_logfc, ...
    'chip_names',chip_names, ...
    'gene_names',gene_names, ...
    'method','milp', ...
    'weighting', 'log');
