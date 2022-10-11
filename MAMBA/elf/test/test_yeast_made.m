
% COBRA toolbox model for iND750
load cobra_yeast_model.mat;
elf = make_elf_model(cobra);

elf.lb = 10000*elf.lb;
elf.ub = 10000*elf.ub;

% expression data from Roberts and Hudson, 2001
load yeast_exp_data.mat;

% there are 4 transitions in the dataset
N = 4;

% obj flux must be 30% of optimal value
frac = 0.3;

% run MADE
[gene_states,genes,sol,models] = made(elf,exp_data(:,1:N),p_vals(:,1:N),frac,'gene_names',gene_names,'method','milp');

