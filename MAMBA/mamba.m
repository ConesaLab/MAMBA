function [mip, gene_states,genes,sol, chip_out] = ...
                               mamba(model,fold_change,pvals, ...
                               mets_ms, mets_names, frac, ...
                               chip_ana, ...
                               varargin)
% MAMBA Metabolic adjusment via multiomic blocks aggregation
%
%   [MIP,GENE_STATES,GENES,SOL,CHIP_OUT] = 
%       MAMBA(MODEL,FOLD_CHANGE,PVALS,METS_MS,METS_NAMES,FRAC,
%             CHIP_ANA,parameters...)
%
%   Inputs
%   MODEL       COBRA or ELF model.  COBRA models will be converted 
%               to ELF models
%               with a warning.
%   FOLD_CHANGE Measured fold change from expression data.  Columns
%               correspond to transitions beteween conditions, rows 
%               correspond to genes.
%   PVALS       P-values for changes.  Format is the same as for
%               FOLD_CHANGE.
%   METS_MS     Measured fold change from metabolomic data. Same format
%               as FOLD_CHANGE but containing metabolites in rows instead.
%   METS_NAMES  The name of the metabolites. Note that the names in this
%               variable must match with the ones included in the metabolic
%               model. Rows in METS_MS must have the same order as in
%               METS_NAMES.
%   FRAC        (optional) Fraction of metabolic objective required
%               in the resulting model (v_obj >= frac*v_obj_max).
%               Default is 0.3.  FRAC can also be a vector giving a
%               separate fraction for each condition.
%   CHIP_ANA    Logic. Whether to perform the concordance analysis with
%               ChIP-seq data or not. 
%
%   Parameters
%   'chip_logfc' Measured fold changes from ChIP-seq data at reaction level.
%                Columns correspond to transitions between conditions and
%                rows correspond to reactions. 
%                (only required if CHIP_ANA = TRUE)
%   'chip_names' Names of the reactions in chip_logfc. 
%                (only required if CHIP_ANA = TRUE)
%   'gene_names' Cell array of names for genes in expression dataset.
%                These correspond to the rows in FOLD_CHANGE and PVALS.
%                If none is given, the rows correspond to ELF.genes.
%   'method'     Method for formulating the MIP problem.  Choices are:
%                    'milp'  Linear programming (default)
%                    'miqp'  Quadratic programming
%                The method changes how constant variables are measured.
%                Linear method uses delta variables (delta = x XOR y).
%                Quadratic method uses (x - y)^2.
%   'weighting'  Method to convert PVALS to weights.  Options include:
%                   'log'     w(p) = -log(p)
%                   'linear'  w(p) = 1 - p
%                   'unit'    w(p) = 1
%   'bounds'     Structure array of condition-specific bounds.  For
%                example, BOUNDS{i}.lb and BOUNDS{i}.ub are the lower
%                and upper bounds for the ith condition.  If not
%                specified, the bounds from MODEL are copied to each
%                condition.
%   'objs'       Cell array of condition-specific objectives.  If not
%                specified, the objective MODEL.c is copied to each
%                condition.
%   'p_thresh'   Threshold above which a P-value is not considered
%                significant.  If a gene increases with a P-value
%                p > P_THRESH, then it will be held constant with P-value
%                1 - p.  The default is 0.5.
%   'p_eps'      P-values below P_EPS are considered equal to P_EPS.  This
%                is used with log weighting to avoid taking the logarithm
%                of very small P-values.  The default is 1e-10.
%   'interaction_matrix'   A matrix (I) describing the interaction between
%                conditions.  If I(i,j) = k, then the kth column in
%                FOLD_CHANGE describes the change in expression between
%                condition i to condition j.  If not given, FOLD_CHANGE
%                assumes 1 -> 2, 2 -> 3, ..., n-1 -> n.
%   'remove_rev' Remove reversibility constraints from the model before
%                converting to a MIP.  May improve performance if activity
%                cycles are not a concern.  (Default = false)
%
%   Outputs
%   MIP         The Mixed Integer Linear Problem formulation.
%   GENE_STATES Binary expression states calculated following MADE rules.  Columns
%               correspond to conditions, rows are genes in both
%               the model and the expression data.
%   GENES       Cell of gene names corresponding to the rows in
%               GENE_STATES.
%   SOL         Solution structure with details from the MILP solver.
%   CHIP_OUT    A vector of (0,1) for each reaciton indicating whether the
%               ChIP-seq information associated to the reaction matches
%               with its flux along time.


% test inputs
assert(nargin >= 3, 'MAMBA requires at least three inputs.');

% make sure the input was an ELF model
model = assert_elf_model(model);

ntrans = size(fold_change,2);   % number of transitions
ncond  = ntrans + 1;            % number of conditions

save_models = (nargout >= 4);

assert(all(size(pvals) == size(fold_change)), ...
       'FOLD_CHANGE and PVALS must have the same dimensions');

% default frac setting
if nargin < 4 || isempty(frac),  frac = 0.3; end
% set frac for each condition
if length(frac) == 1
    fracs = repmat(frac,1,ncond);
end

% parse input param/value pairs
p = inputParser;

valid_methods = {'miqp','milp'};
p.addParamValue('method','milp',@(x) ismember(x,valid_methods));

valid_weights = {'log','linear','unit'};
p.addParamValue('weighting','log',@(x) ismember(x,valid_weights));

p.addParamValue('bounds',[]);
p.addParamValue('objs',[]);

pvalidate = @(x) validateattributes(x,'numeric', ...
                                    {'scalar','real','positive','<=',1});
p.addParamValue('p_thresh',0.5,pvalidate);
p.addParamValue('p_eps',1e-10,pvalidate);

%p.addParamValue('gene_names',[],@(x) validatestring(x,model.genes));
p.addParamValue('gene_names',[]);

p.addParamValue('interaction_matrix',[]);

p.addParamValue('remove_rev',false);

p.addParamValue('chip_logfc',[]);
p.addParamValue('chip_names',[]);

p.parse(varargin{:});

% find the genes to match
if isempty(p.Results.gene_names)
    genes = model.genes;
    igenes = 1:length(genes);
else
    [genes,igenes,igene_names] = intersect(model.genes, ...
                                           p.Results.gene_names);
    fold_change = fold_change(igene_names,:);
    pvals = pvals(igene_names,:);
end
ngenes = length(genes);

p_thresh = p.Results.p_thresh;
p_eps = p.Results.p_eps;

d = zeros(ngenes,ntrans);
P = pvals;

% convert significant fold changes to differences
d(fold_change <  -0.5 & pvals <= p_thresh) = -1;
d(fold_change >= 0.5 & pvals <= p_thresh) =  1;

% shift p-values above p_thresh onto [0, 1-p_thresh]
P(pvals > p_thresh) = 1 - pvals(pvals > p_thresh);

% convert p-values to weights
switch lower(p.Results.weighting)
    case 'log'
        P(P < p_eps) = p_eps;  % avoid taking log of zero
        w = -log(P);
    case 'linear'
        w = 1 - P;
    case 'unit'
        w = ones(size(P));
end

bounds = p.Results.bounds;
objs = p.Results.objs;

I = p.Results.interaction_matrix;
if isempty(I)
    % build sequential interaction matrix (1 -> 2, 2 -> 3, ...)
    I = zeros(ncond);
    for i = 1 : ntrans
        I(i,i+1) = i;
    end
end

% create the MILP
if p.Results.remove_rev
    milp = elf_to_milp(model,true,false);
else
    milp = elf_to_milp(model,true,true);
end

% var indexes are offset by nrxns
vars = igenes + size(model.S,2);

method = p.Results.method;

opt_vals = [];
for i = 1:ntrans
    opt_vals = [opt_vals opt_mets(mets_ms(:,i))];
end
mets_ms(:,1)
% run without growth to determine theoretical matches
[ unc_states, ~ ,unc_error] ...
    = diffadj(milp,vars,d,w,I,bounds,objs,0,method,mets_names, mets_ms, opt_vals);
% run actual MAMBA
[mip, gene_states,sol,con_error,milps] ...
    = diffadj(milp,vars,d,w,I,bounds,objs,frac,method,mets_names, mets_ms, opt_vals);
%unc_error
%con_error
%if unc_error
%    fprintf('Error:  Unconstrained model was infeasible.\n\n');
%elseif con_error
%    fprintf('Error:  The model was infeasible.\n\n');
%else
%verify solutions
%verification_error = any(~sol.verified);
verification_error = 0;
if save_models
    models = cell(1,ncond);
    for i = 1 : ncond
       models{i} = model;
       models{i}.lb = milps{i}.lb;
       models{i}.ub = milps{i}.ub;
       models{i}.c  = milps{i}.c;
     end
end
    
ratio = sol.adj_vals ./ sol.obj_vals;

if chip_ana == 'true'
  chip_out = chip_reg(milp, sol, p.Results.chip_logfc, p.Results.chip_names);
else 
    chip_out = [];
end
    
count = @(x) length(find(x));

%T = round( unc_states(:,2:end) -  unc_states(:,1:end-1));
D = round(gene_states(:,2:end) - gene_states(:,1:end-1));
    
fprintf('\nGene counts:\n');
fprintf('           |  Increasing    Decreasing      Constant\n');
fprintf('Transition | Data /  Fit   Data /  Fit   Data /  Fit\n');
for i = 1 : ntrans
        incs = count( d(:,i) ==  1 );
        incs_matched = count( d(:,i) ==  1 & D(:,i) >=  1 );
        
        decs = count( d(:,i) == -1 );
        decs_matched = count( d(:,i) == -1 & D(:,i) <= -1 );
        
        cons = count( d(:,i) ==  0 );
        cons_matched = count( d(:,i) ==  0 & D(:,i) ==  0 );
        
        [cond1,cond2] = ind2sub(size(I),find(I(:) == i));
        fprintf(' %2i ->%2i   | %4i / %4i   %4i / %4i   %4i / %4i\n', ...
                cond1, cond2, incs, incs_matched, ...
                              decs, decs_matched, ...
                              cons, cons_matched );
end
    
matches = count( D - d == 0 );
all_matches = numel(D);
%theor_matches = count( T - d == 0 );
fprintf('\nTotal match:     %i / %i (%3.1f%%)\n', ...
        matches, all_matches, matches/all_matches*100 );
%fprintf('Adjusted match:  %i / %i (%3.1f%%)\n\n', ...
%            matches, theor_matches, matches/theor_matches*100);

function [opt_vec] = opt_mets (vector)
    in_lim = -1*log10(min(pvals(pvals ~= 0)));
    opt_vec = [];
    if length(vector) == 1
        val = 1;
        opt_vec = [opt_vec; 10.^(in_lim-val)];
    else
        for v = 1:length(vector)
            val = (abs(vector(v)) - min(abs(vector))) * 1 / (max(abs(vector)) - min(abs(vector)));
            opt_vec = [opt_vec; 10.^(in_lim-val)];
        end
    end
end

end

 %set output on error
%if unc_error || con_error
%    gene_states = [];
%    models = {};
%end
    
    
