function elf = make_elf_model(cobra)
% MAKE_ELF_MODEL  Create an ELF model from a COBRA toolbox model
%
%   [ELF] = MAKE_ELF_MODEL(COBRA)
%
%   Create an Enzyme-Linked FBA (ELF) model from a COBRA toolbox model.
%   For details on ELF, see [cite:elf].
%
%   Inputs
%   COBRA   Model structure from the COBRA toolbox.
%
%   Outputs
%   ELF     ELF model structure.  Fields added are:
%               A      ELF constraint matrix coefficients.
%               b      RHS for Ax, resized to accomodate extra constraints
%                      in A (as opposed to S).
%               param  Structure of parameters for ELF methods:
%                         binary     Logical vector, true if a variable
%                                    is binary
%                         ctype      CMPI constraint types for rows in A.
%                         nmets      Number of metabolites in S.
%                         nrxns      Number of reaction fluxes in S.
%                         nrgbs      Number of RGB variables in A.
%                         rev        Logical vector, true if a reaction in
%                                    S is reversible
%                         has_gpr    Logical vector, true if a reaction in
%                                    S has a nonempty GPR
%                         n_has_gpr  Number of reactions in S with
%                                    nonempty GPRs.
%                         n_fr_pairs Number of f/r binary variable pairs.
%
%

S = cobra.S;
gpr = cobra.grRules;
genes = cobra.genes;

ngenes = length(genes);
[nmets,nrxns] = size(S);

has_gpr = cellfun(@(x) ~isempty(x), gpr);
n_has_gpr = sum(has_gpr);

% create lists for RGBs
exprs = cell(size(gpr));
gene_assocs = cell(size(gpr));
for i = 1 : nrxns
    if has_gpr(i)
        exprs{i} = parse_string(gpr{i});
        lists = make_dnf(exprs{i});
        gene_assocs{i} = cellfun(@(x) genes_to_idxs(x), lists, ...
                                 'UniformOutput', false);
    else
        exprs{i} = [];
        gene_assocs{i} = [];
    end
end

rev = (cobra.lb < 0);

% total number of RGBs added
nrgbs = sum((1 + rev) .* cellfun(@(x) length(x), gene_assocs));

% start with all the COBRA fields
elf = cobra;

A = sparse(nmets+2*ngenes+n_has_gpr,nrxns+2*ngenes+nrgbs);
elf.lb = zeros(nrxns + 2*ngenes + nrgbs,1);
elf.ub = zeros(nrxns + 2*ngenes + nrgbs,1);
elf.lb(1:nrxns) = cobra.lb;
elf.ub(1:nrxns) = cobra.ub;
elf.ub(nrxns+1:nrxns+ngenes) = 1;
elf.b = zeros(nmets+2*ngenes+n_has_gpr,1);%

% set parameters
elf.param.binary = [false(1,nrxns) true(1,ngenes) false(1,ngenes+nrgbs)];
elf.param.ctype  = [repmat('=',nmets+ngenes+n_has_gpr,1); ...
                    repmat('<',ngenes,1)];
elf.param.nmets = nmets;
elf.param.nrxns = nrxns;
elf.param.n_has_gpr = n_has_gpr;
elf.param.nrgbs = nrgbs;

coff = nrxns + 2*ngenes;
roff = nmets + ngenes;
roffplus = nmets;
coffplus = nrxns;
fcoff = nrxns + 2*ngenes + nrgbs;
froff = nmets + 2*ngenes + n_has_gpr;
for r = 1 : nrxns
    isozymes = gene_assocs{r};
    if has_gpr(r)
        roff = roff + 1;
        %roffplus = roffplus + 1;
        % balance sum of RGBs with original reaction
        A(roff,r) = -1;
        
        % create RGBs
        for i = 1 : length(isozymes)
            coff = coff + 1;
            coffplus = coffplus + 1;
            A(1:nmets,coff) = S(:,r);       % metabolite coefficients
            A(isozymes{i},coff) = -1;       % gene coefficients
            A(roff,coff) = 1;               % v_orig = sum(v_rgbs)
            %A(roffplus,coffplus) = 1;       % GPRs constraints
            elf.ub(coff) = cobra.ub(r);
            if rev(r)
                coff = coff + 1;
                A(1:nmets,coff) = -S(:,r);
                A(isozymes{i},coff) = -1;
                A(roff,coff) = -1;
                elf.ub(coff) = abs(cobra.lb(r));
            end
        end
        
        % create f and r constraints
        if rev(r)
            disp(rev(r))
            niso = length(isozymes);
            % sum forward reactions
            froff = froff + 1;
            fcoff = fcoff + 1;
            A(froff,coff-2*niso+1:coff-niso) = 1;
            A(froff,fcoff) = -niso*cobra.ub(r);
            % sum reverse reactions
            froff = froff + 1;
            fcoff = fcoff + 1;
            A(froff,coff-niso+1:coff) = 1;
            A(froff,fcoff) = -niso*abs(cobra.lb(r));
            % NOT constraint on f and r
            froff = froff + 1;
            A(froff,[fcoff-1 fcoff]) = 1;
            b(froff) = 1;
            % update parameters
            elf.param.binary(end+1:end+2) = true;
            elf.param.ctype(end+1:end+2) = '<';
            elf.param.ctype(end+1) = '=';
            elf.lb(fcoff-1:fcoff) = 0;
            elf.ub(fcoff-1:fcoff) = 1;
        end
    else
        % leave reaction in S matrix
        A(1:nmets,r) = S(:,r);
    end
end

% add gene exchanges
row_range = nmets+1:nmets+ngenes;
col_range = nrxns+ngenes+1:nrxns+2*ngenes;
A(row_range,col_range) = eye(ngenes);

% add binding constraints
row_range = nmets+ngenes+n_has_gpr+1:nmets+ngenes+n_has_gpr+ngenes;
col_range = nrxns+1:nrxns+2*ngenes;
% compute maximum upper bound for each gene
if ~isfield(cobra,'rxnGeneMat')
    rxnGeneMat = make_rxnGeneMat(cobra,exprs);
else
    rxnGeneMat = cobra.rxnGeneMat;
end
max_up = max([(rxnGeneMat' * cobra.ub) ...
              (rxnGeneMat' * (abs(cobra.lb) .* rev))], [],2);
bounds = repmat(max_up(:)',ngenes,1) .* eye(ngenes);
A(row_range,col_range) = [-bounds eye(ngenes)];
elf.ub(nrxns+ngenes+1:nrxns+2*ngenes) = max_up;
elf.A = A;
elf.c = zeros(size(A,2),1);
elf.c(1:nrxns) = cobra.c;
%elf.b = b;
elf.b(1:nmets) = cobra.b;

% more params
elf.param.has_gpr = has_gpr;
elf.param.rev = rev;
elf.param.n_fr_pairs = sum(rev);
                
function [loc] = genes_to_idxs(gene_list)
    [~,loc] = ismember(gene_list,genes);
    loc = loc + nmets;
end

end


    
