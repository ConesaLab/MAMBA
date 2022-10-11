function [rgm] = make_rxnGeneMat(cobra,exprs)
% MAKE_RXNGENEMAT  Build a rxnGeneMat for cobra models
%
%   [RGM] = MAKE_RXNGENEMAT(COBRA,EXPRS)
%
%   Create a reaction-gene matrix (rxnGeneMat) for a COBRA toolbox model.
%   The (i,j)th entry indicates whether the ith gene appears is the GPR
%   for the jth reaction.
%
%   Inputs
%   COBRA   A COBRA toolbox model
%   EXPRS   (optional) A set of expr objects parsed from the GPR.  This
%           avoids the need to parse the GPR again.
%
%   Outputs
%   RGM     The rxnGeneMat in the form of the COBRA toolbox field of the
%           same name.

nrxns = size(cobra.S,2);
ngenes = length(cobra.genes);

has_gpr = cellfun(@(x) ~isempty(x), cobra.grRules);

if nargin < 2
    exprs = cell(1,nrxns);
    for i = 1 : nrxns
        if has_gpr(i)
            exprs{i} = parse_string(cobra.grRules{i});
        else
            exprs{i} = [];
        end
    end
end

rgm = false(nrxns,ngenes);
for i = 1 : nrxns
    if has_gpr(i)
        rgm(i,:) = ismember(cobra.genes, ...
                            exprs{i}.get_atoms());
    end
end

