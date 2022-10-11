function [tf,violators] = check_bounded_rxns(model,sol,tol)
% This function is deprecated from the iCot package.

if nargin < 3,  tol = 1e-8; end

bounded_rxns = find(model.rxn_idxs);

v = sol(1:model.nrxns);
x = sol(model.nrxns+1:model.ntotal);

violators = bounded_rxns( abs(v(bounded_rxns)) > tol ...
                            & x(model.rxn_idxs(bounded_rxns)) < 0.5 );

tf = isempty(violators);

