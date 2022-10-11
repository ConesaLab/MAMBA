function out = chip_reg(milp, sol, chip_logfc, chip_names)
%
%   [OUT] = CHIP_REG(MILP,SOL,CHIP_LOGFC,CHIP_NAMES)
%
%   Analyze concordande between reaction fluxes and ChIP-seq data.
%
%   Inputs
%   MILP   MILP model
%   SOL    MAMBA solution of the MILP model
%   CHIP_LOGFC Measured fold changes from ChIP-seq data at reaction level.
%                Columns correspond to transitions between conditions and
%                rows correspond to reactions
%   CHIP_NAMES Names of the reactions in chip_logfc
%
%   Outputs
%   OUT A vector of (0,1) for each reaciton indicating whether the
%       ChIP-seq information associated to the reaction matches
%       with its flux along time


[nvars,ntrans] = size(chip_logfc);
ncond = ntrans + 1;

d = zeros(nvars,ntrans);
% convert significant fold changes to differences
d(chip_logfc <  -1.0) = -1;
d(chip_logfc >= 1.0) =  1;      

% Check chip-seq concordance

[rowsA,colsA] = size(milp.A);
out = [];

for m = 1:length(chip_names)
    index = find(not(cellfun('isempty', regexp(milp.rxns, chip_names{m}))));
    sel = [];
    
    for c = 1:ncond
        sel = [sel; sol.x(index+(colsA*(c-1)),:)];
    end

    status = [];
    for t = 1:ntrans
        res = sign (sel(t+1)-sel(t)) == sign(d(m,t));
        status = [status, res];
    end

    out = [out; status];
end





end % function