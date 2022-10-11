function [rownames,colnames] = get_elf_names(elf)
% GET_ELF_NAMES  Returns row and column names for an ELF model
%
%   [ROWNAME,COLNAMES] = GET_ELF_NAMES(ELF)
%
%   Creates cell arrays of names describing each constraint or variable
%   in an ELF model.
%
%   Inputs
%   ELF     ELF model structure
%
%   Outputs
%   ROWNAMES    Names for each constraint:
%                   BALi    The ith balance constraint:
%                               v_i = sum(v_rgb,i).
%                   ExB_x   Binding constraint for gene x:
%                               a_x <= a_x,max * I_x
%                   FRi     Forward-reversible constraints.
%   COLNAMES    Names for each variable:
%                   I_x     Boolean indicator variable for gene x.
%                   a_x     Activity (exchange flux) for enzyme x.
%                   rgbi    The ith RGB reaction flux.
%                   fi      Boolean indicator if the reversible reaction
%                           i carried positive flux.
%                   ri      Boolean indicator if the reversible reaction
%                           i carries negative flux.
%

[nmets,nrxns] = size(elf.S);
ngenes = length(elf.genes);
[rA,cA] = size(elf.A);
nrgbs = elf.param.nrgbs;
n_has_gpr = elf.param.n_has_gpr;

n_fr_pairs = (cA - (nrxns + 2*ngenes + nrgbs)) / 2;

catwith = @(s,array) cellfun(@(x) [s x], array, 'Uniform', false);
catwithn = @(s,n) catwith(s, arrayfun(@(x) num2str(x), ...
                                      1:n, 'Uniform', false));

colnames = [elf.rxnNames' ...
            catwith('I_', elf.genes') ...
            catwith('a_', elf.genes') ...
            catwithn('rgb', nrgbs) ...
            make_fr(n_fr_pairs)];

rownames = [elf.metNames' ...
            elf.genes' ...
            catwithn('BAL', n_has_gpr) ...
            catwith('ExB_', elf.genes') ...
            catwithn('FR', rA - (nmets + 2*ngenes + n_has_gpr))];
      
function [frs] = make_fr(n)
    frs = cell(1,2*n);
    for i = 1 : n
        frs{2*(i-1)+1} = ['f' num2str(i)];
        frs{2*(i-1)+2} = ['r' num2str(i)];
    end
    