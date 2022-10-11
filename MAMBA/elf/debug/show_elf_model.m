function show_elf_model(elf)
% SHOW_ELF_MODEL  Show equations in an ELF model for debugging.
%
%   SHOW_ELF_MODEL(ELF)
%
%   Prints readable equations for an ELF model.  Useful for debugging.
%
%   Inputs
%   ELF     ELF model structure

[nmets,nrxns] = size(elf.S);
ngenes = length(elf.genes);
A = full(elf.A);
nrgbs = elf.param.nrgbs;

[rnames,cnames] = get_elf_names(elf);
          
make_str = @(coef,name) [sprintf('%+g*',coef) name];
% show rgbs
fprintf('----- RGBs -----\n');
for i = nrxns+2*ngenes+1:nrxns+2*ngenes+nrgbs
    fprintf('%s:  ',cnames{i});
    nzs = find(A(1:nmets+ngenes,i) < 0);
    for j = 1 : length(nzs)
        fprintf('%s ', make_str(abs(A(nzs(j),i)),rnames{nzs(j)}));
    end
    if elf.lb(i) < 0
        fprintf('<--> ');
    else
        fprintf('--> ');
    end
    nzs = find(A(1:nmets+ngenes,i) > 0);
    for j = 1 : length(nzs)
        fprintf('%s ', make_str(A(nzs(j),i),rnames{nzs(j)}));
    end
    fprintf('\n');
end

cmpi.show_mip(elf_to_milp(elf),[],[],rnames,cnames,true);
    
