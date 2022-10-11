function [sol] = fba(elf)
% FBA  Run FBA on an ELF model
%
%   [SOL] = FBA(ELF)
%
%   Run Flux Balance Analysis on an ELF model.
%
%   Inputs
%   ELF     ELF model structure
%
%   Outputs
%   SOL     CMPI solution structure

sol = cmpi.solve_milp(elf_to_milp(elf));

