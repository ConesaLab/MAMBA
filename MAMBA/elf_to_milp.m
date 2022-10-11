function [milp] = elf_to_milp(elf,indicators,reverse)
% ELF_TO_MILP  Convert and ELF model to a MILP
%
%   [MILP] = ELF_TO_MILP(ELF,INDICATORS,REVERSE)
%
%   Converts an ELF model to a CMPI MILP structure.
%
%   Inputs
%   ELF  ELF model structure.
%   INDICATORS  If true, binary indicator contraints (I_g) for each
%               gene are included.  If false, the constraints and 
%               variables are zeroed for removal during presolve.
%               (Default = true)
%   REVERSE     If true, binary reversibility constraints (f,r) are
%               added to prevent activity cycles.  If false, the
%               variables and constraints are zeroed for removal
%               during presolve.  (Default = true)
%

if nargin < 3 || isempty(reverse),  reverse = true; end
if nargin < 2 || isempty(indicators),  indicators = true; end

% find rows and columns to be zeroed
czeros = false(1,size(elf.A,2));
rzeros = false(1,size(elf.A,1));
nrxns = size(elf.S,2);
ngenes = length(elf.genes);
nrgbs = elf.param.nrgbs;
n_fr_pairs = elf.param.n_fr_pairs;
if ~indicators
    czeros(nrxns+1:nrxns+ngenes) = true;
end
if ~reverse
    czeros(nrxns+2*ngenes+nrgbs+1:nrxns+2*ngenes+nrgbs+n_fr_pairs) = true;
end
rzeros(any(elf.A(:,czeros) ~= 0,2)) = true;

milp.c = elf.c;
milp.A = elf.A;
milp.b = elf.b;
milp.lb = elf.lb;
milp.ub = elf.ub;

% zero corresponding rows and columns
milp.A(:,czeros) = 0;
milp.lb(czeros) = 0;
milp.ub(czeros) = 0;
milp.A(rzeros,:) = 0;
milp.b(rzeros) = 0;
milp.mets = elf.mets;
milp.S = elf.S;
milp.rxns = elf.rxns;
milp.vartypes = repmat('C',1,size(elf.A,2));
milp.vartypes(elf.param.binary) = 'B';
             
milp.ctypes = elf.param.ctype;
milp.sense = elf.param.ctype;
milp.obj = milp.c;
milp.rhs = milp.b;
milp.vtype = milp.vartypes;
milp.modelsense = 'min';
