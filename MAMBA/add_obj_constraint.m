function [milp] = add_obj_constraint(milp,frac)
% ADD_OBJ_CONSTRAINT  Add minimal objective constraint to a MILP
%
%   [MILP] = ADD_OBJ_CONSTRAINT(MILP,FRAC)
%
%   Adds a constraint the c'x >= FRAC*fmax, where
%       fmax = max c'x s.t. Ax <= b
%
%   Inputs
%   MILP    CMPI MILP problem structure
%   FRAC    Fraction of maximum objective value that must be
%           reached (optional, default = 1.0)
%
%   Outputs
%   MILP    CMPI MILP problem structure with objective constraint

if nargin < 2,  frac = 1.0; end

sol = cmpi.solve_milp(milp);

milp.A = [milp.A; -milp.c(:)'];
milp.b(end+1) = -frac*sol.val;
milp.ctypes(end+1) = '<';

