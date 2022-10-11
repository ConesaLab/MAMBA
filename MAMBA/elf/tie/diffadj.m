function [states,sol,mip_error,models] = ...
                        diffadj(milp,vars,d,w,I,bounds,objs,fracs,method)
% DIFFADJ  Formulate and solve the differential adjustment problem
%
%   [STATES,SOL] = DIFFADJ(MILP,VARS,W,D,I,BOUNDS,OBJS,FRACS,METHOD)
%
%   Attempts to match changes between different conditions while
%   retaining functional models.
%
%   Inputs
%   MILP    MILP model
%   VARS    List of indices for variables to adjust
%   D       List of differences for each transition:
%                1 -> x(t) < x(t')
%                0 -> x(t) = x(t')
%               -1 -> x(t) > x(t')
%           Columns correspond to conditions.  Each row 'i' is a variable
%           with index VARS(i).
%   W       List of weights for each transition.  Structure is the same
%           as for D.
%   I       Interaction matrix.  If I(i,j) = k, then the kth column in
%           D and W describe the transition between conditions i and j.
%           The number of nonzero entries in I should equal the number
%           of columns in D and W.
%   BOUNDS  Structure array of bounds for each condition:
%               BOUNDS{i}.lb
%               BOUNDS{i}.ub
%           If empty, the bounds in MILP are used.
%   OBJS    Cell array of objective coefficients for each condition.  If
%           empty, MILP.c is used.
%   FRACS   Fractions of optimal growth that must be possible in each
%           condition:
%               c'x >= FRACS(i)*c'x_max
%           If only a single number is given, it is used for every
%           condition.  If empty, the default is 0.3.
%   METHOD  Method used for solving the problem.  Options are:
%               'qp'   Mixed-Integer Quadratic (default)
%               'milp' Mixed-Integer Linear
%
%   Outputs
%   STATES     A |variables| x |conditions| matrix of states (levels) for
%              each variable in each condition.
%   SOL        Solution object from CMPI.  Contains the additional fields
%                  'obj_vals' Objective values through original models
%                  'adj_vals' Objective values through adjusted models
%                  'verified' True if obj_val(i) >= FRAC(i)*c'x_max
%   MIP_ERROR  True if an error was encourtered during the optimization.
%   MODELS     Models with STATES applied to VARS.


if nargin < 6 || isempty(bounds)
    bounds = [];
end

if nargin < 7 || isempty(objs)
    objs = [];
end

if nargin < 8 || isempty(fracs)
    fracs = 0.3;
end

if nargin < 9 || isempty(method)
    method = 'qp';
end

switch lower(method)
    case {'qp','miqp'}
        qp = true;
    case {'lp','milp'}
        qp = false;
    otherwise
        error('unknown method: %s',method);
end

[nvars,ntrans] = size(w);
ncond = ntrans + 1;

% check interaction matrix
assert(all(size(I) == [ncond,ncond]), 'I not square or wrong size');
assert(length(find(I)) == ntrans, 'I does not match w and d');
assert(all(ismember(1:ntrans,I(:))), 'I is missing transition indices');

% fill out fracs
if length(fracs) == 1
    fracs = repmat(fracs,1,ncond);
end

% fill out bounds
if length(bounds) <= 1
    for i = 1 : ncond
        if isempty(bounds)
            bounds{i}.lb = milp.lb;
            bounds{i}.ub = milp.ub;
        else
            % one provided; replicate it
            bounds{i}.lb = bounds{1}.lb;
            bounds{i}.ub = bounds{1}.ub;
        end
    end
end

% fill out objs
if length(objs) <= 1
    for i = 1 : ncond
        if isempty(objs)
            objs{i} = milp.c;
        else
            % one provided; replicate it
            objs{i} = objs{1};
        end
    end
end     
    
% make MILPs for each condition
milps = cell(1,ncond);
obj_vals = zeros(1,ncond);
for i = 1 : ncond
    milps{i} = milp;
    milps{i}.c(1:length(objs{i})) = objs{i};
    milps{i}.lb = bounds{i}.lb;
    milps{i}.ub = bounds{i}.ub;
    
    sol = cmpi.solve_milp(milps{i});
    if ~cmpi.is_acceptable_exit(sol)
        error('Error optimizing condition %i',i);
    end
    obj_vals(i) = sol.val;
    
    % add growth constraint
    milps{i}.A(end+1,:) = milps{i}.c(:)';
    milps{i}.b(end+1) = fracs(i)*obj_vals(i);
    milps{i}.ctypes(end+1) = '>';
end

% number of elements to hold constant
n_con = length(find(d(:) ==  0));
% number of variables and constraints added to create delta variables
if qp
    n_cols_per_con = 1;
    n_rows_per_con = 1;
else
    n_cols_per_con = 3;
    n_rows_per_con = 6;
end

% initialize the MIP structure
[rowsA,colsA] = size(milps{1}.A);
n_miprows = ncond*rowsA + n_rows_per_con*n_con;
n_mipcols = ncond*colsA + n_cols_per_con*n_con;
mip.A = sparse(n_miprows,n_mipcols);
mip.lb = zeros(n_mipcols,1);
mip.ub = zeros(n_mipcols,1);
mip.b = zeros(n_miprows,1);
mip.c = zeros(n_mipcols,1);
mip.vartypes = repmat('C',1,n_mipcols);
mip.ctypes = repmat('=',1,n_miprows);
mip.sense = 1; % minimize
if qp
    mip.Q = sparse(n_mipcols,n_mipcols);
end

% populate the MIPs for each condition
coff = 0;
roff = 0;
for i = 1 : ncond
    mip.A(roff+1:roff+rowsA,coff+1:coff+colsA) = milps{i}.A;
    mip.lb(coff+1:coff+colsA) = milps{i}.lb;
    mip.ub(coff+1:coff+colsA) = milps{i}.ub;
    mip.b(roff+1:roff+rowsA) = milps{i}.b;
    mip.vartypes(coff+1:coff+colsA) = milps{i}.vartypes;
    mip.ctypes(roff+1:roff+rowsA) = milps{i}.ctypes;
    
    roff = roff + rowsA;
    coff = coff + colsA;
end

% compute the index of variable V in condition C
idxof = @(v,c) colsA*(c-1) + v;

% create the objective function
for t = 1 : ntrans
    % find the interacting conditions
    [cond1,cond2] = ind2sub(size(I),find(I(:) == t));
    for v = 1 : nvars
        idx1 = idxof(vars(v),cond1);
        idx2 = idxof(vars(v),cond2);
        switch d(v,t)
            case {1}
                % increasing
                mip.c(idx1) = mip.c(idx1) + w(v,t);
                mip.c(idx2) = mip.c(idx2) - w(v,t);
            case {-1}
                % decreasing
                mip.c(idx1) = mip.c(idx1) - w(v,t);
                mip.c(idx2) = mip.c(idx2) + w(v,t);
            case {0}
                % constant
                diff_idx = coff + 1;
                if qp
                    mip.Q(diff_idx,diff_idx) = w(v,t);
                    add_qp_diff_cons(idx1,idx2,diff_idx);
                else
                    mip.c(diff_idx) = w(v,t);
                    add_ip_diff_cons(idx1,idx2,diff_idx);
                end
                coff = coff + n_cols_per_con;
                roff = roff + n_rows_per_con;
        end
    end
end

if qp
    sol = cmpi.solve_miqp(mip);
else
    sol = cmpi.solve_milp(mip);
end
sol.mip = mip;

mip_error = ~cmpi.is_acceptable_exit(sol);
if ~mip_error
    states = zeros(nvars,ncond);
    for c = 1 : ncond
        for v = 1 : nvars
            states(v,c) = sol.x(idxof(vars(v),c));
        end
    end

    % round binary variables
    [~,bins] = intersect(vars,find(mip.vartypes == 'B'));
    states(bins,:) = round(states(bins,:));
    
    % verify solutions
    sol.verified = false(1,ncond);
    sol.adj_vals = zeros(1,ncond);
    sol.obj_vals = obj_vals;
    for i = 1 : ncond
        % copy models and apply gene states
        models{i} = milps{i};
        models{i}.lb(vars) = states(:,i);
        models{i}.ub(vars) = states(:,i);
        
        % run FBA to verify the solution
        kosol = cmpi.solve_milp(models{i});
        sol.verified(i) = cmpi.is_acceptable_exit(kosol);
        if ~isempty(kosol.val)
            sol.adj_vals(i) = kosol.val;
        end
    end
else
    % solver did not return a solution
    states = [];
    models = {};
end

% ------------------------------------------------------------------------

function add_qp_diff_cons(i1,i2,di)
    % add constraints for di = x(i1) - x(i2)
    % TODO:  use off-diagonal entries in Q to eliminate constraint
    mip.A(roff+1,[i1 i2 di]) = [1 -1 -1];
    mip.b(roff+1) = 0;
    mip.ctypes(roff+1) = '=';
    mip.vartypes(coff+1) = mip.vartypes(i1);
    switch mip.vartypes(i1)
        case {'B'}
            mip.ub(coff+1) = 1;
        case {'C'}
            max_flux = 2*max( abs([mip.ub([i1 i2]) mip.lb([i1 i2])]) );
            mip.ub(coff+1) =  max_flux;
            mip.lb(coff+1) = -max_flux;
    end
end

function add_ip_diff_cons(i1,i2,di)
    % add constraints for di = x(i1) XOR x(i2)
    % TODO:  enable MILPs for continuous differences
    I1 = di + 1; % indicator variable 1
    I2 = di + 2; % indicator variable 2
    mip.A(roff+1,[I1 I2 di]) = [2 2 -4];
    mip.b(roff+1) = 3;
    mip.A(roff+2,[I1 I2 di]) = [-2 -2 4];
    mip.b(roff+2) = 1;
    mip.A(roff+3,[i1 i2 I1]) = [-1 -1 -3];
    mip.b(roff+3) = -2;
    mip.A(roff+4,[i1 i2 I1]) = [1 1 3];
    mip.b(roff+4) = 4;
    mip.A(roff+5,[i1 i2 I2]) = [1 1 -3];
    mip.b(roff+5) = 0;
    mip.A(roff+6,[i1 i2 I2]) = [-1 -1 3];
    mip.b(roff+6) = 2;
    
    mip.vartypes(coff+1:coff+3) = 'BBB';
    mip.ctypes(roff+1:roff+6) = '<<<<<<';
    mip.ub(coff+1:coff+3) = 1;
end

end % function
                






    