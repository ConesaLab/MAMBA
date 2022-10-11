function [sol] = solve_milp(milp,solver)
% SOLVE_MILP  Solve a Mixed Integer Linear Programming problem
%
%   [SOL] = SOLVE_MILP(MILP,SOLVER)
%
%   Solves the Mixed Integer Linear Programming problem
%       min sense*c*x
%       subject to
%           A*x (<=/=/>=) b
%           lb <= x <= ub
%
%   Inputs
%   MILP    Problem structure.  Fields include:
%               c         Objective coefficients
%               sense     Direction of optimization:
%                              1    Minimization
%                             -1    Maximization
%               A         Coefficient matrix for constraints.
%               b         Right-hand side of constraints.
%               ctypes    Char array denoting the type of each constraint:
%                             '<'   a*x <= b
%                             '='   a*x == b
%                             '>'   a*x >= b
%               lb        Lower bound for each variable.
%               ub        Upper bound for each variable.
%               vartypes  Char array denoting the type of each variable:
%                             'C'   Continuous
%                             'B'   Binary
%                             'I'   Integer
%               options   Solver options (see below).  If not given, the
%                         options in the global variable CMPI_OPTIONS are
%                         used (if defined).
%                             MaxTime     Maximum solution time (seconds)
%                             MaxIter     Maximum simplex iterations
%                             MaxNodes    Maximum number of nodes
%                             Display     Turn reporting to the screen
%                                         'on' or 'off' (default)
%                             FeasTol     Feasibility tolerance
%                             IntFeasTol  Integer feasibility tolerance
%                             OptTol      Optimality tolerance
%   SOLVER  Solver to use.  If not given, the value of the global
%           variable CMPI_SOLVER is used.
%
%   Outputs
%   SOL     Solution structure with fields:
%               x       Optimal vector
%               val     Objective values (sense*c*x)
%               flag    Exit flag.  Possible values are:
%                           1   Not started
%                           2   Optimal
%                           3   Infeasible
%                           4   Infeasible or unbounded
%                           5   Unbounded
%                           6   Objective worse than user cutoff
%                           7   Iteration limit reached
%                           8   Node limit reached
%                           9   Time limit reached
%                           10  Solution limit reached
%                           11  User interruption
%                           12  Numerical difficulties
%                           13  Suboptimal solution
%               output  Other solver-specific output


if nargin < 2,  solver = cmpi.get_solver(); end

if ~issparse(milp.A)
    milp.A = sparse(milp.A);
end

if ~isfield(milp,'options')
    milp.options = cmpi.get_options();
end

solver = 'gurobi';
switch solver
    case 'gurobi'
        sol = ...
            gurobi(milp);
        if sol.status == 'OPTIMAL'
            sol.flag = 2;
        else
            sol.flag = 3;
        end
        sol.val = sol.objval;
        
    case 'ilog_cplex'        
        opts = cmpi.set_cplex_opts(milp.options);
        
        milp.b = milp.b(:);

        Aineq = [ milp.A(milp.ctypes == '<',:); 
                 -milp.A(milp.ctypes == '>',:)];
        bineq = [ milp.b(milp.ctypes == '<');
                 -milp.b(milp.ctypes == '>')];
             
        Aeq = milp.A(milp.ctypes == '=',:);
        beq = milp.b(milp.ctypes == '=');
        
        [sol.x,sol.val,~,sol.output] = ...
            cplexmilp(milp.sense*milp.c(:), ...
                      Aineq, bineq, ...
                      Aeq, beq, ...
                      [], [], [], ...
                      milp.lb(:), milp.ub(:), ...
                      milp.vartypes, ...
                      [], ...
                      opts);

        sol.val = milp.sense*sol.val;
                   
        sol.flag = cmpi.get_cplex_flag(sol.output.cplexstatus);

    case 'cplex_multi'
        cplex = Cplex();
        if milp.sense == 1
           cplex.Model.sense = 'minimize';
        else
           cplex.Model.sense = 'maximize';
        end
        cplex.Model.obj = milp.c(:);
        milp.A(milp.ctypes == '>',:) = -milp.A(milp.ctypes == '>',:);
        milp.b(milp.ctypes == '>') = -milp.b(milp.ctypes == '>');
        cplex.Model.A = milp.A;
        cplex.Model.ctype = milp.vartypes;
        cplex.Model.lb = milp.lb(:);
        cplex.Model.ub = milp.ub(:);
        cplex.Model.rhs = milp.b(:);
        lhs = -inf*ones(size(milp.b(:)));
        %lhs(:) = -inf;
        lhs(milp.ctypes == '=') = milp.b(milp.ctypes == '=');
        cplex.Model.lhs = lhs;
        
        cplex.writeModel('prob.lp');
        
        cplex.solve();
        sol.x = check_field('x', cplex.Solution);
        sol.flag = get_cplex_flag(cplex.Solution.status);
        sol.val = check_field('objval', cplex.Solution);
        sol.output = cplex.Solution.statusstring;
        
        cplex.Param.mip.pool.relgap.Cur = 1000.0;
        cplex.populate();
        sol.pool = check_field('solution', cplex.Solution.pool);
        
    case 'gams_cplex'
        milp.ctype(milp.ctype == '<') = 'u';
        milp.ctype(milp.ctype == '=') = 's';
        milp.ctype(milp.ctype == '>') = 'l';
        
        [sol.f, sol.x, sol.flag] = ...
            gams_milp(milp.c, ...
                      milp.A, ...
                      milp.b, ...
                      milp.lb, ...
                      milp.ub, ...
                      milp.ctypes, ...
                      milp.vartypes, ...
                      milp.sense);
         sol.val = milp.sense*sol.f;
    case 'glpk'
        [sol.x, sol.val, sol.flag, sol.extra] = ...
            glpk(milp.c, ...
                      milp.A, ...
                      milp.b, ...
                      milp.lb, ...
                      milp.ub, ...
                      milp.ctypes, ...
                      milp.vartypes, ...
                      1);
        sol.obj_val = sol.val;

    otherwise
        error('Unrecognized solver: %s',solver);
end

end



