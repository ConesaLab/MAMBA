classdef cmpi
% CMPI  Common Mathematical Programming Interface
%
%   CMPI defines a common interface for mathematical programming solvers.
%   Its methods help construct and solve LP/QP/MILP/MIQP problems.
%
%   CMPI represents a mathematical program as a structure with various
%   fields.  A description of these fields can be found in the 
%   documentation for the SOLVE_MILP function (help solve_milp).
    

properties (Constant)
    % Initial (default) value for CMPI_SOLVER
    init_SOLVER = 'ilog_cplex';
end

methods (Static)
    function [ctype] = make_ctype(leq,geq,eq)
        % MAKE_CTYPE  Make contraint type identifiers
        %
        %   [CTYPE] = MAKE_CTYPE(LEQ,GEQ,EQ)
        %
        %   Returns a char array defining the constraint types for
        %   LEQ '<=' inequalities, GEQ '<=' inequalities, and 
        %   EQ equalities.
        
        if nargin < 3,   eq = 0; end
        if nargin < 2,  geq = 0; end

        ctype = [repmat('<', 1, leq) ...
                 repmat('>', 1, geq) ...
                 repmat('=', 1, eq)];
    end

    function [vartype] = make_vartype(c,b,i,s,n)
        % MAKE_VARTYPE  Make variable type identifiers
        %
        %   [VARTYPE] = MAKE_VARTYPE(C,B,I,S,N)
        %
        %   Returns a char array defining the variable types for
        %   C continuous, B binary, I integer, S semi-continuous,
        %   and N semi-integer-continuous variables.
        
        if nargin < 5,  n = 0; end
        if nargin < 4,  s = 0; end
        if nargin < 3,  i = 0; end
        if nargin < 2,  b = 0; end

        vartype = [repmat('C', 1, c) ...
                   repmat('B', 1, b) ...
                   repmat('I', 1, i) ...
                   repmat('S', 1, s) ...
                   repmat('N', 1, n)];
    end

    function [sense] = make_sense(s)
        % MAKE_SENSE  Convert optimization sense to an integer
        %
        %   Converts the strings 'min', 'minimize', 'max', and 'maximize'
        %   into the corresponding optimization directions (1 and -1).
        
        switch lower(s)
            case {'min','minimize'}
                sense =  1;
            case {'max','maximize'}
                sense = -1;
        end
    end

    function [tf] = is_acceptable_exit(sol)
        % IS_ACCEPTABLE_EXIT  Determine if a solution was found
        %
        %   Determines if a CMPI solution structure contains an error
        %   code that returns a valid, feasible solution.  Acceptable
        %   exits are:
        %       Optimal
        %       Iteration limit reached
        %       Node limit reached
        %       Time limit reached
        %       Solution limit reached
        
        tf = ( (sol.flag == 2) || (sol.flag == 7) ...
              || (sol.flag == 8) || (sol.flag == 9) ...
              || (sol.flag == 10) ) && ~isempty(sol.x);
    end

    function set_option(option,val)
        % SET_OPTIONS  Set a default solver option
        %
        %   Sets an option in the default option structure.  For a
        %   description of the solver options, see the documentation
        %   for SOLVE_MILP.
        
        global CMPI_OPTIONS
        CMPI_OPTIONS.(option) = val;
    end

    function clear_options()
        % CLEAR_OPTIONS  Clear the default solver options
        global CMPI_OPTIONS
        CMPI_OPTIONS = [];
    end

    function [opts] = get_options()
        % GET_OPTIONS  Return the default solver options
        global CMPI_OPTIONS
        opts = CMPI_OPTIONS;
    end

    function set_solver(solver)
        % SET_SOLVER  Set the default solver name
        global CMPI_SOLVER
        CMPI_SOLVER = solver;
    end

    function [solver] = get_solver()
        % GET_SOLVER  Get the default solver name
        global CMPI_SOLVER
        solver = CMPI_SOLVER;
    end

    function init()
        % INIT  Initialze CMPI
        %
        %   Creates the global variables and default settings for the
        %   solver name and options.
        
        global CMPI_OPTIONS
        global CMPI_SOLVER

        CMPI_SOLVER = cmpi.init_SOLVER;
    end

    % external static methods
    [sol] = solve_milp(milp,solver)
    [tf] = verify_sol(milp,sol,tol)
    [milp] = add_obj_constraint(milp,frac)
    [tf,violators] = check_bounded_rxns(model,sol,tol)
    show_mip(mip,rowidxs,colidxs,rownames,colnames,showvars)
    [opts] = set_cplex_opts(options)
    [flag] = get_cplex_flag(status)
end

end % class
        
        
    
