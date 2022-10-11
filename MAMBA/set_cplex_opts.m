function [opts] = set_cplex_opts(options)
% SET_CPLEX_OPTS  Convert a CMPI options structure into CPLEX opts

    opts = cplexoptimset;
    setifdef('MaxTime','MaxTime');
    setifdef('MaxIter','MaxIter');
    setifdef('MaxNodes','MaxNodes');
    setifdef('FeasTol','EpRHS');
    setifdef('IntFeasTol','TolXInteger');
    setifdef('OptTol','TolFun');
    setifdef('Display','Diagnostics');
    if isfield(options,'Display') && strcmpi(options.Display,'on')
        opts.Display = 'iter';
    end
    
    function setifdef(cmpifield,cplexfield)
        if isfield(options,cmpifield)
            opts.(cplexfield) = options.(cmpifield);
        end
    end

end