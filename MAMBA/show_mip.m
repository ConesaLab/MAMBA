function show_mip(mip,rowidxs,colidxs,rownames,colnames,showvars)
% SHOW_MIP  Show equations for a MIP structure
%
%   SHOW_MIP(MIP,ROWIDXS,COLIDXS,ROWNAMES,COLNAMES,SHOWVARS)
%
%   Print algebraic equations for a MIP structure.
%
%   Inputs
%   MIP      CMPI model structure
%   ROWIDXS  Indexes of rows (constraints) to show.  (Default = all rows)
%   COLIDXS  Indexes of columns (variables) to show.  (Default = all rows)
%   ROWNAMES Cell array of names for each row.  Empty values are replaced
%            with 'ROW_x'.  ROWNAMES may be shorter than the number of
%            rows; extra 'ROW_x' names are added automatically.
%   COLNAMES Cell array of name for each column.  Format is the same as
%            for ROWNAMES.
%   SHOWVARS If true, show bounds on each variable.  (Default = false)

[nrows,ncols] = size(mip.A);

if nargin < 6
    showvars = false;
end

default_colnames = arrayfun(@(x) ['x(' num2str(x) ')'], 1:ncols, ...
                            'Uniform', false);
default_rownames = arrayfun(@(x) ['ROW_' num2str(x)], 1:nrows, ...
                            'Uniform', false);
                        
% default ranges and names
if nargin < 5 || isempty(colnames)
    colnames = default_colnames;
end
if nargin < 4 || isempty(rownames)
    rownames = default_rownames;
end
if nargin < 3 || isempty(colidxs)
    colidxs = 1 : ncols;
elseif isa(colidxs,'logical')
    colidxs = find(colidxs);
end
if nargin < 2 || isempty(rowidxs)
    rowidxs = 1 : nrows;
elseif isa(rowidxs,'logical')
    rowidxs = find(rowidxs);
end

% expand the names if an incomplete list was given
colnames = zip_names(default_colnames,colnames);
rownames = zip_names(default_rownames,rownames);

% show objective
fprintf('\n\n----- Objective -----\n');
fprintf('z = ');
show_coef_list(mip.c);

% show constraints
fprintf('\n\n----- Constraints -----\n');
for r = 1 : length(rowidxs)
    row = rowidxs(r);
    fprintf('%s:  ',rownames{row});
    show_coef_list(mip.A(row,:));
    switch mip.ctypes(r)
        case '<'
            fprintf(' <= ');
        case '>'
            fprintf(' >= ');
        case '='
            fprintf(' = ');
    end
    fprintf('%g\n',mip.b(r));
end

if showvars
    % show variable bounds
    maxlength = max( cellfun(@(x) length(x), colnames) );
    fmt = ['%' num2str(maxlength + 2) 's'];
    fprintf('\n\n----- Variable bounds -----\n');
    for i = 1 : length(colidxs)
        fprintf(fmt,[colnames{colidxs(i)} ':']);
        fprintf('  %s',mip.vartypes(colidxs(i)));
        fprintf('  [%g,%g]\n', [mip.lb(colidxs(i)) mip.ub(colidxs(i))]);
    end
end


function [zipped] = zip_names(default,given)
    zipped = default;
    for j = 1 : length(given)
        if ~isempty(given{j})
            zipped{j} = given{j};
        end
    end
end
    
function show_coef_list(constraint)
    nonzeros = find(constraint(colidxs));
    for j = 1 : length(nonzeros)
        col = nonzeros(j);
        coef = constraint(col);
        if coef < 0
            if j > 1
                coef_str = ' - ';
            else
                coef_str = '-';
            end
        else
            if j > 1
                coef_str = ' + ';
            else
                coef_str = '';
            end
        end
        if abs(coef) ~= 1
            coef_str = [coef_str num2str(abs(coef)) '*'];
        end
        fprintf('%s%s',coef_str,colnames{col});
    end
end

end
    