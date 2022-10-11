function [and_lists] = make_dnf(ex)

ex_list = {ex};

while true
    new_list = {};
    N = length(ex_list);
    for i = 1 : N
        exi = ex_list{i};
        if has_or(exi)
            new_list{end+1} = split_by_or(exi,'left');
            new_list{end+1} = split_by_or(exi,'right');
        else
            new_list{end+1} = exi;
        end
    end
    ex_list = new_list;
    
    if ~any(cellfun(@(x) has_or(x), ex_list))
        break;
    end
end

and_lists = cellfun(@(x) x.get_atoms, ex_list, 'UniformOutput', false);


function [tf] = has_or(ex)
    tf = ex.is_junction ...
            && (ex.OR || has_or(ex.lexpr) || has_or(ex.rexpr));

        
function [sub_ex] = split_by_or(ex,dir)
    if nargin < 2,  dir = 'left'; end
    
    if ~has_or(ex)
        sub_ex = ex;
    elseif ex.OR
        switch dir
            case 'left'
                sub_ex = ex.lexpr;
            case 'right'
                sub_ex = ex.rexpr;
        end
    else  % AND
        sub_ex = ex;
        if has_or(ex.lexpr)
            switch dir
                case 'left'
                    sub_ex.lexpr = split_by_or(ex.lexpr,'left');
                case 'right'
                    sub_ex.lexpr = split_by_or(ex.lexpr,'right');
            end
        else
            switch dir
                case 'left'
                    sub_ex.rexpr = split_by_or(ex.rexpr,'left');
                case 'right'
                    sub_ex.rexpr = split_by_or(ex.rexpr,'right');
            end
        end
    end



                
