classdef expr
    properties (Dependent)
        is_junction
        is_atom
    end
    
    properties
        AND = false
        OR = false
        lexpr
        rexpr
        
        id = ''
        
        NULL = false
    end
    
    methods
        function [tf] = get.is_junction(obj)
            tf = (obj.AND | obj.OR);
        end
        
        function [tf] = get.is_atom(obj)
            tf = ~obj.is_junction;
        end
        
        function [atoms] = get_atoms(obj)
            if obj.is_junction
                atoms = [obj.lexpr.get_atoms obj.rexpr.get_atoms];
            else
                atoms = {obj.id};
            end
        end
        
        function display(obj,indent)
            BASE_INDENT = '   ';
            if nargin < 2,  indent = ''; end
            
            if obj.is_junction
                if obj.AND
                    fprintf('%sAND--\n', indent);
                elseif obj.OR
                    fprintf('%sOR--\n', indent);
                end
                obj.lexpr.display([indent BASE_INDENT]);
                obj.rexpr.display([indent BASE_INDENT]);
            else
                fprintf('%s%s\n', indent, obj.id);
            end
        end
    end
end

                