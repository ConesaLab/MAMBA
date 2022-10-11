classdef stack < handle
    properties (Dependent)
        is_another
        is_empty
        length
    end
    properties (Dependent,Hidden)
        N
    end
    
    properties (SetAccess = private)
        values
    end
    
    methods
        function [obj] = stack(array)
            if nargin == 0
                obj.values = {};
            elseif isa(array,'cell')
                obj.values = array;
            else
                N = length(array);
                obj.values = cell(1,N);
                for i = 1 : N
                    obj.values{i} = array(i);
                end
            end
        end
        
        function [item] = pop(obj)
            assert( obj.is_another, 'stack is empty' );
            item = obj.values{end};
            obj.values = obj.values(1:end-1);
        end
        
        function push(obj,item)
            obj.values{end+1} = item;
        end
        
        function [item] = peek(obj)
            assert( obj.is_another, 'stack is empty' );
            item = obj.values{end};
        end
        
        function [N] = get.length(obj)
            N = length(obj.values);
        end
        function [N] = get.N(obj)
            N = obj.length;
        end
        function [tf] = get.is_another(obj)
            tf = obj.N > 0;
        end
        function [tf] = get.is_empty(obj)
            tf = ~obj.is_another;
        end
        
        function reverse(obj)
            obj.values = fliplr(obj.values);
        end
    end
end
            