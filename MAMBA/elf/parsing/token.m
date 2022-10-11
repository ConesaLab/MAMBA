classdef token
    properties (Dependent)
        is_lparen
        is_rparen
        is_id
        is_and
        is_or
    end
    
    properties
        value
    end
    
    methods
        function [obj] = token(str)
            obj.value = str;
        end
        
        function [tf] = get.is_lparen(obj)
            tf = strcmp(obj.value,'(');
        end
        function [tf] = get.is_rparen(obj)
            tf = strcmp(obj.value,')');
        end
        function [tf] = get.is_and(obj)
            tf = strcmp(obj.value,'&');
        end
        function [tf] = get.is_or(obj)
            tf = strcmp(obj.value,'|');
        end
        function [tf] = get.is_id(obj)
            tf = ~any([obj.is_lparen ...
                       obj.is_rparen ...
                       obj.is_and ...
                       obj.is_or]);
        end
    end
end

