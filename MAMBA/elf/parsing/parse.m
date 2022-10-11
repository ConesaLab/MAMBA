function [expression] = parse(tokens,weak)

if nargin < 2,  weak = false; end

previous = [];
full_expr = false;
while tokens.is_another && ~(full_expr && weak)
    next = tokens.pop();
    if next.is_lparen
        previous = parse(tokens);
        after = tokens.peek();
        assert( after.is_rparen, 'syntax error: ) expected' );
        tokens.pop();
        full_expr = true;
    elseif next.is_rparen
        tokens.push(next);
        break;
    elseif next.is_and || next.is_or
        assert( ~isempty(previous), 'syntax error: nothing on stack' );
        exp = expr();
        if next.is_and
            exp.AND = true;
        else
            exp.OR = true;
        end
        exp.lexpr = previous;
        exp.rexpr = parse(tokens,true);
        previous = exp;
        full_expr = true;
    else
        exp = expr();
        exp.id = next.value;
        previous = exp;
        full_expr = true;
    end
end

expression = previous;
