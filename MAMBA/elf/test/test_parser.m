
str = 'a AND (b OR c) | "abc"';
disp(str);
t1 = parse_string(str)
fprintf('\n\n');

str = '(x && y | z) AND (b OR c) | "abc"';
disp(str);
t2 = parse_string(str)
fprintf('\n\n');

