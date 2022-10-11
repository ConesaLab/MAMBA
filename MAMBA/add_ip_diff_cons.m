function add_ip_diff_cons(i1,i2,di, roff, coff)
    % add constraints for di = x(i1) XOR x(i2)
    % TODO:  enable MILPs for continuous differences
    I1 = di + 1; % indicator variable 1
    I2 = di + 2; % indicator variable 2
    mip.A(roff+1,[I1 I2 di]) = [2 2 -4];
    mip.b(roff+1) = 3;
    mip.A(roff+2,[I1 I2 di]) = [-2 -2 4];
    mip.b(roff+2) = 1;
    mip.A(roff+3,[i1 i2 I1]) = [-1 -1 -3];
    mip.b(roff+3) = -2;
    mip.A(roff+4,[i1 i2 I1]) = [1 1 3];
    mip.b(roff+4) = 4;
    mip.A(roff+5,[i1 i2 I2]) = [1 1 -3];
    mip.b(roff+5) = 0;
    mip.A(roff+6,[i1 i2 I2]) = [-1 -1 3];
    mip.b(roff+6) = 2;
    
    mip.vartypes(coff+1:coff+3) = 'BBB';
    mip.ctypes(roff+1:roff+6) = '<<<<<<';
    mip.ub(coff+1:coff+3) = 1;
end