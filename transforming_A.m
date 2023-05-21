function [A_new] = transforming_A(A,C_0,istrans)
    if istrans
    
        A_new = 2*C_0/pi*atan(A);
    else
        A_new = A;
    end
end