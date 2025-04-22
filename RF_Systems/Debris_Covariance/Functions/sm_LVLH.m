function [A] = sm_LVLH(a)
    
    mu = 3.986e5;

    n = sqrt(mu/a^3);
    A_21 = [(3*n^2) 0 0;
            0 0 0;
            0 0 -n^2];
    A_22 = [0 2*n 0;
         -2*n 0 0;
            0 0 0];
    A = [zeros(3,3), eye(3);
         A_21,  A_22];
end