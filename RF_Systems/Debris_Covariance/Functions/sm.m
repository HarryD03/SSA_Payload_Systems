function [A] = sm(x)
%Purpose: Define the state matrix of the 2 body problem it is ECI centred 

    mu = 3.986e5;
    
    %Extract the state
    rx = x(1);
    ry = x(2);
    rz = x(3);

    r0_norm = norm(x(1:3));
    coef1 = 3*mu/(2*r0_norm^5);
    coef2 = mu/(2*r0_norm^3);
    
    %Form 2BP dynamics matrix
    A_21 = [coef1*rx^2 - coef2, coef1*rx*ry, coef1*rx*rz;
            coef1*rx*ry, coef1*ry^2 - coef2, coef1*ry*rz;
            coef1*rx*rz, coef1*ry*rz, coef1*rz^2-coef1];

    A = [zeros(3), eye(3); A_21, zeros(3)];

end