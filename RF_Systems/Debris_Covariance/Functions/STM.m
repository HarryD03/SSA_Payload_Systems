function [dSTM_dt_flat] = STM(t,x,A)
    mu = 3.986e14;


    STM = reshape(x(1:end), 6,6);           %Get the STM from the inputted vector
    dSTM_dt = A*STM;                        %
    dSTM_dt_flat = reshape(dSTM_dt, 36, 1); 
 

end

