%% Set the SNR threshold by the Neymann-Person Crtierion 
% Ref: Richards, M. A. Fundamentals of Radar Signal Processing. New York: McGraw-Hill, 2005, pp 298–336.

SNRdB = [3 6 9 12];

[Pd,Pfa] = rocsnr(SNRdB,SignalType="NonfluctuatingCoherent");

semilogx(Pfa,Pd)
grid on
xlabel("P_{fa}")
ylabel("P_d")
legend("SNR "+SNRdB+" dB",Location="northwest")
title("Receiver Operating Characteristic (ROC) Curves")
