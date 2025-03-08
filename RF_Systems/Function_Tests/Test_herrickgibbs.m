%% Unit test for herrick gibbs function 

%Test completed using p.467 Vallado
r1 = [3419.85564,6019.82602,2784.60022];
r2 = [2935.91195,6326.18324,2660.59584];
r3 = [2434.95202,6597.38674,2521.523112];

t1 = 0;
t2 = 60+16.48;
t3 = 120+33.04;
mue = 3.986e5;
[v2] = herrickgibbs(r1, r2, r3, t1, t2, t3,mue);