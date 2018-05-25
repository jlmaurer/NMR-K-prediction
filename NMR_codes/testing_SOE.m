% test SOE

name = 'dpnmr_larned_lwph'; 
% [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp] = loadnmrdata(name); 
[d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(name); 

lse = log10(SumEch); 
lse3 = log10(SumEch_3s); 
lsetw = log10(SumEch_twm); 
lsetw3 = log10(SumEch_twm_3s); 

% X = [lp,lt, lse, lse3, lsetw, lsetw3]; 
X = [lt, lse, lsetw3]; 

LM_all = stepwiselm(X, kk)

% [U, S, V] = svd(X); 
% figure; plot(diag(S))
% 
% beta = X\kk
% 

