function [SSEs, corvals, Betti0s,SSE_bs] = topo_mapping(F,S);

% [SSEs, corvals,Betti0s,SSE_bs] = topo_mapping(F,S);
%
% F, S: functional and structural square matrices
% SSEs: sum of square errors as a function of terms in Taylor series, excluidng [the
%   1st term corresponding to SSE between original signals], 
% corvals: correlation between function connectivity and its surroagte, the
%   1st term corresponding to corvals between original signals
% Betti0s: (N-1) x K, 65x10
% SSE_bs: SSE based on Bettis curves for each k, 1 x K
%
% Hualou Liang at Drexel University, 2015
%

N=size(S, 1);
% structure to function maping: F = f(S)
Isubdiag = find(tril(ones(N),-1)); 

corvals=[];
SSEs=[];
etas=[];

% barcode for target funct matrix
Bett0_target=barcode(1-F);

Betti0s=[]; aucs=[]; SSE_bs=[];
pLen = 10; % pLen - set path length of 10 arbitrarily
Fhats = zeros(N,N, pLen); 
for k=1:pLen
    [Fhat,beta,SSE] = matrix_mapping(F, S, k);  % beta
    SSEs = [SSEs SSE];

    % computing barcode based on persistent homology 
    Fhats(:,:,k) = Fhat;
    Betti0=barcode(1-Fhat); 
  
    % compute SSE using area under curve of Betti
    dBetti = (Betti0 - Bett0_target).^2;
    SSE_b = trapz([0; dBetti],  [N:-1:1]);
    Betti0s = [Betti0s Betti0];
    SSE_bs = [SSE_bs SSE_b];

    corval = corr(F(Isubdiag),Fhat(Isubdiag));
    corval_org = corr(F(Isubdiag),S(Isubdiag));
    

    fprintf('\n Number of terms, direct_corr, estimated_corrv: %d  %6.3f  %6.3f',k, corval_org, corval);
    corvals = [corvals corval];
end
fprintf('\n');

SSE_bs = SSE_bs ./N^2;

return



