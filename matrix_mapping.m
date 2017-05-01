function [Fhat,beta,SSE] = matrix_mapping(F, S, K);

% Estimate structure to function mapping: F = f(S) = C + C_0*S^0 + C_1*S^1 + 
% f can be estmated by either pinv or nlinfit, both have a problem with
% rank deficiency, use Tikhonov regularization
%
% Inputs: 
%   F: function matrix,  NxN
%   S: structure matrix, NxN
%   K: Number of terms in Taylor series, with maximum of N-1;
%      K=1, just the first 3 terms, i.e., only direct edges used 
% Outputs:
%   Fhat: estimated F from S, NxN
%   beta: coefficients of polynominal, [C C_0 C_1 C_2 ...C_K], (K+2)x1
%   SSE: sum of square errors
% 
% [Fhat,beta,SSE] = matrix_mapping(F,S,5);
%
% Hualou Liang at Drexel University, 2015
%

N=size(S,1);
Isubdiag = find(tril(ones(N))); 

Y = F(Isubdiag); 
X = ones(size(Isubdiag)); % all-one matrix to represent common input
Xorg = X;

normfact = zeros(K+1,1); 
for i=0:K 
    spow = S^i;
    
    Xorg = [Xorg spow(Isubdiag)]; 
    
    normfact(i+1) = max(abs(spow(:)));
    spow = spow/normfact(i+1);
    X = [X spow(Isubdiag)];
end

%%%  gcv via Tikhonov regularization 
[Ua,s,Va] = csvd(X); 
lambda_gcv = gcv(Ua,s,Y);
beta = tikhonov(Ua,s,Va,Y,lambda_gcv);

beta(2:end) = beta(2:end)./normfact;

% construct the estimated matrix 
Ybar = zeros(N,N); 
Ydat = Xorg*beta; 
Ybar(Isubdiag) = Ydat;
Fhat = Ybar + Ybar' - diag(diag(Ybar)); 

% get estimation quality of Y
dof = N*(N-1)/2 + N - 1;
SSE_fro = norm(F - Fhat,'fro'); % Frobenius norm 
SSE = SSE_fro^2./dof;
return




