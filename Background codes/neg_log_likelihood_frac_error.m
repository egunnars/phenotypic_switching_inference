function value = neg_log_likelihood_frac_error(param_error,scale,C,Ntot,alpha,ratios,T)

K = size(alpha,2);
param_error = param_error.*scale;
param = reshape(param_error(1:K*(K+1)),K,K+1);
error = param_error(K*(K+1)+1);

ninitial = size(alpha,1);
ntimes = size(ratios,3);
nreplicates = size(ratios,4);

[lambda_rel,d,lambda_1,A_rel] = params_extract_frac(param);

value = 0;
B = eye(K,K-1);
I = eye(K,K);
for i=1:ninitial
    for l=1:ntimes
        mean_vector = alpha(i,:)*expm(T(l)*A_rel)/sum(alpha(i,:)*expm(T(l)*A_rel));
        SM = C{i}'*B'*(Ntot(i)^(-1)*covmatrix_frac(alpha(i,:),A_rel,lambda_rel,d,lambda_1,T(l))+error^2*I)*B*C{i};
        invSM = inv(SM);
        for r=1:nreplicates
            dminusmu = (ratios(i,1:end-1,l,r) - mean_vector(1:end-1))*C{i};
            value = value + dminusmu*invSM*dminusmu';
        end
        value = value + nreplicates*log(det(SM));
    end
end

end