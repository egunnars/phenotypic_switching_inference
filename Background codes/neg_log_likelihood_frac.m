function value = neg_log_likelihood_frac(param,scale,C,Ntot,alpha,ratios,T)

param = param.*scale;

ninitial = size(ratios,1);
K = size(ratios,2);
ntimes = size(ratios,3);
nreplicates = size(ratios,4);

[lambda_rel,d,lambda_1,A_rel] = params_extract_frac(param);

value = 0;
B = eye(K,K-1);
for i=1:ninitial
    for l=1:ntimes
        mean_vector = alpha(i,:)*expm(T(l)*A_rel)/sum(alpha(i,:)*expm(T(l)*A_rel));
        SM = Ntot(i)^(-1)*C{i}'*B'*covmatrix_frac(alpha(i,:),A_rel,lambda_rel,d,lambda_1,T(l))*B*C{i};
        invSM = inv(SM);
        for r=1:nreplicates
            dminusmu = (ratios(i,1:end-1,l,r) - mean_vector(1:end-1))*C{i};
            value = value + dminusmu*invSM*dminusmu';
        end
        value = value + nreplicates*log(det(SM));
    end
end

end