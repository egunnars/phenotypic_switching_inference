function [lambda_rel,d,lambda_1,A_rel] = params_extract_frac(param)

K = size(param,1);
d = param(:,1);
lambda_rel = param(2:end,2);
lambda_1 = param(1,2);
nu = param(:,3:end);

nu_augm = zeros(K,K);
for i=1:K
    nu_augm(i,[1:i-1,i+1:K]) = nu(i,:);
end
nu_augm = nu_augm - diag(sum(nu_augm,2));
A_rel = diag([0;lambda_rel]) + nu_augm;

end