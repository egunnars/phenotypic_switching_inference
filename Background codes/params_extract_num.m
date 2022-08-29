function [b,A] = params_extract_num(param,dead)

K = size(param,1);
b = param(:,1)';
lambda = param(:,2)';
nu = param(:,3:end);

nu_augm = zeros(K,K);
for i=1:K
    nu_augm(i,[1:i-1,i+1:K]) = nu(i,:);
end
nu_augm = nu_augm - diag(sum(nu_augm,2));
A = diag(lambda) + nu_augm;

if ~isempty(dead)
   A = [[A,b'-lambda'];zeros(1,K+1)];
end

end