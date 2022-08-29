function SM = covmatrix_frac(alpha,A_rel,lambda_rel,d,lambda_1,t)

K = size(alpha,2);
SM = zeros(K,K);
I = eye(K,K);

for j=1:K
     if alpha(j)>0
         e_j = I(j,:);
         mean_vector = e_j*expm(t*A_rel);
         f = @(s) 2*exp(-lambda_1*s)*expm((t-s).*A_rel)'*diag((d'+lambda_1*ones(1,K)+[0;lambda_rel]').*(e_j*expm(s.*A_rel)))*expm((t-s).*A_rel);
         SM = SM + alpha(j)*(integral(f,0,t,'ArrayValued',true) + exp(-lambda_1*t)*diag(mean_vector) - mean_vector'*mean_vector);
     end
end
M_alpha = alpha*expm(t*A_rel)*ones(K,1);
Q_alpha = I-(1/M_alpha)*kron(ones(K,1),alpha)*expm(t*A_rel);
SM = (1/M_alpha^2)*Q_alpha'*SM*Q_alpha;

end