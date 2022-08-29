function SM = covmatrix_num(n,A,b,t,dead)

K = size(A,1);
SM = zeros(K,K);
I = eye(K,K);

if ~isempty(dead)
    b = [b,0];
end

for j=1:K
     if n(j)>0
         e_j = I(j,:);
         mean_vector = e_j*expm(t*A);
         f = @(s) 2*expm((t-s).*A)'*diag(b.*(e_j*expm(s.*A)))*expm((t-s).*A);
         SM = SM + n(j)*(integral(f,0,t,'ArrayValued',true)+diag(mean_vector)-mean_vector'*mean_vector);
     end
end

end