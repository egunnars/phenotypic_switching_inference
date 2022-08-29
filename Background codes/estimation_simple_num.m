function params = estimation_simple_num(data,N,T)

K = size(data,2);
ntimes = size(data,3);
nreplicates = size(data,4);

params = zeros(K,K+1);

v = 0;
invN = inv(N'*N);
for l=1:ntimes
    for r=1:nreplicates
        D = data(:,:,l,r);
        if sum(eig(invN*N'*D)<0)>0
        else
            Aest = logm(invN*N'*D)/T(l);
            params(:,2) = params(:,2) + diag(Aest) + sum(Aest-diag(diag(Aest)),2);
            for i=1:K
               params(i,3:end) = params(i,3:end) + Aest(i,[1:i-1,i+1:K]);
            end
            v = v+1;
        end
    end
end
params = params/v;

end