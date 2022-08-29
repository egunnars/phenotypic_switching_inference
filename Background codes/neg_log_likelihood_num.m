function value = neg_log_likelihood_num(param,scale,C,N,data,dead,T)

param = param.*scale;

ninitial = size(data,1);
ntimes = size(data,3);
nreplicates = size(data,4);

[b,A] = params_extract_num(param,dead);

if ~isempty(dead)
   for i=1:ninitial
       C{i} = [[C{i},zeros(size(C{i},1),1)];[zeros(1,size(C{i},2)),1]];
   end
   N = [N,zeros(size(N,1),1)];
end

value = 0;
for i=1:ninitial
    for l=1:ntimes
        mean_vector = N(i,:)*expm(T(l)*A);
        SM = C{i}'*covmatrix_num(N(i,:),A,b,T(l),dead)*C{i};
        invSM = inv(SM);
        for r=1:nreplicates
            if ~isempty(dead)
                dminusmu = ([data(i,:,l,r),dead(i,1,l,r)] - mean_vector)*C{i};
            else
                dminusmu = (data(i,:,l,r) - mean_vector)*C{i};
            end
            value = value + dminusmu*invSM*dminusmu';
        end
        value = value + nreplicates*log(det(SM));
    end
end
    
end