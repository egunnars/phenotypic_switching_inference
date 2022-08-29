function [params_simple, optvalue] = estimation_simple_frac(ratios,alpha,Ntot,T,Aineq,bineq,Aeq,beq,lb,ub)

K = size(alpha,2);

if isequal(alpha,eye(K))
    nu0 = (ratios(:,:,1,1)-eye(K))/T(1);
    nu1 = zeros(K,K-1);
    for i=1:K
        nu1(i,:) = nu0(i,[1:i-1,i+1:K]);
        for j=1:K-1
           if nu1(i,j) == 0
              nu1(i,j) = 1/Ntot(i)/T(1);
           end
        end
    end
    x0 = [zeros(K,2), nu1];
    scale = 10.^floor(log10(abs(x0)));
    scale(find(scale == 0)) = 1;
else
    x0 = initialguess_frac([],lb,ub);
    scale = 10.^floor(log10(abs(x0)));
    scale(find(scale == 0)) = 1;
end

[Aineq,bineq,Aeq,beq,lb,ub] = addscaleconstraints_frac(Aineq,bineq,Aeq,beq,lb,ub,scale);
x0 = x0./scale;

options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',50000,'MaxIter',10000,'ScaleProblem',true);
[params_simple, optvalue] = fmincon(@(params) costfcn(params,ratios,alpha,T,scale),x0,Aineq,bineq,Aeq,beq,lb,ub,[],options);
params_simple = params_simple.*scale;

for i=1:K
    for j=1:K-1
       if params_simple(i,j+2) == 0
          params_simple(i,j+2) = 1/Ntot(i)/T(1);
       end
    end
end

end

function value = costfcn(params,ratios,alpha,T,scale)
    ninitial = size(alpha,1);
    K = size(alpha,2);
    ntimes = size(ratios,3);
    nreplicates = size(ratios,4);
    
    params = params.*scale;
    [~,~,~,A_rel] = params_extract_frac(params);

    value = 0;
    e = ones(K,1);
    for i=1:ninitial
        for l=1:ntimes
            for r=1:nreplicates
                value = value + norm(alpha(i,:)*expm(T(l)*A_rel)/(alpha(i,:)*expm(T(l)*A_rel)*e)-ratios(i,:,l,r))^2;
            end
        end
    end
end