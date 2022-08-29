function [x0, scale] = initialguess_num(dead,N,T,x0_in,lb,ub)

ninitial = size(N,1);
K = size(N,2);

if ~isempty(x0_in)
    x0 = x0_in;
    if ~isempty(dead)
        D = zeros(ninitial,K);
        c = zeros(ninitial,1);
        for i=1:ninitial
           for j=1:K
               D(i,j) = N(i,j)*(1/x0_in(j,2))*(exp(x0_in(j,2)*T(1))-1);
           end
           c(i) = mean(dead(i,1,1,:));
        end
        for j=1:K
        x0(:,1) = x0_in(:,2) + inv(D'*D)*D'*c;
        end
    else
        x0(:,1) = max(min(abs(x0(:,2)./rand(K,1)),ub(:,1)),lb(:,1));
    end
else
   [~,x0] = params_generation(K,1);
end

scale = 10.^floor(log10(abs(x0)));
scale(find(scale == 0)) = 1;

end