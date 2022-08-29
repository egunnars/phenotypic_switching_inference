function [x0, scale] = initialguess_frac(x0_in,lb,ub)

K = size(lb,1);
if ~isempty(x0_in)
    x0 = x0_in;
    x0(1,1) = lb(1,1)+rand*(ub(1,1)-lb(1,1));
    x0(1,2) = min(max(lb(1,2),-x0(1,1)),ub(1,2)) + (ub(1,2)-min(max(lb(1,2),-x0(1,1)),ub(1,2)))*rand;
    for k=2:K
       x0(k,1) = min(max(-x0(1,2)-x0(k,2),lb(k,1)),ub(k,1)) + (ub(k,1)-min(max(-x0(1,2)-x0(k,2),lb(k,1)),ub(k,1)))*rand;
    end
else
   [~,~,x0] = params_generation(K,1);
end

scale = 10.^floor(log10(abs(x0)));
scale(find(scale == 0)) = 1;

end
