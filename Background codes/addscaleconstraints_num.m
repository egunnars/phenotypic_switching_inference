function [Aineq_ret,bineq_ret,Aeq_ret,beq_ret,lb_ret,ub_ret] = addscaleconstraints_num(Aineq,bineq,Aeq,beq,lb,ub,scale)

K = size(scale,1);

Aineq_def = [-eye(K,K),eye(K,K),zeros(K,K^2-K)]; 
bineq_def = zeros(K,1);
if ~isempty(Aineq)
    for n=1:size(Aineq,3)
        Aineq_def = [Aineq_def; reshape(Aineq(:,:,n),1,K*(K+1))];
        bineq_def = [bineq_def; bineq(n)];
    end
end
Aineq_ret = Aineq_def;
Aineq_ret = Aineq_ret.*reshape(scale,1,K*(K+1));
bineq_ret = bineq_def;

Aeq_def = [];
beq_def = [];
if ~isempty(Aeq)
    for n=1:size(Aeq,3)
        Aeq_def = [Aeq_def; reshape(Aeq(:,:,n),1,K*(K+1))];
        beq_def = [beq_def; beq(n)];
    end
    Aeq_def = Aeq_def.*reshape(scale,1,K*(K+1));
end
Aeq_ret = Aeq_def;
beq_ret = beq_def;

lb_def = [zeros(K,1),-Inf*ones(K,1),zeros(K,K-1)];
if isempty(lb)
    lb_ret = lb_def./scale;
else
    lb_ret = max(lb_def,lb)./scale;
end
ub_ret = ub./scale;

end