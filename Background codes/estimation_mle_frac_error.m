function [params, error, scale_mle, scale_error, optvalue, feasible] = estimation_mle_frac_error(ratios,C,Ntot,alpha,T,x0,scale,Aineq,bineq,Aeq,beq,lb,ub,tol)

K = size(alpha,2);
[Aineq,bineq,Aeq,beq,lb,ub] = addscaleconstraints_frac(Aineq,bineq,Aeq,beq,lb,ub,scale);
x0 = x0./scale;

Aineq_new = [Aineq,zeros(size(Aineq,1),1)];
if ~isempty(Aeq)
   Aeq_new = [Aeq,zeros(size(Aeq,1),1)];
else
   Aeq_new = Aeq;
end
scale_new = [reshape(scale,1,K*(K+1)),0.01];
lb_new = [reshape(lb,1,K*(K+1)),0];
ub_new = [reshape(ub,1,K*(K+1)),1/scale_new(K*(K+1)+1)];
x0_new = [reshape(x0,1,K*(K+1)),0.01/scale_new(K*(K+1)+1)];

options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',50000,'MaxIter',10000,'ScaleProblem',true);
[params_error, optvalue] = fmincon(@(param) neg_log_likelihood_frac_error(param,scale_new,C,Ntot,alpha,ratios,T),x0_new,Aineq_new,bineq,Aeq_new,beq,lb_new,ub_new,[],options);

params = reshape(params_error(1:K*(K+1)),K,K+1);

feasible = 0;
feasible = feasible + sum(Aineq*reshape(params,K*(K+1),1) - bineq > tol);
if ~isempty(Aeq)
    feasible = feasible + sum(abs(Aeq*reshape(params,K*(K+1),1) - beq) > tol);
end
feasible = feasible + sum(sum(lb > params)) + sum(sum(params > ub));
if feasible == 0
    feasible = 1;
else
    feasible = 0;
end

params = params.*scale;
error = params_error(K*(K+1)+1)*scale_new(K*(K+1)+1);

scale_mle = 10.^floor(log10(abs(params)));
scale_mle(find(scale_mle == 0)) = 1;
scale_error = 10.^floor(log10(abs(error)));
scale_error(find(scale_error == 0)) = 1;

end