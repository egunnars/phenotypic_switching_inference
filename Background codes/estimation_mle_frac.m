function [params, scale_mle, optvalue, feasible] = estimation_mle_frac(ratios,C,Ntot,alpha,T,x0,scale,Aineq,bineq,Aeq,beq,lb,ub,tol)

[Aineq,bineq,Aeq,beq,lb,ub] = addscaleconstraints_frac(Aineq,bineq,Aeq,beq,lb,ub,scale);
x0 = x0./scale;

options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',50000,'MaxIter',10000,'ScaleProblem',true);
[params, optvalue] = fmincon(@(param) neg_log_likelihood_frac(param,scale,C,Ntot,alpha,ratios,T),x0,Aineq,bineq,Aeq,beq,lb,ub,[],options);

K = size(alpha,2);
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
scale_mle = 10.^floor(log10(abs(params)));
scale_mle(find(scale_mle == 0)) = 1;

end