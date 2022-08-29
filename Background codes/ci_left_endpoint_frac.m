function [ci_lb,ci_yvalues_lb,ci_lb_check] = ci_left_endpoint_frac(data,C,Ntot,alpha,T,mle,scale,optlikelihood,Aineq,bineq,Aeq,beq,lb,ub,i,j,quant,tol)

ci_lb_check = 0;

[Aineq,bineq,Aeq,beq,lb,ub] = addscaleconstraints_frac(Aineq,bineq,Aeq,beq,lb,ub,scale);
mle = mle./scale;

options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',50000,'MaxIter',10000,'ScaleProblem',true);
sol = fmincon(@(param) param(i,j),mle,Aineq,bineq,Aeq,beq,lb,ub,@(param) constraints_frac(param,scale,C,Ntot,alpha,data,T,optlikelihood,quant),options);

ci_yvalues_lb = neg_log_likelihood_frac(sol,scale,C,Ntot,alpha,data,T);
sol = sol.*scale;
ci_lb = sol(i,j);

if (abs(ci_yvalues_lb - optlikelihood - quant)/(optlikelihood+quant)<tol) || (ci_lb == lb(i,j) && ci_yvalues_lb < optlikelihood + quant + tol)
    ci_lb_check = 1;
end

end

function [c,ceq] = constraints_frac(param,scale,C,Ntot,alpha,data,T,optvalue,quant)
    c = neg_log_likelihood_frac(param,scale,C,Ntot,alpha,data,T) - optvalue - quant;
    ceq = 0;
end