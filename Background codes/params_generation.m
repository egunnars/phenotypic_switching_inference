function [params_true, params_true_num, params_true_frac, scale_true_num, scale_true_frac, N, Ntot, alpha] = params_generation(K,nparamregimes)

params_true = zeros(K,K+1,nparamregimes);
params_true_num = zeros(K,K+1,nparamregimes);
params_true_frac = zeros(K,K+1,nparamregimes);

scale_true_num = zeros(K,K+1,nparamregimes);
scale_true_frac = zeros(K,K+1,nparamregimes);

N = zeros(K,K,nparamregimes);
Ntot = zeros(1,K,nparamregimes);
alpha = zeros(K,K,nparamregimes);

I = eye(K,K);

for reg=1:nparamregimes
    cond = 0;
    while cond == 0
        params_true(:,:,reg) = [rand(K,2),10.^(-3+2*rand(K,K-1))];
        params_true_num(:,:,reg) = params_true(:,:,reg);
        params_true_num(:,2,reg) = params_true(:,1,reg) - params_true(:,2,reg);

        if sum(sum(params_true_num(:,2,reg)<0)) == K || sum(sum(abs(params_true_num(:,1:2,reg))<0.01)) > 0
        else
            cond = 1;
            
            params_true_frac(:,:,reg) = params_true_num(:,:,reg);
            params_true_frac(2:end,2,reg) = params_true_frac(2:end,2,reg)-params_true_frac(1,2,reg);
            params_true_frac(:,1,reg) = params_true(:,2,reg);

            scale_true_num = 10.^floor(log10(abs(params_true_num)));
            scale_true_num(find(scale_true_num == 0)) = 1;   
            
            scale_true_frac = 10.^floor(log10(abs(params_true_frac)));
            scale_true_frac(find(scale_true_frac == 0)) = 1;
        end
    end
    if min(min(log10(params_true(:,:,reg)))) >= -2
        N(:,:,reg) = 10^3*I;
    else if min(min(log10(params_true(:,:,reg)))) >= -3
        N(:,:,reg) = 10^4*I;
        else
            N(:,:,reg) = 10^5*I;
        end
    end
    Ntot(1,:,reg) = sum(N(:,:,reg),2);
    alpha(:,:,reg) = N(:,:,reg)./Ntot(1,:,reg)';
end

end