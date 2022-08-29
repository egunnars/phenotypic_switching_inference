function [data, ratios, dead] = data_generation(T,N,K,ninitial,ntimes,nreplicates,niter,params_true_num,time_lapse)

nparamregimes = size(params_true_num,3);
data = zeros(ninitial,K,ntimes,nreplicates,nparamregimes,niter);
ratios = zeros(ninitial,K,ntimes,nreplicates,nparamregimes,niter);
dead = zeros(ninitial,1,ntimes,nreplicates,nparamregimes,niter);

for reg=1:nparamregimes
    for k=1:niter
        for i=1:ninitial
            for r=1:nreplicates
                if time_lapse == 0
                    for l=1:ntimes
                        [num, num_dead] = Samplepath_new(params_true_num(:,:,reg),N(i,:,reg),T(l));
                        data(i,:,l,r,reg,k) = num;
                        ratios(i,:,l,r,reg,k) = data(i,:,l,r,reg,k)/sum(data(i,:,l,r,reg,k));
                        dead(i,:,l,r,reg,k) = num_dead;
                    end
                else
                    [num, num_dead] = Samplepath_new_timelapse(params_true_num(:,:,reg),N(i,:,reg),T);
                    data(i,:,:,r,reg,k) = num;
                    for l=1:ntimes
                        ratios(i,:,l,r,reg,k) = data(i,:,l,r,reg,k)/sum(data(i,:,l,r,reg,k));
                    end
                    dead(i,:,:,r,reg,k) = num_dead;
                end
            end
        end
    end
end

end