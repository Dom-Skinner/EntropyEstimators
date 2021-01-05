% A tiny script to correct some of the values that have not found their
% true minimum.
clear
load nvars_comb

ntrials = 4;
hess=true;
gs = false;
rand_init=true;
ct = 1;


for i = 3:size(sig_est,2)
    for j = 1:size(sig_est,1)
        if sig_est(j,i) > sig_est(j,i-1)*(1+4e-3)
            sig_est = run_trials(ntrials,sig_est,i,j,hess,rand_init);
            ct = ct+1;
            if ct == 3
                clearvars -except sig_est
                save('nvars_comb')
                disp(sum(sum(sig_est(:,3:end) > sig_est(:,2:end-1)*(1+1e-6))))
                return
            end
        end
    end
end
