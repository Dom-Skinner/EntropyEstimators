% A tiny script to correct some of the values that have not found their
% true minimum.
clear
load nvars_comb

ntrials = 4;
hess=true;
gs = false;
rand_init=true;
def_idx = [];

for i = 3:size(sig_est,2)
    for j = 1:size(sig_est,1)
        if sig_est(j,i) > sig_est(j,i-1)*(1+4e-3)
            def_idx(end+1,1) = i;
            def_idx(end,2) = j;
        end
    end
end

for i = 1:size(def_idx,1)
            sig_est = run_trials(ntrials,sig_est,def_idx(i,1),def_idx(i,2),hess,rand_init,gs);
end
clearvars -except sig_est
save('nvars_comb')
disp(sum(sum(sig_est(:,3:end) > sig_est(:,2:end-1)*(1+1e-6))))
    
