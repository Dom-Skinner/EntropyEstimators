sig_vals_fine = linspace(0.01,5,100);
Wm = zeros(size(sig_vals_fine));
ent_rate = @(sig,wm) (1-wm) .*abs(log(abs(1 ./ wm))) - sig;
for i = 1:length(sig_vals_fine)
    f = @(wm) ent_rate(sig_vals_fine(i),wm);
    Wm(i) = fzero(f,[eps 1]);
end




TUR_med = zeros(size(Wm));
TUR_97 = zeros(size(Wm));
TUR_3 = zeros(size(Wm));
OPTIM_med = zeros(size(Wm));
OPTIM_97 = zeros(size(Wm));
OPTIM_3 = zeros(size(Wm));



T = 4000;
ntrials = 1000;
nexp = 200;

for i = 1:length(Wm)

qA_var_est = zeros(ntrials,1);
qA_mean= zeros(ntrials,1);
N_emp_mean = zeros(ntrials,1);
t2_norm = zeros(ntrials,1);
tau = zeros(ntrials,1);

for j  = 1:ntrials    
    [q_A,N_emp,TA_emp,TB_emp,TA2_emp] = traj_sample(T,Wm(i),nexp);
    qA_var_est(j) = var(q_A);
    qA_mean(j) = mean(q_A);
    N_emp_mean(j) = mean(N_emp);
    t2_norm(j) = mean(TA2_emp)/mean(TA_emp)^2;
    tau(j) = mean( TA_emp+TB_emp)/2;
end

TUR_ests = 8* (qA_mean.^2) .* (1 - qA_mean).^2 ./ (T.*qA_var_est)  - 2*N_emp_mean./T;
OPTIM_ests = gamma_precomp(t2_norm,true)./tau;

OPTIM_med(i) = median(OPTIM_ests);
OPTIM_3(i) = quantile(OPTIM_ests,0.025);
OPTIM_97(i) = quantile(OPTIM_ests,0.975);

TUR_med(i) = median(TUR_ests);
TUR_3(i) = quantile(TUR_ests,0.025);
TUR_97(i) = quantile(TUR_ests,0.975);

end

TUR_med(TUR_med < 0) = 0;
TUR_3(TUR_3 < 0) = 0;
TUR_97(TUR_97 < 0) = 0;
clearvars -except Wm OPTIM_med OPTIM_3 OPTIM_97 TUR_med TUR_3 TUR_97
save sensor_empirical

