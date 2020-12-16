clear

sig_est = zeros(100,11);
sig_est(:,1) = linspace(0.1,40,size(sig_est,1));
load nvars_comb


plot(sig_est(:,1),sig_est(:,11),'-o');

ntrials = 15;
hess = true;
rand_init = true;
for nvars = 20:20
    for i = 16:16%1:size(sig_est,1)
    sig_est = run_trials(ntrials,sig_est,nvars,i,hess,rand_init);
    end
end

for i = 2:11
    hold on
    plot(sig_est(:,1),sig_est(:,i),'-');
    %s = linspace(5,40,20);
    %plot(s,1 + 1/i + (1 + 1/i - 2/i^2)*2^(2/(i+1)) * exp(-2*s/(i+1)),'o');

end

%s11 = sig_est(:,11);
%load nvars_comb
%sig_est(:,11) = s11;
clearvars -except sig_est
save('nvars_comb')

s = linspace(0.1,40,400);
vq = gamma_precomp(s);
plot(s,vq)
%ctol = 1e-7;
%otol = 1e-7;
%[x,fval,exitflag] = min_via_fmin(3, 15,ctol,otol,x);