clear
load nvars_comb


ntrials = 10;
hess = true;
rand_init = true;
for nvars = 10:10
    for i = 10:14%1:size(sig_est,1)
    sig_est = run_trials(ntrials,sig_est,nvars,i,hess,rand_init);
    end
end

hold on
s = sig_est(:,1);
plot(s,sig_est(10:14,nvars),'o-')
plot(s,gamma_precomp(s,false))



for i = 2:18
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)],'-'); 
    %s = linspace(5,40,50);
    %plot(s,1 + 1/i + (1 + 1/i - 2/i^2)*2^(2/(i+1)) * exp(-2*s/(i+1)),'o');

end
xlim([0,5])

s = zeros(size(sig_est,1),1);
for k = 1:size(sig_est,1)
    s(k) = s_extrap(sig_est(k,13:19),13:19);
end
sig = sig_est(:,1);
plot(sig,s,sig,1+(5/3) ./sig)
plot(sig,s,sig,1+2 ./sig)
clearvars -except sig_est
save('nvars_comb')

s = linspace(2.1,40,400);
vq = gamma_precomp(s,false);
plot(s,vq)
hold on
plot(s,0.54*log(s)./s+ 1 ,'o')

plot(s,(vq-1).*s./log(s))
hold on
plot([1.2;1.4],[1.2;1.4]);
%ctol = 1e-7;
%otol = 1e-7;
%[x,fval,exitflag] = min_via_fmin(3, 15,ctol,otol,x);