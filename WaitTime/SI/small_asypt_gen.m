% Check the small value asymptotics
sig_est = zeros(20,5);
sig_est(:,1) = linspace(0.1,2,20);

ntrials = 20; % may need more than 20 trials for convergence.
hess=true;
gs = false;
rand_init=true;
rng('default')

for nvars = 2:5
    for i = 1:size(sig_est,1)
    sig_est = run_trials(ntrials,sig_est,nvars,i,hess,rand_init,gs);
    end
end

for i = 2:5
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)],'-o');
end
plot([0,2],2 - [0,2]/4,'k--') % This is the asymptotic function
xlabel('Entropy production rate')
ylabel('<t_A^2>')

save('small_asymp','sig_est')