% This script produces some of the asymptotic figures from the
% supplementary. First we check the convergence to the asymptotic
% expression for fixed N, sig-> infinity, then we check the convergence for
% fixed sig as N -> infinity.

%% Check asymptotics as sig->infinity 
clear
load large_s
for i = 2:6
    hold on
    plot(sig_est(:,1),(sig_est(:,i)-1 - 1/i),'o');    
    plot(sig_est(:,1),exp(-2*sig_est(:,1)/(i+1))*(i^2+i-2)/(i^2) * 2^(2/(i+1)),'-');

end
xlim([10,40])
ylim([1e-6,0.1])
set(gca, 'YScale', 'log')
xlabel('Entropy production rate')
ylabel('Difference between t_2 and asymptote,  1 + 1/N')
saveas(gcf,'large_sig.pdf')

%% Check asymptotics as N -> infinity
clear
load n_conv
for k = 1:size(sig_est,1)
s = s_extrap(sig_est(k,12:18),12:18);
hold on
plot(2:18,abs(sig_est(k,2:18)-s))
end
set(gca, 'YScale', 'log')
xlabel('Number of internal states')
ylabel('Estimated difference between finite N and extrapolated value')
saveas(gcf,'converge_n.pdf')