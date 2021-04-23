% Creates the plot for figure 2(a)
clear
load nvars_comb

% Show the first 8 for convergence
for i = 2:8
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)],'-'); 
end
% Show 18 for the 'converged' value, could show extrapolated value, but
% they amount to essentially the same thing for this scale.
hold on
plot([0;sig_est(:,1)],[2;sig_est(:,18)],'k-'); 
xlim([0,5])

saveas(gcf,'fig2.pdf')
