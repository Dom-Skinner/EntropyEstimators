clear
load nvars_comb

for i = 2:10
    hold on
    plot([0;sig_est(:,1)],[2;sig_est(:,i)],'-');
end
plot([0,10],2 - [0,10]/4)
xlim([0,4])
saveas(gcf,'asympt_small.pdf')


for i = 2:6
    hold on
    plot(sig_est(:,1),(sig_est(:,i)-1 - 1/i),'-');    
    plot(sig_est(:,1),exp(-2*sig_est(:,1)/(i+1))*(i^2)/(i^2+i-2) * 2^(2/(i+1)),'o');

end
xlim([10,30])
set(gca, 'YScale', 'log')


clear
load n_conv
for k = 1:size(sig_est,1)
s = s_extrap(sig_est(k,12:18),12:18);
hold on
plot(2:18,abs(sig_est(k,2:18)-s))
end

set(gca, 'YScale', 'log')
saveas(gcf,'converge_n.pdf')