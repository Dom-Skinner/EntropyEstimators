% This script makes plots for the active sensor example

%% First set up problem
WA_fun = @(wp,wm) [-(wp+wm), wp, 0, 0; ... 
                  wm, -(wp+wm), wp, 0; ... 
                  0, wm, -(wp+wm), wp; ...
                0, 0, wm, -(wm+wp)];
n = 4;
vec1 = ones(n,1);
pi_A = ones(n,1)/(n+1);
ent_rate = @(sig,wm) (1-wm) .*abs(log(abs(1 ./ wm))) - sig;

%% make plots of time spent bound distributions, i.e. fig 3b
sig_vals = [0.,1,2,3,4,5];
Wm = zeros(size(sig_vals));
for i = 1:length(sig_vals)
    f = @(wm) ent_rate(sig_vals(i),wm);
    Wm(i) = fzero(f,[eps 1]);
end

t = linspace(0,80,3200);
f = zeros(length(t),length(Wm));

for j = 1:length(Wm)
    %WA = WA_fun(exp(E_vals(j)/5),1);
    WA = WA_fun(1,Wm(j));
    
    for i = 1:length(t)
        f(i,j) = (pi_A')*(WA*WA* expm(WA*t(i))) * vec1;
    end
    K = - (pi_A')*WA * vec1;
    f(:,j) = f(:,j)/K;
end

plot(t,f(:,1),t,f(:,2),t,f(:,3),t,f(:,4),t,f(:,5))
xlim([0,10])
ylim([0,0.5])
%saveas(gcf,'ActiveSensor.pdf')

%% Make plots 3c. Load the data for the empirical testing as it takes a while to generate.
% To generate the data run script ...

sig_vals_fine = linspace(0.01,5,100);
Wm = zeros(size(sig_vals_fine));
for i = 1:length(sig_vals_fine)
    f = @(wm) ent_rate(sig_vals_fine(i),wm);
    Wm(i) = fzero(f,[eps 1]);
end

tA = zeros(length(Wm),1);
tB = zeros(length(Wm),1);
tA2 = zeros(length(Wm),1);

for j = 1:length(Wm)
    
    WA = WA_fun(1,Wm(j));
    K = - (pi_A')*WA * vec1;
    tA(j) = (pi_A')*vec1 / K;
    tA2(j) = -2* (pi_A')*(WA \ vec1) /K;
    tB(j) = 1/(1 + Wm(j));
end


t2_norm = tA2./tA.^2;
vq = gamma_precomp(t2_norm,true);


plot(sig_vals_fine,2*vq./(tA+tB),sig_vals_fine,sig_vals_fine)
load sensor_empirical

hold on
plot(sig_vals_fine,TUR_med)
s2 = [sig_vals_fine, fliplr(sig_vals_fine)];
inBetween = [TUR_3, fliplr(TUR_97)];
fill(s2, inBetween, 'g');


plot(sig_vals_fine,OPTIM_med)
inBetween = [OPTIM_3, fliplr(OPTIM_97)];
fill(s2, inBetween, 'r');
saveas(gcf,'ActiveSensorQuant.pdf')