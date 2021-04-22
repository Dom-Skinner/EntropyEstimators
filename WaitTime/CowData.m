% This script contains everything needed to recreate fig 2c and compute
% entropy estimates
%% Cow experiment 1
probA = readtable('~/Dropbox (MIT)/NonEqDetection/Data/lying_exp_1.txt').Var1;
% probA extracted from fig 1 (b) of Tolkamp et al.
probB = readtable('~/Dropbox (MIT)/NonEqDetection/Data/stand_exp_1.txt').Var1;
% probA extracted from fig 2 (c) of Tolkamp et al. (hence log conversion
% later)
histogram('BinEdges',linspace(0,4,length(probA)+1),'BinCounts',probA,...
    'Normalization','pdf','DisplayStyle', 'stairs')


X_A = linspace(0,4,length(probA)+1);
X_A = 0.5*(X_A(2:end) + X_A(1:end-1));
tA = X_A*probA/sum(probA);
tA2 = (X_A.^2)*probA/sum(probA);

X_B = linspace(0,12,length(probB)+1);
X_B = 0.5*(X_B(2:end) + X_B(1:end-1));
tB = exp(X_B/60)*probB/sum(probB);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For experiment 1, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)


%% Cow experiment 2
probA = readtable('~/Dropbox (MIT)/NonEqDetection/Data/lying_exp_2.txt').Var1;
probB = readtable('~/Dropbox (MIT)/NonEqDetection/Data/stand_exp_2.txt').Var1;
hold on
histogram('BinEdges',linspace(0,4,length(probA)+1),'BinCounts',probA,...
    'Normalization','pdf','DisplayStyle', 'stairs')

X_A = linspace(0,4,length(probA)+1);
X_A = 0.5*(X_A(2:end) + X_A(1:end-1));
tA = X_A*probA/sum(probA);
tA2 = (X_A.^2)*probA/sum(probA);

X_B = linspace(0,12,length(probB)+1);
X_B = 0.5*(X_B(2:end) + X_B(1:end-1));
tB = exp(X_B/60)*probB/sum(probB);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For experiment 2, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)

%% Cow experiment 2
probA = readtable('~/Dropbox (MIT)/NonEqDetection/Data/lying_exp_3.txt').Var1;
probB = readtable('~/Dropbox (MIT)/NonEqDetection/Data/stand_exp_3.txt').Var1;
hold on
histogram('BinEdges',linspace(0,4,length(probA)+1),'BinCounts',probA,...
    'Normalization','pdf','DisplayStyle', 'stairs')
xlabel('Time spent lying (h)')
ylabel('Probability density')
ylim([0,1.6])

X_A = linspace(0,4,length(probA)+1);
X_A = 0.5*(X_A(2:end) + X_A(1:end-1));
tA = X_A*probA/sum(probA);
tA2 = (X_A.^2)*probA/sum(probA);

X_B = linspace(0,12,length(probB)+1);
X_B = 0.5*(X_B(2:end) + X_B(1:end-1));
tB = exp(X_B/60)*probB/sum(probB);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For experiment 3, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)


%saveas(gcf,'lying_hist_' + string(i) + '.pdf')