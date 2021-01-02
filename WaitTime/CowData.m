% For Bmal1
i = 3;
X = readtable('~/Dropbox (MIT)/NonEqDetection/Data/lying_exp_' + string(i) + '.txt');
Y = linspace(0,4,length(X.Var1)+1);
Y = 0.5*(Y(2:end) + Y(1:end-1));
tA = Y*X.Var1/sum(X.Var1);
tA2 = (Y.^2)*X.Var1/sum(X.Var1);

X = readtable('~/Dropbox (MIT)/NonEqDetection/Data/stand_exp_' + string(i) + '.txt');
Y = linspace(0,12,length(X.Var1)+1);
Y = 0.5*(Y(2:end) + Y(1:end-1));
tB = exp(Y/60)*X.Var1/sum(X.Var1);

sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
disp(sig)




i= 3;
X = readtable('~/Dropbox (MIT)/NonEqDetection/Data/lying_exp_' + string(i) + '.txt');
histogram('BinEdges',linspace(0,4,length(X.Var1)+1),'BinCounts',X.Var1,...
    'Normalization','pdf','DisplayStyle', 'stairs')
ylim([0,1.6])
saveas(gcf,'lying_hist_' + string(i) + '.pdf')