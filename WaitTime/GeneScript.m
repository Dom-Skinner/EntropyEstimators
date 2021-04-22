% Script to recreate fig 2b and entropy estimates for the gene data
%% For Bmal1
X = readtable('~/Dropbox (MIT)/NonEqDetection/Data/Bmal1a_off.txt').Var1;
Y = linspace(0,10,length(X)+1);
Y = 0.5*(Y(2:end) + Y(1:end-1));
histogram('BinEdges',linspace(0,10,length(X)+1),'BinCounts',X,...
    'Normalization','pdf','DisplayStyle', 'stairs')
tA = Y*X/sum(X);
tA2 = (Y.^2)*X/sum(X);
tB = 5/60; % from median on figure - estimate but not so important since it is so small.
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For Bmal1 gene, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)

%% For GTGlutaminas
X = readtable('~/Dropbox (MIT)/NonEqDetection/Data/GTGlutaminase.txt').Var1;
Y = linspace(0,10,length(X)+1);
Y = 0.5*(Y(2:end) + Y(1:end-1));
hold on
histogram('BinEdges',linspace(0,10,length(X)+1),'BinCounts',X,...
    'Normalization','pdf','DisplayStyle', 'stairs')
legend({'Bmal1a','GT:Glutaminase'})
xlabel('Time spent in off state (h)')
ylabel('Probability density')
tA = Y*X/sum(X);
tA2 = (Y.^2)*X/sum(X);
tB = 5/60; % from median on figure - estimate but not so important since it is so small.
sig = 2/(tA + tB) * gamma_precomp(tA2/tA^2,true);
fprintf('For GTGlutaminas gene, <tA^2>/<tA>^2 = %4.3f, and sig_T = %4.3f k_B /h\n',tA2/tA^2,sig)