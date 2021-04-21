% This script contains all that is needed to recreate the heartbeat figures
% from fig 4, as well as the heartbeat entropy production estimates. 

%% Get human data
A = readtable('../../../Downloads/MIT-BIH Long-Term ECG Database/ecg.txt');

% We do a very rough estimate of the time between beats by thresholding and
% saying X > 2 is a beat. This is only used visualization purposes not for 
% the entropy calculation, which uses the more precise calculation from
% Umetani et al.
T = A.Var1(1:end);
X = A.Var3(1:end);
X = X > 2;

t_A = [];
t0 = 0;

for i = 2:length(T)
    if X(i) ~= X(i-1)
        t = T(i-1) - t0;
        t0 = T(i-1);
        if ~X(i-1)
            t_A(1+end) = t;
        end
    end
end

t_human_scaled = t_A/mean(t_A); % For histograms
T_human = T(520:776) - T(520); % For example ECG trace
X_human = A.Var3(520:776);

%% Make plots of Beat histogram
figure
subplot(1,3,1)
histogram(t_human_scaled,'Normalization','pdf','BinWidth',0.08)
title('Human')
xlim([0,1.7])
ylim([0,4.5])


i = 1;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/dog/Dog_%02d/peaks_Dog_%02d.mat'; 
load(sprintf(formatSpec,i,i),'Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
subplot(1,3,2)
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
title('Dog')
xlim([0,1.7])
ylim([0,4.5])



i= 20;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/mouse/Mouse_%02d/peaks_Mouse_%02d.mat';
load(sprintf(formatSpec,i,i),'Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
subplot(1,3,3)
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
title('Mouse')
xlim([0,1.7])
ylim([0,4.5])


saveas(gcf,'EcgHist.pdf')

%% Make plots of ECG trace

figure

subplot(1,3,1)
plot(T_human,X_human)
title('Human')


i = 1;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/dog/Dog_%02d/elctrography_Dog_%02d.mat';    
load(sprintf(formatSpec,i,i),'Data','Fs')
subplot(1,3,2)
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))
title('Dog')

i= 19;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/mouse/Mouse_%02d/elctrography_Mouse_%02d.mat';
load(sprintf(formatSpec,i,i),'Data','Fs')
subplot(1,3,3)
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))
title('Mouse')

saveas(gcf,'EcgTrace.pdf')

%% Calculate entropy estimates
% Here we use the data directly as measured in Umetani et al. and Behar et
% al. (see main text and SI).

clear

s_asy = @(t2) 1/t2 + 4*log(1/t2);

% Human data directly from table 2 of Umetani et al.
% For 10-19 year old humans
var_A = 81^2; % square of stdev
tA = 1000*60/80; % convert from beats per min to ms
var_A_scaled = var_A/tA^2;
sig = (1/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for 10-19 year old human\n',sig)

% For 40-49 year old humans
var_A = 60^2; % square of stdev
tA = 1000*60/78; % convert from beats per min to ms
var_A_scaled = var_A/tA^2;
sig = (1/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for 40-49 year old human\n',sig)

% For dogs
var_A = 69.20^2;
tA = 482.63;
var_A_scaled = var_A/tA^2;
sig = (1/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for dogs\n',sig)

% For mice
var_A = 10.39^2;
tA = 108.46;
var_A_scaled = var_A/tA^2;
sig = (1/tA) * s_asy(var_A_scaled)*1000;
fprintf('Entropy production ~%5.1f k_B/s for mice\n',sig)
