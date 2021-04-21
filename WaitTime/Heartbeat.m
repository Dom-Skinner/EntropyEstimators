s_to_n_dog = zeros(17,1);
tau_dog = zeros(size(s_to_n_dog));
for i = 1:length(s_to_n_dog)
    formatSpec = '~/Downloads/PZ_Databases/matlab_format/dog/Dog_%02d/peaks_Dog_%02d.mat';
    
   
    load(sprintf(formatSpec,i,i),'Data','Fs')
    wt = (Data(2:end) - Data(1:end-1))/Fs;
    wt = wt(wt < 5); % remove outliers
    histogram(wt/mean(wt))
    s_to_n_dog(i) = mean(wt.^2)/mean(wt)^2;
    tau_dog(i) = mean(wt);
end

s_to_n_mouse = zeros(26-19+1,1);
tau_mouse = zeros(size(s_to_n_mouse));
for i = 19:26
    formatSpec = '~/Downloads/PZ_Databases/matlab_format/mouse/Mouse_%02d/peaks_Mouse_%02d.mat';
    load(sprintf(formatSpec,i,i),'Data','Fs')
    wt = (Data(2:end) - Data(1:end-1))/Fs;
    wt = wt(wt < 5); % remove outliers
    %histogram(wt)
    s_to_n_mouse(i-18) = mean(wt.^2)/mean(wt)^2;
    tau_mouse(i-18) = mean(wt);
end


s_to_n_rabbit = zeros(4*3,1);
tau_rabbit = zeros(size(s_to_n_rabbit));
for i = 1:4
    for j = 1:3
        formatSpec = '~/Downloads/PZ_Databases/matlab_format/rabbit/Rabbit_%02d_part_%01d/peaks_Rabbit_%02d_part_%01d.mat';
        load(sprintf(formatSpec,i,j,i,j),'Data','Fs')
        wt = (Data(2:end) - Data(1:end-1))/Fs;
        wt = wt(wt < 5); % remove outliers
        %histogram(wt)
        s_to_n_rabbit(i + 4*(j-1)) = mean(wt.^2)/mean(wt)^2;
        tau_rabbit(i+4*(j-1)) = mean(wt);
    end
end


x = [s_to_n_dog; s_to_n_mouse; s_to_n_rabbit];
g1 = repmat({'Dog'},length(s_to_n_dog),1);
g2 = repmat({'Mouse'},length(s_to_n_mouse),1);
g3 = repmat({'Rabbit'},length(s_to_n_rabbit),1);
g = [g1; g2; g3];
boxplot(x,g)


%% Make plots of Beat histogram
clear
i = 1;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/dog/Dog_%02d/peaks_Dog_%02d.mat'; 
load(sprintf(formatSpec,i,i),'Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
%wt = wt(wt < 5); % remove outliers
subplot(1,3,1)
title('Dog')
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
xlim([0,1.7])
ylim([0,4.5])


j = 2;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/rabbit/Rabbit_%02d_part_%01d/peaks_Rabbit_%02d_part_%01d.mat';
load(sprintf(formatSpec,i,j,i,j),'Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
%wt = wt(wt < 5); % remove outliers
subplot(1,3,2)
title('Rabbit')
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
xlim([0,1.7])
ylim([0,4.5])


i= 20;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/mouse/Mouse_%02d/peaks_Mouse_%02d.mat';
load(sprintf(formatSpec,i,i),'Data','Fs')
wt = (Data(2:end) - Data(1:end-1))/Fs;
%wt = wt(wt < 5); % remove outliers
subplot(1,3,3)
title('Mouse')
histogram(wt/mean(wt),'Normalization','pdf','BinWidth',0.08)
xlim([0,1.7])
ylim([0,4.5])



saveas(gcf,'EcgHist.pdf')

%% Make plots of ECG trace
clear
i = 1;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/dog/Dog_%02d/elctrography_Dog_%02d.mat';    
load(sprintf(formatSpec,i,i),'Data','Fs')
subplot(1,3,1)
title('Dog')
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))

j = 2;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/rabbit/Rabbit_%02d_part_%01d/elctrography_Rabbit_%02d_part_%01d.mat';
load(sprintf(formatSpec,i,j,i,j),'Data','Fs')
subplot(1,3,2)
title('Rabbit')
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))

i= 19;
formatSpec = '~/Downloads/PZ_Databases/matlab_format/mouse/Mouse_%02d/elctrography_Mouse_%02d.mat';
load(sprintf(formatSpec,i,i),'Data','Fs')
subplot(1,3,3)
title('Mouse')
plot((1:2*Fs)/(Fs),Data(500+(1:2*Fs)))

saveas(gcf,'EcgTrace.pdf')