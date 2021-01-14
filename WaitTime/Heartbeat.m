s_to_n_dog = zeros(17,1);
tau_dog = zeros(size(s_to_n_dog));
for i = 1:length(s_to_n_dog)
    formatSpec = '~/Downloads/PZ_Databases/matlab_format/dog/Dog_%02d/peaks_Dog_%02d.mat';
    load(sprintf(formatSpec,i,i),'Data','Fs')
    wt = (Data(2:end) - Data(1:end-1))/Fs;
    wt = wt(wt < 5); % remove outliers
    %histogram(wt)
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
