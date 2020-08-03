function [trans,cond_trans] = gen_stats_bc(Path,Times,k,n)
% Essentially the same code as gen_stats, but for the brownian clock case,
% contains coarse graining in time as well as space. TODO - integrate with
% gen stats.
nbin = 3;
Path = Path(1:k:end); % coarse grain in time
Times = Times(end);
Path = mod(floor(nbin*n*Path/2/pi),nbin); % Coarse grain in space.
Path(Path==0) = nbin;

not_repeat = Path(2:end) ~= Path(1:end-1);
Path = Path([true;not_repeat]);

trans = zeros(nbin,nbin);
for i = 1:nbin
    for j = 1:nbin
        if i ~= j
            trans(i,j) = sum((Path(1:end-1) == i) & (Path(2:end) == j))/Times;
            if trans(i,j) == 0
                trans(i,j) = 1/Times; % If none observed assume 1 was observed (bad stats either way)
                disp("Warning: some transitions not observed")
            end
        end
    end
end
    
cond_trans = zeros(nbin,1);
cond_trans(1) = sum((Path(1:end-2) == nbin) & (Path(2:end-1) == 1) & (Path(3:end) == 2))/Times;
cond_trans(nbin) = sum((Path(1:end-2) == nbin-1) & (Path(2:end-1) == nbin) & (Path(3:end) == 1))/Times;
for i = 2:nbin-1
    cond_trans(i) = sum((Path(1:end-2) == i-1) & (Path(2:end-1) == i) & (Path(3:end) == i+1))/Times;
end
if any(cond_trans == 0)
    cond_trans(cond_trans == 0) = 1/Times; % If none observed assume 1 was observed (bad stats either way)
    disp("Warning: some transitions not observed")
end
end
