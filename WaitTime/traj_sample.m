function [q_A,N_emp,TA_emp,TB_emp,TA2_emp] = traj_sample(T,Wm,ntrials)
% Simulate the waiting times of the active sensor for a window T. return
% the empirically observed qA.

n = 4;
vec1 = ones(n,1);
pi_A = ones(n,1)/(n+1);

WA_fun = @(wp,wm) [-(wp+wm), wp, 0, 0; ... 
                  wm, -(wp+wm), wp, 0; ... 
                  0, wm, -(wp+wm), wp; ...
                0, 0, wm, -(wm+wp)];
            
t = linspace(0,40,2000);
f = zeros(length(t),1);
f_cs = zeros(length(t),1);

% Define the partial transition rate matrix WA, and calculate the
% cumalitive probability density from analytic formula
WA = WA_fun(1,Wm);
K = - (pi_A')*WA * vec1;

for i = 1:length(t)
    f(i) = (pi_A')*(WA*WA* expm(WA*t(i))) * vec1/K;
    f_cs(i) = (pi_A')*(WA*expm(WA*t(i)) - WA) * vec1/K;
end

% sample from wait time distribution by sampling uniform dist and transform
f_cs_inv = @(x_in) interp1(f_cs,t,x_in,'pchip'); 

q_A = zeros(ntrials,1);
N_emp = zeros(ntrials,1);
TA_emp = zeros(ntrials,1);
TA2_emp = zeros(ntrials,1);
TB_emp = zeros(ntrials,1);
for n = 1:ntrials
    X = rand(round(T) + 1e4,1);
    TA_rand = f_cs_inv(X);
    TB_rand =  exprnd(1/(1 + Wm),size(X));

    TA_glob = zeros(2*length(TA_rand),1);
    TA_glob(1:2:end) = TA_rand;
    TB_glob = zeros(2*length(TB_rand),1);
    TB_glob(2:2:end) = TB_rand;
    T_tot = cumsum(TA_glob+TB_glob);



    T_start = 1e4;
    T_end = T_start + T;
    idx_start = find(T_tot > T_start,1,'first');
    idx_end = find(T_tot > T_end,1,'first');

    if isempty(idx_end)
        error("Trajectory window too long")
    end

    q_A(n) = (sum(TA_glob(idx_start+1:idx_end)) + ...
        (TA_glob(idx_start) > 0) * (T_tot(idx_start) - T_start) -...
        (TA_glob(idx_end) >0) * (T_tot(idx_end) - T_end)  )/T;
    N_emp(n) = idx_end - idx_start;
    TA_emp(n) = 2*mean(TA_glob(idx_start+1: idx_end-1)); % factor of 2 because every other entry is 0.
    TB_emp(n) = 2*mean(TB_glob(idx_start+1: idx_end-1)); 
    TA2_emp(n) = 2*mean(TA_glob(idx_start+1: idx_end-1).^2); 
end
end