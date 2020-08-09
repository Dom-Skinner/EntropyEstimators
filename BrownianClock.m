% Calculate the entropy of a Brownian clock: a biased Brownian motion that
% takes place on a periodic domain.

%% Visualize typical trajectories
dt = 0.001;
n = 25/dt;
[Path1,~,~] = simulate(bm(0.5,1),n,'DeltaTime', dt);
[Path2,~,~] = simulate(bm(1,1),n,'DeltaTime', dt);
[Path3,Times,~] = simulate(bm(2,1),n,'DeltaTime', dt);

plot(Times,Path1 - floor(Path1/2/pi)*2*pi,'.', ... % periodic with period 2 pi
    Times,Path2- floor(Path2/2/pi)*2*pi,'.', ...
    Times, Path3- floor(Path3/2/pi)*2*pi,'.')
xlabel('time')
%writematrix([Times,Path1,Path2,Path3],'BrownianClockTraj.txt')

%% Calculate entropies for a range of Brownian motion parameters. Compare
% with the naive estimator as well as the analytic solution

mu = linspace(0.,3,15)'; % Simulate for these mu's with sigma =  1

% Dimensional analysis tells us how to coarse grain given the parameters of
% Brownian motion. However, we have to find the constant numerically.
c1 = [0.005, 0.01, 0.015];
c2 = [1, 2, 4];

logfn = @(x,y) (x-y)*log(x/y); % useful for the naive estimator
sig2  = zeros(length(mu),length(c1),length(c2)); % for our estimator
sig1  = zeros(length(mu),length(c1),length(c2)); % for the naive estimator
sig_t = 2*mu.^2; % The true rate of entropy production
for i = 1:length(mu)
    dt = 0.001;
    n = 10000/dt; 
    [Path,Times,~] = simulate(bm(mu(i),1),n,'DeltaTime', dt);
    for j = 1:length(c1)
        for k = 1:length(c2)
            deltat = round(c1(j)*1/mu(i)^2  / dt);   % Determine the spatial
            n_bins = round(c2(k)*2*pi / 3 / mu(i) ); % and temporal coarse
            n_bins = max([n_bins;1]);                % graining.
            deltat = max([deltat;1]);
            deltat = min([deltat;1e4]);
            
            X =  floor(3*n_bins*Path(1:deltat:end)/2/pi);
            [trans,cond_trans] = gen_stats(X,Times(end),3);
            for kk = 1:3
                s = circshift(1:3,kk);
                [~,v] = EntropyEst(trans(s(2),s(3)), trans(s(1),s(2)), trans(s(2),s(1)), cond_trans(s(2)));
                sig2(i,j,k) = sig2(i,j,k) + 0.5*v;
                sig1(i,j,k) = sig1(i,j,k) + logfn(trans(s(1),s(2)),trans(s(2),s(1)));
            end
        end
    end
end


%% Plot the results
% too many lines to label, but we see that using the improved estimator, we
% can get a reasonable estimate of the entropy production rate, beyond that
% of the naive estimator
plot(mu,squeeze(sig2(:,1,1)),'o-',mu,squeeze(sig2(:,2,1)),'o-',mu,squeeze(sig2(:,3,1)),'o-', ...
    mu,squeeze(sig2(:,1,2)),'s-',mu,squeeze(sig2(:,2,2)),'s-',mu,squeeze(sig2(:,3,2)),'s-', ...
    mu,squeeze(sig2(:,1,3)),'d-',mu,squeeze(sig2(:,2,3)),'d-',mu,squeeze(sig2(:,3,3)),'d-', ...
   mu,squeeze(sig1(:,1,1)),'o',mu,squeeze(sig1(:,2,1)),'o',mu,squeeze(sig1(:,3,1)),'o', ...
   mu,squeeze(sig1(:,1,2)),'s',mu,squeeze(sig1(:,2,2)),'s',mu,squeeze(sig1(:,3,2)),'s', ...
   mu,squeeze(sig1(:,1,3)),'d',mu,squeeze(sig1(:,2,3)),'d',mu,squeeze(sig1(:,3,3)),'d', ...
    mu,sig_t,'k')
xlabel('mu')
ylabel('entropy production rate')

%% Save results to text file
combined_res = zeros(length(sig2(:)),5);
n = 1;
for i = 1:length(mu)
    for j = 1:length(c1)
        for k = 1:length(c2)
            combined_res(n,1) = mu(i);
            combined_res(n,2) = c1(j);
            combined_res(n,3) = c2(k);
            combined_res(n,4) = sig1(i,j,k);
            combined_res(n,5) = sig2(i,j,k);
            n = n+1;
        end
    end
end
writematrix(combined_res,'BrownianClock.txt')