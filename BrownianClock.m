%% For Plotting trajectories
dt = 0.001;
n = 500/dt;
[Path1,~,~] = simulate(bm(0.5,1),n,'DeltaTime', dt);
[Path2,~,~] = simulate(bm(1,1),n,'DeltaTime', dt);
[Path3,Times,~] = simulate(bm(2,1),n,'DeltaTime', dt);

writematrix([Times,Path1,Path2,Path3],'BrownianClockTraj.txt')

%% For calculating entropies
mu = linspace(0,3,15)';
deltat = [1;5;25;100];
n_bins = [1;3;6];
logfn = @(x,y) (x-y)*log(x/y);
sig2  = zeros(length(mu),length(deltat),length(n_bins));
sig1  = zeros(length(mu),length(deltat),length(n_bins));
sig_t = 2*mu.^2;
for i = 1:length(mu)
    dt = 0.001;
    n = 1000/dt;
    [Path,Times,~] = simulate(bm(mu(i),1),n,'DeltaTime', dt);
    for j = 1:length(deltat)
        for k = 1:length(n_bins)
            [trans,cond_trans] = gen_stats_bc(Path,Times,deltat(j),n_bins(k));
            for kk = 1:3
                s = circshift(1:3,kk);
                [~,v] = OptimGS(trans(s(2),s(3)), trans(s(1),s(2)), trans(s(2),s(1)), cond_trans(s(2)));
                sig2(i,j,k) = sig2(i,j) + 0.5*v;
                sig1(i,j,k) = sig1(i,j) + logfn(trans(s(1),s(2)),trans(s(2),s(1)));
            end
        end
    end
end


plot(mu,squeeze(sig2(:,1,1)),'o-',mu,squeeze(sig2(:,2,1)),'o-',mu,squeeze(sig2(:,3,1)),'o-', ...
    mu,squeeze(sig2(:,1,2)),'s-',mu,squeeze(sig2(:,2,2)),'s-',mu,squeeze(sig2(:,3,2)),'s-', ...
    mu,squeeze(sig2(:,1,3)),'d-',mu,squeeze(sig2(:,2,3)),'d-',mu,squeeze(sig2(:,3,3)),'d-', ...
   mu,squeeze(sig1(:,1,1)),'o',mu,squeeze(sig1(:,2,1)),'o',mu,squeeze(sig1(:,3,1)),'o', ...
   mu,squeeze(sig1(:,1,2)),'s',mu,squeeze(sig1(:,2,2)),'s',mu,squeeze(sig1(:,3,2)),'s', ...
   mu,squeeze(sig1(:,1,3)),'d',mu,squeeze(sig1(:,2,3)),'d',mu,squeeze(sig1(:,3,3)),'d', ...
    mu,sig_t,'k')

combined_res = zeros(length(sig2(:)),5);
n = 1;
for i = 1:length(mu)
    for j = 1:length(deltat)
        for k = 1:length(n_bins)
            combined_res(n,1) = mu(i);
            combined_res(n,2) = deltat(j);
            combined_res(n,3) = n_bins(k);
            combined_res(n,4) = sig1(i,j,k);
            combined_res(n,5) = sig2(i,j,k);
            n = n+1;
        end
    end
end
writematrix(combined_res,'BrownianClock.txt')