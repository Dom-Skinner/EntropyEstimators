function vq = gamma_precomp(x_in,invert)
load('nvars_comb.mat','sig_est')
sig = sig_est(:,1);
f = min(sig_est(:,2:end),[],2);

if invert
    x_in(x_in > 2) = 2;
    vq = interp1([2;f],[0;sig],x_in,'pchip');
else
    vq = interp1([0;sig],[2;f],x_in,'pchip');
end
end