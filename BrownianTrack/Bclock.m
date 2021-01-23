N = 40;

real_space = @(fk,x) 2*real(sum( fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x))) + real(fk(1));

real_space_f = @(fk,x) real(sum( fk .* exp( 1j * 2* pi* (-(N/2):(N/2))' * x)));
make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
reduce = @(fk) fk((numel(fk)-1)/2+1 :end);

ck = zeros(N/2+1,1);
ck(1) = 1;
ck(2) = 0.25*1j;

ak = zeros(N/2+1,1);
ak(1) = 1;
ak(2) = 0.25;

Fk = zeros(N/2+1,1);
Dk = zeros(N/2+1,1);
Dk(1) = 0.1;
Fk(1) = 0.0;

L = L_matrix(Fk,Dk,ck);

ds = zeros(N/2+1,1);
ds(1) = 1;
T = 2*reduce(L\(L \ make_full(ds)));
p = -0.5*reduce((L') \ make_full(ak));
T2 = convolve(T,ak);
sig = entropy_rate(ak,ck,p,Fk,Dk);
%x = linspace(0,1,400);
%hold on
%plot(x,real_space(ak,x),x,real_space_p(ak,x),x,real_space_p2(ak,x))

t2 = linspace(1.1,2,12);
sig = zeros(size(t2));
N = 40;
M = 7;
for i = 1:length(t2)
    [a,c,Dr,Fr,Tf,pf,fval,~] = min_cont(N,M, t2(i),1e-6,1e-6);
    sig(i) = fval;
end
%xs = linspace(0,1,400);
%real_space = @(fk,x) 2*real(sum( fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x))) + real(fk(1));
[a,c,Dr,Fr,Tf,pf,fval,~] = min_cont(N,M, 1.05,2e-6,2e-6);
plot(xs,real_space(a,xs),xs,real_space(c,xs),xs,real_space(Dk,xs),xs,real_space(pf,xs))

entropy_rate(a,c,pf,Fr,Dr)
plot(sig,t2)