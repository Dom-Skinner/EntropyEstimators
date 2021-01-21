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
