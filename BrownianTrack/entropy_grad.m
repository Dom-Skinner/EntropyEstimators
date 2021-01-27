function [P,Q,R] = entropy_grad(a,b,p)
N = 2*length(a)-2;
offset = N/2+1;

P = zeros(N+1,1);
Q = zeros(N+1,1);
R = zeros(N+1,1);


x = linspace(0,1,400);
x = x(1:end-1);

rs = @(fk) 2*real(sum( fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x))) + real(fk(1));

a = rs(a);
b = rs(b);
p = rs(p);

for l = -N/2:N/2
    P(l+offset) = sum( exp(2*pi*1j*l*x).* b .* log(b./a))/length(x);
    Q(l+offset) = sum( exp(2*pi*1j*l*x).* (-2*p.*b./a + (log(a./b) + 1) ) )/length(x);
    R(l+offset) = sum( exp(2*pi*1j*l*x).* (2*p.*(log(b./a) + 1) - a./b))/length(x);
end

end