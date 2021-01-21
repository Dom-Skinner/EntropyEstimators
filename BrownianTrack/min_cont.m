function [x,fval,exitflag] = min_cont(N,M, t2,ctol,otol)
% Calculate a lower bound for entropy production given that we have seen
% nVW transitions V -> W, nVU transitions V -> U, nUV transitions U -> V,
% and nVUPUWV transitions W -> V -> U.
% We do this by optimizing over a system with 4 hidden states within the V
% macrostate.

% Sanity checks
assert(isscalar(N+M+t2),'Error: Non-scalar input')
assert( floor(N) == N && (t2 > 0) && (N > 0) && (M>0))

ObjectiveFunction = @(x) obj(x,N,M);
nvars = 4*M + 3;
UB = (1e4)*ones(nvars,1);
LB = -(1e4)*ones(nvars,1);
LB(end-1:end) = 0;

% Set up the constraints for the problem. Constraints come from probability
% consv. and specifics about the entropy

ConstraintFunction = @(x) simple_constraint(x,N,M,t2); % the one non-linear constraint

% Set up the problem with an initial guess and solve
ar = zeros(M,1);
ai = zeros(M,1);
cr = zeros(M+1,1);
ci = zeros(M,1);
D = 0.01;
F = 1;
ar(2) = 0.1;
ci(2) = 0.1;
cr(1) = 1;
x0 = [ar;ai;cr;ci;D;F];

A = zeros(M,length(x0));
M0 = 4*M;
for r = 1:M0
A(r,1:M) = -2*cos(2*pi*(1:M)/M0*r);
A(r,M+1:2*M) = 2*sin(2*pi*(1:M)/M0*r);
end
%A(2,2*M + 1:4*M+1) = 1;
b = 1*ones(M0,1);

options = optimoptions('fmincon','Display','iter','ConstraintTolerance',ctol,...
    'OptimalityTolerance',otol,'MaxFunctionEvaluations', 8e+03);

[x,fval,exitflag,~]  = fmincon(ObjectiveFunction,x0,A,b,[],[],LB,UB,ConstraintFunction,options);

function [a,c,Dk,Fk] = unpack(x,N,M)
    a_r = x(1:M);
    a_i = x(M+1:2*M);
    c_r = x(2*M+1:3*M+1);
    c_i = x(3*M+2:4*M+1);
    D_ = x(4*M+2);
    F_ = x(4*M+3);

    c = zeros(N/2+1,1);
    c(1) = c_r(1);
    c(2:M+1) = c_r(2:end) + 1j*c_i;

    a = zeros(N/2+1,1);
    a(1) = 1;
    a(2:M+1) = a_r + 1j*a_i;

    Dk = zeros(N/2+1,1);
    Dk(1) = D_;

    Fk = zeros(N/2+1,1);
    Fk(1) = F_;
end

function [c, ceq] = simple_constraint(x,N,M,t2)
    make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
    reduce = @(fk) fk((numel(fk)-1)/2+1 :end);
    [ak,ck,Dk,Fk] = unpack(x,N,M);
    
    L = L_matrix(Fk,Dk,ck);
    ds = zeros(N/2+1,1);
    ds(1) = 1;
    T = 2*reduce(L\(L \ make_full(ds)));
    p = -0.5*reduce((L') \ make_full(ak));
    T2 = convolve(T,ak);

    ceq = [T2 - t2; real(p(1))-0.5];
    c = [];
end

function f = obj(x,N,M)
    make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
    reduce = @(fk) fk((numel(fk)-1)/2+1 :end);
    [ak,ck,Dk,Fk] = unpack(x,N,M);
    L = L_matrix(Fk,Dk,ck);
    p = -0.5*reduce((L') \ make_full(ak));
    f = entropy_rate(ak,ck,p,Fk,Dk);
    
    
    
    xs = linspace(0,1,400);
    real_space = @(fk,x) 2*real(sum( fk(2:end) .* exp( 1j * 2* pi* (1:length(fk)-1)' * x))) + real(fk(1));
    plot(xs,real_space(ak,xs),xs,real_space(ck,xs),xs,real_space(p,xs))
end

end