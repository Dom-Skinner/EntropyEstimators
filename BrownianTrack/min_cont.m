function [a,c,Dr,Fr,Tf,pf,fval,exitflag] = min_cont(N,M, t2,ctol,otol)
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
ar(2) = 0.2*rand();
ci(2) = 0.2*rand();
cr(1) = 1;
x0 = [ar;ai;cr;ci;D;F];

M0 = 24*M;
A = zeros(2*M0,length(x0));
for r = 1:M0
A(r,1:M) = -2*cos(2*pi*(1:M)/M0*r);
A(r,M+1:2*M) = 2*sin(2*pi*(1:M)/M0*r);

A(r+M0,2*M+1) = -1;
A(r+M0,(2*M+2):(3*M+1)) = -2*cos(2*pi*(1:M)/M0*r);
A(r+M0,(3*M+2):(4*M+1)) = 2*sin(2*pi*(1:M)/M0*r);

end
%A(2,2*M + 1:4*M+1) = 1;
b = [ones(M0,1);zeros(M0,1)];

options = optimoptions('fmincon','Display','iter','ConstraintTolerance',ctol,...
    'OptimalityTolerance',otol,'MaxFunctionEvaluations', 8e+04);

%[x,fval,exitflag,~]  = fmincon(ObjectiveFunction,x0,A,b,[],[],LB,UB,ConstraintFunction,options);

problem = createOptimProblem('fmincon','x0',x0,'objective',ObjectiveFunction,...
        'lb',LB,'ub',UB, 'Aineq',A,'bineq',b,'nonlcon',ConstraintFunction,'options',options);
gs = GlobalSearch;
[x,fval,exitflag,~] = run(gs,problem);

[a,c,Dr,Fr] = unpack(x,N,M);
Lf = L_matrix(Fr,Dr,c);
dsf = zeros(N/2+1,1);
dsf(1) = 1;
Tf = 2*reduce(Lf\(Lf \ make_full(dsf)));
pf = -0.5*reduce((Lf') \ make_full(a));

function [c, ceq] = simple_constraint(x,N,M,t2)
    %make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
    %reduce = @(fk) fk((numel(fk)-1)/2+1 :end);
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
    %make_full = @(fk) [conj(fk(end:-1:1)); fk(2:end)];
    %reduce = @(fk) fk((numel(fk)-1)/2+1 :end);
    [ak,ck,Dk,Fk] = unpack(x,N,M);
    L = L_matrix(Fk,Dk,ck);
    p = -0.5*reduce((L') \ make_full(ak));
    f = entropy_rate(ak,ck,p,Fk,Dk);
end
function f = make_full(fk)
    f = [conj(fk(end:-1:1)); fk(2:end)];
end

function f = reduce(fk)
    f = fk((numel(fk)-1)/2+1 :end);
end

end