function [x,fval,exitflag] = min_via_fmin_alt(n_int, t2,ctol,otol,x0)
% Calculate a lower bound for entropy production given that we have seen
% nVW transitions V -> W, nVU transitions V -> U, nUV transitions U -> V,
% and nVUPUWV transitions W -> V -> U.
% We do this by optimizing over a system with 4 hidden states within the V
% macrostate.

% Sanity checks
assert(isscalar(n_int+t2),'Error: Non-scalar input')
assert( floor(n_int) == n_int && (t2 > 0) && (n_int > 0))

ObjectiveFunction = @(x) obj(x,n_int);
nvars = n_int^2 + n_int;    % Number of variables

UB_h = 1e5;
UB = zeros(nvars,1);
UB_D = diag(UB_h*ones(n_int,1));
UB(1:n_int^2) = reshape(UB_D,n_int^2,1);
UB(n_int^2:end) = UB_h;

LB = zeros(nvars,1);   % Lower bound
LB(UB==0) = -UB_h;

% Set up the constraints for the problem. Constraints come from probability
% consv. and specifics about the entropy
A = zeros(2*n_int,nvars);
for i = 1:n_int
    rows = zeros(n_int,n_int);
    rows(i,:) = 1;
    A(i,1:n_int^2) = -reshape(rows,n_int^2,1); % constraint on row sum
    A(i+n_int,1:n_int^2) = -reshape(rows',n_int^2,1); % constraint on column sum
end
b = zeros(2*n_int,1);

Aeq = zeros(2,nvars);
Aeq(1,n_int^2 + 1: end) = 1; % constraint on sum(x)
Aeq(2,1:n_int^2) = 1; % constraint on sum(B)
beq = [1;1];

ConstraintFunction = @(x) simple_constraint(x,n_int,0.5*t2); % the one non-linear constraint

% Set up the problem with an initial guess and solve
if islogical(x0)
    B0 = -rand(n_int,n_int);
    for i= 1:size(B0,1)
        B0(i,i) = B0(i,i) - sum(B0(i,:)) - sum(B0(:,i))+ rand();
    end
    x0 = ones(nvars,1);
    x0(1:n_int^2) = reshape(B0,n_int^2,1);
    x0(n_int^2 + 1 : end) = 1/n_int;
end
options = optimoptions('fmincon','Display','iter','ConstraintTolerance',ctol,...
    'OptimalityTolerance',otol,'MaxFunctionEvaluations', 8e+03,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

[x,fval,exitflag,~]  = fmincon(ObjectiveFunction,x0,A,b,Aeq,beq,LB,UB,ConstraintFunction,options);
fval = 0.5*fval;

function [c, ceq,DC,DCeq] = simple_constraint(x,n_int,t2)
    % The variance of wait time.
    B = reshape(x(1:n_int^2),n_int,n_int);
    xs = x(n_int^2+1:end);
    q0 = (B\ xs);
    q1 = ((B')\ xs);
    c = xs'* q0- t2;
    DC = [reshape(-q1*q0',n_int^2,1); (q0+q1)];
    ceq = [];
    DCeq = [];
end

function [f, gradf] = obj(x,n_int)
    % The one non-linear equality constraint
    B = reshape(x(1:n_int^2),n_int,n_int);
    sig = 0;
    g = @(x,y) (x-y)*log( (abs(x)+eps)/(abs(y) + eps));
    h = @(x,y) (x-y)/(x + sign(x)*eps) +log( (abs(x)+eps)/(abs(y) + eps));
    for k = 1:(n_int-1)
        for l = (k+1):n_int
               % if (abs(B(k,l)) > eps) && (abs(B(l,k)) > eps)
                    sig = sig - g(B(k,l),B(l,k));%(B(k,l) - B(l,k))*log((abs(B(k,l))+eps)/(abs(B(l,k))+eps));
                %end
        end
    end
    
    c_ = zeros(n_int,1);
    d_ = zeros(n_int,1);
    for k = 1:n_int
        c_(k) = sum(B(k,:));
        d_(k) = sum(B(:,k));
        %if (c_k > eps) && (d_k > eps)
                    sig = sig + g(c_(k),d_(k));%(c_k - d_k)*log((c_k+eps)/(d_k+eps));
        %end
    end
    
    f = sig;
    M = zeros(n_int,n_int);
    for k = 1:n_int
        for l = 1:n_int
            M(k,l) = - h(B(k,l),B(l,k)) + h(c_(k),d_(k)) + h(d_(l),c_(l));
        end
    end
    gradf = [reshape(M,n_int^2,1);zeros(n_int,1)];
end
end