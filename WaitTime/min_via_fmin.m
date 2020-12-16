function [x,fval,exitflag] = min_via_fmin(n_int, sig_f,ctol,otol,x0,hess)
% Calculate a lower bound for entropy production given that we have seen
% nVW transitions V -> W, nVU transitions V -> U, nUV transitions U -> V,
% and nVUPUWV transitions W -> V -> U.
% We do this by optimizing over a system with 4 hidden states within the V
% macrostate.

% Sanity checks
assert(isscalar(n_int+sig_f),'Error: Non-scalar input')
assert( floor(n_int) == n_int && (sig_f > 0) && (n_int > 0))

ObjectiveFunction = @(x) time_var(x,n_int);
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

ConstraintFunction = @(x) simple_constraint(x,n_int,2*sig_f); % the one non-linear constraint
hessfcn = @(x,lambda) hessianfcn(x,lambda,n_int);
% Set up the problem with an initial guess and solve
if islogical(x0)
    B0 = -rand(n_int,n_int);
    for i= 1:size(B0,1)
        B0(i,i) = B0(i,i) - sum(B0(i,:)) - sum(B0(:,i))+ rand();
    end
    %B0 = B0/sum(B0(:));
    x0 = ones(nvars,1);
    x0(1:n_int^2) = reshape(B0,n_int^2,1);
    x0(n_int^2 + 1 : end) = 1/n_int;
end
options = optimoptions('fmincon','Display','iter','ConstraintTolerance',ctol,...
    'OptimalityTolerance',otol,'MaxFunctionEvaluations', 8e+03,...
    'MaxIterations',8e3,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
if hess
    options.HessianFcn = hessfcn;
    options.DiffMaxChange = 1;
end
[x,fval,exitflag,~]  = fmincon(ObjectiveFunction,x0,A,b,Aeq,beq,LB,UB,ConstraintFunction,options);
fval = 2*fval;

function [f,gradf] = time_var(x,n_int)
    % The variance of wait time.
    B = reshape(x(1:n_int^2),n_int,n_int);
    xs = x(n_int^2+1:end);
    q0 = (B\ xs);
    q1 = ((B')\ xs);
    f = xs'* q0;
    gradf = [reshape(-q1*q0',n_int^2,1); (q0+q1)];
end

function [c, ceq,DC,DCeq] = simple_constraint(x,n_int,sig_f)
    % The one non-linear equality constraint
    B = reshape(x(1:n_int^2),n_int,n_int);
    sig = 0;
    g = @(x,y) (x-y)*log( (x+eps)/(y + eps));
    h = @(x,y) (x-y)/(x + eps) +log( (x+eps)/(y + eps));
    for k = 1:(n_int-1)
        for l = (k+1):n_int
               % if (abs(B(k,l)) > eps) && (abs(B(l,k)) > eps)
               if l ~= k
                    sig = sig + g(-B(k,l),-B(l,k));%(B(k,l) - B(l,k))*log((abs(B(k,l))+eps)/(abs(B(l,k))+eps));
               end
                %end
        end
    end
    
    c_ = zeros(n_int,1);
    d_ = zeros(n_int,1);
    for k = 1:n_int
        c_(k) = -sum(B(k,:));
        d_(k) = -sum(B(:,k));
        %if (c_k > eps) && (d_k > eps)
                    sig = sig + g(-c_(k),-d_(k));%(c_k - d_k)*log((c_k+eps)/(d_k+eps));
        %end
    end
    ceq = [];
    c = sig-sig_f;
    M = zeros(n_int,n_int);
    for k = 1:n_int
        for l = 1:n_int
            if l ~=k
                M(k,l) = - h(-B(k,l),-B(l,k)) + h(-c_(k),-d_(k)) + h(-d_(l),-c_(l));
            else
                M(k,l) =  h(-c_(k),-d_(k)) + h(-d_(l),-c_(l));
            end
        end
    end
    DC = [reshape(M,n_int^2,1);zeros(n_int,1)];
    DCeq = [];
end

function [p,q] = index_ref(k,n_int)
    p = mod(k,n_int) + n_int*(mod(k,n_int) == 0);
    q = floor((k-1)/n_int)+1;
end
        
    function Hout = hessianfcn(x,lambda,n_int)
        lam = lambda.ineqnonlin;
        B = reshape(x(1:n_int^2),n_int,n_int);
        Binv = inv(B);
        xs = x(n_int^2+1:end);
        y = (B\ xs);
        z = ((B')\ xs);
        Hout = zeros(length(x),length(x));
        for k = 1:n_int^2
            [p,q] = index_ref(k,n_int);
            for l = 1:n_int^2
                [r,s] = index_ref(l,n_int);
                Hout(k,l) = z(r)*Binv(s,p)*y(q) + z(p)*Binv(q,r)*y(s);
                %Hout(l,k) = Hout(k,l);
            end
            
            for l = 1:n_int
                Hout(k,l+n_int^2) = -Binv(l,p)*y(q) - z(p)*Binv(q,l);
                Hout(l+n_int^2,k) = Hout(k,l+n_int^2);
            end
            
        end
        for p = 1:n_int
            for q = 1:n_int
                Hout(p + n_int^2,q + n_int^2) = Binv(p,q) + Binv(q,p);
                %Hout(q + n_int^2,p + n_int^2) = Binv(p,q) + Binv(q,p);
            end
        end
        ha = @(a,b) (a+b+2*eps)./(a+eps).^2;
        hb = @(a,b) -(a+b+2*eps)./((a+eps).*(b+eps));
        for k = 1:n_int^2
            [p,q] = index_ref(k,n_int);
            for l = 1:n_int^2
                [i1,j1] = index_ref(l,n_int);
                
                if (p == i1) && (q ==j1)
                    Hout(k,l) = Hout(k,l) + lam*ha(-B(i1,j1),-B(j1,i1));
                    %Hout(l,k) = Hout(k,l);
                end
                
                if (q == i1) && (p == j1)
                    Hout(k,l) = Hout(k,l) + lam*hb(-B(i1,j1),-B(j1,i1));
                    %Hout(l,k) = Hout(k,l);
                end
                
                if p == i1
                    Hout(k,l) = Hout(k,l) + lam*ha(sum(B(i1,:)),sum(B(:,i1)));
                    %Hout(l,k) = Hout(k,l);
                end
                
                if q == i1
                    Hout(k,l) = Hout(k,l) + lam*hb(sum(B(i1,:)),sum(B(:,i1)));
                    %Hout(l,k) = Hout(k,l);
                end
                
                
                if q == j1
                    Hout(k,l) = Hout(k,l) + lam*ha(sum(B(:,j1)),sum(B(j1,:)));
                   % Hout(l,k) = Hout(k,l);
                end
                
                if p == j1
                    Hout(k,l) = Hout(k,l) + lam*hb(sum(B(:,j1)),sum(B(j1,:)));
                    %Hout(l,k) = Hout(k,l);
                end
            end
        end
    end


end