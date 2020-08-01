function [xgs,v] = EntropyEst(nVW, nVU, nUV, nVUPUWV)
% Calculate a lower bound for entropy production given that we have seen
% nVW transitions V -> W, nVU transitions V -> U, nUV transitions U -> V,
% and nVUPUWV transitions W -> V -> U.
% We do this by optimizing over a system with 4 hidden states within the V
% macrostate.

% Sanity checks
assert(isscalar(nVW + nVU + nUV + nVUPUWV),'Error: Non-scalar input')
assert( (nVW >= 0) && (nVU >= 0) && (nUV >= 0) && (nVUPUWV >= 0))

min_size = 1e-10; % can't have any input too small for numerical purposes
nVW     = max(nVW,min_size); 
nVU     = max(nVU,min_size);
nUV     = max(nUV,min_size);
nVUPUWV = max(nVUPUWV,min_size);
if min([nVW,nVU,nUV,nVUPUWV]) == min_size % warn user if inputs changed
    warning('Some inputs too small, replaced by default tolerance')
end

n_int = 4; % number of internal states
ObjectiveFunction = @(x) sys_entropy(x,n_int);
nvars = 3*n_int;    % Number of variables
LB = zeros(nvars,1);   % Lower bound
UB = [nUV*ones(n_int,1); nVU*ones(n_int,1); nVW*ones(n_int,1)];  % Upper bound

% Set up the constraints for the problem. Constraints come from probability
% consv. and specifics about the entropy
A = zeros(2 + n_int,3*n_int);
A(1,1:n_int) = 1;
A(2,1+2*n_int:3*n_int) = 1;
for i = 1:n_int
    A(2+i,i) = -1;
    A(2+i,i+n_int) = 1;
    A(2+i,i+2*n_int) = -1;
end
b = [nUV;nVW;zeros(n_int,1)];
Aeq = [ones(1,n_int),-1*ones(1,n_int),zeros(1,n_int)];
beq = nUV-nVU;
ConstraintFunction = @(x) simple_constraint(x,nVUPUWV,n_int); % the one non-linear constraint

% Set up the problem with an initial guess and solve
x0 = [nUV*ones(n_int,1)/n_int;nVU*ones(n_int,1)/n_int;nVW*ones(n_int,1)/n_int];
ConTol = 1e-4;
opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',800,...
                    'MaxFunctionEvaluations',180*nvars,...
                    'FunctionTolerance',1e-5,'ConstraintTolerance',ConTol);
problem = createOptimProblem('fmincon','objective',ObjectiveFunction,...
            'x0',x0,'lb',LB,'ub',UB,'nonlcon',ConstraintFunction,...
            'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'options',opts);
gs = GlobalSearch('BasinRadiusFactor',0.1);
[xgs,v,EXITFLAG] = run(gs,problem);


% Sometimes optimization fails to converge, due to numerical issues. In
% this case, we rerun with a larger constraint tolerance. Running with a
% large constraint tolerance will only find solutions with a lower
% objective value than the constrained minimum, which is what is needed for
% a lower bound. We warn the user in any case.
while EXITFLAG < 1
    ConTol = 4*ConTol;
    opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',800,...
                    'MaxFunctionEvaluations',180*nvars,...
                    'FunctionTolerance',1e-5,'ConstraintTolerance',ConTol);
    problem = createOptimProblem('fmincon','objective',ObjectiveFunction,...
                    'x0',x0,'lb',LB,'ub',UB,'nonlcon',ConstraintFunction,...
                    'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'options',opts);
    [xgs,v,EXITFLAG] = run(gs,problem);
    if (ConTol > 1e-2) && (EXITFLAG < 0)
        warning('Large constraint tolerance had to be used')
    end
end


function sig = sys_entropy(x,n_int)
    % the objective function, in this case the entropy for the system
    niV = x(1:n_int);
    nVi = x(1+n_int:2*n_int);
    niW = x(1+2*n_int:3*n_int);
    nWi = niW + niV - nVi;
    sig = real(sum( (niV - nVi).* log(niV./(nVi))) + sum( (niW - nWi).* log(niW./(nWi))));
end

function [c, ceq] = simple_constraint(x,nVUPUWV,n_int)
    % The one non-linear equality constraint
    niV = x(1:n_int);
    nVi = x(1+n_int:2*n_int);
    niW = x(1+2*n_int:3*n_int);
    c = [];
    ceq = -sum(nVi.*niW ./(niV+niW + 1e-10)) + nVUPUWV;
end
end