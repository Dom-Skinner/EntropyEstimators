clear

sig_est = zeros(20,20);


sig_est(:,1) = linspace(0.1,40,size(sig_est,1));
%{
for nvars = 13:13
    for i = 1:size(sig_est,1)
        ctol = 1e-6;
        otol = 1e-6;
        fvals = 20*ones(10,1);
        for j = 1:length(fvals)
            hess = true;
            [x,fval,exitflag] = min_via_fmin(nvars, sig_est(i,1),ctol,otol,false,hess);
            if entropy_r(reshape(x(1:nvars^2),nvars,nvars)) < sig_est(i) && exitflag >= 0
                fvals(j) = fval;
            else
                disp('g')
            end

        end
        %[x,fval,exitflag] = min_via_fmin(nvars, sig(i),0.1*ctol,0.1*otol,x);
        sig_est(i,nvars)  = min(fvals);
    end
end
%}

%%{
ctol = 1e-6;
otol = 1e-6;
for n_int = 2:13
    for i = 1:size(sig_est,1)


    ent_rate = @(wp) (n_int+1)*(wp-1) .*log(wp) - 0.1*sig_est(i,1);
    f = @(wp) ent_rate(wp);
    wpvals = linspace(1.01,10,5);
    for j = 1:length(wpvals)
        wp = wpvals(j);%fzero(f,[1 1e5]);

        B0 = -ones(n_int,n_int);
        for k= 1:size(B0,1)
            B0(k,k) = n_int+wp-1;
            if k < n_int
                B0(k,k+1) = -wp;
            end
        end
        B0 = B0/sum(B0(:));
        x0 = ones(n_int^2 + n_int,1);
        x0(1:n_int^2) = reshape(B0,n_int^2,1);
        x0(n_int^2 + 1 : end) = 1/n_int;
        hess = true;
        [x,fval,exitflag] = min_via_fmin(n_int, sig_est(i,1),ctol,otol,x0,hess);
        %hess = true;
        %[x,fval,exitflag] = min_via_fmin(n_int, sig_est(i,1),ctol,otol,x,hess);

        if (fval < sig_est(i,n_int)) || (0 == sig_est(i,n_int))
            sig_est(i,n_int) = fval;
        end
    end
    end
end
%}
for i = 2:13
    hold on
    plot(sig_est(:,1),sig_est(:,i),'+-');
end


%clearvars -except sig v nvars
%save('nvars_' + string(nvars))

%ctol = 1e-7;
%otol = 1e-7;
%[x,fval,exitflag] = min_via_fmin(3, 15,ctol,otol,x);