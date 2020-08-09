function [sig,sig_1] = EntropyEstWrap(X,T)
    % I always end up using this same block of code, so it is convient to
    % wrap it in a tiny function.
    [trans,cond_trans] = gen_stats(X,T,3);

    [~,v1] = EntropyEst(trans(3,1), trans(1,2), trans(2,1), cond_trans(2));
    [~,v2] = EntropyEst(trans(1,2), trans(2,3), trans(3,2), cond_trans(3));
    [~,v3] = EntropyEst(trans(2,3), trans(3,1), trans(1,3), cond_trans(1));
    sig = 0.5*(v1+v2+v3);
    
    % Now do the naive estimator
    logfn = @(x,y) (x-y)*log(x/y);
    sig_1 = logfn(trans(1,2),trans(2,1));
    sig_1 = sig_1 + logfn(trans(2,3),trans(3,2));
    sig_1 = sig_1 + logfn(trans(3,1),trans(1,3));
end