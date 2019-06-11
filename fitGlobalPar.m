%
% optimizes lysine pool parameters [a,b,r] for given dataset by minimizing
% squared sum of residuals between measured and predicted labeling for all given
% proteins simultaneously
%
% g = fitGlobalPar(d, p, dispFun)
%
% input:
%
% d                 : structure for experimental data:
%  .t   (T, 1)      : times after heavy-isotope feeding ("pulse time")
%  .y   (N, T)      : corresponding heavy-isotope amplitudes at times t for N proteins
%
% p                 : starting global parameters structure:
%  .gPar            : parameters for global pool equations, namely
%         (1)       : (a) protein degradation rate
%         (2)       : (b) lysine feeding / excretion rate
%         (3)       : (r) ratio of protein / lysine pool size
%   tau0            : fitting starting value for protein time constant
%
% dispFun           : function handle for custom display at each iteration
%                     (may be empty for default output)
%
% output:
% 
% g                 : struct containing fit result:
%     .gPar         (1,3)   [a,b,r], a=protein synthesis/degregation rate,
%                           b=lysine uptake/removal rate, r=bound/free pool
%                           size ratio
%     .tau0                 initial time constant for each fit
%     .RNrm                 total squared sum of residuals from all proteins
%     .R            (P,N)   residuals of P valid proteins in N measurements
%     .conf         (3,2)   [a,b,r] confidence intervals (lower and upper)
%     .J            (P,3)   Jacobian from fit for P proteins (for parameters [a,b,r])
%
% Mihai Alevra, 2019


function g = fitGlobalPar(d, p, dispFun)

    function s = locFit (gp)

        n = n+1;
        p.gPar = gp;
        o = fitSingle(d, p);        % do individual fits with next global pars

        % residuals: fill missing / bad values with average residuals, because lsqnonlin doesn't like NaNs
        s = abs(o.R); is = isfinite(s) & repmat(o.tau,[1 size(s,2)]) > 0.25; % unrealistic fits: very short tau
        s(~is) = nanmean(s(is));
     
        dispFun(o, s); 
    end

    function dispFit (o, s)
        if n==1
            fprintf('time       fit [     a      |     b      |     r      ]   tot. residuals\n');
        end
        
        fprintf('%s % 4d  [ % .3e | % .3e | % .3e ]   Err %.5e \n', timeMsg, n, o.a, o.b, o.r, sum(s(:)));
    end
    

    if nargin < 3
        dispFun = @dispFit;
    end

    n = 0;  % fit counter
    
%     opts = optimset('Display','off');
%     opts = optimoptions('patternsearch','Cache','on');
    opts = optimset('FinDiffRelStep',1e-3,'Display','off');
    
    [g.gPar, g.RNrm, g.R,~,~,~, g.J] = lsqnonlin (@locFit, p.gPar, [1e-3 1e-3 1e-3], [inf inf inf], opts);
    g.conf = nlparci(g.gPar, g.R, 'jacobian', g.J,'alpha',.05);
    g.tau0 = p.tau0;

end
