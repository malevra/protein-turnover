% function o = fitSingle (d,p)
%
%
% input:
%
% d                 : structure for experimental data:
%  .t   (T, 1)      : times after heavy-isotope feeding ("pulse time")
%  .y   (N, T)      : corresponding heavy-isotope amplitudes at times t for N proteins
%
% p                 : global parameters structure:
%  .gPar            : parameters for global pool equations, namely
%         (1)       : (a) protein degradation rate
%         (2)       : (b) lysine feeding / excretion rate
%         (3)       : (r) ratio of protein / lysine pool size
%   tau0            : fitting starting value for protein time constant
%
% ouput:
%
% o                 : result structure:
%  .tau    (N, 1)   : individual protein time constants
%  .fY     (N, T)   : predicted heavy-protein amplitudes
%  .ci_tau (N, 2)   : confidence intervals (95% confidence) for time constants
%  .ci_fY  (N, 2)   : confidence intervals (95% confidence) for predicted amplitudes
%  .RNrm   (N, 1)   : normalized residual sum
%  .R      (N, 1)   : raw residuals
%  .a,b,r,C,t1,t2,A : global derived parameters
%
% Mihai Alevra, 2019


function o = fitSingle (d,p)

    % predicted protein-of-interest evolution with time constant tp at times t
    % (and global parameters t1,t2,A)
    function yy = amp_pulse(tp,t)
        yy = -exp(-t./tp).*(-(A.*t1)./(t1-tp)+(t2.*(A-1.0))./(t2-tp)+1.0)+exp(-t./tp).*(exp(t./tp)-(A.*t1.*exp(-t./t1+t./tp))./(t1-tp)+(t2.*exp(-t./t2+t./tp).*(A-1.0))./(t2-tp));
    end

    % residuals between prediction and data
    function ys = dif_pulse(p)
        ys = amp_pulse(p,ti) - yi;
    end

    % set / calculate parameters
  
    a = p.gPar(1);
    b = p.gPar(2);
    r = p.gPar(3);
    C = sqrt(-4*a*b + (a+b+a*r)^2);
    t1 = 2/(a+b+a*r+C);
    t2 = 2/(a+b+a*r-C);
    A  = -(a-b+a*r-C)/(2*C);
    N = size(d.y,1); T = numel(d.t);

    % initialize result arrays
    
    o.tau     = NaN(N,1);           % time constant for proteins-of-interest
    o.fY      = NaN(N,T);           % fitted amplitudes
    o.RNrm    = NaN(N,1);           % normalized summed residuals
    o.R       = NaN(N,T);
    
    o.ci_tau  = NaN(N,2);           % confidence interval (95%) for time constants
    o.ci_fY   = NaN(N,T,2);         % confidence interval (95%) for predicted amplitudes
    
    o.a = a;
    o.b = b;
    o.r = r;
    o.C = C;
    o.t1 = t1;
    o.t2 = t2;
    o.A = A;
    
    opts = optimset('Display','off');
    
    % loop over all proteins
    for n=1:N
        iOk = isfinite(d.y(n,:));       % select only finite measurements
        yi  = d.y(n,iOk);
        ti  = d.t(iOk);
        
        if numel(unique(ti))>0          % need at least 3 different pulse time points for fit

            [o.tau(n),o.RNrm(n), o.R(n,iOk),~,~,~, J] = lsqnonlin (@dif_pulse, p.tau0, 0.25, 1000, opts);
            
            % predicted amplitudes and confidence intervals
            
            % confidence bounds for single fit parameter tau
            o.ci_tau(n,:) = nlparci(o.tau(n), o.R(n,iOk), 'jacobian', J,'alpha',.05);

            % predicted amplitude 
            o.fY(n,iOk) = amp_pulse(o.tau(n),ti);
            
            % predicted amplitude at lower/upper confidence bounds
            o.ci_fY(n,iOk,1) = amp_pulse(o.ci_tau(n,1),ti);
            o.ci_fY(n,iOk,2) = amp_pulse(o.ci_tau(n,2),ti);
            
        end
    end
end