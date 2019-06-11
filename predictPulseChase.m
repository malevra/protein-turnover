%
% predicts labeling of proteins for pulse and chase experiments
%
% [yP, yC] = predictPulseChase(t, tStop, tau, gPar)
% 
% input:
%   t           (1,T)   measurement times [days]
%   tStop               end of pulse [days]
%   tau         (1,P)   time constants for P proteins
%   gPar        (1,3)   [a,b,r] lysine pool parameters
%
% output:
%   yP          (T,P)   labeling (heavy/total) for pure pulse experiment
%   yC          (T,P)   labeling (heavy/total) for pulse and chase experiment 
%                       (stopped heavy-isotope feeding after tStop days)
%
% Mihai Alevra, 2019


function [yP, yC] = predictPulseChase(t, tStop, tau, gPar)

    % amplitudes during pulse at times t and POI time constant p
    function yy = amp_pulse(p,t)
        yy = -exp(-t./p).*(-(A.*t1)./(t1-p)+(t2.*(A-1.0))./(t2-p)+1.0)+exp(-t./p).*(exp(t./p)-(A.*t1.*exp(-t./t1+t./p))./(t1-p)+(t2.*exp(-t./t2+t./p).*(A-1.0))./(t2-p));
    end

    % amplitudes during chase at times t and POI time constant p
    function yy = amp_chase(p,t)
        hp0 = amp_pulse(p,tStop);
        tN = t - tStop;
        yy = exp(-tN./p).*((A1.*t1.*exp(-tN./t1).*exp(tN./p))./(t1-p)+(A2.*t2.*exp(-tN./t2).*exp(tN./p))./(t2-p))-exp(-tN./p).*(-hp0+(A1.*t1)./(t1-p)+(A2.*t2)./(t2-p));
    end

    % solved labeled lysine for pulsed experiment
    function yy = Hs (t)
        yy = exp(-t/t2)*(A - 1) - A*exp(-t/t1) + 1;
    end

    % protein-bound labeled lysine for pulsed experiment
    function yy = Hp (t)
        yy = exp(-t/t1)/(C*t2) - exp(-t/t2)/(C*t1) + 1;
    end


    % set / calculate pool parameters  
    a = gPar(1);
    b = gPar(2);
    r = gPar(3);
    C = sqrt(-4*a*b + (a+b+a*r)^2);
    t1 = 2/(a+b+a*r+C);
    t2 = 2/(a+b+a*r-C);
    A  = -(a-b+a*r-C)/(2*C);
    A1 = ((C-a+b+a*r)*Hs(tStop) - 2*a*r*Hp(tStop)) / (2*C);
    A2 = ((C+a-b-a*r)*Hs(tStop) + 2*a*r*Hp(tStop)) / (2*C);
    
    T = numel(t);
    N = numel(tau);
    
    yP = NaN(T,N);
    yC = NaN(T,N);
    iP = t <= tStop;
    
    for n=1:N
        yP( : ,n)   = amp_pulse(tau(n),t     );
        yC( iP,n)   = amp_pulse(tau(n),t( iP));
        yC(~iP,n)   = amp_chase(tau(n),t(~iP));
    end

end