%
% predict labeling amplitudes for proteins with incomplete cleavage
%
% light:    no  18-O labeling
% medium:   one 18-O labeling
% heavy:    two 18-O labeling
%
% in:
%
%   t        (1,T)      measurement times [days]
%   tau      (1,P)      time constant(s) for P proteins
%   gPar     (1,3)      [a,b,r] lysine pool parameters
%   
% out:
%   yL       (T,P)      predicted labeling (light/total ) for P proteins at T times
%   yM       (T,P)      predicted labeling (medium/total) for P proteins at T times
%   yH       (T,P)      predicted labeling (heavy/total ) for P proteins at T times
%
% Mihai Alevra, 2019

function [yL, yM, yH] = predictDoubleLabel(t, tau, gPar)


    % prediction function for single-labeling
    function y = fM (tpoi)
        y = -exp(-t./tpoi).*((t2.*2.0)./(t2-tpoi)-(t2.*2.0)./(t2-tpoi.*2.0)+(A.*t1.*2.0)./(t1-tpoi)-(A.*t2.*2.0)./(t2-tpoi)+(A.*t2.*4.0)./(t2-tpoi.*2.0)-(A.^2.*t1.*2.0)./(t1-tpoi.*2.0)-(A.^2.*t2.*2.0)./(t2-tpoi.*2.0)-(t2.*exp(-t./t2+t./tpoi).*2.0)./(t2-tpoi)+(t2.*exp((t.*-2.0)./t2+t./tpoi).*2.0)./(t2-tpoi.*2.0)+(A.*t1.*t2.*4.0)./(-t1.*t2+t1.*tpoi+t2.*tpoi)-(A.*t1.*exp(-t./t1+t./tpoi).*2.0)./(t1-tpoi)+(A.*t2.*exp(-t./t2+t./tpoi).*2.0)./(t2-tpoi)-(A.*t2.*exp((t.*-2.0)./t2+t./tpoi).*4.0)./(t2-tpoi.*2.0)-(A.^2.*t1.*t2.*4.0)./(-t1.*t2+t1.*tpoi+t2.*tpoi)+(A.^2.*t1.*exp((t.*-2.0)./t1+t./tpoi).*2.0)./(t1-tpoi.*2.0)+(A.^2.*t2.*exp((t.*-2.0)./t2+t./tpoi).*2.0)./(t2-tpoi.*2.0)+(A.^2.*t1.*t2.*exp(-t./t1-t./t2+t./tpoi).*4.0)./(-t1.*t2+t1.*tpoi+t2.*tpoi)-(A.*t1.*t2.*exp(-t./t1-t./t2+t./tpoi).*4.0)./(-t1.*t2+t1.*tpoi+t2.*tpoi));
    end

    % prediction function for double-labeling
    function y = fH (tpoi) 
        y = exp(-t./tpoi).*(exp(t./tpoi)+(t2.*exp((t.*(t2-tpoi))./(t2.*tpoi)).*(A-1.0).*2.0)./(t2-tpoi)+(A.^2.*t1.*exp((t.*(t1-tpoi.*2.0))./(t1.*tpoi)))./(t1-tpoi.*2.0)+(t2.*exp((t.*(t2-tpoi.*2.0))./(t2.*tpoi)).*(A-1.0).^2)./(t2-tpoi.*2.0)-(A.*t1.*exp((t.*(t1-tpoi))./(t1.*tpoi)).*2.0)./(t1-tpoi)+(tpoi.^2.*(t1.^2.*tpoi.^2.*3.0+t1.^3.*t2-t1.*tpoi.^3.*2.0-t1.^3.*tpoi-t2.*tpoi.^3.*2.0+t1.*t2.*tpoi.^2.*5.0-t1.^2.*t2.*tpoi.*4.0+A.*t1.^2.*t2.^2.*2.0+A.^2.*t1.*t2.^3+A.^2.*t1.^3.*t2-A.*t1.^2.*tpoi.^2.*4.0+A.*t2.^2.*tpoi.^2.*2.0-A.^2.*t1.^3.*tpoi-A.^2.*t2.^3.*tpoi-A.^2.*t1.^2.*t2.^2.*2.0+A.^2.*t1.^2.*tpoi.^2+A.^2.*t2.^2.*tpoi.^2-A.*t1.^3.*t2.*2.0+A.*t1.^3.*tpoi.*2.0+A.*t1.*t2.*tpoi.^2.*2.0-A.*t1.*t2.^2.*tpoi.*5.0+A.*t1.^2.*t2.*tpoi.*3.0-A.^2.*t1.*t2.*tpoi.^2.*2.0+A.^2.*t1.*t2.^2.*tpoi+A.^2.*t1.^2.*t2.*tpoi).*2.0)./((-t1.*t2+t1.*tpoi+t2.*tpoi).*(t1.*tpoi.*-3.0+t1.^2+tpoi.^2.*2.0).*(t2.*tpoi.*-3.0+t2.^2+tpoi.^2.*2.0))+(A.*t1.*t2.*exp(-(t.*(-t1.*t2+t1.*tpoi+t2.*tpoi))./(t1.*t2.*tpoi)).*(A-1.0).*2.0)./(-t1.*t2+t1.*tpoi+t2.*tpoi));
    end

    % set / calculate pool parameters  
    a = gPar(1);
    b = gPar(2);
    r = gPar(3);
    C = sqrt(-4*a*b + (a+b+a*r)^2);
    t1 = 2/(a+b+a*r+C);
    t2 = 2/(a+b+a*r-C);
    A  = -(a-b+a*r-C)/(2*C);

    % get data info
    T = numel(t);
    N = numel(tau);
    
    % init vars
    yM = NaN(T,N);
    yH = NaN(T,N);

    % calculate light / medium / heavy labeling for each tau
    for n=1:N
        yM(:,n) = fM(tau(n));
        yH(:,n) = fH(tau(n));
    end
    
    yL = 1-yM-yH;
    
end