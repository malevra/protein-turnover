% 
% estimates pool parameters gPar [a,b,r] from direct measurements of free lysine
% pool
%
% in:
%       t           time points for measurements
%       y           measured amplitudes (as heavy/total)
%       gPar0       starting values for [a,b,r]
%
% out:
%       gPar        fitted pool parameters [a,b,r]
%                   or as pool structure:
%   g.
%       gPar        [a,b,r] fit result
%       conf        confidence intervals for gPar
%       tau0        default tau0 (for compability with other scripts)


function [gPar, g] = fitDirectPool(t,y,gPar0)

    % pool function
    function yp = poolfun(g)
        A = g(1);
        B = g(2);
        R = g(3);
        C = sqrt(-4*A*B + (A+B+A*R)^2);
        T1 = 2/(A+B+A*R+C);        
        T2 = 2/(A+B+A*R-C);
        Amp  = -(A-B+A*R-C)/(2*C);
        yp = 1 - Amp*exp(-t ./ T1) - (1-Amp)*exp(-t ./ T2);
    end

    % residuals
    function df = difPool(g)
        df = y - poolfun(g);
    end


    % entry point
    opts = optimset('Display','off');
    [g.gPar, g.RNrm, g.R,~,~,~, g.J] = lsqnonlin (@difPool, gPar0, [1e-3 1e-3 1e-3], [inf inf inf], opts);    
    
    g.conf = nlparci(g.gPar, g.R, 'jacobian', g.J,'alpha',.05);
    gPar = g.gPar;

end