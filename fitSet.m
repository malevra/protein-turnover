% function r = fitSet(d,p)
%
% Fits a whole dataset of labeled isotopes after pulsing to obtain
% individual protein lifetimes. Machine replicates are condensed into
% average data before fitting. Different treatments/conditions are fitted
% individually.
%
% input:
%
% d                 : structure for experimental data:
%  .cnd (1, C)      : cell array of C unique conditions
%  .cndI(1, T)      : indices of conditions (.cnd) of T datasets
%  .t   (1, T)      : times after heavy-isotope feeding ("pulse time" in days)
%  .sNr (1, T)      : biological replicate number (sample number)
%  .pID (N, 1)      : cell array of N protein IDs
%  .y   (N, T)      : measured labeled fraction for N proteins in T datasets
%
% p                 : global parameters structure:
%  .gPar            : parameters for global
%         (1)       : (a) protein degradation rate
%         (2)       : (b) lysine feeding / excretion rate
%         (3)       : (r) ratio of protein / lysine pool size
%  .tau0            : starting value for fitting of protein time constants
% 
%
% output:
%
% r (1,C)           : fit result structure (one for each unique condition)
%  .cnd             : name of current condition
%  .t   (1, T)      : times after heavy-isotope feeding ("pulse time" in days)
%  .sNr (1, T)      : biological replicate number (sample number)
%  .y   (N, T)      : measured labeled fraction for N proteins in T datasets
%  .fit             : fit data structure (see file fitSingle.m)
%
% Mihai Alevra, 2019

function r = fitSet(d,p)

    % find machine replicates
    [uD,~,uI] = unique([d.cndI' d.t' d.sNr'],'rows','stable');
    
    % condense data by averaging over all machine replicates
    dA.cnd  = d.cnd;
    dA.cndI = uD(:,1)';
    dA.t    = uD(:,2)';
    dA.sNr  = uD(:,3)';
    dA.pID  = d.pID;
    dA.y    = NaN(numel(d.pID),size(uD,1));
    
    for n=1:size(uD,1)
        dA.y(:,n) = nanmean(d.y(:,uI==n),2);
    end

    % fitting loop, separately for each condition/treatment etc.

    for c = 1:numel(dA.cnd)
        
        iSel = dA.cndI==c;                  % indices for this condition
        
        r(c).cnd = dA.cnd(c);               % collect corresponding data
        r(c).t = dA.t(iSel);
        r(c).sNr = dA.sNr(iSel);
        r(c).y = dA.y(:,iSel);
        r(c).fit = fitSingle(r(c),p);       % fit protein lifetimes
        
    end

end