% exports fit result of a dataset to Excel file
%
% input:
%
% r (1,C)           : fit result structure (one for each unique condition)
%  .cnd             : name of current condition
%  .t   (1, T)      : times after heavy-isotope feeding ("pulse time" in days)
%  .sNr (1, T)      : biological replicate number (sample number)
%  .y   (N, T)      : measured labeled fraction for N proteins in T datasets
%  .fit             : fit data structure (see file fitSingle.m)
% filename          : name of .xls / .xlsx file
% tau2half          : true if time constants should be written as half-lives (*log(2))

function fitSetWriteXLS (r, pID, filename, tau2half)

    nPrt = size(r(1).y,1);          % number of proteins (assumed identical for all sets!)
    nCnd  = numel(r);               % number of conditions
    nCol = 4;                       % number of columns per condition
    A = cell(2+nPrt,1+nCnd*nCol);   % empty cell array for output
    
    if nargin<4
        tau2half = false;
    end
        
    if tau2half
        timeStr = 'tHalf';
        H = log(2);
    else
        timeStr = 'tau';
        H = 1;
    end
    
    A(3:nPrt+2,1)   = pID;          % protein IDs
    
    for c=1:nCnd
        o = (c-1)*nCol+1;             % offset        
        % headers
        A{1,1}          = 'condition';
        A{2,1}          = 'variable';
        A(1,o+(1:nCol)) = repmat(r(c).cnd,[1 nCol]);
        A{2,o+1}        = timeStr;
        A{2,o+2}        = [timeStr '_c1'];
        A{2,o+3}        = [timeStr '_c2'];
        A{2,o+4}        = 'RNorm';
        % data
        A(3:nPrt+2,o+1) = num2cell(H * r(c).fit.tau);                   % time constand / half-life
        A(3:nPrt+2,o+2) = num2cell(H * r(c).fit.ci_tau(:,1));           % lower confidence bound
        A(3:nPrt+2,o+3) = num2cell(H * r(c).fit.ci_tau(:,2));           % upper confidence bound
        A(3:nPrt+2,o+4) = num2cell(r(c).fit.RNrm);                  % squared sum of residuals
    end
    
    T = cell2table(A);    
    writetable(T,filename,'WriteVariableNames',false);

end