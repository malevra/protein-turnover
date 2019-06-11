% compare residuals for different global pool parameters (from saved
% sessions and optionally from extra pool control experiment)
%
% in: 
%
%       ses1{1,P}       cell array of P sessions exported from from turnoverGUI, in 
%                       which the lysine pool was fitted using different datasets
%                       to create such an array, use: {ses1,ses2, ...} for
%                       a set of individually exported sessions
%       gPar_extra      (optional) additional pool parameter [a,b,r]
%                       obtained externally (e.g. from control pool
%                       experiment, fitted by fitDirectPool)
%       
% 
% out:
%       R  (P,N)        matrix of normalized residuals for N datasets,
%                       using P given pool parameters, compared to the 
%                       residuals using their "own" pool parameters. Use 
%                       these values to decide whether the residuals are 
%                       only slightly larger (1 < R <~ 2) so that a common 
%                       pool can be used, or significantly larger (R >~ 2) 
%                       where it is advisable to use individual pool 
%                       parameters for different conditions
%       r  (P,N)        matrix of N fitted datasets using P different pools
%
%
% example:
%
%   [R,r] = comparePools({ses1,ses2,ses3});
%
%                       compares data residuals with lysine pools from ses1,
%                       ses2, ses3. results in (3,3) matrix if pools are
%                       obtained from different data sets
%
%   [R,r] = comparePools({ses1,ses2},[0.08 0.3 8]);
%                                            
%                       compares data residuals with lysine pools from ses1,
%                       ses2 and from extra pool parameters. results in
%                       (3,2) matrix (3 pools for 2 datasets)
%
% Mihai Alevra, 2019

function [R, r] = comparePools (sessions, gPar_extra)

    % note: sessions is a cell because user might combine sessions at
    % different stages, e.g. with 

    % get conditions used for global fits
    for n=1:numel(sessions)
        cndI(n) = sessions{n}.gfInfo.cond;
        p(n) = sessions{n}.pool;
    end
    
    % add external pool parameter if given
    if nargin > 1
        p(end+1)=p(1);  % copy field structure
        p(end).gPar = gPar_extra;  % set gPar from given directPool
    end
   
    d = sessions{1}.data;   % assumed to be identical for all sessions    
    cndI = unique(cndI);
    iSel = ismember(d.cndI,cndI);   % conditions to check
    
    d.cnd = d.cnd(cndI);    
    d.t   = d.t(iSel);
    d.sNr = d.sNr(iSel);
    d.mNr = d.mNr(iSel);
    d.y   = d.y(:,iSel);
    d.cndI= d.cndI(iSel);
    [~,d.cndI] = ismember(d.cndI,cndI);
    
    for n=1:numel(p)
        fprintf('fitting subset with pool par %d/%d [a=%6.2f b=%6.2f r=%6.2f]\n',...
                 n,numel(p),p(n).gPar(1), p(n).gPar(2), p(n).gPar(3));
        r(n,:) = fitSet(d,p(n));
    end
    
    for np=1:numel(p)
        for nd=1:numel(cndI)
            R(np,nd) = nansum(r(np,nd).fit.RNrm) / nansum(r(nd,nd).fit.RNrm);            
            txtDat{nd} =  ['data_' num2str(nd)];
        end
        txtPar{np} =  ['par_'  num2str(np)];
    end
    
    % display datasets that were compared    
    fprintf('\ndatasets used:\n\n');
    disp([txtDat' d.cnd']);
    
    fprintf('\nnorm. residuals:\n\n');
                
    T = array2table(R,'VariableNames',txtDat,'Rownames',txtPar);
    disp(T);
    
end