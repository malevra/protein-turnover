% generate header-only XLS file based on user input to fill in measurements

function genTemplateXLS (filename, sC, Times, nB, nM)
    
    nC = numel(sC);
    nT = numel(Times);

    C  = repmat(reshape(sC,1,[]),   [nT*nB*nM, 1]); C=reshape(C,[],1);
    T  = repmat(num2cell(Times),    [nB*nM,nC]);    T=T(:);
    B  = repmat(num2cell(1:nB),     [nM,nC*nT]);    B=B(:);
    M  = repmat(num2cell(1:nM)',    [1,nC*nT*nB]);  M=M(:);
    
    N = numel(C);
    
    Tc = cell([4 N+1]);
    Tc(1:4,1) = {'condition';'pulsing time';'biol. replicate';'mach. replicate'};
    Tc(1:4,2:N+1) = [C T B M]';

    Tt = cell2table(Tc);
    
    % (xlswrite has some problems in Linux)
    writetable(Tt,filename,'WriteVariableNames',false);

end