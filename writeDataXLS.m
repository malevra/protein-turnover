function writeDataXLS(d, filename)

    N = numel(d.cndI);
    P = numel(d.pID);
    
    Tc = cell([4+P N+1]);
    Tc(1:4,1) = {'condition';'pulsing time';'biol. replicate';'mach. replicate'};
    Tc(4+(1:P),1) = d.pID;
    
    Tc(1,1+(1:N)) = d.cnd(d.cndI);
    Tc(2,1+(1:N)) = num2cell(d.t);
    Tc(3,1+(1:N)) = num2cell(d.sNr);
    Tc(4,1+(1:N)) = num2cell(d.mNr);
    
    Tc(4+(1:P),1+(1:N)) = num2cell(d.y ./ (1-d.y));

    Tt = cell2table(Tc);
    
    % (xlswrite has some problems in Linux)
    writetable(Tt,filename,'WriteVariableNames',false);

end