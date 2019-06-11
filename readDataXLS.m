% read xls / xlsx file of pulsed data experiment into a structure
function d = readDataXLS(file)

    [~,~,raw]           = xlsread(file);

    cnd                 = raw(1,2:end);
    cnd = cnd(cellfun('isclass', cnd, 'char'));

    N                   = numel(cnd)+1;
    [d.cnd, ~, d.cndI]  = unique(cnd,'stable'); d.cndI=d.cndI'; % unique conditions and corresponding data indices
    d.t                 = cell2mat(raw(2,2:N));                 % pulse times        
    d.sNr               = cell2mat(raw(3,2:N));                 % biological replicate
    d.mNr               = cell2mat(raw(4,2:N));                 % machine replicate
    d.pID               = raw(5:end,1);                         % protein IDs
    Y                   = raw(4+(1:numel(d.pID)),2:N);          % data range
    Yc                  = cellfun(@(x) ~isa(x,'double'),Y);     % check for possible 0x0 char reads
    Y(Yc)               = {NaN};                                % set those to NaNs
    d.y                 = cell2mat(Y);                          % heavy/light ratios,
    if ~isempty(d.y)
        d.y                 = d.y./(d.y+1);                     % converted to fraction of heavy
    end
    
end
