% simDataset
%
%
% script that generates artificial data for turnoverGUI to test its
% functions.
%
% it generates data for four conditions,
%
%   1. poolA_ctrl: pool parameters A, taus are gamma-distributed (~ 1..100 days)
%   2. poolA_mod : pool parameters A, modified distribution 
%                  (1/3 identical, 2/3 increased or decreased by factor of 2)
%   3. poolB_ctrl: pool parameters B, same taus as in 1.
%   4. poolB_mod : pool parameters B, same taus as in 2.

clear

% two distinct pools
p(1).tau0 = 10;
p(1).gPar = [.04 .4 8];
p(2).tau0 = 10;
p(2).gPar = [.08 .3 8];

% default set of lifetimes
tau = gamrnd(4,5,[1 300]);                  % gamma-distributed taus (~1-100, somewhat similar to experimental distribution)

% modified set of lifetimes
tau(2,:) = tau(1,:);
tau(2,2:3:end) = tau(2,2:3:end) * 2;        % modified: 1/3 doubled tau
tau(2,3:3:end) = tau(2,3:3:end) * 0.5;      % modified: 1/3 halved  tau
tau=tau';

% dummy protein names
pID=cell(size(tau,1),1);
for n=1:size(tau,1)
    pID{n} = sprintf('id_%03d',n);
end
pID=pID';

% text blocks for conditions and pools
ctxt = {'ctrl','mod'};                      % conditions
ptxt = {'poolA','poolB'};                   % pools

% columns of data
R = numel(ctxt)*numel(ptxt);                % total conditions
S = 2;                                      % biol. replicates
M = 3;                                      % mach. replicates

TT = [5 14 21];                             % time points

t   = repmat(TT, [S*M,         1]);         % total time vector
T   = numel(TT)*S*M;
sNr = repmat(1:S,[  M, numel(TT)]);
mNr = mod(0:T-1,M)+1;

t=t(:)'; sNr=sNr(:)';

d.cnd = {};
d.cndI = NaN(1,R*T);
d.t = NaN(1,R*T);
d.sNr = NaN(1,R*T);
d.mNr = NaN(1,R*T);
d.y = NaN(size(tau,1),R*T);

for m=1:numel(ptxt)                         % pools
    for n=1:numel(ctxt)                     % conditions
        r = (m-1)*numel(ptxt)+n;            % dataset index
        ni = (r-1)*T + (1:T);

        d.cnd{r} = [ctxt{n} '_' ptxt{m}];   % dataset text
        d.t(ni) = t;
        d.sNr(ni) = sNr;
        d.mNr(ni) = mNr;
        d.cndI(ni) = r;
        d.pID = pID;
        
        d.y(:,ni) = predictPulseChase(t,14,tau(:,n),p(m).gPar)';
    end
end

% make some multiplicative gaussian noise
d.y = d.y .* normrnd(1,.1,size(d.y));

% create basic session with this data
ts.data = d;
ts.pool.gPar = [.05 .5 10];
ts.pool.tau0 = 10;
ts.pool.t    = linspace(0,25,200);

clearvars -except ts d tau p