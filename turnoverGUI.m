% graphical user interface for protein turnover interpretation.
%
% Use without parameter for new session, or continue previously saved
% session, e.g. " turnoverGUI(turnoverSession) "
%
%
% saving the current session (rightmost tab) saves session variable to
% workspace, using possibly following fields:
%
% variable name    (size)   explanation
%
% turnoverSession           session variable (struct)
%
%   .pool                   lysine pool parameters:
%     .gPar         (1,3)   [a,b,r], a=protein synthesis/degregation rate,
%                           b=lysine uptake/removal rate, r=bound/free pool
%                           size ratio
%     .tau0                 initial time constant for each fit
%     .t          (1,200)   example time points for plotting free lysine pool curve
%     .y          (1,200)   heavy/total time evolution for free lysine pool
%     .yC         (1,200)   confidence bound for free lysine pool curve
%     .iter                 number of global pool fit iterations till convergence
%     .RNrm                 total squared sum of residuals from all proteins
%     .R            (P,N)   residuals of P valid proteins in N measurements
%     .conf         (3,2)   [a,b,r] confidence intervals (lower and upper)
%     .J            (P,3)   Jacobian from fit for P proteins (for parameters [a,b,r])
%     
%   .datafile               name of imported .xlsx file
%   .data                   imported H/total data, containing:
%     .cnd          {1xC}   cell array of C strings describing different
%                           experimental conditions
%     .cndI         (1,M)   index assining all M measurements the correct
%                           condition (.cnd{i})
%     .t            (1,M)   pulsing time for all M measurements
%     .sNr          (1,M)   biol. replicate number for all M measurements
%     .pID          {P,1}   protein names (P_ID) for all P proteins
%     .y            (P,M)   measured heavy/total ratio for P proteins at M
%                           measurements
%
%   .gfInfo                 struct containing info about global fit selection
%     .cond                 index of condition used for global fit
%     .minN                 minimum valid measurements for protein to be 
%                           used in global fit
%     .pOk          (p,1)   indices of p proteins selected by criterium above
%
%   .dsFit          (1,C)   array structs (1 for each condition) containing
%                           turnover data & fit results:
%     .cnd                  name of condition
%     .t            (1,m)   pulse times of each measurement m
%     .sNr          (1,m)   biological replicate number for each m
%     .y            (P,m)   measured heavy/total labeling for each P and m
%     .fit                  struct containing fit results:
%       .tau        (P,1)   turnover time constants for each P
%       .fY         (P,m)   predicted heavy/total labeling by fit
%       .R          (P,m)   residuals between measurements and fit
%       .RNrm       (P,1)   squared sum of residuals (all m) for each P
%       .ci_tau     (P,2)   confidence bounds for .tau
%       .ci_fY      (P,m,2) lower/upper prediction interval
%       .a .b .r .C .t1 .t2 .A      pool (derived) parameters used
%   .dsMaxT                 upper tau limit for plotting results
%
%   .st                     struct containing statistical dataset comparison:
%     .condI        (1,2)   indices of the two conditions compared
%     .cond         {1,2}   names of the two conditions compared
%     .tau1         (P,1)   time constants for P and condition 1
%     .tau2         (P,1)   time constants for P and condition 2
%     .P            (P,1)   p-value for significance test
%     .F            (P,1)   fold-increase between condition 1 and 2
%     .thP                  p-val threshold to select significant pairs
%     .thF                  fold-increase threshold to select significant pairs
%     .iSort        (I,1)   indices of proteins selected by thresholds and 
%                           sorted by fold-increase
%   .stMaxP                 threshold for p-value selection
%   .stMaxF                 threshold for fold-increase selection
%
% for detailed usage, please refer to the published protocol
% 
% Mihai Alevra, 2019

function turnoverGUI(s)
   

    % entry point: check matlab/toolbox, input, setup variables & GUI
    
    % check for matlab version / toolboxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if verLessThan('Matlab','8.4')
        warndlg('expected Matlab version > 8.4 (R2014b), this version is not tested!');
    end
    
    try 
        if verLessThan('optim','7.1')
            warndlg('expected Optimization Toolbox ver. > 7.1, this version is not tested!');
        end
        
    catch
        errordlg('Optimization Toolbox not found! aborting...');
    end
    
    try 
        if verLessThan('stats','9.1')
            warndlg('expected Optimization Toolbox ver. > 9.1, this version is not tested!');
        end
        
    catch
        errordlg('Optimization Toolbox not found! aborting...');
    end
        
    % end check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    g = []; % global gui structure
    initUI; % create GUI window / elements
    
    
    % check if a session structure is supplied: update UI with data
    if exist('s','var') && ~isempty(s)
        updateUI;
    else
        s = []; % session structure
        % default pool parameters from Fornasiero et al. 2018
        s.pool.gPar = [0.0342767364, 0.4448654275, 11.8365728939];
        s.pool.tau0 = 10; % default start value for POI time constant
    end

    % general UI setup
    function initUI
        
        % depends on OS! 'FixedWidth' Font does not work on Linux.
        if isunix
            FixWidthFont = 'Monospaced';
        elseif ismac
            FixWidthFont = 'FixedWidth';
        elseif ispc
            FixWidthFont = 'FixedWidth';
        end
        
        % generate UI figure and tabs
        g.fig = figure('units','pixels','position',[100 100 800 400],...
            'menubar','none','name','turnover GUI','numbertitle','off',...
            'WindowButtonMotionFcn',@mouseMove,'WindowButtonUpFcn',@mouseUp);
        g.mouseDown = false;
        g.mouseObj = 0;
        g.tg  = uitabgroup(g.fig,'Position',[0 0 1 1]);
        g.tab(1) = uitab(g.tg,'Title','template');  %,'backgroundcolor',[1 1 1]);
        g.tab(2) = uitab(g.tg,'Title','import');
        g.tab(3) = uitab(g.tg,'Title','pool fit');
        g.tab(4) = uitab(g.tg,'Title','dataset fit');
        g.tab(5) = uitab(g.tg,'Title','statistics');
        g.tab(6) = uitab(g.tg,'Title','session');
        
        g.sesSave = uicontrol('parent',g.tab(6),'units','norm',...
            'position',[.825 .05 .125 .1],'style','pushbutton',...
            'string','save session','Callback',@saveSession);
        g.sesName = uicontrol('parent',g.tab(6),'units','norm',...
            'position',[.550 .05 .250 .1],'style','edit','String','turnoverSession');
        g.sesNameL= uicontrol('parent',g.tab(6),'units','norm',...
            'position',[.550 .15 .250 .05],'style','text','String','variable name');
        
        g.tg.SelectedTab = g.tab(2);
        
        % elements for xls generation tab
        g.xlsCondL      = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.05 .7 .25 .1],'style','text','String','Condition(s)');
        g.xlsCond       = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.05 .4 .25 .3],'style','edit','Max',2);
        g.xlsDaysL      = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.35 .7 .25 .1],'style','text','String','Pulse times [days]');
        g.xlsDays       = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.35 .4 .25 .3],'style','edit','Max',2);
        g.xlsRepL       = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.65 .7 .30 .1],'style','text','String','Number of Replicates');
        g.xlsBRepL      = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.65 .5 .125 .1],'style','text','String','biological');
        g.xlsMRepL      = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.825 .5 .125 .1],'style','text','String','machine');
        g.xlsBRep       = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.65 .4 .125 .1],'style','edit','String','1');        
        g.xlsMRep       = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.825 .4 .125 .1],'style','edit','String','1');
        g.xlsGenP       = uicontrol('parent',g.tab(1),'units','norm',...
            'position',[.825 .05 .125 .1],'style','pushbutton',...
            'string','generate','Callback',@xlsGen);
        
        % elements for import tab
        g.impLoad       = uicontrol('parent',g.tab(2),'units','norm',...
            'position',[.825 .05 .125 .1],'style','pushbutton',...
            'string','import XLS','Callback',@impLoad);
        g.impInfoL      = uicontrol('parent',g.tab(2),'units','norm',...
            'position',[.25 .7 .5 .1],'style','text','String','summary of imported data');
        g.impInfo       = uicontrol('parent',g.tab(2),'units','norm',...
            'position',[.25 .2 .5 .5],'style','text','String','','backgroundcolor',[1 1 1],...
            'HorizontalAlignment','left','FontName',FixWidthFont);
        
        % elements for global fit tab
        %       filtering of proteins
        g.gfDSelL       = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.025 .85 .25 .05],'style','text','String','data selection');
        g.gfDSelCondL   = uicontrol('parent',g.tab(3),'units','norm','enable','inactive',...
            'position',[.025 .75 .25 .05],'style','text','String','condition');
        g.gfDSelCond    = uicontrol('parent',g.tab(3),'units','norm','callback',@gfDSel,...
            'position',[.025 .70 .25 .05],'style','popupmenu','String',{''});
        g.gfDSelProtL   = uicontrol('parent',g.tab(3),'units','norm','enable','inactive',...
            'position',[.025 .60 .25 .05],'style','text','String','min. datapoints [# proteins]');
        g.gfDSelProt    = uicontrol('parent',g.tab(3),'units','norm','callback',@gfDSel,...
            'position',[.025 .55 .25 .05],'style','popupmenu','String',{'all'});
        %       parameter initial values
        g.gfParsL = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.325 .85 .25 .05],'style','text','String','pool parameters');
        g.gfParsAL      = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.325 .75 .12 .05],'style','text','String','a');
        g.gfParsA       = uicontrol('parent',g.tab(3),'units','norm','callback',@gfPars,...
            'position',[.325 .70 .12 .05],'style','edit','String','0.034277');
        g.gfParsBL      = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.325 .60 .12 .05],'style','text','String','b');
        g.gfParsB       = uicontrol('parent',g.tab(3),'units','norm','callback',@gfPars,...
            'position',[.325 .55 .12 .05],'style','edit','String','0.444865');
        g.gfParsRL      = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.325 .45 .12 .05],'style','text','String','r');
        g.gfParsR       = uicontrol('parent',g.tab(3),'units','norm','callback',@gfPars,...
            'position',[.325 .40 .12 .05],'style','edit','String','11.836573');
        %       derived parameters
        g.gfParsT1L     = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.455 .75 .12 .05],'style','text','String','tau 1');
        g.gfParsT1      = uicontrol('parent',g.tab(3),'units','norm','enable','inactive',...
            'position',[.455 .70 .12 .05],'style','edit','String','1.000','backgroundcolor',[.8 .8 .8]);
        g.gfParsT2L     = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.455 .60 .12 .05],'style','text','String','tau 2');
        g.gfParsT2      = uicontrol('parent',g.tab(3),'units','norm','enable','inactive',...
            'position',[.455 .55 .12 .05],'style','edit','String','1.000','backgroundcolor',[.8 .8 .8]);
        g.gfParsAmpL    = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.455 .45 .12 .05],'style','text','String','A');
        g.gfParsAmp     = uicontrol('parent',g.tab(3),'units','norm','enable','inactive',...
            'position',[.455 .40 .12 .05],'style','edit','String','1.000','backgroundcolor',[.8 .8 .8]);
        %       pool plot for current parameters
        g.gfPlotL       = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.625 .85 .35 .05],'style','text','String','pool plot');
        g.gfPlotAx      = axes('parent',g.tab(3),'box','on',...
            'position',[.650 .30 .30 .50],'ylim',[0 1]);
        xlabel(g.gfPlotAx,'time [d]'); ylabel(g.gfPlotAx,'amplitude'); hold(g.gfPlotAx,'on');
        g.gfPlotPlC1    = plot(g.gfPlotAx,NaN,NaN,'LineWidth',1,'color',[.8 .5 .5 .5]); 
        g.gfPlotPlC2    = plot(g.gfPlotAx,NaN,NaN,'LineWidth',1,'color',[.8 .5 .5 .5]);
        g.gfPlotPl      = plot(g.gfPlotAx,NaN,NaN,'LineWidth',2,'color',[.1 .1 .5]);
        %       message log/table for fitting
        g.gfLogL        = uicontrol('parent',g.tab(3),'units','norm','HorizontalAlignment','left',...
            'position',[.025 .325 .10 .05],'style','text','String','fit iterations');

        g.gfLogT        = uitable('parent',g.tab(3),'units','norm','ColumnName',{'time','a','b','r','residuals'},...
            'position',[.025 .025 .55 .30], 'data',cell(0,5),'ColumnWidth',{'auto','auto','auto','auto'},...
            'ColumnFormat',{'char', 'numeric', 'numeric', 'numeric', 'numeric'});
        
        %       POI initial time constant for fitting
        g.gfParsT0L     = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.625 .10 .12 .05],'style','text','String','POI tau0');
        g.gfParsT0      = uicontrol('parent',g.tab(3),'units','norm','callback',@gfPars,...
            'position',[.625 .05 .12 .05],'style','edit','String','10.0');
        %       fit start button
        g.gfStart       = uicontrol('parent',g.tab(3),'units','norm',...
            'position',[.825 .05 .125 .1],'style','pushbutton',...
            'string','global fit','Callback',@gfStart);
        
        % elements for dataset fit
        
        g.dsExport      = uicontrol('parent',g.tab(4),'units','norm',...
            'position',[.825 .05 .125 .1],'style','pushbutton',...
            'string','export fit','Callback',@dsExport);
        g.dsFit         = uicontrol('parent',g.tab(4),'units','norm',...
            'position',[.675 .05 .125 .1],'style','pushbutton',...
            'string','start fit','Callback',@dsStart);
        g.dsViewSel     = uicontrol('parent',g.tab(4),'units','norm',...
            'position',[.35 .05 .2 .05],'style','popupmenu',...
            'string',{''},'Callback',@dsView);
        g.dsViewTauOrHL = uicontrol('parent',g.tab(4),'units','norm',...
            'position',[.35 .1 .1 .05],'style','checkbox','Value',1,...
            'string','half-life','Callback',@dsView);
        g.dsViewXLimL    = uicontrol('parent',g.tab(4),'units','norm',...
            'position',[.225 .10 .10 .05],'style','text',...
            'string','x axis limit');
        g.dsViewXLim    = uicontrol('parent',g.tab(4),'units','norm',...
            'position',[.225 .05 .10 .05],'style','edit',...
            'string','50','Callback',@dsView);
        g.dsFig1        = axes('parent',g.tab(4),'units','norm',...
            'position',[.1 .30 .35 .60],'box','on');
        xlabel(g.dsFig1,'half-life [d]'); ylabel(g.dsFig1,'amplitude');
        g.dsFig2        = axes('parent',g.tab(4),'units','norm',...
            'position',[.55 .30 .35 .60],'box','on');
        xlabel(g.dsFig2,'half-life [d]'); ylabel(g.dsFig2,'count');
        
        % elements for statistics (comparison of two conditions)
        
        g.stExport      = uicontrol('parent',g.tab(5),'units','norm',...
            'position',[.825 .05 .125 .1],'style','pushbutton',...
            'string','export stats','Callback',@stExport);
        g.stExportSel   = uicontrol('parent',g.tab(5),'units','norm',...
            'position',[.825 .15 .125 .05],'style','checkbox',...
            'string','threshold','Value',0);
        g.stCalc        = uicontrol('parent',g.tab(5),'units','norm',...
            'position',[.675 .05 .125 .1],'style','pushbutton',...
            'string','compare','Callback',@stCompare);
        g.stCond1       = uicontrol('parent',g.tab(5),'units','norm',...
            'position',[.500 .225 .225 .05],'style','popupmenu','string',{''});
        g.stCond2       = uicontrol('parent',g.tab(5),'units','norm',...
            'position',[.725 .225 .225 .05],'style','popupmenu','string',{''});
        g.stProtT        = uitable('parent',g.tab(5),'units','norm','ColumnName',{'pID','lt','log2(f)','-log10(p)'},...
            'position',[.5 .325 .45 .575], 'data',cell(0,4),'ColumnWidth',{'auto','auto','auto','auto'},...
            'ColumnFormat',{'char', 'numeric', 'numeric', 'numeric'});
        
        g.stAx          = axes('parent',g.tab(5),'units','norm',...
            'position',[.075 .15 .400 .75],'box','on','ButtonDownFcn',@mouseDown); hold(g.stAx,'on'); 
        xlabel(g.stAx,'log2(f)'); ylabel(g.stAx,'-log10(p)');
        g.stLineF     = plot(g.stAx,NaN,NaN,'linewidth',3,'color',[0  0 .5 .5],'hittest','off');
        g.stLineP     = plot(g.stAx,NaN,NaN,'linewidth',3,'color',[0 .5  0 .5],'hittest','off');
        g.stSelPos    = plot(g.stAx,NaN,NaN,'o','color',[ 0 .5  0 ],'hittest','off');
        g.stSelNeg    = plot(g.stAx,NaN,NaN,'o','color',[.5  0  0 ],'hittest','off');        
        g.stPoints    = plot(g.stAx,NaN,NaN,'.','color',[ 0  0  0 ],'MarkerSize',.5,'hittest','off');
        
        g.stMxL         = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.525 .15 .05 .05],'style','text','String','sel');
        g.stThL         = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.600 .15 .05 .05],'style','text','String','max');
        g.stThFoldL      = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.500 .10 .025 .05],'style','text','String','f');
        g.stThFold      = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.525 .10 .050 .05],'style','edit','String','0.5');
        g.stMxFold      = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.600 .10 .050 .05],'style','edit','String','1.5');
        g.stThPValL     = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.500 .05 .025 .05],'style','text','String','p');
        g.stThPVal      = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.525 .05 .05 .05],'style','edit','String','10.0');
        g.stMxPVal      = uicontrol('parent',g.tab(5),'units','norm','callback',@stCompare,...
            'position',[.600 .05 .05 .05],'style','edit','String','50.0');        
    end

    % updates UI elements using entries in session structure s
    function updateUI
        if isfield(s,'data') && ~isempty(s.data)            
            g.impInfo.String = ...
                {['conditions :   ' sprintf('%s ',s.data.cnd{:})],...
                 ['times      :   ' sprintf('%d ',unique(s.data.t))],...
                 ['proteins   :   ' sprintf('%d ',numel(s.data.pID))],...
                 '','data successfully imported.'};
            g.gfDSelCond.String = s.data.cnd;
            if isfield(s,'gfInfo')
                g.gfDSelCond.Value  = s.gfInfo.cond;
            end
            g.gfDSelProt.Enable = 'on';
            if isfield(s,'gfInfo')
                g.gfDSelProt.String = 1:sum(s.data.cndI==s.gfInfo.cond);
                g.gfDSelProt.Value  = s.gfInfo.minN;
            else
                g.gfDSelProt.String = 1:sum(s.data.cndI==1);
                g.gfDSelProt.Value  = sum(s.data.cndI==1);
            end
            
            g.gfDSelCond.Enable = 'on';
            if isfield(s,'pool')
                g.gfParsA.String    = num2str(s.pool.gPar(1));
                g.gfParsB.String    = num2str(s.pool.gPar(2));
                g.gfParsR.String    = num2str(s.pool.gPar(3));
                g.gfParsT0.String   = num2str(s.pool.tau0);
            end            
            
            g.dsViewSel.String = s.data.cnd;
            g.stCond1.String = s.data.cnd;
            g.stCond2.String = s.data.cnd;            
            
            gfDSel;
            updatePoolPlot;
        end
    end

    % selection of condition or proteins filtering changed
    function gfDSel(varargin)
        
        % update conditions / protein selection thresholds if data loaded
        if isfield(s,'data') && ~isempty(s.data)
    
            % if no global fit done yet, init variables
            if ~isfield(s,'gfInfo')
                s.gfInfo.cond = 0;  % updated below
                s.gfInfo.minN = 0;  % updated below
                s.gfInfo.pOk = [];  % updated below
            end
            
            % if selected condition changed (or not initialized), refresh
            % possible protein selection thresholds
            if s.gfInfo.cond ~= g.gfDSelCond.Value      % selected codition changed
                s.gfInfo.cond = g.gfDSelCond.Value;     % save new selection
                g.gfDSelProt.String = 1:sum(s.data.cndI==s.gfInfo.cond); % update possible protein number thresholds
                g.gfDSelProt.Value  =   sum(s.data.cndI==s.gfInfo.cond); % select maximum (tightest threshold)  
            end
            
            % in any case, refresh protein selection threshold
            s.gfInfo.minN = g.gfDSelProt.Value;     % selected minimum number of measurements per protein
            
            % filter only proteins that are measured at least minN times
            s.gfInfo.pOk  = find(sum(isfinite(s.data.y(:,s.data.cndI==s.gfInfo.cond)),2)>=s.gfInfo.minN);
            % show in label how many proteins were actually selected using
            % the current threshold
            g.gfDSelProtL.String = sprintf('min. data [%d prot]',numel(s.gfInfo.pOk));
        end
    end

    % update pool parameters from GUI
    function gfPars(varargin)
        s.pool.gPar = [str2double(g.gfParsA.String) ...
                       str2double(g.gfParsB.String) ...
                       str2double(g.gfParsR.String)];
        s.pool.tau0 =  str2double(g.gfParsT0.String);
        updatePoolPlot;
    end

    % compare two conditions
    function stCompare(varargin)
        C1 = g.stCond1.Value;                   % first condition selected
        C2 = g.stCond2.Value;                   % second condition selected
        
        if C1>0 && C2>0 && C1~=C2               % if valid and not identical
            title(g.stAx,sprintf('%s | %s',s.data.cnd{C1},s.data.cnd{C2}),'interpreter','none','FontWeight','normal');
        end
        
        if isfield(s,'dsFit')                   % only compare if datasets already have fit results
            tau1 = s.dsFit(C1).fit.tau;                             % time constant for condition 1
            tau2 = s.dsFit(C2).fit.tau;                             % time constant for condition 2
            se1  = diff(s.dsFit(C1).fit.ci_tau,[],2) / (2*1.96);    % standart deviation 1 from confidence interval 
            se2  = diff(s.dsFit(C2).fit.ci_tau,[],2) / (2*1.96);    % standart deviation 2 from confidence interval 
            se12 = sqrt(se1.^2+se2.^2);                             % combined standard deviation
            
            Z    = abs(tau1 - tau2)./se12;                          % z statistic from absolute difference in time constant
            P    = exp(-0.717.*Z-0.416.*Z.^2);                      % p-val from z statistic
            
            % store details about this comparison
            s.st.cndI = [C1 C2];
            s.st.cond = s.data.cnd([C1 C2]);
            s.st.tau1 = tau1;
            s.st.tau2 = tau2;
            s.st.P    = P;
            
            % update F and P thresholds from GUI
            s.st.thF = str2double(g.stThFold.String);
            s.st.thP = str2double(g.stThPVal.String);            
            
            % calculate fold-increase
            s.st.F   = s.st.tau2./s.st.tau1;
            
            % select by thresholding in F and P
            iSelP = find(log2(s.st.F) >  s.st.thF & -log10(s.st.P) > s.st.thP);
            iSelN = find(log2(s.st.F) < -s.st.thF & -log10(s.st.P) > s.st.thP);
                     
            % sort by fold-increase
            [~,iSortP] = sort(log2(s.st.F(iSelP)),'descend');
            [~,iSortN] = sort(-log10(-s.st.F(iSelN)),'descend');
            
            % store selected and sorted indices (first positive, then
            % negative fold-increase)
            s.st.iSort = [iSelP(iSortP); iSelN(iSortN)];
            
            % fill table
            data = cell(numel(s.st.iSort),4);
            if numel(s.st.iSort)>0
                data(:,1) = s.data.pID(s.st.iSort);
                data(:,2) = mat2cell(s.st.tau1(s.st.iSort),ones(1,numel(s.st.iSort)));
                data(:,3) = mat2cell(log2(s.st.F   (s.st.iSort)),ones(1,numel(s.st.iSort)));
                data(:,4) = mat2cell(-log10(s.st.P   (s.st.iSort)),ones(1,numel(s.st.iSort)));
            end
            g.stProtT.Data = data;
                     
            % volcano plot
            s.stMaxP = str2double(g.stMxPVal.String);
            s.stMaxF = str2double(g.stMxFold.String);

            set(g.stLineP,'XData',[-1 1]*s.stMaxF,'YData',[1 1]*s.st.thP);
            set(g.stLineF,'XData',[-1 -1 NaN 1 1]*s.st.thF,'YData',[0 1 NaN 0 1]*s.stMaxP);

            % update data points/selection markers
            set(g.stPoints,'XData',log2(s.st.F),'YData',-log10(s.st.P));
            set(g.stSelPos,'XData',log2(s.st.F(iSelP)),'YData',-log10(s.st.P(iSelP)));
            set(g.stSelNeg,'XData',log2(s.st.F(iSelN)),'YData',-log10(s.st.P(iSelN)));
                         
            g.stAx.YLim = [0 s.stMaxP];
            g.stAx.XLim = [-s.stMaxF s.stMaxF];
            
        end
    end

    % handle mouse click (to move threshold lines)
    function mouseDown(h,ev)
        if h==g.stAx                % only do something in comparison plot
            g.mouseDown = true;     % mouse drag has started
            g.mouseObj = h;         % remember axis in which drag started
            mouseMove(h,ev);        % interpred mouse position
        else                        % if clicked somewhere else
            g.mouseDown = true;     % mouse down, but unimportant
            g.mouseObj = 0;         % remember wrong axes
        end
    end

    % handle mouse button release
    function mouseUp(~,~)
        g.mouseDown = false;        % possible mouse dragging stopped
        if g.mouseObj==g.stAx       % if drag was started in plot axes
            g.mouseObj = 0;         % reset dragged axes
            stCompare;              % update statistical comparison / selection
        end 
    end

    % handle mouse drag
    function mouseMove(~,~)
        if g.mouseDown && g.mouseObj~=0         % if dragging in right plot
            cp = get(g.stAx,'CurrentPoint');    % get mouse pos

            s.st.thP = abs(cp(1,2));            % update p-val threshold
            s.st.thF = abs(cp(1,1));            % update f-val threshold
            g.stThPVal.String = s.st.thP;       % update GUI with number
            g.stThFold.String = s.st.thF;       % update GUI with number

            % update threshold line positions
            set(g.stLineP,'XData',[-1 1]*s.stMaxF,'YData',[1 1]*s.st.thP);
            set(g.stLineF,'XData',[-1 -1 NaN 1 1]*s.st.thF,'YData',[0 1 NaN 0 1]*s.stMaxP);
            drawnow;                    
        end
    end

    % export results of comparison
    function stExport(varargin)
        if isfield(s,'st')
            if g.stExportSel.Value > 0      % if only thresholded values should be exported
                I = s.st.iSort;
            else                            % otherwise: export all proteins
                I = 1:numel(s.data.pID);
            end
            
            stExpFile = uiputfile({'*.xls;*.xlsx'},'save comparison data...','comparison.xlsx');
            if stExpFile ~= 0
                nPrt = numel(I);

                if g.dsViewTauOrHL.Value
                    tScale = log(2);
                    tScaleT= 'half_life';                
                else
                    tScale = 1;
                    tScaleT= 'time_constant';
                end                    
                
                A = cell(1+nPrt, 5);
                
                A{1,1} = 'pID';
                A{1,2} = [tScaleT ' (' s.st.cond{1} ')'];
                A{1,3} = [tScaleT ' (' s.st.cond{2} ')'];
                A{1,4} = 'log2(fold increase)';
                A{1,5} = '-log10(p-value)';
                
                A(1+(1:nPrt),1) = s.data.pID(I);
                A(1+(1:nPrt),2) = num2cell(s.st.tau1(I)*tScale);
                A(1+(1:nPrt),3) = num2cell(s.st.tau2(I)*tScale);
                A(1+(1:nPrt),4) = num2cell( log2 (s.st.F(I)));
                A(1+(1:nPrt),5) = num2cell(-log10(s.st.P(I)));
                
                T = cell2table(A);
                writetable(T,stExpFile,'WriteVariableNames',false);                
            end
        end
    end

    % fit whole dataset (all conditions, all proteins)
    function dsStart(varargin)
        if isfield(s,'data')
            g.dsFit.String = 'fitting...';      % update button text while fitting
            drawnow;
            s.dsFit = fitSet(s.data,s.pool);    % fit all turnover time constants (external function)
            g.dsFit.String = 'start fit';       % change back button text
            g.dsViewSel.String = s.data.cnd;
        end

        % use xlim entered in GUI
        s.dsMaxT = str2double(g.dsViewXLim.String);
        
        % plot labeling vs. tau
        dsView;
    end
    
    % change / update view of current condition, if fit data available
    function dsView(varargin)
        if isfield(s,'dsFit')
            n = g.dsViewSel.Value;
            tMax = str2double(g.dsViewXLim.String);
            if g.dsViewTauOrHL.Value
                tScale = log(2);
                tScaleT= 'half-life [d]';
            else
                tScale = 1;
                tScaleT= 'time constant [d]';
            end
            
            xlabel(g.dsFig1,tScaleT); ylabel(g.dsFig1,'rel. ampl.');
            xlabel(g.dsFig2,tScaleT); ylabel(g.dsFig2,'count');            
            
            T = unique(s.dsFit(n).t);
            colM = hsv(numel(T)*2);
                        
            cla(g.dsFig1); hold(g.dsFig1,'on');   
            cla(g.dsFig2); hold(g.dsFig2,'on');
            for t=1:numel(T)
                
                tSel = find(s.dsFit(n).t == T(t));
                tPlt = repmat(s.dsFit(n).fit.tau*tScale,[1 numel(tSel)]);
                YEx  = s.dsFit(n).fit.fY(:,tSel);
                YTh  = s.dsFit(n).y(:,tSel);
                
                plot(g.dsFig1, tPlt(:), YTh(:), '.','MarkerSize',1,'color',colM(t,:)*.8);
                plot(g.dsFig1, tPlt(:), YEx(:), '.','MarkerSize',1,'color',colM(t,:)*.4);
                
            end          
            
            histogram(g.dsFig2, s.dsFit(n).fit.tau*tScale, linspace(0,tMax,100),'edgecolor','none','facecolor',[.1 .1 .4]);
            g.dsFig1.XLim = [0 tMax]; g.dsFig1.YLim = [0 1];
            g.dsFig2.XLim = [0 tMax]; g.dsFig1.YLim = [0 1];
            
        end
    end

    % export dataset fits to table
    function dsExport(varargin)
        dsExpFile = uiputfile({'*.xls;*.xlsx'},'save dataset fit...','dataset.xlsx');
        if dsExpFile ~= 0            
            fitSetWriteXLS (s.dsFit, s.data.pID, dsExpFile, g.dsViewTauOrHL.Value);            
        end        
        disp(['exported dataset results to file ' dsExpFile]);
    end

    % save current progress to main workspace variable "turoverSession"
    function saveSession(varargin)
        assignin('base',g.sesName.String,s);
    end
    

    % calculate pool curve and derived parameters
    function [y,A,B,R,T1,T2,Amp] = poolFun(gPar,t)
        A = gPar(1);
        B = gPar(2);
        R = gPar(3);
        C = sqrt(-4*A*B + (A+B+A*R)^2);
        T1 = 2/(A+B+A*R+C);        
        T2 = 2/(A+B+A*R-C);
        Amp  = -(A-B+A*R-C)/(2*C);
        y = 1 - Amp*exp(-t ./ T1) - (1-Amp)*exp(-t ./ T2);
    end

    % plot free lysine pool curve using current pool parameters
    function updatePoolPlot(varargin)
        
        if isfield(s,'data') && ~isempty(s.data)
            s.pool.t = linspace(0,max(s.data.t)*1.5,200);
        else
            s.pool.t = linspace(0,30,200);
        end

        % calculate pool curve / derived parameters
        [s.pool.y,A,B,R,T1,T2,Amp] = poolFun(s.pool.gPar,s.pool.t);
        
        % update GUI
        g.gfParsA.String    = num2str(A);
        g.gfParsB.String    = num2str(B);
        g.gfParsR.String    = num2str(R);        
        g.gfParsT1.String   = num2str(T1);
        g.gfParsT2.String   = num2str(T2);
        g.gfParsAmp.String  = num2str(Amp);

        set(g.gfPlotPl,'XData',s.pool.t,'YData',s.pool.y);
        set(g.gfPlotAx,'XLim',[0 max(s.pool.t)]);
        
        % in addition, plot pool confidence interval if existent
        if isfield(s.pool,'yC') && any(isfinite(s.pool.yC(:)))
            set(g.gfPlotPlC1,'XData',s.pool.t, 'YData',s.pool.yC(1,:));
            set(g.gfPlotPlC2,'XData',s.pool.t, 'YData',s.pool.yC(2,:));
        else % no confidence bounds
            set(g.gfPlotPlC1,'XData',NaN, 'YData',NaN);
            set(g.gfPlotPlC2,'XData',NaN, 'YData',NaN);
        end        
    end    

    % global pool fit
    % fit all valid proteins, vary pool parameters until parameters with
    % minimal total residuals found
    function gfStart(varargin)
        
        if isfield(s,'data') && ~isempty(s.data)
            
            % select data to fit
            CondSel = s.data.cndI==s.gfInfo.cond;
            D.t = s.data.t (CondSel);
            D.y = s.data.y(s.gfInfo.pOk,CondSel);
            D.pID = s.data.pID (s.gfInfo.pOk);
            
            s.pool.iter = 0;    % reset current iteration
            g.gfLogT.Data = cell(0,5); % reset GUI table
            G = fitGlobalPar(D,s.pool,@fitStep); % call external global fit
            
            % store results
            s.pool.gPar = G.gPar;
            s.pool.tau0 = G.tau0;
            s.pool.RNrm = G.RNrm;
            s.pool.R    = G.R;
            s.pool.conf = G.conf;
            s.pool.J    = G.J;
                        
            % calculate prediction interval for pool curve using Jacobian
            % of fit
            [Ypred,delta] = nlpredci(@poolFun,s.pool.t,s.pool.gPar,s.pool.R,'Jacobian',s.pool.J);
            s.pool.yC   = [Ypred + delta; Ypred - delta];   % confidence bound for example pool curve
            
            updateUI;
            msgbox('optimization of pool parameters complete','pool fit');
            
        else
            error('no Data loaded for global fit');
        end
        
    end

    % callback for external global fit function, to update current fit
    % progress (called after each global fit step)
    function fitStep(o,residuals)
        s.pool.gPar = [o.a, o.b, o.r];
        s.pool.iter = s.pool.iter+1;
        if mod(s.pool.iter,4)==1 % "real" new fit steps come every 4 calls, others are for jacobian estimation...
            g.gfLogT.Data = cat(1, {timeMsg, o.a, o.b, o.r, sum(residuals(:).^2)}, g.gfLogT.Data);
            updatePoolPlot;
            drawnow;
        end
    end

    % generates XLS template for user to fill in datasets
    function xlsGen(varargin)
        % parse input fields (only simple error checking!)
        s.xlsGen.cond = cellstr(g.xlsCond.String);
        for n=1:size(g.xlsDays.String,1)
            s.xlsGen.time(n) = str2double(g.xlsDays.String(n,:));
        end
        s.xlsGen.bRep = str2double(g.xlsBRep.String);
        s.xlsGen.mRep = str2double(g.xlsMRep.String);
        
        if numel(s.xlsGen.cond)>0 && numel(s.xlsGen.time)>0 && all(isfinite(s.xlsGen.time)) ...
           && isfinite(s.xlsGen.bRep) && s.xlsGen.bRep > 0 ...
           && isfinite(s.xlsGen.mRep) && s.xlsGen.mRep > 0
       
            s.xlsGen.file = uiputfile({'*.xls;*.xlsx'},'save template...','template.xlsx');
            if s.xlsGen.file ~= 0
                genTemplateXLS(s.xlsGen.file,...
                               s.xlsGen.cond,...
                               s.xlsGen.time,...
                               s.xlsGen.bRep,...
                               s.xlsGen.mRep);
            end
        end
    end

    % imports data from XLS into current session
    function impLoad(varargin)
        [fname,fpath] = uigetfile({'*.xls;*.xlsx'},'load data...');
        if fname ~= 0
            s.datafile = [fpath fname];
            s.data = readDataXLS(s.datafile); 
            updateUI;
        end
    end

end