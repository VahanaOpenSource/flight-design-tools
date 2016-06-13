function [polMatName, pol] = runPolars(AFIL,varargin)
%function [polMatName, pol] = runPolars(AFIL,varargin)
% runPolars computes viscous 2D airfoil polars by writing/executing an Xfoil batch
% file and mining the data from polar accumulation files.
% Output polars are intended for use in 1,2 and 3D lookup tables that can
% have variables of AOA, Re#, and Mach

%% input parser ----------------------
ip = inputParser;
 addRequired(ip,'AFIL'       ,                           @ischar);
addParameter(ip,'aSweep'     , -6:1:18,                  @isnumeric);
addParameter(ip,'rSweep'     , linspace(.2,3,5)*1e6,     @isnumeric);
addParameter(ip,'Mach'       , 0,                        @isnumeric);
addParameter(ip,'nCrit'      , 9,                        @isnumeric);
addParameter(ip,'xTripTop'   , 1,                        @isnumeric);
addParameter(ip,'xTripBot'   , 1,                        @isnumeric);
addParameter(ip,'polDir'     , [pwd filesep 'polars'],   @ischar);
addParameter(ip,'xfoilDir'   , pwd,                      @ischar);
addParameter(ip,'keepBatch'  , false,                    @islogical);
addParameter(ip,'keepPolars' , false,                    @islogical);
addParameter(ip,'saveMatFile', true,                     @islogical);
addParameter(ip,'noConFlag'  , 'interp',                 @ischar);
parse(ip,AFIL,varargin{:});

%% Initialize the "pol" struct ----------------------
pol = struct;

pol.alpha = ip.Results.aSweep;
pol.Re    = ip.Results.rSweep;
pol.Mach  = ip.Results.Mach;
pol.nCrit = ip.Results.nCrit;
pol.xTrip = [ip.Results.xTripTop ip.Results.xTripBot];

pol.dim = [length(pol.alpha) length(pol.Re) length(pol.Mach)];
pol.d1  = 'alpha';
pol.d2  = 'Re';
pol.d3  = 'Mach';
pol.CD  = zeros(pol.dim); 
pol.CM  = pol.CD; 
pol.CL  = pol.CD;
pol.topXtr   = pol.CD;
pol.botXtr   = pol.CD;
pol.conFlag  = pol.CD;
pol.foilName = ip.Results.AFIL;
pol.pFile    = cell(pol.dim(2),pol.dim(3));

%% generate a filename for this polar
[~, foilName, ~] = fileparts(ip.Results.AFIL);
pol.name = [foilName '_topXtr' num2str(pol.xTrip(1)) '_botXtr' num2str(pol.xTrip(2)) ...
    '_nCrit' num2str(pol.nCrit)];

%% Set up the XFOIL Run File ----------------------
XFinfname = [ip.Results.xfoilDir filesep pol.name '.XFin.txt']; %xfoil run script
XFin = fopen(XFinfname,'w');
%make sure the polar output directory exists
if exist(ip.Results.polDir,'dir') ~= 7
    mkdir(ip.Results.polDir);
end

%% XFOIL setup and viscous parameters ----------------------
    fprintf(XFin,'plop\nG\n\n'); %graphics off
    if numel(ip.Results.AFIL)==4 && all(isstrprop(ip.Results.AFIL, 'digit')) %detect NACA naming convention
        fprintf(XFin,'naca\n%s\n',ip.Results.AFIL);
    elseif numel(ip.Results.AFIL)==5 && all(isstrprop(ip.Results.AFIL, 'digit'))
        fprintf(XFin,'naca\n%s\n',ip.Results.AFIL);
    else        
        fprintf(XFin,'%s\n',['load ' ip.Results.AFIL]);
    end
    
    %Improve paneling
    fprintf(XFin, 'gdes\n');
    fprintf(XFin, 'unit\n'); % make sure airfoil is unit length in chord
    fprintf(XFin, 'cadd\n'); % add corner points exceeding angle theshhold
    fprintf(XFin, '5\n');    % 5 corner point to add
    fprintf(XFin, '1\n');    % angle threshhold of 1 degree
    fprintf(XFin, '-0.1 1.1\n\n');
    
    fprintf(XFin, 'operv\n');
    fprintf(XFin, 'vpar\n');
    fprintf(XFin, 'n\n%2.0d\n',ip.Results.nCrit);
    fprintf(XFin, 'xtr\n%2.2f\n%2.2f\n\n', ip.Results.xTripTop, ip.Results.xTripBot);

for mc = 1:pol.dim(3) %Mach count
    for rc = 1:pol.dim(2) %Reynolds count
        %% polar file naming and checking ----------------------        
        % make up unique xFoil polar file name that uses less than 24 characters
        % use up to 16 characters of the airfoil name, then random number for uniqueness
        randNumString = num2str(round(rand(24,1)*9)'); randNumString(randNumString == ' ') = [];
        runName = randNumString;
        foilNameLen = length(foilName);
        if foilNameLen > 16; foilName = foilName(1:16); end
        runName(1:foilNameLen+1) = [foilName '_']; %overwrite up to 16 characters of the random string
        pol.pFile{rc,mc} = [ip.Results.polDir filesep runName '.pol'];
        
        %% AOA sweep parameters ----------------------
        if mc == 1 && rc == 1
            fprintf(XFin, ['visc ' num2str(pol.Re(rc)) '\n']);
        else
            fprintf(XFin, ['re ' num2str(pol.Re(rc)) '\n']);
        end
        fprintf(XFin,'Mach\n%2.2f\n',pol.Mach(mc));
        fprintf(XFin, 'pacc\n');
        fprintf(XFin,'%s\n\n',pol.pFile{rc,mc});
        
        %find the center of the alpha sweep - start sweeps here for better
        %convergence
        [~,mai]=min(abs(pol.alpha)); %mid-alpha index
        for ac = mai:pol.dim(1) %the positive alphas
            fprintf(XFin,['ALFA ' num2str(pol.alpha(ac)) '\n']);
        end
        fprintf(XFin, 'init\n');
        for ac = mai-1:-1:1 %the negative alphas
            fprintf(XFin,['ALFA ' num2str(pol.alpha(ac)) '\n']);
        end
        fprintf(XFin, 'pacc\n');        
        fprintf(XFin, 'init\n');
    end % Reynolds Number
    fprintf(XFin, 'pdel 0\n');
end %Mach Number
fprintf(XFin, '\n\n\nquit\n');
fclose(XFin);

%% Run XFOIL ----------------------
if ispc
    executable = 'xfoil';  
else
    executable = 'Xfoil.app/Contents/Resources/./xfoil';
end
disp(['launching airfoil: ' foilName ' in XFOIL for ' num2str(pol.dim(2)*pol.dim(3),'%3.0f') ' Mach/Re runs, each consisting of ' ...
    num2str(pol.dim(1),'%3.0f') ' angles of attack']); tic
%system command
[~,~] = system([ip.Results.xfoilDir filesep executable ' < ' '"' XFinfname '"']); toc
%delete batch file
if ~ip.Results.keepBatch 
    try delete(XFinfname);catch;end
end

%% Parse .Pol Files ----------------------
for mc = 1:pol.dim(3)
    for rc = 1:pol.dim(2)
        %import the xfoil polar
        clear p
        p = importPolars(pol.pFile{rc,mc}); %polar for this reynolds number, mach
        p.conFlag = p.alpha*0 + 1; %flag all of these points as converged (conFlag == 1)
        %check for unconverged points, not written to output file -
        %detected by length of output polar points vs. # of requested
        %alphas
        if length(p.alpha) < length(pol.alpha)                
                aint = setxor(round(p.alpha,2),round(pol.alpha,2)); %alphas to interpolate
                unConv         = table; %initialize table to store unconverged points
                unConv.alpha   = aint;
                if strcmpi(ip.Results.noConFlag,'interp') && ~isempty(p.alpha)
                    unConv.CL      = interp1(p.alpha,p.CL,aint,'linear','extrap');
                    unConv.CD      = interp1(p.alpha,p.CD,aint,'linear','extrap');
                    unConv.CDp     = interp1(p.alpha,p.CDp,aint,'linear','extrap');
                    unConv.CM      = interp1(p.alpha,p.CM,aint,'linear','extrap');
                    unConv.topXtr  = interp1(p.alpha,p.topXtr,aint,'linear','extrap');
                    unConv.botXtr  = interp1(p.alpha,p.botXtr,aint,'linear','extrap');
                else
                    unConv.CL      = aint*nan;
                    unConv.CD      = aint*nan;
                    unConv.CDp     = aint*nan;
                    unConv.CM      = aint*nan;
                    unConv.topXtr  = aint*nan;
                    unConv.botXtr  = aint*nan;
                end
            unConv.conFlag = aint*0; %mark these points as unconverged
            %mix the interpolated values into the alpha sweep
            if length(aint) + length(p.alpha) ~= length(pol.alpha); keyboard; end
            p = [p; unConv]; 
            p = sortrows(p,1); %re-sort with new interpolated values that did not converge
        elseif length(p.alpha) > length(pol.alpha)
            disp(['polar file ' pol.pFile{rc,mc} ' has more points than called for']);keyboard;
        end %alpha points needing interpolating
        pol.CL(:,rc,mc)      = p.CL;
        pol.CD(:,rc,mc)      = p.CD;
        pol.CM(:,rc,mc)      = p.CM;
        pol.topXtr(:,rc,mc)  = p.topXtr;
        pol.botXtr(:,rc,mc)  = p.botXtr;
        pol.conFlag(:,rc,mc) = p.conFlag;
    end
end

% warn user that some points did not converge
numNonConv = nnz(~pol.conFlag);
warning([num2str(numNonConv) ' XFOIL alphas did not converge. see "Pol.conFlag" = 0 for unconverged points']);

%% save polar structure to .mat file ----------------------
if ip.Results.saveMatFile
    polMatName = [ip.Results.polDir filesep pol.name '.mat'];
    save(polMatName,'pol');
    disp(['saved polars to: ' polMatName]);
else
    polMatName = [];
end

%% clean up polars directory ----------------------
if ~ip.Results.keepPolars
    for ii = 1:numel(pol.pFile)
        delete(pol.pFile{ii})
    end
end

