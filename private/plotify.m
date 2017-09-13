function h=plotify(varargin)
%PLOTIFY Create axes for subplots. 
%
% Inspired by:
% - subplot_tight by Nikolay S. (http://vision.technion.ac.il/~kolian1/).
% - tight_plot by William T. Adler

% Author:   Luigi Acerbi
% Email:    luigi.acerbi@gmail.com
% Date:     09/15/2016

% Return help
if nargin == 0 && nargout == 0
    help plotify;
    return;
end

% Options
options.gutter = [];
options.margins = [];
options.ticklength = [];
options.title = [];
options.labels = [];
options.square = [];
options.position = [];
options.fontsize = [];
options.fontname = [];

if isstruct(varargin{1})
    
elseif nargin == 1 || ischar(varargin{2})   % Input SUBGRID
    subgrid = varargin{1};
    [m,n] = size(subgrid);
    argidx0 = 2;
else                                        % Input M,N
    m = varargin{1};
    n = varargin{2};
    if ~isscalar(m) || ~isscalar(n) || ~(m == round(m)) || ~(n == round(n)) || m < 1 || n < 1
        error('plotify:WrongSizes', 'M and N need to be nonnegative integers.');
    end
    subgrid = reshape(1:(m*n),n,m)';
    argidx0 = 3;
end

% Number of subgraphs
ng = max(subgrid(:));

% Parse inputs
if numel(varargin) >= argidx0
    % varargin(argidx0:end)
    options = parseoptions(options,varargin{argidx0:end});
end

gutter = options.gutter;
margins = options.margins;
ticklength = options.ticklength;
titlestr = options.title;
labels = options.labels;
square = options.square;
position = options.position;
fontsize = options.fontsize;
fontname = options.fontname;

if isempty(gutter)  % Total gutter space
    gutter = [.05, .05]; %horizontal, vertical
%    gutter = [.002, .002]; %horizontal, vertical
end

if numel(gutter) == 1
    gutter(2)=gutter;
elseif numel(gutter) > 2
    error('plotify:GutterWrongSize','GUTTER must be of length 1 or 2.')
end

if isempty(margins)
    % margins = [.06 .01 .04 .04]; % L R B T
    margins = [.05 .02 .05 .05]; % L R B T
end

% Specify margins
if numel(margins) == 1
    margins = margins*ones(1,4);
elseif numel(margins) == 2
    margins = [margins(1),margins(1),margins(2),margins(2)];
end

if isempty(ticklength)
    ticklength = 0.005; % Tick length in total figure size
end

if isempty(square)
    square = false(1,ng);
else
    square = logical([square(:)',zeros(1,ng-numel(square))]);
end

if isempty(fontsize)
    fontsize = 14;
end

if isempty(fontname)
    fontname = 'Arial';
end

% Set default font name
set(0,'defaultUicontrolFontName',fontname);
set(0,'defaultUitableFontName',fontname);
set(0,'defaultAxesFontName',fontname);
set(0,'defaultTextFontName',fontname);
set(0,'defaultUipanelFontName',fontname);

if ~isnumeric(subgrid) || any(subgrid(:) < 0)
    error('plotify:SubGridWrong','The SUBGRID matrix must contain nonnegative integers.');
end

grect = zeros(ng,4);

% Check that all panels exist and that they are rectangles
for g = 1:ng
    if all(subgrid ~= g)
        error('plotify:MissingPanel', ['Panel #' num2str(g) ' missing in SUBGRID matrix.']);
    end
    % Find the upper left corner
    for ii = 1:m
        for jj = 1:n
            if subgrid(ii, jj) == g && grect(g, 1) == 0; grect(g, 1:2) = [ii jj]; end
            if subgrid(ii, jj) == g; grect(g, 3:4) = [ii jj]; end
        end
    end    
    if any((grect(g, 3:4) - grect(g, 1:2)) < 0)
        error('plotify:MisshapenPanels', ...
            'The panels in SUBGRID matrix have inconsistent shapes.');
    end
    if square(g) && ( (grect(g,3)-grect(g,1)) ~= (grect(g,4)-grect(g,2)) )
        warning(['Panel #' num2str(g) ' is required to be SQUARE but it is not square in SUBGRID matrix.']);
    end    
end

% Define lengths
Lmargin = margins(1);
Rmargin = margins(2);
Bmargin = margins(3);
Tmargin = margins(4);

% Gutter
% gutter = gutter ./ [n-1,m-1];
gutter(~isfinite(gutter)) = 0;

effheight = 1-Bmargin-Tmargin;              % available height
effwidth = 1-Lmargin-Rmargin;               % available width
sqheight = (effheight - (m-1)*gutter(2))/m; % height of each square
sqwidth = (effwidth - (n-1)*gutter(1))/n;   % width of each square

h0 = axes('Position', [0,0,1,1]);   % Whole figure
box off;
axis off;

for g = 1:ng
    r = grect(g,:);
    
    height = (1+r(3)-r(1))*(sqheight + gutter(2)) - gutter(2);    % panel height
    width = (1+r(4)-r(2))*(sqwidth + gutter(1)) - gutter(1);      % panel width
    bottom = (m-r(3))*(sqheight+gutter(2))+Bmargin;           % bottom pos
    left   = (r(2)-1)*(sqwidth +gutter(1))+Lmargin;            % left pos
    
    pos_vec = [left bottom width height];
    
    h(g)=axes('Position', pos_vec);
end

drawnow;

for g = 1:ng
    base_pos(g,:) = get(h(g),'Position');
end

for g = 1:ng
    if square(g); axes(h(g)); axis square; end
end

% Resize figure
if ~isempty(position)
    if ischar(position)
        if strncmpi(position,'screen',6)
            position = get(0,'ScreenSize');
            set(gcf,'Position',position);
        elseif strncmpi(position,'max',3)
            warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            try
                jFrame = get(handle(gcf),'JavaFrame');
                jFrame.setMaximized(true);
            catch
                set(gcf,'units','normalized','OuterPosition',[0 0 1 1]);
            end
        end
    else
        set(gcf,'Position',position);
    end
end

h(end+1) = h0;
axes(h0);

% Figure title
if ~isempty(titlestr)
    text(0.5,1-Tmargin/2,titlestr,'HorizontalAlignment','center','FontSize',fontsize,'FontWeight','bold');
end

drawnow;
pause(0.05);

% Set tick length to the same size
equalticklength(h,ticklength);

for g = 1:ng
    
    set(h(g),'TickDir','out');
    
    rect = get(h(g),'Position');
    height = rect(4);
    bottom = rect(2);
    left = rect(1);    
    
    if ~isempty(labels) && g <= numel(labels) && ~isempty(labels{g})
        text(left-gutter(1)*0.75,bottom+height,labels{g}, ...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',fontsize);
    end
end


set(gcf,'Color','w');

end

%--------------------------------------------------------------------------

function options = parseoptions(options,varargin)
%PARSEOPTIONS Parse options either as struct or variable arguments in name/value format.
%   OPTIONS = PARSEOPTIONS(OPTIONS,'PROPERTY1',VALUE1,'PROPERTY2',VALUE2,...)
%   sets the fields propertyX in default OPTIONS structure to valueX. Input 
%   field names are not case sensitive. 
%
%   OPTIONS = PARSEOPTIONS(OPTIONS,NEWOPTS) assigns fields in struct NEWOPTS
%   to OPTIONS. Matching fields are not case sensitive for OPTIONS.
%
%   OPTIONS = PARSEOPTIONS(OPTIONS,NEWOPTS,'PROPERTY1',VALUE1,'PROPERTY2',VALUE2,...) 
%   first assigns values from struct NEWOPTS, and then name/value pairs.
%

%   Author: Luigi Acerbi
%   Email:  luigi.acerbi@gmail.com
%   Date:   Sep/08/2016

if nargin < 1; help parseoptions; return; end

if isempty(options)
    error('parseOptions:emptyOptions','Default OPTIONS struct should be nonempty.');
end

if isempty(varargin)
    return;
end

deff = fields(options)';

if isstruct(varargin{1})                    % Input as NEWOPTS struct
    newopts = varargin{1};
    for f = fields(newopts)'
        idx = find(strcmpi(f{:}, deff),1);
        if isempty(idx)
            error('parseOptions:unknownProperty', ...
                ['Unknown property ''' f{:} ''' in NEWOPTS.']);
        else
            options.(deff{idx}) = newopts.(f{:});
        end
    end
    varargin(1) = [];
end

if ~isempty(varargin)
    if ischar(varargin{1})                      % Input in name/value format
        % check for correct number of inputs
        if mod(numel(varargin),2) == 1
            error('parseOptions:wrongInputFormat', ...
                'Name and value input arguments must come in pairs.');
        end

        % parse arguments
        for i = 1:2:numel(varargin)
            if ischar(varargin{i})
                idx = find(strcmpi(varargin{i}, deff),1);
                if isempty(idx)
                    error('parseOptions:unknownProperty', ...
                        ['Unknown property name ''' varargin{i} '''.']);
                else
                    options.(deff{idx}) = varargin{i+1};
                end            
            else
                error('parseOptions:wrongInputFormat', ...
                    'Name and value input arguments must come in pairs.');
            end
        end
    else
        error('parseOptions:wrongInputFormat', ...
                'Input should come as a NEWOPTS struct and/or with name/value pairs.');
    end
end

end

%--------------------------------------------------------------------------

function equalticklength(h,ticklength)
%EQUALTICKLENGTH Make all tick lengths equal in all figure panels.


if nargin < 2 || isempty(ticklength); ticklength = 0.0035; end

ng = numel(h);

for g = 1:ng
    rect = get(h(g),'Position');
    height = rect(4);
    width = rect(3);
    bottom = rect(2);
    left = rect(1);
    
    %if square(g)
    %    axislen = min([height,width]);        
    %    ticklen = ticklength./axislen^2;        
    %    axislen = height;
    %else
    
    hfig = get(h(g),'Parent');
    rectfig = get(hfig, 'Position');
    width_px = rect(3) * rectfig(3);
    height_px = rect(4) * rectfig(4);
        
    %set(0,'units','pixels');
    %set(h(g),'Units','normalized');
            
        axislen = max([height_px,width_px]);
        ticklen = ticklength/axislen*max(rectfig(3:4));
    %end
    
    set(h(g),'TickLength',ticklen*[1 2.5]);    
end

end