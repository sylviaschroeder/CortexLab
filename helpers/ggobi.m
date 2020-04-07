% [colors glyphs] = ggobi(data, colors, glyphs, point_names, var_names)
%
% initiates an ggobi process from the array data
% colors is an optional array of colors.
% The entries of colors should be numbers in the
% range 1 to 10, which are then translated to the
% color names used by xgobi.
%
% glyphs is an optional array of glyph numbers -
% to find out what number means what glyph look at
% the 'Glyph' menu on the 'Brush' screen of ggobi.
% "Points" is 31.
%
% Any Infs or NaNs will cause removal of that data line
%
% point_names gives labels for each points.  This should be a
% cell array or character array.
%
% to get labels for columns, either pass a structure array as data
% so you get the structure fields as labels, or use var_names
%
% optional output arguments give final brushing - also not yet
% implemented.

function ggobi(Input, colors, glyphs, labels, FieldNames)
% create xml file
docNode = com.mathworks.xml.XMLUtils.createDocument('ggobidata');
ggobidata = docNode.getDocumentElement;

data = docNode.createElement('data');

% is data a structure array?
if isstruct(Input) | istable(Input)
    % convert data to a numeric array
    FieldNames = fieldnames(Input);
    if istable(Input); FieldNames = setdiff(FieldNames, 'Properties'); end
    for i=1:length(FieldNames)
        d = getfield(Input, FieldNames{i});
        if isnumeric(d)
            d2(:,i) = d(:);
        elseif iscellstr(d)
            [c, ia, ic] = unique(d);
            if length(c)>1
                fprintf('Converting %s to numeric: ', FieldNames{i});
                d2(:,i) = ic;
                for j=1:length(c);
                    fprintf('%s -> %d, ', c{j}, ia(j));
                end
                fprintf('\n');
            else
                fprintf('Dropping %s since all entries the same\n', FieldNames{i});
            end
        else
            fprintf('Didn''t know what to do with %s\n', FieldNames{i});
        end
    end
    Input = d2;
elseif nargin<5 % i.e. FieldNames not provided
    % make up variable names
    FieldNames = cell(size(Input,2),1);
    for i=1:size(Input,2)
        FieldNames{i} = sprintf('v %d', i);
    end
end

% write variable names
variables = docNode.createElement('variables');
data.appendChild(variables);
variables.setAttribute('count', num2str(size(Input,2)));

for v=1:length(FieldNames);
    var = docNode.createElement('realvariable');
    var.setAttribute('name', FieldNames{v});
    variables.appendChild(var);
end

nPoints = size(Input,1);

% Look for missing data
% delete rows with NaNs or infs, but remember original row number
Infinities = ~isfinite(Input);
OriginalRow = 1:size(Input,1);
if (any(Infinities(:)))
    % find rows with infinities
    InfinityRows = find(any(Infinities,2));
    % delete those rows
    Input(InfinityRows, :) = [];
    OriginalRow(InfinityRows) = [];
    if exist('colors'), if length(colors)>1, colors(InfinityRows) = []; end; end;
    if exist('glyphs'), if length(glyphs)>1, glyphs(InifinityRows) = []; end; end;
    if exist('labels'), if length(labels)>1, labels(InifinityRows) = []; end; end;
    disp([num2str(length(InfinityRows)), ' rows containing infinities were deleted.']);
end;

nPoints = size(Input,1);
nColors = 9;
nGlyphs = 5;

% now write main data
records = docNode.createElement('records');
data.appendChild(records);
records.setAttribute('count', sprintf('%d', nPoints));
for r=1:nPoints
    rec = docNode.createElement('record');
    records.appendChild(rec);
    
    if exist('labels') &&~isempty(labels)
        rec.setAttribute('label', labels{r});
    end
    
    if exist('colors')&&~isempty(colors)
        rec.setAttribute('color', sprintf('%d', mod(colors(r), nColors)));
    end
    if exist('glyphs')&&~isempty(glyphs)
        rec.setAttribute('glyphType', sprintf('%d', mod(glyphs(r), nGlyphs)));
        rec.setAttribute('glyphSize', '1');
    end
    rec.appendChild(docNode.createTextNode(sprintf('%f ', Input(r,:))));
end

ggobidata.appendChild(data);
xmlwrite('gobidata.xml', docNode);
system('start ggobi gobidata.xml');