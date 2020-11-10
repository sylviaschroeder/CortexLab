i = 0;

i = i+1;
db(i).mouse_name    = 'M150121_SS038';
db(i).date          = '2015-02-17';
db(i).expts         = 1:2; % which experiments to cluster together
db(i).planesToProcess = 1;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[256 256 256 128 128],[512 512]};

i = i+1;
db(i).mouse_name    = 'M150305_SS041';
db(i).date          = '2015-04-23';
db(i).expts         = 1:4; % which experiments to cluster together
db(i).planesToProcess = 2:5;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[128 128 128 128],512};

i = i+1;
db(i).mouse_name    = 'M150410_SS044';
db(i).date          = '2015-04-28';
db(i).expts         = [2 3 4 5]; % which experiments to cluster together
db(i).planesToProcess = 2:5;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';

i = i+1;
db(i).mouse_name    = 'M150410_SS044';
db(i).date          = '2015-05-29';
db(i).expts         = 1:4; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[128 256 128],512};

i = i+1;
db(i).mouse_name    = 'M150410_SS045';
db(i).date          = '2015-05-04';
db(i).expts         = [1 2 3 4]; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[128 256 128],512};

i = i+1;
db(i).mouse_name    = 'M150410_SS045';
db(i).date          = '2015-05-05';
db(i).expts         = 1:4; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[170 170 172],[256 256]};

i = i+1;
db(i).mouse_name    = 'M150610_SS047';
db(i).date          = '2015-11-23';
db(i).expts         = 1:5; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[170 170 172],[256 256]};

i = i+1;
db(i).mouse_name    = 'M150610_SS047';
db(i).date          = '2015-12-03';
db(i).expts         = 1:5; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[170 170 172],512};

i = i+1;
db(i).mouse_name    = 'M150611_SS048';
db(i).date          = '2015-11-09';
db(i).expts         = 1:5; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[170 170 172],[256 256]};

i = i+1;
db(i).mouse_name    = 'M150611_SS048';
db(i).date          = '2015-12-02';
db(i).expts         = [1 2 3 4 6]; % which experiments to cluster together
db(i).planesToProcess = 2:4;
db(i).AlignToRedChannel = 1;
db(i).microID       = 'b';
db(i).splitBlocks   = {[170 170 172],[256 256]};

%% Don't consider those experiments for now

% i = i+1;
% db(i).mouse_name    = 'M150114_SS035';
% db(i).date          = '2015-02-10';
% db(i).expts         = 2; % which experiments to cluster together
% db(i).planesToProcess = 2:4;
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';
% db(i).splitBlocks   = {[128 128 128 128],[512]};

% i = i+1;
% db(i).mouse_name    = 'M150305_SS041';
% db(i).date          = '2015-03-30';
% db(i).expts         = 1:2; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';
% 
% i = i+1;
% db(i).mouse_name    = 'M150305_SS041';
% db(i).date          = '2015-04-11';
% db(i).expts         = 1:3; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';
% 
% i = i+1;
% db(i).mouse_name    = 'M150305_SS041';
% db(i).date          = '2015-04-15';
% db(i).expts         = 1:5; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';

% i = i+1;
% db(i).mouse_name    = 'M150323_SS042';
% db(i).date          = '2015-04-14';
% db(i).expts         = 1:4; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';

% i = i+1;
% db(i).mouse_name    = 'M150323_SS042';
% db(i).date          = '2015-04-22';
% db(i).expts         = 1:4; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';

% i = i+1;
% db(i).mouse_name    = 'M150410_SS044';
% db(i).date          = '2015-04-30';
% db(i).expts         = 1:4; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';

% i = i+1;
% db(i).mouse_name    = 'M150410_SS044';
% db(i).date          = '2015-05-15';
% db(i).expts         = [2:5 7]; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';

% i = i+1;
% db(i).mouse_name    = 'M150410_SS045';
% db(i).date          = '2015-05-28';
% db(i).expts         = [5 6]; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';

% i = i+1;
% db(i).mouse_name    = 'M150610_SS047';
% db(i).date          = '2015-11-04';
% db(i).expts         = 1:5; % which experiments to cluster together
% db(i).AlignToRedChannel = 1;
% db(i).microID       = 'b';