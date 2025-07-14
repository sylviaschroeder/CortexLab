function db = db_columns_vane

% Experiments used in addition to data for Neuron paper:
% M150305_SS041, 2015-04-11, 1-3 (Vane pre-processed)
% M150305_SS041, 2015-04-23, 2
% M150410_SS044, 2015-04-28, 2
% M150410_SS044, 2015-04-30, 2-4 (Vane pre-processed)
% M150410_SS044, 2015-05-15, 3-5 + 7 (Vane pre-processed)
% M150410_SS045, 2015-05-04, 1
% M150410_SS045, 2015-05-05, 1
% M150610_SS047, 2015-12-03, 2-3
% M150611_SS048, 2015-11-09, 3-4
% 
% Experiments not used here but in the Neuron paper:
% any experiments in darkness or gray screens
% M150121_SS038, 2015-02-17, 1+2

k = 0;

%% Retinal Boutons
% large Delta-F-over-F values have been corrected already for these data for Neuron paper
k = k + 1;
db(k).subject = 'M160706_SS066';
db(k).date = '2016-09-01';
db(k).expGratings = 2;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 7:14;
db(k).newSuite2p = false;
db(k).depth = 55;

k=k+1; %2
db(k).subject = 'M160923_SS069';
db(k).date = '2016-10-11';
db(k).expGratings = 3;
db(k).expNoise = [];
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 7:14;
db(k).newSuite2p = false;
db(k).depth = 29;

k=k+1; %3
db(k).subject = 'M160923_SS069';
db(k).date = '2016-10-13';
db(k).expGratings = 1;
db(k).expNoise = 4;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = [9 10 12 14 15];
db(k).newSuite2p = false;
db(k).depth = 28;

k=k+1; %4
db(k).subject = 'M160923_SS069';
db(k).date = '2016-10-21';
db(k).expGratings = 1;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 39; % above plane 4

k=k+1; %5
db(k).subject = 'M160923_SS070';
db(k).date = '2016-10-18';
db(k).expGratings = 1;
db(k).expNoise = 4;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 7:12;
db(k).newSuite2p = false;
db(k).depth = 14; % above plane 6

k=k+1; %6
db(k).subject = 'M160923_SS071';
db(k).date = '2016-10-18';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 10:14;
db(k).newSuite2p = false;
db(k).depth = 32; % above plane 6

k=k+1; %7
db(k).subject = 'M170821_SS075';
db(k).date = '2017-09-13';
db(k).expGratings = 3;
db(k).expNoise = 2;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 5;
db(k).newSuite2p = false;
db(k).depth = 30; % above plane 4

k=k+1; %8
db(k).subject = 'SS075';
db(k).date = '2017-09-27';
db(k).expGratings = 3;
db(k).expNoise = 2;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 5;
db(k).newSuite2p = false;
db(k).depth = 24.5; % above plane 4

k=k+1; %9
db(k).subject = 'M170821_SS076';
db(k).date = '2017-09-12';
db(k).expGratings = [];
db(k).expNoise = 2;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 7;
db(k).newSuite2p = false;
db(k).depth = 20; % above plane 4

k=k+1; %10
db(k).subject = 'M170821_SS076';
db(k).date = '2017-09-27';
db(k).expGratings = 2;
db(k).expNoise = [];
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 50; % above plane 4

k=k+1; %11
db(k).subject = 'SS076';
db(k).date = '2017-10-02';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 38;

k=k+1; %12
db(k).subject = 'SS076';
db(k).date = '2017-10-04';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 85; % above plane 4

k=k+1; %13
db(k).subject = 'SS077';
db(k).date = '2017-09-28';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 37; % above plane 4

k=k+1; %14
db(k).subject = 'SS077';
db(k).date = '2017-10-03';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 30; % above plane 4

k=k+1; %15
db(k).subject = 'SS077';
db(k).date = '2017-10-05';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 6;
db(k).newSuite2p = false;
db(k).depth = 50; % above plane 5

k=k+1; %16
db(k).subject = 'SS078';
db(k).date = '2017-09-28';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 5;
db(k).newSuite2p = false;
db(k).depth = 44; % above plane 4

k=k+1; %17
db(k).subject = 'SS078';
db(k).date = '2017-10-04';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 5;
db(k).newSuite2p = false;
db(k).depth = 49; % above plane 4

k=k+1; %18
db(k).subject = 'SS078';
db(k).date = '2017-10-05';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 5;
db(k).newSuite2p = false;
db(k).depth = 30; % above plane 4

%% SC Neurons

k=k+1;
db(k).subject = 'M150305_SS041'; % Vane
db(k).date = '2015-04-11';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = 1;
db(k).planes = 1:4; % set to be 8 um apart
db(k).newSuite2p = true;
db(k).depth = 0; % top-most neurons; dirst seen 103 microns above recording depth

k=k+1;
db(k).subject = 'M150305_SS041';
db(k).date = '2015-04-23';
db(k).expGratings = 3;
db(k).expNoise = 4;
db(k).expStatic = 2;
db(k).expBars = [];
db(k).planes = 2:5;
db(k).newSuite2p = false;
db(k).depth = 10.8;

k=k+1;
db(k).subject = 'M150410_SS044';
db(k).date = '2015-04-28';
db(k).expGratings = 3;
db(k).expNoise = 4;
db(k).expStatic = [];
db(k).expBars = 2;
db(k).planes = 2:5;
db(k).newSuite2p = false;
db(k).depth = 5.8;

k=k+1;
db(k).subject = 'M150410_SS044'; % Vane
db(k).date = '2015-04-30';
db(k).expGratings = 3;
db(k).expNoise = 4;
db(k).expStatic = 2;
db(k).expBars = [];
db(k).planes = 1:3; % set to be 10 um apart
db(k).newSuite2p = true;
db(k).depth = 8.2;

k=k+1;
db(k).subject = 'M150410_SS044'; % Vane
db(k).date = '2015-05-15';
db(k).expGratings = 4;
db(k).expNoise = 7;
db(k).expStatic = 5;
db(k).expBars = 3;
db(k).planes = 1:3; % set to be 10 um apart
db(k).newSuite2p = true;
db(k).depth = 6.1;

k=k+1;
db(k).subject = 'M150410_SS044'; % not on Vane's list; dataset added now (08/2021)
db(k).date = '2015-05-29';
db(k).expGratings = 4;
db(k).expNoise = 1;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 2:4;
db(k).newSuite2p = true;
db(k).depth = 3.6; % depth information per ROI looks fine

k=k+1;
db(k).subject = 'M150410_SS045';
db(k).date = '2015-05-04';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = [];
db(k).expBars = 1;
db(k).planes = 2:4;
db(k).newSuite2p = false;
db(k).depth = 13.3;

k=k+1;
db(k).subject = 'M150410_SS045';
db(k).date = '2015-05-05';
db(k).expGratings = 2;
db(k).expNoise = 3;
db(k).expStatic = 1;
db(k).expBars = [];
db(k).planes = 2:4;
db(k).newSuite2p = false;
db(k).depth = 7.5;

k=k+1;
db(k).subject = 'M150610_SS047'; % re-processed by Vane, but using my old data
db(k).date = '2015-11-23';
db(k).expGratings = 1;
db(k).expNoise = 2;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 2:4;
db(k).newSuite2p = false;
db(k).depth = 17;

k=k+1;
db(k).subject = 'M150610_SS047';
db(k).date = '2015-12-03';
db(k).expGratings = 1;
db(k).expNoise = 4;
db(k).expStatic = 2;
db(k).expBars = 3;
db(k).planes = 2:4;
db(k).newSuite2p = false;
db(k).depth = 28;

k=k+1;
db(k).subject = 'M150611_SS048';
db(k).date = '2015-11-09';
db(k).expGratings = 1;
db(k).expNoise = 5;
db(k).expStatic = 3;
db(k).expBars = 4;
db(k).planes = 2:4;
db(k).newSuite2p = false;
db(k).depth = 30;

k=k+1;
db(k).subject = 'M150611_SS048'; % re-processed by Vane, but using my old data
db(k).date = '2015-12-02';
db(k).expGratings = 1;
db(k).expNoise = 4;
db(k).expStatic = [];
db(k).expBars = [];
db(k).planes = 2:4;
db(k).newSuite2p = false;
db(k).depth = 42;

%% Not used
% k=k+1;
% db(k).subject = 'M150114_SS035'; % Vane; registered movie bad quality
% db(k).date = '2015-02-10';
% db(k).expGratings = 2;
% db(k).expGrayScreen = [];
% db(k).expNoise = 1;
% db(k).expDark = [];
% db(k).expStatic = [];
% db(k).expBars = [];
% db(k).planes = 2:4;

% k=k+1;
% db(k).subject = 'M150305_SS041'; % Vane; registered movie bad quality
% db(k).date = '2015-03-30';
% db(k).expGratings = 2;
% db(k).expGrayScreen = [];
% db(k).expNoise = 1;
% db(k).expDark = [];
% db(k).expStatic = [];
% db(k).expBars = [];
% db(k).planes = 2:3;

% k=k+1;
% db(k).subject = 'M150305_SS041'; % Vane; registered movie bad quality
% db(k).date = '2015-04-15';
% db(k).expGratings = 2;
% db(k).expGrayScreen = 5;
% db(k).expNoise = 3;
% db(k).expDark = [];
% db(k).expStatic = [];
% db(k).expBars = 1;
% db(k).planes = 2:5;

% k=k+1;
% db(k).subject = 'M150323_SS042'; % Vane; registered movie bad quality
% db(k).date = '2015-04-14';
% db(k).expGratings = 2;
% db(k).expGrayScreen = 4;
% db(k).expNoise = 3;
% db(k).expDark = [];
% db(k).expStatic = [];
% db(k).expBars = 1;
% db(k).planes = 2:5;

% k=k+1;
% db(k).subject = 'M150323_SS042'; % Vane; registered movie bad quality
% db(k).date = '2015-04-22';
% db(k).expGratings = 2;
% db(k).expGrayScreen = 4;
% db(k).expNoise = 3;
% db(k).expDark = [];
% db(k).expStatic = [];
% db(k).expBars = 1;
% db(k).planes = 2:5;

% k=k+1;
% db(k).subject = 'M150410_SS045'; % Vane; registered movie bad quality
% db(k).date = '2015-05-28';
% db(k).expGratings = [];
% db(k).expGrayScreen = [];
% db(k).expNoise = 6;
% db(k).expDark = [];
% db(k).expStatic = [];
% db(k).expBars = [];
% db(k).planes = 2:4;