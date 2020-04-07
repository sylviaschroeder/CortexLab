k = 0;

% Path should be set to ZSERVER\DATA\Subjects!

% k = k + 1;
% db(k).subject = 'SS079';
% db(k).date = '2017-09-01';
% db(k).expOri = 3;
% db(k).expNoise = 2;
% db(k).darkTime = [50 3100];
% db(k).OTdepth = [0 600; 1650 2000];

% k = k + 1;
% db(k).subject = 'SS079';
% db(k).date = '2017-09-02';
% db(k).expOri = 3;
% db(k).expNoise = 2;
% db(k).darkTime = [80 2000];
% db(k).OTdepth = [1300 1700; 2900 3200];
% db(k).OTunits = []; %[732 17];
% db(k).OTtimes = {[], []};
% 
% k = k + 1;
% db(k).subject = 'SS080';
% db(k).date = '2017-09-05';
% db(k).expOri = 3;
% db(k).expNoise = 2;
% db(k).darkTime = [50 2190;4050 5540];
% db(k).OTdepth = [150 700; 1500 1800];
% 
% k = k + 1;
% db(k).subject = 'SS080';
% db(k).date = '2017-09-06';
% db(k).expOri = 3;
% db(k).expNoise = 2;
% db(k).darkTime = [50 3060;4850 5840];
% db(k).OTdepth = [600 900; 1750 2000];
 
k = k + 1; % 1
db(k).subject = 'SS081';
db(k).date = '2017-09-14';
db(k).expOri = 3;
db(k).expNoise = 2;
db(k).darkTime = [50 2340;3950 4700];
db(k).OTdepth = [0 500];

k = k + 1; % 2 
db(k).subject = 'SS081';
db(k).date = '2017-09-15';
db(k).expOri = 2;
db(k).expNoise = 3;
db(k).darkTime = [50 2500];
db(k).OTdepth = [0 800];

% From here on, path should be set to ZUBJECTS\Subjects!!!
% 
k = k + 1; % 3
db(k).subject = 'SS082';
db(k).date = '2017-11-23';
db(k).expOri = 2;
db(k).expNoise = 3;
db(k).darkTime = [30 1800];
db(k).OTdepth = [500 1400];
% db(k).OTunits = [2487 815 1012];
% db(k).OTtimes = {[], [], [1000 1800]};

k = k + 1; % 4
db(k).subject = 'SS082';
db(k).date = '2017-11-24';
db(k).expOri = 2;
db(k).expNoise = 3;
db(k).darkTime = [40 1900];
db(k).OTdepth = [700 1700];
% db(k).OTunits = [920];
% db(k).OTtimes = {[300 1700]};
% 
% k = k + 1;
% db(k).subject = 'SS082';
% db(k).date = '2017-11-25';
% db(k).expOri = 2;
% db(k).expNoise = 3;
% db(k).darkTime = [50 2320];
% db(k).OTdepth = [000 1300];
% 
% k = k + 1;
% db(k).subject = 'SS083';
% db(k).date = '2017-11-28';
% db(k).expOri = 2;
% db(k).expNoise = 3;
% db(k).darkTime = [40 2440];
% db(k).OTdepth = [300 1100];
% 
% k = k + 1;
% db(k).subject = 'SS083';
% db(k).date = '2017-11-29';
% db(k).expOri = 2;
% db(k).expNoise = 3;
%  
% k = k + 1;
% db(k).subject = 'SS083';
% db(k).date = '2017-11-30';
% db(k).expOri = 2;
% db(k).expNoise = 3;
% db(k).darkTime = [40 2040];
% db(k).OTdepth = [300 1100];
% 
% k = k + 1;
% db(k).subject = 'SS084';
% db(k).date = '2017-12-01';
% db(k).expOri = 2;
% db(k).expNoise = 3;
% db(k).darkTime = [30 2060];
% db(k).OTdepth = [0 700];
% 
% k = k + 1;
% db(k).subject = 'SS084';
% db(k).date = '2017-12-02';
% db(k).expOri = 2;
% db(k).expNoise = 3;
% 
k = k + 1; % 5
db(k).subject = 'SS084';
db(k).date = '2017-12-03';
db(k).expOri = 2;
db(k).expNoise = 3;
db(k).darkTime = [50 2100];
db(k).OTdepth = [0 800];
db(k).OTunits = [21 22];
db(k).OTgood = [1 1];
db(k).OTtimes = {[]};
db(k).cutAmplitudes = {[60 2100 3500]}; %each row: [amplitude threshold, start time, end time]

% ------------------------------------------------------------------------

k = k + 1; % 6
db(k).subject = 'SS096';
db(k).date = '2018-03-05';
db(k).expFlicker = 2;
db(k).expOri = 6;
db(k).expNoise = 5;
db(k).darkTime = [60 2500];
db(k).OTdepth = {[],[0 800]};

k = k + 1; % 7
db(k).subject = 'SS096';
db(k).date = '2018-03-06';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise =3;
db(k).darkTime = [50 2470];
% db(k).OTunits = {[],[10]};
% db(k).OTtimes = {{[]},{[]}};
db(k).OTdepth = {[],[0 1300]};

k = k + 1; % 8
db(k).subject = 'SS096';
db(k).date = '2018-03-07';
db(k).expFlicker = 3;
db(k).expOri = 5;
db(k).expOriTropic = 6;
db(k).expNoise = 4;
db(k).darkTime = [50 2480];
db(k).darkTime_IR_off = [4266 4591];
db(k).OTunits = {[],[20 48 49 52]};
db(k).OTgood = {[],[0 0 0 1]};
db(k).OTtimes = {{[]},{[],[],[],[1000 5700]}};
db(k).OTdepth = {[],[0 1000]};

k = k + 1; % 9
db(k).subject = 'SS096';
db(k).date = '2018-03-08';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [40 2620];
db(k).darkTime_IR_off = [4097 4519];
db(k).OTunits = {[],[19 14 17 16 21  22 23 24 25 30  31 32 34 35 36  70 37 44]};
db(k).OTgood = {[], [1  1  0  0  1   1  0  0  1  1   0  0  1  1  0   1  1  1]};
db(k).OTdepth = {[],[0 800]};
db(k).OTtimes = {{[]},{[],[],[],[],[15 600;1810 4050], ...
    [0 500;1800 2240;2625 5500],[],[],[500 5500],[450 1810; 2621 5500], ...
    [],[],[0 600;770 4050;4780 5500],[],[], ...
    [],[500 5500],[]}};
db(k).OTampSTDs = {[], [0 0 0 0 10  5 0 0 0 10  0 0 10 0 0  0 0 0]};

k = k + 1; % 10
db(k).subject = 'SS096';
db(k).date = '2018-03-09';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [30 2620];
db(k).darkTime_IR_off = [4031 4450];
db(k).OTunits = {[],[0 6 2]};
db(k).OTgood = {[], [1  1  1]};
db(k).OTtimes = {{[]},{[790 5500],[240 700;2280 2520;2622 5500], ...
    [190 1300;2622 5500]}};
db(k).OTampSTDs = {[], [0 10 5]};
db(k).OTdepth = {[],[0 800]};

% k = k + 1;
% db(k).subject = 'SS097';
% db(k).date = '2018-03-05';
% db(k).expFlicker = 2;
% db(k).expOri = 4;
% db(k).expNoise = 3;
% db(k).darkTime = [40 2380];
% db(k).OTdepth = {[],[]};

k = k + 1; % 11
db(k).subject = 'SS097';
db(k).date = '2018-03-06';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [40 2430];
db(k).OTdepth = {[],[0 1000]};

k = k + 1; % 12
db(k).subject = 'SS097';
db(k).date = '2018-03-07';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [40 3050];
db(k).OTdepth = {[0 1000],[0 1000]};

k = k + 1; % 13
db(k).subject = 'SS097';
db(k).date = '2018-03-08';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [30 2490];
db(k).OTdepth = {[],[0 800]};

k = k + 1; % 14
db(k).subject = 'SS097';
db(k).date = '2018-03-09';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [50 2380];
db(k).OTdepth = {[0 800],[0 800]};

k = k + 1; % 15
db(k).subject = 'SS098';
db(k).date = '2018-03-12';
db(k).expFlicker = 2;
db(k).expOri = 5;
db(k).expNoise = 4;
db(k).darkTime = [50 2440];
db(k).OTdepth = {[0 1000],[0 1000]};

k = k + 1; % 16
db(k).subject = 'SS098';
db(k).date = '2018-03-13';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [50 2440];
db(k).darkTime2 = [4300 4650];
db(k).OTdepth = {[0 800],[0 800]};

k = k + 1; % 17
db(k).subject = 'SS098';
db(k).date = '2018-03-14';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [40 2450];
db(k).darkTime2 = [3906 4236];
db(k).OTdepth = {[0 800],[0 800]};

k = k + 1; % 18
db(k).subject = 'SS098';
db(k).date = '2018-03-15';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [50 2600];
db(k).darkTime2 = [3972 4295];
% db(k).OTdepth = {[],[0 500]};
db(k).OTdepth = {[500 1300],[0 800]};

k = k + 1; % 19
db(k).subject = 'SS098';
db(k).date = '2018-03-16';
db(k).expFlicker = 4;
db(k).expOri = 6;
db(k).expOriTropic = 7;
db(k).expNoise = 5;
db(k).darkTime = [80 2540];
db(k).darkTime_IR_off = [4262 4597];
db(k).OTunits = {[52 51 57 59 65  66 67 69 70 71  73 72 74],[]};
db(k).OTgood = {[ 1  1  1  1  1   1  1  1  0  0   1  1  1], []};
db(k).OTtimes = {{[],[],[],[],[], ...
    [160 495;1845 2430;2600 4250;4600 5500],[], ...
        [160 495;1400 1840;1596 2431;2620 5500],[],[], ...
    [],[150 4250;4750 5500],[160 5500]}, ...
    {[]}};
db(k).OTampSTDs = {[0 0 0 0 5  5 0 0 0 0  0 0 7], []};
% db(k).OTdepth = {[],[0 750]};
db(k).OTdepth = {[0 1200],[0 800]};

k = k + 1; % 20
db(k).subject = 'SS099';
db(k).date = '2018-03-12';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [30 2430];
db(k).darkTime2 = [3916 4295];
db(k).OTdepth = {[0 1000],[0 800]};
% 
k = k + 1; % 21
db(k).subject = 'SS099';
db(k).date = '2018-03-13';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [40 2650];
db(k).darkTime_IR_off = [4054 4377];
db(k).OTunits = {[42 44 45 48],[]};
db(k).OTgood = {[ 1  1  1  1 ],[]};
db(k).OTtimes = {{[],[],[],[]},{[]}};
db(k).OTampSTDs = {[0 0 0 0],[]};
db(k).OTdepth = {[500 1300],[0 800]};

k = k + 1; % 22
db(k).subject = 'SS099';
db(k).date = '2018-03-14';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [40 2470];
db(k).darkTime2 = [4020 4542];
db(k).OTdepth = {[0 1000],[500 1300]};

k = k + 1; % 23
db(k).subject = 'SS099';
db(k).date = '2018-03-15';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expNoise = 3;
db(k).darkTime = [30 2440];
db(k).darkTime2 = [3815 4172];
db(k).OTdepth = {[0 800],[600 1300]};

k = k + 1; % 24
db(k).subject = 'SS099';
db(k).date = '2018-03-16';
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [30 2420];
db(k).darkTime_IR_off = [3864 4199];
db(k).OTunits = {[1 3],[]};
db(k).OTgood =  {[0 1],[]};
db(k).OTtimes = {{[],[]},{[]}};
db(k).OTampSTDs = {[0 7.5],[]};
db(k).OTdepth = {[0 1000],[0 1000]};