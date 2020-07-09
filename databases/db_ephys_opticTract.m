k = 0;
% 
k = k + 1; % 1 (2)
db(k).subject = 'SS084';
db(k).date = '2017-12-03';
db(k).probeNames = {'K1'};
db(k).hemispheres = {'left'};
db(k).expOri = 2;
db(k).expNoise = 3;
db(k).darkTime = [50 2100];
db(k).OTdepth = {[0 800]};
db(k).OTprobe = 1;
db(k).OTunits = {[21 22]};
db(k).OTgood = {[1 1]};
db(k).OTtimes = {{[],[]}};
db(k).OTampSTDs = {[20 0]};

% ------------------------------------------------------------------------

k = k + 1; % 2 (1)
db(k).subject = 'SS096';
db(k).date = '2018-03-07';
db(k).probeNames = {'K1', 'ZO'};
db(k).hemispheres = {'left','right'};
db(k).expFlicker = 3;
db(k).expOri = 5;
db(k).expOriTropic = 6;
db(k).expNoise = 4;
db(k).darkTime = [50 2480];
db(k).darkTime_IR_off = [4266 4591];
db(k).OTprobe = 2;
db(k).OTunits = {[20 48 49 52]};
db(k).OTgood = {[0 0 0 1]};
db(k).OTtimes = {{[],[],[],[1000 5700]}};
db(k).OTdepth = {[0 1000]};

k = k + 1; % 3 (9)
db(k).subject = 'SS096';
db(k).date = '2018-03-08';
db(k).probeNames = {'K1', 'ZO'};
db(k).hemispheres = {'left','right'};
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [40 2620];
db(k).darkTime_IR_off = [4097 4519];
db(k).OTprobe = 2;
db(k).OTunits = {[19 14 17 16 21  22 23 24 25 30  31 32 34 35 36  70 37 44]};
db(k).OTgood =  {[1  1  0  0  1   0  0  0  1  1   0  0  0  1  0   1  1  1]};
db(k).OTdepth = {[0 800]};
db(k).OTtimes = {{[],[0 1450;1810 5500],[],[],[15 600;1810 4050], ...
    [0 500;1800 2240;2625 5500],[],[],[500 5500],[450 1810; 2621 5500], ...
    [],[],[0 600;770 4050;4780 5500],[],[], ...
    [],[500 5500],[]}};
db(k).OTampSTDs = {[0 20 0 0 10  0 0 0 0 10  0 0 0 0 0  0 0 0]};

k = k + 1; % 4 (3)
db(k).subject = 'SS096';
db(k).date = '2018-03-09';
db(k).probeNames = {'K1', 'ZO'};
db(k).hemispheres = {'left','right'};
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [30 2620];
db(k).darkTime_IR_off = [4031 4450];
db(k).OTprobe = 2;
db(k).OTunits = {[0 6 2]};
db(k).OTgood = {[1  1  1]};
db(k).OTtimes = {{[790 5500],[240 700;2280 2520;2622 5500], ...
    [190 1300;2622 5500]}};
db(k).OTampSTDs = {[0 10 5]};
db(k).OTdepth = {[0 800]};

k = k + 1; % 5 (9) unit 66 was excluded because there is no running
db(k).subject = 'SS098';
db(k).date = '2018-03-16';
db(k).probeNames = {'K1', 'ZO'};
db(k).hemispheres = {'left','right'};
db(k).expFlicker = 4;
db(k).expOri = 6;
db(k).expOriTropic = 7;
db(k).expNoise = 5;
db(k).darkTime = [80 2540];
db(k).darkTime_IR_off = [4262 4597];
db(k).OTprobe = 1;
db(k).OTunits = {[52 51 57 59 65  66 67 69 70 71  73 72 74]};
db(k).OTgood = {[ 1  1  1  1  1   0  1  0  0  0   1  1  1]};
db(k).OTtimes = {{[],[500 5500],[],[],[], ...
    [160 495;1845 2430;2600 4250;4600 5500],[], ...
        [160 495;1400 1596;1840 2431;2620 5500],[],[], ...
    [],[150 4250;4750 5500],[160 5500]}};
db(k).OTampSTDs = {[0 0 0 0 5  5 0 0 0 0  0 0 7]};
db(k).OTdepth = {[0 1200]};

k = k + 1; % 6 (1)
db(k).subject = 'SS099';
db(k).date = '2018-03-16';
db(k).probeNames = {'K1', 'ZO'};
db(k).hemispheres = {'left','right'};
db(k).expFlicker = 2;
db(k).expOri = 4;
db(k).expOriTropic = 5;
db(k).expNoise = 3;
db(k).darkTime = [30 2420];
db(k).darkTime_IR_off = [3864 4199];
db(k).OTprobe = 1;
db(k).OTunits = {[1 3]};
db(k).OTgood =  {[0 1]};
db(k).OTtimes = {{[],[1450 5000]}};
db(k).OTampSTDs = {[0 7.5]};
db(k).OTdepth = {[0 1000]};