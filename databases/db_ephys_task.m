function db = db_ephys_task

k = 0;

k=k+1;
db(k).subject = 'SS087';
db(k).date = '2017-12-12';
db(k).probes = {'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','LA','RP'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS087';
db(k).date = '2017-12-13';
db(k).probes = {'K1', 'K2', 'K3'};
db(k).probeLocations = {'RP','LP','LA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS087';
db(k).date = '2017-12-14';
db(k).probes = {'K1', 'K2', 'K3'};
db(k).probeLocations = {'RP','LP','LA'};
db(k).expTL = 3;
db(k).expTask = 4;
db(k).expNoise = 5;
db(k).expPassive = 6;

k=k+1;
db(k).subject = 'SS087';
db(k).date = '2017-12-15';
db(k).probes = {'K1', 'K2', 'K3'}; % error for K2, K2 and K3 not sorted
db(k).probeLocations = {'RP','LP','LA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS088';
db(k).date = '2018-01-30';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 5;
db(k).expNoise = 6;

k=k+1;
db(k).subject = 'SS088';
db(k).date = '2018-01-31';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS088';
db(k).date = '2018-02-01';
db(k).probes = {'K1', 'K2', 'ZO'};
db(k).probeLocations = {'LP','RP','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS088';
db(k).date = '2018-02-02';
db(k).probes = {'K1', 'K2', 'ZO'};
db(k).probeLocations = {'LP','RP','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS088';
db(k).date = '2018-02-03';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS088';
db(k).date = '2018-02-04';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

% With SC inactivation ----------------------------------------------------
% k=k+1;
% db(k).subject = 'SS088';
% db(k).date = '2018-05-10';
% db(k).probes = {'K1', 'K2'};
% db(k).probeLocations = {'L','R'};
% db(k).expTL = 2;
% db(k).expTaskInact = 4;
% db(k).expNoise = 5;
% 
% k=k+1;
% db(k).subject = 'SS088';
% db(k).date = '2018-05-11';
% db(k).probes = {'K1', 'K2'};
% db(k).probeLocations = {'L','R'};
% db(k).expTL = 3;
% db(k).expTaskInact = 4;
% db(k).expNoise = 5;
% db(k).expPassiveInact = 6;
% -------------------------------------------------------------------------

k=k+1;
db(k).subject = 'SS089';
db(k).date = '2018-02-06';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'}; % ZO not sorted
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 4;
db(k).expNoise = 5;
db(k).expPassive = 6;

k=k+1;
db(k).subject = 'SS089';
db(k).date = '2018-02-07';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

k=k+1;
db(k).subject = 'SS089';
db(k).date = '2018-02-08';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 6;

k=k+1;
db(k).subject = 'SS089';
db(k).date = '2018-02-09';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 5;
db(k).expTask = 6;
db(k).expNoise = 7;
db(k).expPassive = 8;

k=k+1;
db(k).subject = 'SS089';
db(k).date = '2018-02-10';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = [3 4];
db(k).expNoise = 5;
db(k).expPassive = 6;

k=k+1;
db(k).subject = 'SS089';
db(k).date = '2018-02-11';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 5;

% SS090: inactivation of SC during visual stimulation only (no task)

% k=k+1;
% db(k).subject = 'SS091';
% db(k).date = '2018-04-25';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 2;
% db(k).expTask = [3 4];
% db(k).expNoise = 5;
% db(k).expPassive = 6;
% 
% k=k+1;
% db(k).subject = 'SS091';
% db(k).date = '2018-04-26';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 2;
% db(k).expTask = [3 5 6];
% db(k).expNoise = 4;
% 
% k=k+1; % signals was used
% db(k).subject = 'SS091';
% db(k).date = '2018-04-30';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 5;

k=k+1;
db(k).subject = 'SS092';
db(k).date = '2018-03-20';
db(k).probes = {'K1', 'K2', 'ZO'};
db(k).probeLocations = {'LP','RP','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 4;
db(k).expPassive = 6;

k=k+1;
db(k).subject = 'SS092';
db(k).date = '2018-03-21';
db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
db(k).probeLocations = {'LP','RP','LA','RA'};
db(k).expTL = 2;
db(k).expTask = 3;
db(k).expNoise = 5;
db(k).expPassive = 6;

% k=k+1;
% db(k).subject = 'SS092';
% db(k).date = '2018-04-17';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 3;
% db(k).expTask = 4;
% db(k).expNoise = 5;
% db(k).expPassive = 6;
% 
% k=k+1;
% db(k).subject = 'SS092';
% db(k).date = '2018-04-18';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 4;
% db(k).expPassive = 5;
% 
% k=k+1;
% db(k).subject = 'SS092';
% db(k).date = '2018-04-19';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 4;
% db(k).expPassive = 5;
% 
% k=k+1;
% db(k).subject = 'SS092';
% db(k).date = '2018-04-20';
% db(k).probes = {'K1', 'K2', 'K3', 'ZO'};
% db(k).probeLocations = {'LP','RP','LA','RA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 4;
% db(k).expPassive = 5;
% 
% % from here on: signals was used
% k=k+1;
% db(k).subject = 'SS093';
% db(k).date = '2018-05-24';
% db(k).probes = {'K1', 'K2', 'K3'};
% db(k).probeLocations = {'LP','RP','LA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 4;
% db(k).expNoiseGray = 5;
% db(k).expPassive = 6;
% 
% k=k+1;
% db(k).subject = 'SS093';
% db(k).date = '2018-05-25';
% db(k).probes = {'K1', 'K2', 'K3'};
% db(k).probeLocations = {'LP','RP','LA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 4;
% db(k).expNoiseGray = 5;
% db(k).expPassive = 6;
% 
% k=k+1;
% db(k).subject = 'SS093';
% db(k).date = '2018-05-26';
% db(k).probes = {'K1', 'K2', 'K3'};
% db(k).probeLocations = {'LP','RP','LA'};
% db(k).expTL = 2;
% db(k).expTask = 3;
% db(k).expNoise = 4;
% db(k).expNoiseGray = 5;
% db(k).expPassive = 6; % not mentioned in lab book but data look fine
% 
% k=k+1;
% db(k).subject = 'SS093';
% db(k).date = '2018-05-27';
% db(k).probes = {'K1', 'K2', 'K3'};
% db(k).probeLocations = {'LP','RP','LA'};
% db(k).expTL = 3;
% db(k).expTask = 9;
% db(k).expNoise = 10;
% db(k).expNoiseGray = 11;
% db(k).expNoiseLeftOnly = 13;
% db(k).expNoiseRightOnly = 14;
% db(k).expNoiseLeftCovered = 15;
% db(k).expNoiseRightCovered = 16;
% 
% k=k+1;
% db(k).subject = 'SS093';
% db(k).date = '2018-05-28';
% db(k).probes = {'K1', 'K2', 'K3'};
% db(k).probeLocations = {'LP','RP','LA'};
% db(k).expTL = 2;
% db(k).expTask = 4;
% db(k).expNoise = 5;
% db(k).expNoiseGray = 6;
% db(k).expNoiseLeftOnly = 7;
% db(k).expNoiseRightOnly = 10;
% db(k).expNoiseLeftCovered = 8;
% db(k).expNoiseRightCovered = 9;