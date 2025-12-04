function db = db_orientationDataFromSchroederLab

k = 0;

%% contains sparse noise
% FG006 10-11 and 11-11, FG004 09-06 and 10-06 and FG009 31-05

k = k + 1;
db(k).subject = 'FG006';
db(k).date = '2023-11-10';
db(k).gratingExperiments = "all";

k = k + 1;
db(k).subject = 'FG006';
db(k).date = '2023-11-11';
db(k).gratingExperiments = "all";

k = k + 1;
db(k).subject = 'FG004';
db(k).date = '2023-06-09';
db(k).gratingExperiments = "all";

k = k + 1;
db(k).subject = 'FG004';
db(k).date = '2023-06-10';
db(k).gratingExperiments = "all";

k = k + 1;
db(k).subject = 'FG009';
db(k).date = '2024-05-31';
db(k).gratingExperiments = 1; % of 2

% k = k + 1;
% db(k).subject = 'FG004';
% db(k).date = '2023-06-11';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG005';
% db(k).date = '2023-08-15';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG005';
% db(k).date = '2023-08-18';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG005';
% db(k).date = '2023-08-21';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG007';
% db(k).date = '2024-04-17';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG007';
% db(k).date = '2024-04-18';
% db(k).gratingExperiments = 1; % of 2
% 
% k = k + 1;
% db(k).subject = 'FG007';
% db(k).date = '2024-04-19';
% db(k).gratingExperiments = [1 2]; % of 3
% 
% k = k + 1;
% db(k).subject = 'FG007';
% db(k).date = '2024-04-20';
% db(k).gratingExperiments = [1 2]; % of 3
% 
% k = k + 1;
% db(k).subject = 'FG008';
% db(k).date = '2024-05-21';
% db(k).gratingExperiments = [1 2]; % of 4
% 
% k = k + 1;
% db(k).subject = 'FG008';
% db(k).date = '2024-05-22';
% db(k).gratingExperiments = [1 2]; % of 3
% 
% k = k + 1;
% db(k).subject = 'FG008';
% db(k).date = '2024-05-23';
% db(k).gratingExperiments = [1 2]; % of 4
% 
% k = k + 1;
% db(k).subject = 'FG009';
% db(k).date = '2024-06-01';
% db(k).gratingExperiments = [1 2]; % of 3 (after 2nd, there is 'gratingsStep'???)
% 
% k = k + 1;
% db(k).subject = 'FG009';
% db(k).date = '2024-06-02';
% db(k).gratingExperiments = [1 2]; % of 3 (before 1st, there is 'gratingsStep'???)
% 
% k = k + 1;
% db(k).subject = 'FG009';
% db(k).date = '2024-06-03';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG010';
% db(k).date = '2024-10-16';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG010';
% db(k).date = '2024-10-18';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG010';
% db(k).date = '2024-10-19';
% db(k).gratingExperiments = "all";
% 
% k = k + 1;
% db(k).subject = 'FG010';
% db(k).date = '2024-10-20';
% db(k).gratingExperiments = "all";