function db = db_saccadeDataFromSchroederLab

k = 0;

%% contains sparse noise

k = k + 1;
db(k).subject = 'Quille';
db(k).date = '2023-07-24';

k = k + 1; % mismatch between TTL pulses and number of circle stimuli
db(k).subject = 'Quille';
db(k).date = '2023-09-06';

k = k + 1;
db(k).subject = 'Quille';
db(k).date = '2023-09-28';

%contains only circle paradigm (no sparse noise)
k = k + 1;
db(k).subject = 'Giuseppina';
db(k).date = '2023-01-24';

k = k + 1;
db(k).subject = 'Tara';
db(k).date = '2023-12-07';

k = k + 1;
db(k).subject = 'Tara';
db(k).date = '2023-12-12';

k = k + 1;
db(k).subject = 'Tara';
db(k).date = '2023-12-19';

k = k + 1;  % mismatch between TTL pulses and number of circle stimuli
db(k).subject = 'Uma';
db(k).date = '2023-12-06';

k = k + 1; % calcium traces were updated (before mismatch between num of traces and planeIDs)
db(k).subject = 'Uma';
db(k).date = '2023-12-13';

k = k + 1;
db(k).subject = 'Dublin';
db(k).date = '2024-04-10';

k = k + 1;
db(k).subject = 'Ely';
db(k).date = '2024-06-20';