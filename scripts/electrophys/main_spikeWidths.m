%% Folders
folderData = 'C:\STORAGE\OneDriveData\SharePoint\Schroeder, Sylvia - Lab\DATA\electrophys\';

%% Define datasets
db_ephys_driftingGratings

%% plot spikewidths
spikeWidths = [];
for k = 1:length(db)
    data = load(fullfile(folderData, db(k).subject, db(k).date, ...
        sprintf('%02d_data.mat', db(k).exp)));
    spikeWidths = [spikeWidths; data.spikeWidths'];
end

bins = 0:.05:1.5;
figure
hist(spikeWidths, bins)
xlim(bins([1 end]))
xlabel('Spike width (ms)')
ylabel('# Neurons')