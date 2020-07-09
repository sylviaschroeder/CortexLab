subject = 'SS088';
dates = {'2018-04-17', '2018-04-18', '2018-04-20', '2018-04-22', '2018-04-23', ...
    '2018-04-30', '2018-05-01', '2018-05-02', '2018-05-03', ...
    '2018-05-04'};

dataFolder = '\\zubjects.cortexlab.net\Subjects';
resultsFolder = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\behaviour\';

session = 1;
outcome = [];

for d = 1:length(dates)
    exps = dir(fullfile(dataFolder, subject, dates{d}));
    exps = cellfun(@str2num, {exps.name}, 'UniformOutput', false);
    for ex = 1:length(exps)
        if isempty(exps{ex})
            continue
        end
        blockFile = ls(fullfile(dataFolder, subject, dates{d}, ...
            num2str(exps{ex}), '*_Block.mat'));
        if isempty(blockFile)
            continue
        end
        data = load(fullfile(dataFolder, subject, dates{d}, ...
            num2str(exps{ex}), blockFile));
        block = data.block;
        
        resp = block.events.responseValues;
        resp(resp == 1) = 2; % go right
        resp(resp == 0) = 3; % no-go
        resp(resp == -1) = 1; % go left
        numTrials = length(resp);
        if numTrials < 20
            continue
        end
        laserType = double(cat(2,block.paramsValues.laserAmp)>0)';
        if size(laserType,2) < 2
            laserType = repmat(laserType, 1, 2);
        end
        
        out.response = resp';
        out.contrast = cat(2, block.paramsValues.stimulusContrast)';
        out.repeatNum = block.events.repeatNumValues';
        out.laserType = laserType;
        out.feedback = double(block.events.feedbackValues');
        out.session = ones(numTrials,1) .* session;
        out.date = repmat(dates(d), numTrials, 1);
        out.expNum = ones(numTrials,1) .* exps{ex};
        
        if isempty(outcome)
            outcome = out;
        else
            fields = fieldnames(outcome);
            for f = 1:length(fields)
                outcome.(fields{f}) = [outcome.(fields{f}); out.(fields{f})];
            end
%             outcome = [outcome; out];
        end
        
        session = session + 1;
    end
end

save(fullfile(resultsFolder, subject, 'SC_inactivation_data.mat'), 'outcome')