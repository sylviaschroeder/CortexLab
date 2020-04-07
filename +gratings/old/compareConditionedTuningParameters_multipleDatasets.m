function compareConditionedTuningParameters_multipleDatasets

condition = 'pupilSize';

animal(1).name = 'M141002_SS026';
animal(1).expDate(1).name = '2014-10-29';
animal(1).expDate(1).exp(1).ID = 4;
animal(1).expDate(1).exp(1).planes = 2:3;

animal(2).name = 'M141003_SS027';
animal(2).expDate(1).name = '2014-10-23';
animal(2).expDate(1).exp(1).ID = 2;
animal(2).expDate(1).exp(1).planes = 1:4;

animal(2).expDate(2).name = '2014-10-30';
animal(2).expDate(2).exp(1).ID = 3;
animal(2).expDate(2).exp(1).planes = 2:4;

animal(3).name = 'M141007_SS029';
animal(3).expDate(1).name = '2014-11-12';
animal(3).expDate(1).exp(1).ID = 4;
animal(3).expDate(1).exp(1).planes = 2:3;

fitParameters = [];
adjustedRsquared = [];

for an = 1:length(animal)
    for session = 1:length(animal(an).expDate);
        for exp = 1:length(animal(an).expDate(session).exp)
            info = ppbox.infoPopulate(animal(an).name, ...
                animal(an).expDate(session).name, ...
                animal(an).expDate(session).exp(exp).ID);
            [negativeCond, positiveCond, conditionTime, labels] = ...
                ssLocal.classifyData(info, condition, 1);
            for plane = animal(an).expDate(session).exp(exp).planes
                filePath = fullfile(info.folderProcessed, ...
                    sprintf('%s_plane%03d_ROI', info.basename2p, plane));
                % load meta
                load(filePath, 'meta')
                [~, stimSequence, stimMatrix, frameTimes] = ...
                    ssLocal.getStimulusResponseInfo(meta);
                ind = strcmp(meta.ROI.CellClasses, 's');
                [pars, Rs] = ssLocal.compareOrientationBetweenConditions( ...
                    meta.F(:,ind), frameTimes, ...
                    [negativeCond, positiveCond], conditionTime, labels, ...
                    stimMatrix, stimSequence, plane);
                fitParameters = [fitParameters; pars];
                adjustedRsquared = [adjustedRsquared; Rs];
            end
        end
    end
end

ssLocal.compareConditionedTuningParameters(fitParameters, labels);
ssLocal.compareConditionedTuningParameters(fitParameters, labels, adjustedRsquared);