animal = 'M150410_SS044';
date = '2015-04-28';
exp = 2:5;
planes = [2:5];
neurons = {[8 16 17 40 46 47 72 73],[5 27 40 45 113 126 137 141 142 211], ...
    [11 33 34 37 38 75 114 173 208 215],[22 28 76 133]};
% neurons = {[1:100]};

gap = 100;

folderROIData = 'C:\DATA\InfoStructs';
    
folder = [folderROIData filesep animal filesep date filesep];
infos = [];

for iExp = 1:length(exp)
    f = [folder filesep num2str(exp(iExp))];
    fileStart = [date '_' num2str(exp(iExp)) '_' animal];
    for iPlane = 1:length(planes)
        file = sprintf('%s_2P_plane%03d_ROI.mat',fileStart,planes(iPlane));
        data=load(fullfile(f, file));
        if iExp == 1 && iPlane == 1
            infos = data.meta;
        else
            infos(iExp,iPlane) = data.meta;
        end
    end
end

for iPlane = 1:length(planes)
    F = [];
    Npil = [];
    Fcorr = [];
%     slopes = [];
    for iExp = 1:length(exp)
        F = [F; infos(iExp,iPlane).F(:,neurons{iPlane}); ...
            NaN(gap,length(neurons{iPlane}))];
        Npil = [Npil; infos(iExp,iPlane).Npil(:,neurons{iPlane}); ...
            NaN(gap,length(neurons{iPlane}))];
        Fcorr = [Fcorr; infos(iExp,iPlane).Fcorr(:,neurons{iPlane}); ...
            NaN(gap,length(neurons{iPlane}))];
%         slopes = [slopes; infos(iExp,iPlane).NpilSlopes(:,neurons{iPlane}); ...
%             NaN(gap,length(neurons{iPlane}))];
    end
    mins = min([F;Npil;Fcorr], [], 1);
    maxs = max([F;Npil;Fcorr], [], 1);
    for iCell = 1:length(neurons{iPlane})
        figure('Position', [1920 350 1920 750])
%         figure('Position', [4 350 1920 750])
        subplot(2,1,1)
        plot(F(:,iCell))
        hold on
        plot(Npil(:,iCell))
        ylim([mins(iCell) maxs(iCell)])
        title(sprintf('Plane %d, Cell %d', planes(iPlane), neurons{iPlane}(iCell)))
        ylabel('F')
        ax1 = gca;
        subplot(2,1,2)
        plot(Fcorr(:,iCell))
        ylim([mins(iCell) maxs(iCell)])
%         plot(slopes(:,iCell))
        ylabel('Slopes')
        xlabel('Time (in frames)')
        linkaxes([ax1 gca],'x')
    end
end