%% Load database
db_driftingGratings

%% Folders
% loading data
folderROIData = 'C:\STORAGE\OneDriveData\SharePoint\Schroeder, Sylvia - Lab\DATA\InfoStructs';
folderResults = 'C:\STORAGE\OneDriveData\SharePoint\Schroeder, Sylvia - Lab\RESULTS\Preprocessing\greenFiltered_plots';

%% Parameters
smoothStd = 0.25; % in sec
lowCutOff = 0.5;
dist = 8;
shuffles = 200;

%% Filter red and green traces
for k = 1:length(db)
    folder = [folderROIData filesep db(k).subject filesep ...
        db(k).date filesep num2str(db(k).exp)];
    fileStart = [db(k).date '_' num2str(db(k).exp) '_' ...
        db(k).subject];
    file = [fileStart '_2P_plane%03d_ROI.mat'];
    running = [];
    for iPlane = 1:length(db(k).planes)
        if iPlane > 1 && isempty(running)
            continue
        end
        fPlots = fullfile(folderResults, [db(k).subject '_' db(k).date '_' ...
            num2str(db(k).exp)], num2str(db(k).planes(iPlane)));
        if ~exist(fPlots, 'dir')
            mkdir(fPlots)
        end
        % load meta
        data=load(fullfile(folder, sprintf(file,db(k).planes(iPlane))));
        meta = data.meta;
        meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', 'cortexlab.net');
        frameTimes = ppbox.getFrameTimes(meta);
        sr = 1 / median(diff(frameTimes));
        Wn = lowCutOff / (sr/2);
        b = fir1(42, Wn, 'high');
        hamWin = hamming(round(5 * sr));
        hamWin = hamWin ./ sum(hamWin);
        
        if iPlane == 1
            ballData = nonVis.getRunningSpeed(meta);
            if ~isempty(ballData)
                stdSamples = round(smoothStd / median(diff(ballData.t)));
                convWindow = normpdf(-3*stdSamples:3*stdSamples, 0, stdSamples);
                running = conv(ballData.total, convWindow, 'same');
                runningTime = ballData.t;
            else
                running = [];
                continue
            end
        end
        runPlane = interp1(runningTime, running, frameTimes, 'pchip')';
        shifts = randi(length(runPlane), 1, shuffles);
        for iCell = 1:size(meta.R,2)
            if meta.ROI.isDuplicate(iCell) == 1 || meta.ROI.isSwitchOn(iCell) == 1
                continue
            end
            R_f = filtfilt(b, 1, meta.R(2:end,iCell));
            R_c = conv(abs(R_f), hamWin, 'same');
            F_f = filtfilt(b, 1, meta.F_final(2:end,iCell));
            F_c = conv(abs(F_f), hamWin, 'same');
            
            % shuffle test for significance of R
            R_c_sh = NaN(length(R_c), shuffles);
            F_c_sh = NaN(length(R_c), shuffles);
            F_final_sh = NaN(length(R_c), shuffles);
            for sh = 1:shuffles
                R_c_sh(:,sh) = circshift(R_c, shifts(sh));
                F_c_sh(:,sh) = circshift(F_c, shifts(sh));
                F_final_sh(:,sh) = circshift(meta.F_final(2:end,iCell), shifts(sh));
            end
            % (1) filtered and convolved red trace
            r_R_c = corr(runPlane(2:end), R_c_sh);
            r_F_c = corr(runPlane(2:end), F_c_sh);
            r_F_final = corr(runPlane(2:end), F_final_sh);
        
            figure('Position', [1 41 1920 1083])
            ax = [0 0 0];
            
            subplot(8,5,1:4)
            plot(frameTimes(2:end), runPlane(2:end), 'k')
            axis tight
            ylabel('Running')
            title(sprintf('Cell %d', iCell))
            ax(1) = gca;
            
            subplot(8,5,[6:9, 11:14])
            plot(frameTimes(2:end), (meta.R(2:end,iCell) - mean(meta.R(2:end,iCell))) ./ ...
                std(meta.R(2:end,iCell)), 'Color', [0.5 0 0])
            hold on
            plot(frameTimes(2:end), (R_f - mean(R_f)) ./ std(R_f)-dist, 'Color', [.8 .2 0])
            plot(frameTimes(2:end), (R_c - mean(R_c)) ./ std(R_c)-2*dist, ...
                'Color', [.6 .4 0], 'LineWidth', 2)
            axis tight
            legend('R','R filt','R conv')
            ax(2) = gca;
            
            subplot(8,5,setdiff(16:39, 20:5:35))
            plot(frameTimes(2:end), (meta.F(2:end,iCell) - mean(meta.F(2:end,iCell))) ./ ...
                std(meta.F(2:end,iCell)), 'Color', [0 0 .6])
            hold on
            plot(frameTimes(2:end), (meta.F_woNpil(2:end,iCell) - mean(meta.F_woNpil(2:end,iCell))) ./ ...
                std(meta.F_woNpil(2:end,iCell))-dist, 'Color', [0 0.2 .8])
            plot(frameTimes(2:end), (meta.F_delta(2:end,iCell) - mean(meta.F_delta(2:end,iCell))) ./ ...
                std(meta.F_delta(2:end,iCell))-2*dist, 'Color', [0 0.6 .6])
            plot(frameTimes(2:end), (meta.F_final(2:end,iCell) - mean(meta.F_final(2:end,iCell))) ./ ...
                std(meta.F_final(2:end,iCell))-3*dist, 'Color', [0 0.8 .2])
            plot(frameTimes(2:end), (F_f - mean(F_f)) ./ std(F_f)-4*dist, ...
                'Color', [0 0.6 0])
            plot(frameTimes(2:end), (F_c - mean(F_c)) ./ std(F_c)-5*dist, ...
                'Color', [0 0.3 0], 'LineWidth', 2)
            axis tight
            legend('F','wo npil','delta','red corrected','F_{red} filt','F_{red} conv')
            ax(3) = gca;
            linkaxes(ax, 'x')
            
            subplot(3, 5, 5)
            plot(runPlane(2:end), R_c, '.', 'Color', [.6 .4 0])
            axis tight
            title(sprintf('R = %.2f [%.2f  %.2f]', corr(runPlane(2:end), ...
                R_c), prctile(r_R_c,[2.5 97.5])))
            xlabel('Running')
            ylabel('F red corrected')
            set(gca, 'box', 'off')
            
            subplot(3, 5, 10)
            plot(runPlane(2:end), meta.F_woNpil(2:end,iCell), '.', 'Color', [0 0.8 .2])
            axis tight
            title(sprintf('R = %.2f [%.2f  %.2f]', corr(runPlane(2:end), ...
                meta.F_final(2:end,iCell)), prctile(r_F_final,[2.5 97.5])))
            xlabel('Running')
            ylabel('F red corrected')
            set(gca, 'box', 'off')
            
            subplot(3, 5, 15)
            plot(runPlane(2:end), F_c, '.', 'Color', [0 0.3 0])
            axis tight
            title(sprintf('R = %.2f [%.2f  %.2f]', corr(runPlane(2:end), ...
                F_c), prctile(r_F_c,[2.5 97.5])))
            xlabel('Running')
            ylabel('Filtered F')
            set(gca, 'box', 'off')
            
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            print(fullfile(fPlots, sprintf('filteredGreen%03d.jpg', iCell)), ...
                '-djpeg','-r0')
            close gcf
        end
    end
end