%% Folders
% file location depending on PCs
folderPC = 'C:\STORAGE\OneDrive - University College London';
% folderPC = 'C:\Users\skgtchr\OneDrive - University College London'; % ZEN notebook
% folderPC = 'C:\Storage\OneDrive - University College London'; % Lenovo laptop
% data folders
folderROIData = fullfile(folderPC, 'Lab\DATA\InfoStructs');

%% Go through each dataset
make_db_boutons
correct = false(1,length(db));

for k = 1:length(db)
    folder = fullfile(folderROIData, db(k).mouse_name, ...
        db(k).date);
    file = [db(k).date '_%d_' db(k).mouse_name '_2P_plane%03d_ROI.mat'];
    for iPlane = 1 %:length(db(k).planesToProcess)
        F = [];
        F0 = [];
        F_woNpil = [];
        F_delta = [];
        F_final = [];
        Npil = [];
        slopes = [];
        intercepts = [];
%         redSlopes = [];
        
        for iExp = 1:length(db(k).expts)
            % load raw plane data for gratings
%             data = load(fullfile(folder, num2str(db(k).expts(iExp)), ...
%                 sprintf(file,db(k).expts(iExp), db(k).planesToProcess(iPlane))));
            data = load(fullfile(folder, num2str(iExp), ...
                sprintf(file,iExp, db(k).planesToProcess(iPlane))));
            meta = data.meta;
            F = [F; meta.F];
            F0 = [F0; meta.F0];
            F_woNpil = [F_woNpil; meta.F_woNpil];
            F_delta = [F_delta; meta.F_delta];
            F_final = [F_final; meta.F_final];
            Npil = [Npil; meta.Npil];
            slopes = [slopes; meta.NpilSlopes];
            intercepts = [intercepts; meta.NpilIntercepts];
%             redSlopes = [redSlopes; meta.redFiltSlopes];
%             if iExp == 1
%                 slopes = meta.NpilSlopes;
%                 intercepts = meta.NpilIntercepts;
%                 redSlopes = meta.redFiltSlopes;
%             end
        end
        if all(all(slopes == slopes(1,:),1))
            correct(k) = true;
        end
    end
end