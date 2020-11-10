dataset = 10;
planes = 2:4; % 2:5;
exp = 1;
folderSVD = 'C:\STORAGE\OneDriveData\SharePoint\Schroeder, Sylvia - Lab\DATA\F\M150611_SS048\2015-12-02\1_2_3_4_6';
folderSVDRed = 'C:\STORAGE\OneDriveData\SharePoint\Schroeder, Sylvia - Lab\DATA\F\M150611_SS048\2015-12-02\1_2_3_4_6\red';
fileSVD = 'SVD_M150611_SS048_2015-12-02_plane%d.mat';
folderMeta = 'C:\STORAGE\OneDriveData\SharePoint\Schroeder, Sylvia - Lab\DATA\InfoStructs\M150611_SS048\2015-12-02\1'; % 'C:\DATA\InfoStructs\M150305_SS041\2015-04-23\3';
fileMeta = '2015-12-02_1_M150611_SS048_2P_plane%03d_ROI.mat';

smoothing = 3; % sec
filtPoly = 3;
smoothStd = .25; % sec (to divide into running and not running)

data = load(fullfile(folderMeta, sprintf(fileMeta,planes(1))));
meta = data.meta;
meta.folderTL = strrep(meta.folderTL, 'ioo.ucl.ac.uk', 'cortexlab.net');
frameTimes = ppbox.getFrameTimes(meta);
ballData = nonVis.getRunningSpeed(meta);

%% calculate average movie frame during running vs. not running
stdSamples = round(smoothStd / median(diff(ballData.t)));
convWindow = normpdf(-3*stdSamples:3*stdSamples, 0, stdSamples);
running2 = conv(ballData.total, convWindow, 'same');
running2 = running2 / median(diff(ballData.t)) / 53;
running2 = interp1(ballData.t, running2, frameTimes, 'pchip');
running2 = running2 > 5;

for iPlane = planes
    data = load(fullfile(folderSVD, sprintf(fileSVD, iPlane)));
    U = data.U;
    x = size(U,2);
    y = size(U,1);
    U = double(reshape(U,x*y,[]));
    V = data.Vcell;
    V = double(V{exp});
    
    img_running = ( U * sum(V(:,running2==1),2) ) ./ sum(running2);
    img_running = reshape(img_running, y, x);
    img_still =  ( U * sum(V(:,running2==0),2) ) ./ sum(running2==0);
    img_still = reshape(img_still, y, x);
    
    figure
    subplot(1,3,1)
    imagesc(img_running)
    axis image
    subplot(1,3,2)
    imagesc(img_still)
    axis image
    subplot(1,3,3)
    imagesc(img_running-img_still)
    axis image
end

%% calculate correlation maps
filtWindow = ceil(smoothing / median(diff(ballData.t)));
if mod(filtWindow,2) == 0
    filtWindow = filtWindow-1;
end
total = ballData.total ./ median(diff(ballData.t)) ./ 53;
running = sgolayfilt(total, filtPoly, filtWindow);
running = interp1(ballData.t, running, frameTimes, 'pchip')';
T = length(running);

for iPlane = planes
    data = load(fullfile(folderSVDRed, sprintf(fileSVD, iPlane)));
    U_r = data.U;
    x = size(U_r,2);
    y = size(U_r,1);
    U_r = double(reshape(U_r,x*y,[]));
    V_r = data.Vcell;
    V_r = double(V_r{exp});
    
    avg_r = U_r * sum(V_r,2);
    rhos_red_running = (U_r * (V_r * running) - avg_r * mean(running)) ./ ...
        (sqrt(sum(U_r * (V_r*V_r') .* U_r,2) - avg_r.^2./T) .* ...
        sqrt(running'*running - mean(running)^2*T));
    rhos_red_running = reshape(rhos_red_running,y,x,1);
    
    data = load(fullfile(folderSVD, sprintf(fileSVD, iPlane)));
    U_g = data.U;
    x = size(U_g,2);
    y = size(U_g,1);
    U_g = double(reshape(U_g,x*y,[]));
    V_g = data.Vcell;
    V_g = double(V_g{exp});
    
    avg_g = U_g * sum(V_g,2);
    rhos_green_running = (U_g * (V_g * running) - avg_g * mean(running)) ./ ...
        (sqrt(sum(U_g * (V_g*V_g') .* U_g,2) - avg_g.^2./T) .* ...
        sqrt(running'*running - mean(running)^2*T));
    rhos_green_running = reshape(rhos_green_running,y,x,1);
    
    rhos_green_red = (sum(U_g * (V_g * V_r') .* U_r, 2) - ...
        avg_g .* avg_r ./ T) ./ ...
        (sqrt(sum(U_g * (V_g*V_g') .* U_g,2) - avg_g.^2./T) .* ...
        sqrt(sum(U_r * (V_r*V_r') .* U_r,2) - avg_r.^2./T));
    rhos_green_red = reshape(rhos_green_red,y,x,1);
    
    % plot average green fluorescence image
    figure('Position',[2 622 531 494])
    imagesc(reshape(avg_g,y,x)./T)
    greens = gray;
    greens(:,[1 3]) = 0;
    colormap(greens)
    axis image
    title(sprintf('Plane %d: Mean green image', iPlane))
    
    % plot correlation map: green fluorescence and running
    maxi = max(abs(rhos_green_running(:)));
    figure('Position',[536 622 596 494])
    imagesc(rhos_green_running,[-maxi maxi])
    colormap parula
    colorbar
    axis image
    title(sprintf('Plane %d: Green with running',iPlane))
    
    % plot mixture of average green and correlation with running
    g = (avg_g - min(avg_g(:))) ./ (max(avg_g(:)) - min(avg_g(:)));
    g = reshape(g, y, x);
    r = rhos_green_running ./ maxi;
    r(r < 0) = 0;
    b = rhos_green_running ./ maxi;
    b(b > 0) = 0;
    b = -b;
    image = cat(3, r, g, b);
    figure
    imshow(image)
    title(sprintf('Plane %d: Green: avg. green, Red-Blue: pos.-neg. corr. with running', iPlane))
    set(gcf,'Position',[1136 622 531 494])
    
    % plot average red fluorescence image
    figure('Position',[2 42 531 494])
    imagesc(reshape(avg_r,y,x)./T)
    reds = gray;
    reds(:,[2 3]) = 0;
    colormap(reds)
    axis image
    title(sprintf('Plane %d: Mean red image', iPlane))
    
    % plot correlation map: red fluorescence and running
    maxi = max(abs(rhos_red_running(:)));
    figure('Position',[536 42 596 494])
    imagesc(rhos_red_running,[-maxi maxi])
    colormap parula
    colorbar
    axis image
    title(sprintf('Plane %d: Red with running',iPlane))
    
    % plot mixture of average green and correlation with running
    g = (avg_r - min(avg_r(:))) ./ (max(avg_r(:)) - min(avg_r(:)));
    g = reshape(g, y, x);
    r = rhos_red_running ./ maxi;
    r(r < 0) = 0;
    b = rhos_red_running ./ maxi;
    b(b > 0) = 0;
    b = -b;
    image = cat(3, r, g, b);
    figure
    imshow(image)
    title(sprintf('Plane %d: Green: avg. red, Red-Blue: pos.-neg. corr. with running', iPlane))
    set(gcf,'Position',[1136 42 531 494])
    
    % plot correlation map: red and green fluorescence
    maxi = max(abs(rhos_green_red(:)));
    figure('Position',[1924 622 596 494])
    imagesc(rhos_green_red,[-maxi maxi])
    colormap parula
    colorbar
    axis image
    title(sprintf('Plane %d: Green with red',iPlane))
    
    % plot mixture of average green, average red and corr. green-red
    g = (avg_g - min(avg_g(:))) ./ (max(avg_g(:)) - min(avg_g(:)));
    g = reshape(g, y, x);
    r = (avg_r - min(avg_r(:))) ./ (max(avg_r(:)) - min(avg_r(:)));
    r = reshape(r, y, x);
    b = rhos_green_red ./ maxi;
    b(b < 0) = 0;
    image = cat(3, r, g, b);
    figure
    imshow(image)
    title(sprintf('Plane %d: Green: avg. green, Red: avg. red, Blue: green-red corr.', iPlane))
    set(gcf,'Position',[2524 622 531 494])
    
    
    % get semi-partial correlation between green and running, after
    % removing red from green
    % linear regression to predict green from red for each pixel
    var = sum((U_r * (V_r*V_r')) .* U_r,2);
    cross = sum((U_r * (V_r*V_g')) .* U_g,2);
    intercepts = (var .* avg_g - avg_r .* cross) ./ (T .* var - avg_r.^2);
    slopes = (T * cross - avg_r .* avg_g) ./ (T .* var - avg_r.^2);
    
%     residuals = 
end