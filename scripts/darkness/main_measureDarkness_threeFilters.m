%% Load data
load('\\ZSERVER.cortexlab.net\Lab\Share\Sylvia\2019-02_Darkness\zmazeThreeFilterData_20190312.mat')

%% Constants, specs, ...
gain_60dB = 0.75e6; % in V/A
gain_70dB = 2.38e6;
filterLimits = [335 610; 412 569; 275 375];
responsivities = [0.3; 0.25; 0.13]; % in A/W, at 475/420/325 nm
diodeArea = 10000^2; % in um^2
irradiance = [2.4 2.12 1.64]'.*10^18; % in photons/s/m^2, at 475/420/325 nm

%% Convert measures to Watt (and photons) and summarise
% conversion to Watt: data (V) / gain (V/A) / responsivity (A/W)

% use measurements at complete darkness (with cap) as reference -> units of
% final results will be in STDs of this measurement, converted to measures
% at 70 dB (highest sensitivity)
meanZero = mean(cap70dB) / gain_70dB ./ responsivities;
stdZero = std(cap70dB) / gain_70dB ./ responsivities;
meanZero60 = mean(cap60dB) / gain_60dB ./ responsivities;
stdZero60 = std(cap60dB) / gain_60dB ./ responsivities;

% get means and STDs of all measures
meanOff = ([mean([S.darkLeft70dB])' mean([S.darkCenter70dB])' mean([S.darkRight70dB])'] ./ ...
    gain_70dB ./ responsivities) - meanZero;
stdOff = [std([S.darkLeft70dB])' std([S.darkCenter70dB])' std([S.darkRight70dB])'] ./ ...
    gain_70dB ./ responsivities;

meanBlack = ([mean([S.blackLeft70dB])' mean([S.blackCenter70dB])' mean([S.blackRight70dB])'] ./ ...
    gain_70dB ./ responsivities) - meanZero;
stdBlack = [std([S.blackLeft70dB])' std([S.blackCenter70dB])' std([S.blackRight70dB])'] ./ ...
    gain_70dB ./ responsivities;
meanBlack60 = ([mean([S.blackLeft60dB])' mean([S.blackCenter60dB])' mean([S.blackRight60dB])'] ./ ...
    gain_60dB ./ responsivities) - meanZero;
stdBlack60 = [std([S.blackLeft60dB])' std([S.blackCenter60dB])' std([S.blackRight60dB])'] ./ ...
    gain_60dB ./ responsivities;

% determine conversion from 60dB to 70dB
% we have two measurements under identical conditions using the 2
% sensitivities: mean (m) and mean+STD (m+s);
% we assume that measures under different sensitivities are related
% linearly. Therefore: m_2 = a*m_1 + b, and m_2+s_2 = a*(m_1+s_1) + b
% solving for a and b: a = s_2 / s_1, b = m_2 - (s_2/s_1)*m_1
a = stdBlack ./ stdBlack60;
% fprintf('slopes for 70->60dB (left, center, right): %.4f, %.4f, %.4f\n', a)
b = meanBlack - a .* meanBlack60;
% fprintf('intercepts for 70->60dB (left, center, right): %.4f, %.4f, %.4f\n', b)
a = mean(a,2);
b = mean(b,2);

% apply conversion to measurements of gray screen
gray60 = permute(cat(3, [S.grayLeft60dB], [S.grayCenter60dB], [S.grayRight60dB]), [2 3 1]);
meanGray = ((a .* mean(gray60,3) + b) ./ gain_70dB ./ responsivities) - meanZero;
stdGray = std(gray60,0,3) ./ gain_70dB ./ responsivities;
% for comparison: convert measurements measured with 60dB directly to Watt
meanGray60 = (mean(gray60,3) ./ gain_60dB ./ responsivities) - meanZero;
stdGray60 = std(gray60,0,3) ./ gain_60dB ./ responsivities;

% convert to photons/s/um^2
scale = 1/diodeArea .* irradiance;
% scale = 1/diodeDiameter * irradianceAt500 * (pi*0.0015^2); % multiply by pupil area to get photons/s

%% Plot results
x = [.8 1 1.2  1.8 2 2.2  2.8 3 3.2  3.8 4 4.2  4.8 5 5.2];

for f = 1:3
    figure
    hold on
    errorbar([0 x], [0 meanOff(f,:) meanBlack(f,:) meanBlack60(f,:) ...
        meanGray(f,:) meanGray60(f,:)] .* scale(f), ...
        [stdZero(f) stdOff(f,:) stdBlack(f,:) stdBlack60(f,:) ...
        stdGray(f,:) stdGray60(f,:)] .* scale(f), 'ko')
%     plot(0, stdZero(f) * scale(f), 'kv')
    set(gca, 'XTick', 0:5, 'XTickLabel', {'dark','screens off','black', ...
        'black (60dB)','gray','gray (60dB)'}, ...
        'box', 'off')
    xlim([-.5 5.7])
    ylabel('Luminance (in photons/s/um^2)')
    title(sprintf('Filter %s',S(f).Name))
end

%% Convert measures (relative to complete darkness) and summarise

% use measurements at complete darkness (with cap) as reference -> units of
% final results will be in STDs of this measurement, converted to measures
% at 70 dB (highest sensitivity)
meanZero = mean(cap70dB);
stdZero = std(cap70dB);

% get means and STDs of all measures
meanOff = [mean(darkLeft70dB) mean(darkCenter70dB) mean(darkRight70dB)];
stdOff = [std(darkLeft70dB) std(darkCenter70dB) std(darkRight70dB)];

meanBlack = [mean(blackLeft70dB) mean(blackCenter70dB) mean(blackRight70dB)];
stdBlack = [std(blackLeft70dB) std(blackCenter70dB) std(blackRight70dB)];

% determine conversion from 60dB to 70dB
% we have two measurements under identical conditions using the 2
% sensitivities: mean (m) and mean+STD (m+s);
% we assume that measures under different sensitivities are related
% linearly. Therefore: m_2 = a*m_1 + b, and m_2+s_2 = a*(m_1+s_1) + b
% solving for a and b: a = s_2 / s_1, b = m_2 - (s_2/s_1)*m_1
a = [std(blackLeft70dB)/std(blackLeft60dB) ...
    std(blackCenter70dB)/std(blackCenter60dB) ...
    std(blackRight70dB)/std(blackRight60dB)];
fprintf('slopes for 70->60dB (left, center, right): %.4f, %.4f, %.4f\n', a)
b = [mean(blackLeft70dB) - a(1) * mean(blackLeft60dB) ...
    mean(blackCenter70dB) - a(2) * mean(blackCenter60dB) ...
    mean(blackRight70dB) - a(3) * mean(blackRight60dB)];
fprintf('intercepts for 70->60dB (left, center, right): %.4f, %.4f, %.4f\n', b)
a = mean(a);
b = mean(b);

% apply conversion to measurements of gray screen
gray = a .* [grayLeft60dB grayCenter60dB grayRight60dB] + b;
meanGray = mean(gray,1);
stdGray = std(gray,0,1);

% scale all measures to STDs of complete darkness
meanOff = (meanOff - meanZero) ./ stdZero;
stdOff = stdOff ./ stdZero;
meanBlack = (meanBlack - meanZero) ./ stdZero;
stdBlack = stdBlack ./ stdZero;
meanGray = (meanGray - meanZero) ./ stdZero;
stdGray = stdGray ./ stdZero;

%% Plot results

x = [.8 1 1.2 1.8 2 2.2 2.8 3 3.2];

figure
hold on
errorbar(x, [meanOff meanBlack meanGray], [stdOff stdBlack stdGray], 'ko')
plot(0, 1, 'kv')
set(gca, 'XTick', 0:3, 'XTickLabel', {'dark','screens off','black','gray'}, ...
    'box', 'off', 'YScale', 'log')
xlim([-.5 3.7])
ylabel('Luminance (in STDs of dark)')