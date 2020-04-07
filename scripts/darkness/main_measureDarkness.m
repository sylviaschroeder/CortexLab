%% Load data
load('\\ZSERVER.cortexlab.net\Lab\Share\Sylvia\2019-02_Darkness\zmazeThreeFilterData_20190312.mat')

%% Constants, specs, ...
gain_60dB = 0.75e6; % in V/A
gain_70dB = 2.38e6;
filterLimits = [335 610; 412 569; 275 375];
responsivities = [0.3; 0.25; 0.13]; % in A/W, at 475/420/325 nm
diodeArea = 10000^2; % in um^2
irradianceAt500 = [2.4 2.12 2.01]2.52e18; % in photons/s/m^2, at 475/420/325 nm

%% Convert measures to Watt (and photons) and summarise
% conversion to Watt: data (V) / gain (V/A) / responsivity (A/W)

% use measurements at complete darkness (with cap) as reference -> units of
% final results will be in STDs of this measurement, converted to measures
% at 70 dB (highest sensitivity)
meanZero = mean(capLeft70dB) / gain_70dB / respAt500;
stdZero = std(capLeft70dB) / gain_70dB / respAt500;

% get means and STDs of all measures
meanOff = ([mean(darkLeft70dB) mean(darkCenter70dB) mean(darkRight70dB)] ./ ...
    gain_70dB ./ respAt500)  - meanZero;
stdOff = [std(darkLeft70dB) std(darkCenter70dB) std(darkRight70dB)] ./ ...
    gain_70dB ./ respAt500;

meanBlack = ([mean(blackLeft70dB) mean(blackCenter70dB) mean(blackRight70dB)] ./ ...
    gain_70dB ./ respAt500) - meanZero;
stdBlack = [std(blackLeft70dB) std(blackCenter70dB) std(blackRight70dB)] ./ ...
    gain_70dB ./ respAt500;
meanBlack60 = ([mean(blackLeft60dB) mean(blackCenter60dB) mean(blackRight60dB)] ./ ...
    gain_60dB ./ respAt500) - meanZero;
stdBlack60 = [std(blackLeft60dB) std(blackCenter60dB) std(blackRight60dB)] ./ ...
    gain_60dB ./ respAt500;

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
meanGray = (mean(gray,1) ./ gain_70dB ./ respAt500) - meanZero;
stdGray = std(gray,0,1) ./ gain_70dB ./ respAt500;
% for comparison: convert measurements measured with 60dB directly to Watt
meanGray60 = ([mean(grayLeft60dB) mean(grayCenter60dB) mean(grayRight60dB)] ./ ...
    gain_60dB ./ respAt500) - meanZero;
stdGray60 = [std(grayLeft60dB) std(grayCenter60dB) std(grayRight60dB)] ./ ...
    gain_60dB ./ respAt500;

% convert to photons/s/um^2
scale = 1/diodeArea * irradianceAt500;
% scale = 1/diodeDiameter * irradianceAt500 * (pi*0.0015^2); % multiply by pupil area to get photons/s

%% Plot results
x = [.8 1 1.2  1.8 2 2.2  2.8 3 3.2  3.8 4 4.2  4.8 5 5.2];

figure
hold on
errorbar(x, [meanOff meanBlack meanBlack60 meanGray meanGray60] .* scale, ...
    [stdOff stdBlack stdBlack60 stdGray stdGray60] .* scale, 'ko')
plot(0, stdZero * scale, 'kv')
set(gca, 'XTick', 0:5, 'XTickLabel', {'dark','screens off','black', ...
    'black (60dB)','gray','gray (60dB)'}, ...
    'box', 'off', 'YScale', 'log')
xlim([-.5 5.7])
ylabel('Luminance (in photons/s/um^2)')

%% Convert measures (relative to complete darkness) and summarise

% use measurements at complete darkness (with cap) as reference -> units of
% final results will be in STDs of this measurement, converted to measures
% at 70 dB (highest sensitivity)
meanZero = mean(capLeft70dB);
stdZero = std(capLeft70dB);

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