function plotPredictionErrors(errors)

% errors    structure resulting from 
%           gratings.getConditionedOrientationTuning_separated
%           each field: [1 x neuron]
%           fields: .cond_regr, cond_sep, nonCond_regre, nonCond_sep

maxi = max([errors.cond_regr, errors.cond_sep, errors.nonCond_regr, ...
    errors.nonCond_sep]);


figure('Position', [795 680 1070 420])

subplot(1,2,1)
plot(errors.cond_regr, errors.cond_sep, 'k.')
hold on
plot([0 maxi], [0 maxi], 'k:')
axis([0 maxi 0 maxi])
xlabel('Ridge regression')
ylabel('Stimulus-time separated')
title('Prediction errors under classification of non-visual data')

subplot(1,2,2)
plot(errors.nonCond_regr, errors.nonCond_sep, 'k.')
hold on
plot([0 maxi], [0 maxi], 'k:')
axis([0 maxi 0 maxi])
xlabel('Ridge regression')
ylabel('Stimulus-time separated')
title('Prediction errors without classification of non-visual data')


figure('Position', [795 680 1070 420])

subplot(1,2,1)
plot(errors.cond_regr, errors.nonCond_regr, 'k.')
hold on
plot([0 maxi], [0 maxi], 'k:')
axis([0 maxi 0 maxi])
xlabel('With classification')
ylabel('Without classification')
title('Prediction errors of ridge regression filters')

subplot(1,2,2)
plot(errors.cond_sep, errors.nonCond_sep, 'k.')
hold on
plot([0 maxi], [0 maxi], 'k:')
axis([0 maxi 0 maxi])
xlabel('With classification')
ylabel('Without classification')
title('Prediction errors of stimulus-time separated filters')