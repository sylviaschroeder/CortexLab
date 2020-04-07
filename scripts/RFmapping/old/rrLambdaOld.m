function [ lambda, RMSE, meanRMSE, VE ] = rrLambda( predX, cellResps, testLambda, cvType, cvPram, opF, groupVar )
%UNTITLED Summary of this function goes here
%   Finds optimal lambda for ridge regression

% Find optimal lambda parameter with crossvalidation. 
% cvtype sets cv partition type, 1 is kfolds, 2 is holdout, 3 is leave-one-out, 4 is contigous folds. cvPram sets
% number of folds, proportion for data for test set, or number of contigous folds.

rng('default')

if nargin == 7
    if cvType == 1
        cInd = cvpartition(groupVar, 'KFold', cvPram);
    elseif cvType == 2
        cInd = cvpartition(groupVar, 'Holdout', cvPram(1));
    else
        disp('Must use KFold or Holdout')
        return
    end
else
    if cvType == 1
        cInd = cvpartition(size(predX,1), 'KFold', cvPram);
    elseif cvType == 2
%         cInd = cvpartition(size(predX,1), 'Holdout', cvPram(1));
        [trainSetInd, testSetInd, finalSetInd] = dividerand(size(predX,1),1-sum(cvPram), cvPram(1), cvPram(2));
    elseif cvType == 3
        cInd = cvpartition(size(predX,1), 'LeaveOut');
        cvPram = cInd.NumObservations;
    elseif cvType == 4 % Divide data into n contiguous blocks
        [trainSetInd, testSetInd, finalSetInd] = divideblock(size(predX,1), 1-sum(cvPram), cvPram(1), cvPram(2)); 
    elseif cvType == 5 % For when an independent variables is a time series, split data so that training set comes before test set
        fullSet = 1:size(predX,1);
        finalSetInd = fullSet(end-round(cvPram(2)*size(predX,1)):end);
        testSetInd = fullSet(finalSetInd-round(cvPram(1)*size(predX,1)):finalSetInd-1);
        trainSetInd = fullSet(1:testSetInd(1)-1);
        cvType = 4; % Procedure is otherwise identical to type 4
    end
end

if nargin < 6
    opF = 0;
end


if cvType == 1 || cvType == 3
    
    if cvType == 1
        testRMSE = zeros(cvPram,length(testLambda),size(cellResps,2));
    else
        testRMSE = zeros(size(cellResps,2),length(testLambda),size(cellResps,2));
    end
    VE = zeros(size(cellResps,2), length(testLambda));
    
    for l = 1:length(testLambda)
        
         predY = cell(1,cvPram);
         realY = predY;
         
        
        for cv = 1:cvPram
            
            trainSetInd = logical(training(cInd,cv));
            
            testSetInd = logical(test(cInd,cv));
            
            trainX = predX(trainSetInd,:);
            
            trainY = cellResps(trainSetInd,:);
                
            trainX = [ones(size(trainX,1),1) trainX]; % Add intercept
            rVec = eye(size(trainX,2)).*testLambda(l);
            rVec(1,:) = []; % Remove lambda for intercept
            trainX = [trainX; rVec];
            trainY = [trainY; zeros(size(rVec,1),size(trainY,2))];
            
            if opF == 1
            
                disp(['Crossvalidation fold ', num2str(cv), ' lamda ', num2str(l),'.']);
            
            end
                                
            trainX = gpuArray(trainX);
            trainY = gpuArray(trainY);
            
            k = gather(trainX\trainY);
                                 
            testX = predX(testSetInd,:);
            testX = [ones(size(testX,1),1) testX];
            testY = cellResps(testSetInd,:);
            
            predY{cv} = testX*k;
            realY{cv} = testY;
            testRMSE(cv,l,:) = sqrt(nansum((testX*k - testY).^2)/size(testX,1));
            
        end
        
        allY = vertcat(realY{:});
        allPredY = vertcat(predY{:});
        SStot = sum((allY-mean(allY)).^2);
        SSres = sum((allPredY - allY).^2);
        VE(:,l) = 1 - SSres./SStot;
        
    end
    
    
elseif cvType == 2 || cvType == 4
    
    testRMSE = zeros(length(testLambda),size(cellResps,2));
    VE = zeros(size(cellResps,2), length(testLambda));
    
    for l = 1:length(testLambda)
        
%         if cvType == 2
%             
%             trainSetInd = logical(training(cInd));
%             
%             testSetInd = logical(test(cInd));
%             
%         end
        
        trainX = predX(trainSetInd,:);
        trainY = cellResps(trainSetInd,:);
        
        trainX = [ones(size(trainX,1),1) trainX];
        rVec = eye(size(trainX,2)).*testLambda(l);
        rVec(1,:) = [];
        trainX = [trainX; rVec];
        trainY = [trainY; zeros(size(rVec,1),size(trainY,2))];
        
        if opF == 1
            
                disp(['Crossvalidation by ', num2str(cvPram(1)), ' hold out', ' lamda ', num2str(l),'.']);
                
        end
        
        trainX = gpuArray(trainX);
        trainY = gpuArray(trainY);
        
        k = gather(trainX\trainY);
        
        testX = predX(testSetInd,:);
        testX = [ones(size(testX,1),1) testX];
        testY = cellResps(testSetInd,:);
        
        testRMSE(l,:) = sqrt(sum((testX*k - testY).^2)/size(testX,1));
        
       
          
    end
   
    
end


if cvType == 1 || cvType == 3
    meanRMSE = squeeze(nanmean(testRMSE)); % Take mean across CV folds
elseif cvType == 2 || cvType == 4
    meanRMSE = nanmean(testRMSE,2); 
end

if size(testRMSE,3) > 1 && length(testLambda) > 1
    meanRMSE = nanmean(meanRMSE,2); 
    [~,idx] = min(meanRMSE); % Find lambda with minimum RMSE
    
    lambda = testLambda(idx);
    if cvType == 1 || cvType == 3
        VE = VE(:,idx);
    else
        
        finalY = cellResps(finalSetInd,:);
        finalX = predX(finalSetInd,:);
        finalX = [ones(size(finalX,1),1) finalX];
            
        rVec = eye(size(finalX,2)).*lambda;
        rVec(1,:) = [];
        finalX = [finalX; rVec];
        finalY = [finalY; zeros(size(rVec,1),size(finalY,2))];
        
        SStot = sum((finalY-mean(finalY)).^2);
        SSres = sum((finalX*k - finalY).^2);
        
        VE = (1 - SSres./SStot)';
        
    end
    
    RMSE = meanRMSE(idx);
    
elseif length(testLambda) > 1
     [~,idx] = min(meanRMSE);
    
    lambda = testLambda(idx);
    
    if cvType == 1 || cvType == 3
        VE = VE(:,idx);  
    else
        finalY = cellResps(finalSetInd,:);
        finalX = predX(finalSetInd,:);
        finalX = [ones(size(finalX,1),1) finalX];
        
        
        rVec = eye(size(finalX,2)).*lambda;
        rVec(1,:) = [];
        finalX = [finalX; rVec];
        finalY = [finalY; zeros(size(rVec,1),size(finalY,2))];
        
        SStot = sum((finalY-mean(finalY)).^2);
        SSres = sum((finalX*k - finalY).^2);
        
        VE = (1 - SSres./SStot)';
    end
    RMSE = meanRMSE(idx); 
else
    lambda = testLambda;
    RMSE = mean(meanRMSE);    
end


end

