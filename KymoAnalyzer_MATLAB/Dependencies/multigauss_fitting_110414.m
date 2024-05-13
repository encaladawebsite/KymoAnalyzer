function multigauss_fitting(data)

global mu 
global sigmasq
global prob 
global sigma 
global clusterNum
global BIC
global classification

% Script was modified from Laptrack71 code by SN to run multi Gaussian
% Fitting for KymoAnalyzer Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: load the data.
%
% If data size is too large, there will be a memory issue. So some
% downsampling is required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleData = data;
if (length(sampleData) > 5000)
    DS_ratio = round(length(sampleData) / 5000.);
    sampleData = sampleData(1 : DS_ratio : end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: basic normality test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normality_test_with_R_110414(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: launch R connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the MClust here
errFlag = 0;
try
    status = openR;
catch
    disp('Unable to launch R');
    errFlag = 1;
    
end

if ~status
    if ~strcmp(lastwarn,'Already connected to an R server.')
        disp('Unable to launch R');
        errFlag = 1;
    else
        fprintf('Already connected to an R server.\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: load the MCLUST package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[result, runStatus] = evalR('library(mclust)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: clustering parameter control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Commented by Sylvia: removes manual input of Cluster Analysis
% Configuration, instead values defined now to default values

% prompt={'Please set the lower bound for cluster number','Please set the upper bound for cluster number', 'Set a cutoff max value', 'Initial intensity correction coefficient'};
% name='Configure cluster analysis';
% numlines=1;
% defaultanswer={'1', '9', '1e30', '1'};
% answer=inputdlg(prompt,name,numlines,defaultanswer);
% 
% minModeNum = str2num(char(answer{1}));
% maxModeNum = str2num(char(answer{2}));
% cutoffMax = str2double(char(answer{3}));
% intensityCorrectionCoef = str2double(char(answer{4}));

minModeNum=1;
maxModeNum=9;
cutoffMax=1e30;
intensityCorrectionCoef=1;

sampleData = sampleData * intensityCorrectionCoef;
sampleData = sampleData(find(sampleData <= cutoffMax));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: load the data and run the clustering
% In the following, the code is changed to reflect the new functions in
% MCLUST v3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % put observations into R
putRdata('inputData', sampleData);

tempStr = ['clusterBIC <- mclustBIC(inputData,', num2str(minModeNum), ':', num2str(maxModeNum), ')'];
evalR(sprintf(tempStr));
evalR(sprintf('clusterModel <- mclustModel(inputData, clusterBIC)'));
mu = evalR('clusterModel$parameters$mean');
sigmasq = evalR('clusterModel$parameters$variance$sigmasq');
prob = evalR('clusterModel$parameters$pro');
sigma = sqrt(sigmasq);

% take care of possible equal sigma
if length(sigma) == 1
    sigma = repmat(sigma,1,length(mu));
    sigmasq = repmat(sigmasq, 1, length(mu));
end

BIC = evalR('clusterBIC');
classification = evalR('clusterModel$z');

% % find number of Gaussians
numGauss = evalR('clusterModel$G');
clusterNum = numGauss;

end
 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 7: exit from R, plotting results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lowerBound = 0.00;
% sampleData = sampleData(find(sampleData >= lowerBound));
% 
% ButtonName=questdlg('Would you like to set bin centers directly?', 'Config histogram generation',...
%                        'Yes', 'No','No');
% switch ButtonName,
%     case 'Yes',
%         prompt={'Please enter histogram range minimum', 'Please enter histogram range maximum', 'Please enter histogram bin number'};
%         name='Multiple Gauss Fit';
%         numlines=1;
%         defaultanswer={'0', '4', '30'};
%         answer=inputdlg(prompt,name,numlines,defaultanswer);
%         histogramMin = str2num(char(answer{1}));
%         histogramMax = str2num(char(answer{2}));
%         
%         %% histogramBinSize = str2num(char(answer{3}));
%         histogramBinNum = str2num(char(answer{3}));
%         histogramBinCenter = linspace(histogramMin, histogramMax,histogramBinNum);
%         binNum = length(histogramBinCenter);
%         actualDistX = histogramBinCenter;
%         HR_num = 200;  %HR stands for HIGH-RESOLUTION
%         HR_actualDistX = linspace(-0.5, histogramMax, HR_num);
%         %%HR_actualDistX = linspace(-0.05, 2.5, HR_num);
%     case 'No'
%         prompt={'Please enter histogram bin number'};
%         name='Multiple Gauss Fit';
%         numlines=1;
%         defaultanswer={'30'};
%         answer=inputdlg(prompt,name,numlines,defaultanswer);
%         binNum = str2num(char(answer{1}));
%         actualDistX = linspace(-lowerBound, max(sampleData), binNum);
%         HR_num = 200;  %HR stands for HIGH-RESOLUTION
%         HR_actualDistX = linspace(-lowerBound, max(sampleData), HR_num);
%         %% HR_actualDistX = linspace(-0.05, 2.5, HR_num);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 8: calculate separated clusters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one_over_sqrt_2pi = 1 / sqrt(2 * pi);
% separatedCluster = struct('gridY', [], 'gridX', []);
% 
% for i = 1 : clusterNum
%     separatedCluster(i).gridX = HR_actualDistX;
%     for j = 1 : HR_num
%         temp = - 0.5 * (HR_actualDistX(j) - mu(i))^2 / sigmasq(i);
%         separatedCluster(i).gridY = [separatedCluster(i).gridY; prob(i) * one_over_sqrt_2pi / sqrt(sigmasq(i)) * exp(temp)]; 
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 9: calculate the sum of clusters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % normalization later on.
% for i = 1 : HR_num
%     HR_actualDistY(i) = 0;
%     for j = 1 : clusterNum
%         temp = - 0.5 * (HR_actualDistX(i) - mu(j))^2 / sigmasq(j);
%         HR_actualDistY(i) = HR_actualDistY(i) + prob(j) * one_over_sqrt_2pi / sqrt(sigmasq(j)) * exp(temp); 
%     end
%     tempSum = 0;
%     for k = 1 : clusterNum
%         tempSum = tempSum + separatedCluster(k).gridY(i);
%     end
%     
%     if (abs(tempSum - HR_actualDistY(i)) > 1e-6)
%         uiwait(warndlg('Error in fitting calculation'));
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 10: calculate source RAW distribution curve
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% actualDistY = zeros(length(actualDistX), 1);
% for i = 1 : length(actualDistX)
%     actualDistY(i) = 0;
%     for j = 1 : clusterNum
%         temp = - 0.5 * (actualDistX(i) - mu(j))^2 / sigmasq(j);
%         actualDistY(i) = actualDistY(i) + prob(j) * one_over_sqrt_2pi / sqrt(sigmasq(j)) * exp(temp); 
%     end
% end
% % figure;
% % plot(actualDistX, actualDistY, 'r-');
% 
% % % Make a simple integration check
% % stepSize = actualDistX(2) - actualDistX(1);
% % tempIntegration = sum(actualDistY) * stepSize;
% 
% [histN, histCenter] = hist(sampleData, actualDistX);
% histArea = (histCenter(2) - histCenter(1)) * sum(histN);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fit by area matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histH = figure;
% hist(sampleData, actualDistX);
% hold on;
% plot(HR_actualDistX, HR_actualDistY  * histArea, 'r-');
% for i = 1 : clusterNum
%     plot(separatedCluster(i).gridX, separatedCluster(i).gridY  * histArea, 'c-');
% end
% title('Fit by area-matching');
% 
% % [peakY, idx] = max(HR_actualDistY);
% % fprintf('Peak of distribution is at = %f \n', HR_actualDistX(idx));
% % 
% % if strcmp(ButtonName, 'Yes')
% %     xlim([-0.1, histogramMax]);
% % else
% %     xlim([-0.1, max(sampleData)]);
% % end
% 
% % George's Code (can be removed)
% [GECN,GECX] = hist(sampleData, actualDistX);
% histH = figure;
% bar(GECX,100*GECN./sum(GECN),1);
% ylabel('Percentage (%)');
% hold on;
% plot(HR_actualDistX, (HR_actualDistY  * histArea)/sum(GECN)*100, 'r-');
% for i = 1 : clusterNum
%     plot(separatedCluster(i).gridX, (separatedCluster(i).gridY  * histArea)/sum(GECN)*100, 'c-');
% end
% title('Fit by area-matching (%)');
% 
% % [peakY, idx] = max(HR_actualDistY);
% % fprintf('Peak of distribution is at = %f \n', HR_actualDistX(idx));
% % 
% % if strcmp(ButtonName, 'Yes')
% %     xlim([-0.1, histogramMax]);
% % else
% %     xlim([-0.1, max(sampleData)]);
% % end George's Code
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fit by maximum matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histH = figure;
% hist(sampleData, actualDistX);
% hold on;
% plot(HR_actualDistX, HR_actualDistY  / max(HR_actualDistY) * max(histN), 'r-');
% for i = 1 : clusterNum
%     plot(separatedCluster(i).gridX, separatedCluster(i).gridY  / max(actualDistY) * max(histN), 'c-');
% end
% 
% % George's Code (can be removed)
% title('Fit by maximum-matching');
% xlim([-0.1, max(sampleData)]);
% 
% [GECN,GECX] = hist(sampleData, actualDistX);
% histH = figure;
% bar(GECX,100*GECN./sum(GECN),1);
% ylabel('Percentage (%)');
% hold on;
% plot(HR_actualDistX, (HR_actualDistY  / max(HR_actualDistY) * max(histN))/sum(GECN)*100, 'r-');
% for i = 1 : clusterNum
%     plot(separatedCluster(i).gridX, (separatedCluster(i).gridY  / max(actualDistY) * max(histN))/sum(GECN)*100, 'c-');
% end
% 
% 
% title('Fit by maximum-matching (%)');
% %end George's Code
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fit by sum matching
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histH = figure;
% hist(sampleData, actualDistX);
% hold on;
% plot(HR_actualDistX, HR_actualDistY  * sum(histN) / sum(actualDistY), 'r-', 'LineWidth', 1);
% for i = 1 : clusterNum
%     plot(separatedCluster(i).gridX, separatedCluster(i).gridY  * sum(histN) / sum(actualDistY), 'c-', 'LineWidth', 1);
% end
% title('Fit by sum-matching');
% 
% % if strcmp(ButtonName, 'Yes')
% %     xlim([-0.1, histogramMax]);
% % else
% %     xlim([-0.1, max(sampleData)]);
% % end
% 
% % George's Code (can be removed)
% [GECN,GECX] = hist(sampleData, actualDistX);
% histH = figure;
% bar(GECX,100*GECN./sum(GECN),1);
% ylabel('Percentage (%)');
% hold on;
% plot(HR_actualDistX, (HR_actualDistY  * sum(histN) / sum(actualDistY))/sum(GECN)*100, 'r-', 'LineWidth', 1);
% for i = 1 : clusterNum
%     plot(separatedCluster(i).gridX, (separatedCluster(i).gridY  * sum(histN) / sum(actualDistY))/sum(GECN)*100, 'c-', 'LineWidth', 1);
% end
% title('Fit by sum-matching (%)');
% 
% % if strcmp(ButtonName, 'Yes')
% %     xlim([-0.1, histogramMax]);
% % else
% %     xlim([-0.1, max(sampleData)]);
% % end George's Code
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 9: print a summary of the clustering results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('\nMultiple Normal Distribution Fitting for:\n%s\n', dataFileName);
% for i = 1 : clusterNum
%     fprintf('----------------------------- Cluster %d -----------------------------\n', i);
%     fprintf('MEAN = %f, STD =%f, WEIGHT =%f\n',  mu(i), sqrt(sigmasq(i)), prob(i));
% end
% fprintf('------------------------------------------------------------------------\n\n', i);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 10: save the result into a MAT file
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % We need a MAT file to save the results
% prompt={'Please proivide a tag for this analysis'};
% name='Save cluster analysis results';
% numlines=1;
% defaultanswer={'exp'};
% answer=inputdlg(prompt,name,numlines,defaultanswer);
%  
% idx1 = find(dataFileName == '.');
% idx2 = find(dataFileName == filesep);
% dataFileNameBody = dataFileName(idx2(end) + 1 : idx1(end) - 1);
% analysisResultFile = [dataFileName(1 : idx2(end)), char(answer{1}), '_clustAnaRecord.mat'];
% 
% save(analysisResultFile, 'sampleData', 'mu', 'sigma', 'prob', 'numGauss', 'BIC', 'classification');







    