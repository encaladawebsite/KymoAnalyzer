function normality_test_with_R(data)

sampleData = data;



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

[result, runStatus] = evalR('library(nortest)');
% % put observations into R
putRdata('inputData', sampleData);

evalR('testResult <- ad.test(inputData)');
testResult = evalR('testResult');
fprintf('--------------------- Anderson-Darling Test -----------------------\n');
fprintf('p-value is = %s\n', testResult{2});
fprintf('-------------------------------------------------------------------\n');

evalR('testResult <- lillie.test(inputData)');
testResult = evalR('testResult');
fprintf('------------------------- Lilliefors Test -------------------------\n');
fprintf('p-value is = %s\n', testResult{2});
fprintf('-------------------------------------------------------------------\n');


% evalR('testResult <- sf.test(inputData)');
% testResult = evalR('testResult');
% fprintf('-------------------- Shapiro-Francia Test -------------------------\n');
% fprintf('p-value is = %s\n', testResult{2});
% fprintf('-------------------------------------------------------------------\n');

c = 1;

 

    