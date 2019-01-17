function [MSE] = ARX_error(arxOrder, trainData,testData)

order = [max(round(arxOrder), 0)];
opt = arxOptions('Focus','simulation');
trainARX = arx(trainData, order, opt);
testARX = sim(trainARX, testData(:, 2:8));

cost_func = 'MSE';
MSE = goodnessOfFit(testARX,testData(:,1),cost_func);

end

