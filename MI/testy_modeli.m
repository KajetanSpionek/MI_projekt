clear all;
%
load("testData.mat");
load("M_reduced.mat");
%modele zrobione w regressionLearner - wpisac command window: regressionLearner
load("baggedtree.mat");
load("model_neuronowy.mat");
load("Rgaussianqudratic.mat");
load("SVMgaussian.mat");
load("finetree.mat");

error=0;
pocz=1000;
kon=1199;
%%
%siec neuronowa
% test na fragmencie zbioru testowego

% Y_sim = net(testData(:,2:8)');
% for i=pocz:1:kon
%     error=error+(testData(i,1)-Y_sim(i))^2;
% end
% figure(20); hold on;
% title(error);
% plot(testData(pocz:kon,1)');
% plot(Y_sim(pocz:kon));
% legend('process', 'simulation');
% 
% %test na calosci zbioru 
% figure(21); hold on;
% Y_sim = net(M_reduced(1:9999,1:7)');
% plot(M_reduced(1:9999,8)');
% plot(Y_sim);

%%
%inne modele - nauczone na 80 procent danych
%zamienic tylko RGaussianQuadratic na:
% baggedforest
% finetree
% SVMGaussian

Y_sim=RGaussianQuadratic.predictFcn(testData(:,2:8));

for i=pocz:1:kon
    error=error+(testData(i,1)-Y_sim(i))^2;
end
figure(20); hold on;
title(error);
plot(testData(pocz:kon,1));
plot(Y_sim(pocz:kon));
%%
%%inne modele - nauczone na 50 procent danych
% load("testData2.mat");
% load("testData_inne.mat");
% load("baggedtree2.mat");
% load("Rgaussianqudratic2.mat");
% load("SVMgaussian2.mat");
% load("finetree2.mat");
% 
% %zamienic tylko RQuadratic2 na:
% % RandomForest2
% % FineTree2
% % SVMFineGaussian2
% 
% Y_sim=RQuadratic2.predictFcn(testData_inne(:,2:8));
% 
% for i=pocz:1:kon
%     error=error+(testData_inne(i,1)-Y_sim(i))^2;
% end
% figure(20); hold on;
% title(error);
% plot(testData_inne(pocz:kon,1));
% plot(Y_sim(pocz:kon));
