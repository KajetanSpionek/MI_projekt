clear variables;
%wczytanie danych
M = xlsread('piec_dane.xlsx',3);
M(:,23)=M(:,24);
%wypelnienie dziur z interpolacja
M=fillmissing(M,'spline');

n = 50; % median filter smoothing parameter
time=1:10080;
time=time';
%%
%wrzucam dane z macierzy do wektora i normalizuje
%nastepnie wygladzam szumy
%i robie wykresy porownawcze sygnalu przed i po wygladzeniu
%i tak 23 razy
CMill_Flow = featureNormalize((M(:,1)));
CMill_Flow_filtered=medfilt1(CMill_Flow,20);
% figure(1);
% plot(time,CMill_Flow);hold on;
% plot(time,CMill_Flow_filtered)

CMill_Power = featureNormalize(M(:,2));
CMill_Power_filtered=medfilt1(CMill_Power,50);
% figure(2);
% plot(time,CMill_Power);hold on;
% plot(time,CMill_Power_filtered)

CMill_PAir_Flow = featureNormalize(M(:,3));
CMill_PAir_filtered=medfilt1(CMill_PAir_Flow,50);
% figure(3);
% plot(time,CMill_PAir_Flow);hold on;
% plot(time,CMill_PAir_filtered)

CMill_PPreasure=featureNormalize(M(:,4));
CMill_PPreasure_filtered=medfilt1(CMill_PPreasure,50);
% figure(4);
% plot(time,CMill_PPreasure);hold on;
% plot(time,CMill_PPreasure_filtered)

CMill_PAir_Temp=featureNormalize(M(:,5));
CMill_PAir_Temp_filtered=medfilt1(CMill_PAir_Temp,50);
% figure(5);
% plot(time,CMill_PAir_Temp);hold on;
% plot(time,CMill_PAir_Temp_filtered)

CMill_PAir_Preasure=featureNormalize(M(:,6));
CMill_PAir_Preasure_filtered=medfilt1(CMill_PAir_Preasure,50);
% figure(6);
% plot(time,CMill_PAir_Preasure);hold on;
% plot(time,CMill_PAir_Preasure_filtered)

CMill_PAir_HOT=featureNormalize(M(:,7));
CMill_PAir_HOT_filtered=medfilt1(CMill_PAir_HOT,50);
% figure(7);
% plot(time,CMill_PAir_HOT);hold on;
% plot(time,CMill_PAir_HOT_filtered)

CMill_PAir_COLD=featureNormalize(M(:,8));
CMill_PAir_COLD_filtered=medfilt1(CMill_PAir_COLD,75);
% figure(8);
% plot(time,CMill_PAir_COLD);hold on;
% plot(time,CMill_PAir_COLD_filtered)

Main_Steam_Flow= featureNormalize((M(:,9)));
Main_Steam_Flow_filtered=medfilt1(Main_Steam_Flow,75);
% figure(9);
% plot(time,Main_Steam_Flow);hold on;
% plot(time,Main_Steam_Flow_filtered)

Fan1_Power= featureNormalize((M(:,10)));
Fan1_Power_filtered=medfilt1(Fan1_Power,70);
% figure(10);
% plot(time,Fan1_Power);hold on;
% plot(time,Fan1_Power_filtered)

Fan2_Power= featureNormalize((M(:,11)));
Fan2_Power_filtered=medfilt1(Fan2_Power,50);
% figure(11);
% plot(time,Fan2_Power);hold on;
% plot(time,Fan2_Power_filtered)

Fan1_Preassure= featureNormalize((M(:,12)));
Fan_Preassure_filtered=medfilt1(Fan1_Preassure,50);
% figure(12);
% plot(time,Fan1_Preassure);hold on;
% plot(time,Fan_Preassure_filtered)

Fan2_Preasure= featureNormalize((M(:,13)));
Fan2_Preasure_filtered=medfilt1(Fan2_Preasure,50);
% figure(13);
% plot(time,Fan2_Preasure);hold on;
% plot(time,Fan2_Preasure_filtered)

Total_Air_Flow= featureNormalize((M(:,14)));
Total_Air_Flow_filtered=medfilt1(Total_Air_Flow,70);
% figure(14);
% plot(time,Total_Air_Flow);hold on;
% plot(time,Total_Air_Flow_filtered)

Air_before_fan1= featureNormalize((M(:,15)));
Air_before_fan1_filtered=medfilt1(Air_before_fan1,70);
% figure(15);
% plot(time,Air_before_fan1);hold on;
% plot(time,Air_before_fan1_filtered)


Air_before_fan2= featureNormalize((M(:,16)));
Air_before_fan2_filtered=medfilt1(Air_before_fan2,50);
% figure(16);
% plot(time,Air_before_fan2);hold on;
% plot(time,Air_before_fan2_filtered)

Preasure_after_fan1= featureNormalize((M(:,17)));
Preasure_after_fan1_filtered=medfilt1(Preasure_after_fan1,50);
% figure(17);
% plot(time,Preasure_after_fan1);hold on;
% plot(time,Preasure_after_fan1_filtered)

Preasure_after_fan2= featureNormalize((M(:,18)));
Preasure_after_fan2_filtered=medfilt1(Preasure_after_fan2,50);
% figure(18);
% plot(time,Preasure_after_fan2);hold on;
% plot(time,Preasure_after_fan2_filtered)

Air_before_preheater1= featureNormalize((M(:,19)));
Air_before_preheater1_filtered=medfilt1(Air_before_preheater1,10);
% figure(19);
% plot(time,Air_before_preheater1);hold on;
% plot(time,Air_before_preheater1_filtered)

Air_before_preheater2= featureNormalize((M(:,20)));
Air_before_preheater2_filtered=medfilt1(Air_before_preheater2,10);
% figure(20);
% plot(time,Air_before_preheater2);hold on;
% plot(time,Air_before_preheater2_filtered)

Air_after_preheater1= featureNormalize((M(:,21)));
Air_after_preheater1_filtered=medfilt1(Air_after_preheater1,10);
% figure(21);
% plot(time,Air_after_preheater1);hold on;
% plot(time,Air_after_preheater1_filtered)

Air_after_preheater2= featureNormalize((M(:,22)));
Air_after_preheater2_filtered=medfilt1(Air_after_preheater2,50);
% figure(22);
% plot(time,Air_after_preheater2);hold on;
% plot(time,Air_after_preheater2_filtered)

CMill_Primary_coaltemp= featureNormalize((M(:,23)));
CMill_Primary_coaltemp_filtered=medfilt1(CMill_Primary_coaltemp,80);
% figure(23);
% plot(time,CMill_Primary_coaltemp);hold on;
% plot(time,CMill_Primary_coaltemp_filtered)
%%

%teraz wrzucam te wygladzone sygnaly z powrotem do macierzy
M_filtered=[CMill_Flow_filtered,CMill_Power_filtered,CMill_PAir_filtered,CMill_PPreasure_filtered,CMill_PAir_Temp_filtered,CMill_PAir_Preasure_filtered,CMill_PAir_HOT_filtered,CMill_PAir_COLD_filtered,Main_Steam_Flow_filtered,Fan1_Power_filtered,Fan2_Power_filtered,Fan_Preassure_filtered,Fan2_Preasure_filtered,Total_Air_Flow_filtered,Air_before_fan1_filtered,Air_before_fan2_filtered,Preasure_after_fan1_filtered,Preasure_after_fan2_filtered,Air_before_preheater1_filtered,Air_before_preheater2_filtered,Air_after_preheater1_filtered,Air_after_preheater2_filtered,CMill_Primary_coaltemp_filtered];

%tu sobie przygotowuje tablice do zapisania jak sa skorelowane sygnaly
rho_Pearson=zeros(23,23);
p_value_Pearson=zeros(23,23);
%w tej petli sprawdzam, jak kazdy sygnal jest skorelowany z kazdym innym
%najpierw liniowo - korelacja Pearsona
for i=1:1:23
    for j= 1:1:23
        [rho_Pearson(i,j),p_value_Pearson(i,j)] = corr(M_filtered(:,i),M_filtered(:,j),'Type','Pearson');
    end
end
%potem jeszcze korelacja Spearmana, ona moze wykryc zaleznosci nieliniowe
rho_Spearman=zeros(23,23);
p_value_Spearman=zeros(23,23);
for i=1:1:23
    for j= 1:1:23
        [rho_Spearman(i,j),p_value_Spearman(i,j)] = corr(M_filtered(:,i),M_filtered(:,j),'Type','Spearman');
    end
end

%nastepnie odsylam do excela, tam te tabelki ladnie wygladaja,
%zakolorowalem na czerwono duze korelacje, a na zielono te mniej
%skorelowane

%po przyjzeniu sie tym tabelom w excelu stwierdzam, ze z powodzeniem mozna
%wyrzucic wiekszosc wejsc, po sa ze zoba w duzym stopniu skorelowane. Wiec
%bylyby kilkukrotnie te same informacje wprowadzane do modelu.

%zostaja wejscia 1,2,7,8,11,20,21 Robie nowa macierz wejsc tylko z nimi
%stary index                                 nowy index
%1           CMill_Flow_filtered             1
%2           CMill_Power_filtered            2
%7           CMill_PAir_HOT_filtered         3
%8           CMill_PAir_COLD_filtered        4
%11          Fan2_Power_filtered             5
%20          Air_before_preheater2_filtered  6
%21          Air_after_preheater1_filtered   7
M_reduced=[M_filtered(:,1),M_filtered(:,2),M_filtered(:,7),M_filtered(:,8),M_filtered(:,11),M_filtered(:,20),M_filtered(:,21),M_filtered(:,23)];
M_reduced=M_reduced(randperm(size(M_reduced, 1)), :); %% wymieszanie
% U_train=M_reduced(1:8000,:);
% U_test=M_reduced(8001:10080,:);
% 
% %Wyjscie nazywam inaczej, zeby sie juz nie mylilo
% Y_train=M_filtered(1:8000,23);
% Y_test=M_filtered(8001:10080,23);
% 
% trainData = [Y_train U_train];
% testData = [Y_test U_test];
% end
U_train=zeros(4*10080/5,7);
U_test=zeros(10080/5,7);
Y_train=zeros(4*10080/5,1);
Y_test=zeros(10080/5,1);
k=1;
l=1;
for i=1:1:10080
    if mod(i,5)==0
    U_test(l,:)=M_reduced(i,1:7);
    Y_test(l,1)=M_reduced(i,8);
    l=l+1;
    else
    U_train(k,:)=M_reduced(i,1:7);
    Y_train(k,1)=M_reduced(i,8);
    k=k+1;
    end
end

% U_train=zeros(10080/2,7);
% U_test=zeros(10080/2,7);
% Y_train=zeros(10080/2,1);
% Y_test=zeros(10080/2/5,1);
% k=1;
% l=1;
% for i=1:1:10080
%     if mod(i,2)==0
%     U_test(l,:)=M_reduced(i,1:7);
%     Y_test(l,1)=M_reduced(i,8);
%     l=l+1;
%     else
%     U_train(k,:)=M_reduced(i,1:7);
%     Y_train(k,1)=M_reduced(i,8);
%     k=k+1;
%     end
% end



trainData = [Y_train U_train];
testData = [Y_test U_test];
%%2 potega
trainData_nlin=zeros(size(trainData,1),size(trainData,2)*2-1);
testData_nlin=zeros(size(testData,1),size(testData,2)*2-1);

trainData_nlin(:,1)=trainData(:,1);
testData_nlin(:,1)=testData(:,1);
for j=2:1:size(trainData,2)    
            trainData_nlin(:,2*j-2)=trainData(:,j); 
            trainData_nlin(:,2*j-1)=trainData(:,j).^2;
            testData_nlin(:,2*j-2)=testData(:,j); 
            testData_nlin(:,2*j-1)=testData(:,j).^2;
end



%%3 potega
% trainData_nlin=zeros(size(trainData,1),size(trainData,2)*3-2);
% testData_nlin=zeros(size(testData,1),size(testData,2)*3-2);
% 
% trainData_nlin(:,1)=trainData(:,1);
% testData_nlin(:,1)=testData(:,1);
% for j=2:1:size(trainData,2)    
%             trainData_nlin(:,3*j-4)=trainData(:,j); 
%             trainData_nlin(:,3*j-3)=trainData(:,j).^2;
%             trainData_nlin(:,3*j-2)=trainData(:,j).^3;
%             testData_nlin(:,3*j-4)=testData(:,j); 
%             testData_nlin(:,3*j-3)=testData(:,j).^2;
%             testData_nlin(:,3*j-2)=testData(:,j).^3;
% 
% end

%%4 potega
% trainData_nlin=zeros(size(trainData,1),size(trainData,2)*4-3);
% testData_nlin=zeros(size(testData,1),size(testData,2)*4-3);
% trainData_nlin(:,1)=trainData(:,1);
% testData_nlin(:,1)=testData(:,1);
% for j=2:1:size(trainData,2)    
%             trainData_nlin(:,4*j-6)=trainData(:,j); 
%             trainData_nlin(:,4*j-5)=trainData(:,j).^2;
%             trainData_nlin(:,4*j-4)=trainData(:,j).^3;
%             trainData_nlin(:,4*j-3)=trainData(:,j).^4;
%             testData_nlin(:,4*j-6)=testData(:,j); 
%             testData_nlin(:,4*j-5)=testData(:,j).^2;
%             testData_nlin(:,4*j-4)=testData(:,j).^3;
%             testData_nlin(:,4*j-3)=testData(:,j).^4;
% end
%%
%sprawdzam funkcje autokorelacji i korelacji czastkowej dla wyjscia
%select chanel
% ch = CMill_Primary_coaltemp_filtered;
% figure(1)
% autocorr(ch,200);
% figure(2)
% parcorr(ch,50);
%jesli dobrze rozumiem jak to dziala, to wtedy najwiekszy wplyw na wyjscie
%ma wartosc z chwili (t-1). (t-2),(t-3) i (t-4) jeszcze cos robia, a dalej
%juz ten wplyw jest mizerny (z wykresu czastkowej). Wiec na=4

%%identyfikacja opoznienia - pod spodem jakies metody skopiowane, ni chuja
%%nie wiem o co biega xd
%% method 1
% cross coorelation of siganlas
% y = M_filtered(:,23); %to zostawiamy
% u = M_reduced(:,7); %tu mozna zmieniac na inne wejscia
% N = length(y);
% lambda_est = sum(u.^2)/N;
% crossCorr = conv(u,y);
% figure(1);
% plot(crossCorr);
% title('odpowiedz impulsowa');
% xlabel('probka');
% ylabel('amplituda');

%% method 2
% estimate impulse respone 
% Ryu_est = crossCorr/N;
% g_est = Ryu_est/lambda_est;

% figure(2)
% plot(g_est);
% title('Impulse respone g(t)')

% Fft analysy

% sysFft = real(fft(g_est));
% figure(3)
% plot(sysFft)
% title('Real part of impulse respone G_h')

% logOfAbsOfFftIR = log(abs(sysFft));
% figure(4)
% plot(logOfAbsOfFftIR)
% title('Log of Abs value of Real part of FFT of impulse respone')
% 
% invFftSys = ifft(logOfAbsOfFftIR);
% figure(5)
% plot(invFftSys(1:20));
% title('realstrum')
% 
% %% method 4
% % mathlab methods
% impResp = impulseest(iddata(y,u))
% delay_est = delayest(iddata(y,u))
% figure(7);
% plot(impResp);

%% MNK
% nr=1;
% na=0; %to chyba mamy
% 
% nb=[26 21 18 14 22 2 5]; %trzeba wymyslic, jak dobrac
% nk=[23 8 15 20 27 5 2];%trzeba wymyslic, jak dobrac
% 
% % nb=[26 26 21 21 18 18 14 14 22 22 2 2 5 5];
% % nk=[23 23 8 8 15 15 20 20 27 27 5 5 2 2];
% 
% % nb=[26 26 26 21 21 21 18 18 18 14 14 14 22 22 22 2 2 2 5 5 5];
% % nk=[23 23 23 8 8 8 15 15 15 20 20 20 27 27 27 5 5 5 2 2 2];
% % 
% % nb=[26 26 26 26 21 21 21 21 18 18 18 18 14 14 14 14 22 22 22 22 2 2 2 2 5 5 5 5];
% % nk=[23 23 23 23 8 8 8 8 15 15 15 15 20 20 20 20 27 27 27 27 5 5 5 5 2 2 2 2];
% opt = arxOptions('Focus','simulation');
% %ARX MNK
% ARX_MNK_mod = arx(trainData,[na nb nk], opt);
% % ARX_MNK_mod = oe(trainData_nlin,[na nb nk]);
% figure(nr)
% compare(testData, ARX_MNK_mod)
% % compare(testData_nlin, ARX_MNK_mod);
% figure(nr+1)
% compare(trainData, ARX_MNK_mod);
% % compare(trainData_nlin, ARX_MNK_mod);
% % ARX PEM
% ARX_PEM_mod = pem(trainData, ARX_MNK_mod,'Focus','simulation');
% % ARX_PEM_mod = pem(trainData_nlin, ARX_MNK_mod,'Focus','simulation');
% figure(nr+2)
% compare(testData, ARX_PEM_mod);
% % compare(testData_nlin, ARX_PEM_mod);
% nr=nr+3;


%%
%% select ARX order with PSO

%parametry dobrane przez pso
%na=0    nb=[26    21    18    14    22     2     5] 
%nk= [ 23     8    15 20    27     5     2]

% addpath('psopt');
% % len = length(na) + length(nb);
% % len = length(nk);
% len = 1 + 7 + 7;
% % fun = @(x) ARX_error(x, na, nb, trainData,testData);
% % fun = @(x) ARX_error(x, nk, trainData,testData);
% fun = @(x) ARX_error(x, trainData,testData);
% 
% startVec = [0; 30];
% startMatrix = repmat(startVec, [1, len]);
% opts = psooptimset('Generations', 40,'PopInitRange', startMatrix, 'PopulationSize', 20, ...
%                    'SocialAttraction', 1.25, 'CognitiveAttraction', 0.5, 'StallGenLimit', 150);
% [out, val] = pso(fun, len, [], [], [], [], startMatrix(1, :), startMatrix(2, :), [], opts);
% disp(out);
% disp(val);
% 
% params = round(out);
%%
%%
%%siec neuronowa
% czas symulacji 
simStart = 1;
tsim = 50000;

%dane wejsciowe
inputs = U_train';
targets = Y_train';

%utworzenie sieci dopasowujacej
hiddenLayerSize = [10 5 2];
net = fitnet(hiddenLayerSize,'trainlm');

%podzial danych na uczace, testujace i walidujace
net.divideParam.trainRatio = 65/100;
net.divideParam.testRatio = 20/100;
net.divideParam.valRatio = 15/100;

%uczenie sieci
[net, tr] = train(net, inputs, targets);

%testowanie sieci
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);

%wizualizacja sieci
view(net) 


% wykresy
figure, plotperform(tr)
figure, plottrainstate(tr)
figure, plot(targets)
hold on
plot(outputs)
legend('process', 'simulation');
figure, plotregression(targets,outputs)
figure, ploterrhist(errors)

% %test on test data
pocz=1;
kon=150;
random_testData=testData(randperm(size(testData, 1)), :);
load('model_neuronowy.mat')
Y_sim = net(random_testData(:,2:8)');
figure(20); hold on;
plot(random_testData(pocz:kon,1)');
plot(Y_sim(pocz:kon));
legend('process', 'simulation');
%%
figure(21); hold on;
Y_sim = net([M_filtered(:,1),M_filtered(:,2),M_filtered(:,7),M_filtered(:,8),M_filtered(:,11),M_filtered(:,20),M_filtered(:,21)]');
plot(M_filtered(:,23)');
plot(Y_sim);

