clear;clc;close all
load('Data_clean_processed.mat')

%% MAKE A BIG TIME-TABLE DATA

% First experiment: Smoking
% Instruments involved: 
% AeroTrak (PN,T,RH), CO_met (T,RH,P), CPC (PN?)
% Ptrak (), DustTrak (), SidePak ()
% PND, PMSD
% LCS_G1, LCS_G2_01, LCS_G2_02


% AeroTrak_Met
T1 = datetime(AeroTrak_met(:,1:6));
AT_DoY1 = AeroTrak_met(:,7);
AT_T = AeroTrak_met(:,8);
AT_RH = AeroTrak_met(:,9);
AT_met = timetable(T1,AT_DoY1,AT_T,AT_RH);
clear AT_DoY1 AT_T AT_RH 

% AeroTrak_PN
T1 = datetime(AeroTrak_PN(:,1:6));
AT_DoY2 = AeroTrak_PN(:,7);
AT_0p3_0p5   = AeroTrak_PN(:,8);
AT_0p5_1p0   = AeroTrak_PN(:,9);
AT_1p0_2p5   = AeroTrak_PN(:,10);
AT_2p5_5p0   = AeroTrak_PN(:,11);
AT_5p0_10p0  = AeroTrak_PN(:,12);
AT_10p0_45p0 = AeroTrak_PN(:,13);
AT_PN = timetable(T1,AT_DoY2,AT_0p3_0p5, ...
    AT_0p5_1p0,AT_1p0_2p5,AT_2p5_5p0,AT_5p0_10p0,AT_10p0_45p0);
clear AT_0p3_0p5 AT_0p5_1p0 AT_1p0_2p5 ...
    AT_2p5_5p0 AT_5p0_10p0 AT_10p0_45p0

% ClasOhlson_met
T1 = datetime(ClasOhlson_met(:,1:6));
CO_DoY = ClasOhlson_met(:,7);
CO_T   = ClasOhlson_met(:,8);
CO_RH  = ClasOhlson_met(:,9);
CO_P   = ClasOhlson_met(:,10);
CO_met = timetable(T1,CO_DoY,CO_T,CO_RH,CO_P);
clear T1 CO_DoY CO_T CO_RH CO_P

% CPC
T1 = datetime(CPC(:,1:6));
CPC_DoY = CPC(:,7);
CPC_c   = CPC(:,8);
CPC1 = timetable(T1,CPC_DoY,CPC_c);
clear T1 CPC_DoY CPC_c

% Ptrak
T1 = datetime(Ptrak(:,1:6));
Ptrak_DoY = Ptrak(:,7);
Ptrak_c   = Ptrak(:,8);
Ptrak1 = timetable(T1,Ptrak_DoY,Ptrak_c);
clear T1 Ptrak_DoY Ptrak_c

% DustTrak
T1 = datetime(DustTrak(:,1:6));
DustTrak_DoY = DustTrak(:,7);
DustTrak_c   = DustTrak(:,8:12);
DustTrak1 = timetable(T1,DustTrak_DoY,DustTrak_c);
clear T1 DustTrak_DoY DustTrak_c

% SidePak
T1 = datetime(SidePak(:,1:6));
SidePak_DoY = SidePak(:,7);
SidePak_c   = SidePak(:,8);
SidePak1 = timetable(T1,SidePak_DoY,SidePak_c);
clear T1 DustTrak_DoY DustTrak_c


% PND
T1=datetime(PND(:,1:6));
PND_DoY = PND(:,7);
PND_c   = PND(:,8:end);
PND1 = timetable(T1,PND_DoY,PND_c);
PND2 = sortrows(PND1);
PND3 = PND2(1:end-1,:);
%PND2 = unique(PND1);
clear T1 PND_DoY PND_c

% PMSD
T1=datetime(PMSD(:,1:6));
PMSD_DoY = PMSD(:,7);
PMSD_c   = PMSD(:,8:end);
PMSD1 = timetable(T1,PMSD_DoY,PMSD_c);
PMSD2 = sortrows(PMSD1);
PMSD3 = PMSD2(1:end-1,:);
%PMSD2 = unique(PMSD1);
clear T1 PMSD_DoY PMSD_c

% LCS_G1:
T1 = datetime(ISEE_LCS_G1(:,1:6));
LCS_G1 = ISEE_LCS_G1(:,8);
LCS_G1_T = timetable(T1,LCS_G1);
clear T1 LCS_G1

% LCS_G2_01 (PM2.5 and MET)
T1 = datetime(ISEE_LCS_G201(:,1:6));
LCS_G2_01 = ISEE_LCS_G201(:,8); % PM2.5
LCS_G2_01_met = ISEE_LCS_G201_met(:,8:end); % RH, T and P
LCS_G2_01_T = timetable(T1,LCS_G2_01,LCS_G2_01_met);
clear T1 LCS_G2_01 LCS_G2_01_met



% LCS_G2_02 (PM2.5 and MET)
T1 = datetime(ISEE_LCS_G202(:,1:6));
LCS_G2_02 = ISEE_LCS_G202(:,8); % PM2.5
LCS_G2_02_met = ISEE_LCS_G202_met(:,8:end); % RH, T and P
LCS_G2_02_T = timetable(T1,LCS_G2_02,LCS_G2_02_met);
clear T1 LCS_G2_02 LCS_G2_02_met

%DATA_ts = synchronize(AT_met,AT_PN,CO_met,CPC1,Ptrak1,PND1,PMSD1, ...
%    LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));

% https://www.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html

DATA = synchronize(AT_met, AT_PN, CO_met, CPC1, Ptrak1, DustTrak1, SidePak1, ...
    PND3, PMSD3, ... % PND1, PMSD1,...
    LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));


             
%% Met variables\
close all

figure(1); fig=gcf;
plot(DATA.T1,DATA.AT_T,'b','LineWidth',2);hold on;
plot(DATA.T1,DATA.CO_T,'r','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_01_met(:,2),'y','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_02_met(:,2),'c','LineWidth',2);
title('Temperature')
legend('AeroTrak','ClasOhlson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([20 50]); grid on
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

figure(2); fig=gcf;
plot(DATA.T1,DATA.AT_RH,'b','LineWidth',2);hold on;
plot(DATA.T1,DATA.CO_RH,'r','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_01_met(:,1),'y','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_02_met(:,1),'c','LineWidth',2);
title('Relative Humidity')
legend('AeroTrak','ClasOhlson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([20 70]); grid on
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

figure(3); fig=gcf;
%plot(DATA.T1,DATA.AT_P,'b','LineWidth',2);hold on;
plot(DATA.T1,DATA.CO_P,'r','LineWidth',2);hold on;
plot(DATA.T1,DATA.LCS_G2_01_met(:,3),'y','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_02_met(:,3),'c','LineWidth',2);
title('Pressure')
legend('ClasOhlson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([895 905]); grid on
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

%% Aerosol Concentration
figure(4); fig=gcf;
plot(DATA.T1,DATA.AT_1p0_2p5,'b-','LineWidth',2);hold on;
plot(DATA.T1,DATA.AT_2p5_5p0,'b*','LineWidth',2);
plot(DATA.T1,DATA.CPC_c,'r','LineWidth',2);
plot(DATA.T1,DATA.Ptrak_c,'g','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_01(:,1),'y','LineWidth',2);
plot(DATA.T1,DATA.LCS_G2_02(:,1),'c','LineWidth',2);
title('Aerosol particles')
ylabel('PM$_{2.5}$ [$\mu g / m^3$]','interpreter','latex')
xlabel('Time')
legend('AT (PM$_1$ - PM${2.5}$)','AT (PM${2.5}$ - PM$_1$)', ...
    'CPC ($>10nm$)','Ptrak ($>25nm$)', ...
    '$\mathcal{L}_{2a}$ (PM${2.5}$)','$\mathcal{L}_{2b}$ (PM${2.5}$)','interpreter','latex')
ylim([10e-3 10e5]); grid on
set(gca, 'YScale', 'log')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

%% Aerosol concentration (using banana plot)
clc
AT = AeroTrak_PN(:,8:13);
figure(5); fig=gcf;
%pcolor(DATA.PND_c')
%pcolor(DATA.PND_c)
%pcolor(DATA.PMSD_c)
%imagesc(DATA.PND_c)
%image(DATA.PMSD_c')
subplot(411)
imagesc(AT');hold on
title('AeroTrak')
set(gca, 'YScale', 'log')
colorbar
hold off
subplot(412)
plot(DATA.LCS_G1(:,1),'b','LineWidth',2); hold on
plot(DATA.LCS_G2_01(:,1),'y','LineWidth',2); 
plot(DATA.LCS_G2_02(:,1),'c','LineWidth',2);
title('Low-Cost Sensors PM$_{2.5}$','interpreter','latex');grid on
legend('$\mathcal{L}_1$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
set(gca, 'YScale', 'log')
ylim([1e-2 1e4])
%colorbar
hold off
subplot(413)
imagesc(DATA.PND_c');hold on
%C = repmat(colorscale,10,1); pcolor(DATA.T1,DATA.PND_c); hold on
title('DMPS 1')
set(gca, 'YScale', 'log')
colorbar
hold off
subplot(414)
image(DATA.PMSD_c'); hold on
title('DMPS 2');grid on
set(gca, 'YScale', 'log')
colorbar
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

%%

inst_res = linspace(1,9,9);
%inst_res_AT = linspace(1,6,6);
%inst_res_AT = linspace(1,6,6);
inst_res_AT = [0.5,1,2.5,5,10,45];


clc
figure(6);fig=gcf;
subplot(211)
%pcolor(DATA.PND_c')
%pcolor(DATA.PND_DoY(:,1)',inst_res',DATA.PND_c')
h = pcolor(AeroTrak_PN(:,7)',inst_res_AT',AeroTrak_PN(:,8:13)');
set(h, 'EdgeColor', 'none')
xlabel('Day of the year','interpreter','latex')
ylabel('Particle size ($\mu$m)','interpreter','latex')
title('AeroTrak','interpreter','latex')
%yticks(gca,inst_res_AT)
%yticks(gca,[1 2 3 4 5 6 7])
%yticklabels(gca,{'0.3 \mum','0.5 \mum','1 \mum','2.5 \mum','5 \mum','10 \mum','45 \mum','interpreter','latex'})
colorbar
colormap jet
%caxis([50 10^5])
set(gca, 'YScale', 'log','colorscale','log')
subplot(212)
plot(DATA.AT_DoY1,DATA.LCS_G1(:,1),'b','LineWidth',2); hold on
plot(DATA.AT_DoY1,DATA.LCS_G2_01(:,1),'y','LineWidth',2); 
plot(DATA.AT_DoY1,DATA.LCS_G2_02(:,1),'c','LineWidth',2);
xlabel('Day of the year','interpreter','latex')
ylabel('PM$_{2.5}$ concentration','interpreter','latex')
title('Low-Cost Sensors PM$_{2.5}$','interpreter','latex');grid on
legend('$\mathcal{L}_1$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
set(gca, 'YScale', 'log')
ylim([1e-2 1e4])
colorbar
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

%% MATRIX PLOT
close all;clc

DATA_array2 = [AT(:,2), ...  % AT_PN.AT_1p0_2p5, ...
        CPC1.CPC_c ...
        Ptrak1.Ptrak_c, ...
        PND3.PND_c(:,4), ...
        DustTrak(:,8), ...
        SidePak(:,8), ...
        LCS_G1_T.LCS_G1, ...
        LCS_G2_01_T.LCS_G2_01, ...
        LCS_G2_02_T.LCS_G2_02, ...
        AT_met.AT_T, ...
        CO_met.CO_T, ...
        LCS_G2_01_T.LCS_G2_01_met(:,2), ...
        LCS_G2_02_T.LCS_G2_02_met(:,2), ...
        AT_met.AT_RH, ...
        CO_met.CO_RH, ...
        LCS_G2_01_T.LCS_G2_01_met(:,1), ...
        LCS_G2_02_T.LCS_G2_02_met(:,1), ...
        CO_met.CO_P, ...
        LCS_G2_01_T.LCS_G2_01_met(:,3), ...
        LCS_G2_02_T.LCS_G2_02_met(:,3)
        ];


labelX = {'PN$_{1-2.5}$ (AT)', ...
        'PN$_{1-2.5}$ (CPC)', ...
        'PN$_{1-2.5}$ (Ptrak)', ...
        'PN$_{1-2.5}$ (PND)', ...
        'PM$_{2.5}$ (DustTrak)', ...
        'PM$_{2.5}$ (SidePak)', ...
        'PM$_{2.5}$ ($\mathcal{L}_1$)', ...
        'PM$_{2.5}$ ($\mathcal{L}_{2a}$)', ...
        'PM$_{2.5}$ ($\mathcal{L}_{2b}$)', ...
        'Temp (AT)', ... 
        'Temp (CO)', ...
        'Temp ($\mathcal{L}_{2a}$)', ...
        'Temp ($\mathcal{L}_{2b}$)', ...
        'RH (AT)', ...
        'RH (CO)', ...
        'RH ($\mathcal{L}_{2a}$)', ...
        'RH ($\mathcal{L}_{2b}$)', ...
        'P (CO)', ...
        'P ($\mathcal{L}_{2a}$)' ...
        'P ($\mathcal{L}_{2b}$)'
    };

FS=26;
labelY = labelX;


%corrP = corrcoef(DataTot,'Rows','complete');
%CT = 'Pearson';
CT = 'Spearman';
%CT = 'Kendall';

corrP = corr(DATA_array2,'Type',CT,'Rows','pairwise');

clc
figure(2);fig=gcf;

if true
    Rms= abs(corrP);
    [m n]=size(Rms);
    imagesc(Rms,'CdataMapping','scaled');
    colorbar    
    hold on;
    for i = 1:m
        for j = 1:n
            nu = Rms(i,j);
            val = num2str(round(nu,2));
        end
    end
    hold off;
end
ax = gca;
%title('TSI and Low-Cost Sensors of PM_{2.5}', 'FontSize', FS);
% Set colorbar limits
caxis([0 1])
%colorbar;
h = colorbar;
%xlabel(h, 'R')
%colormap(flipud(copper))
%colormap(copper)
set(ax, 'XTick', 1:length(labelX), 'XTickLabel', labelX, 'TickLabelInterpreter', 'latex')
set(ax, 'YTick', 1:length(labelY), 'YTickLabel', labelY, 'TickLabelInterpreter', 'latex')
xtickangle(45)
%xtickangle(-180)
set(findall(fig,'-property','FontSize'),'FontSize',FS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% UPF proxy test:
% https://www.mathworks.com/help/deeplearning/gs/fit-data-with-a-neural-network.html
%X = [DATA_array2(:,[8,12,16,19])]';
X = [DATA_array2(:,[8,12])];
Y = AT(:,1);

tr = 1:roundn(0.7 * size(X,1),0);
te = tr(end)+1 : size(X,1);
Xtr = X(tr,:)';
Ytr = Y(tr,:)';
Xte = X(te,:)';
Yte = Y(te,:)';

% WE NEED TO NORMALIZE
 
% Create a Fitting Network
hiddenLayerSize = 15;
net = fitnet(hiddenLayerSize);

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
 
% Train the Network
[net,training] = train(net,Xtr,Ytr);
 
% Test the Network
%outputs = net(inputs);
Yest = net(Xte);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
 
figure(10); fig = gcf;
minV=0;
maxV=1000;
scatter(Yte,Yest);
xlabel(['PN [$\mu g/m^3$] (Aero Trak)'],'Interpreter','latex')
ylabel(['PN [$\mu g/m^3$] (Virtual Sensor)'],'Interpreter','latex')
hold on
xlim([minV maxV]); ylim([minV maxV])
xLine = linspace(minV,maxV);
yLine = linspace(minV,maxV);
plot(xLine,yLine,'r');
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);












%% KEROSENE heater experiment

clearvars -except DATA

%clear;close all;
load('Data_processed_Heaters/Data_processed_Kerosene_Heaters.mat')

% DATA: 
% AeroTrak (6 colmns), CPC, DustTrak, Ptrak, SidePak
% LCS_G1, LCS_G2_01, LCS_G2_02
% PMD, PMSD, PND, PNSD

figure(1); fig=gcf;
hold on
plot(CPC(:,1),CPC(:,2),'b.');
plot(Ptrak(:,1),Ptrak(:,2),'c.');
plot(DustTrak(:,1),DustTrak(:,2),'c.');
plot(SidePak(:,1),SidePak(:,2),'c*');
plot(ISEE_LCS_G1(:,1),ISEE_LCS_G1(:,2),'g-'); 
plot(ISEE_LCS_G201(:,1),ISEE_LCS_G201(:,2),'r-');
plot(ISEE_LCS_G202(:,1),ISEE_LCS_G202(:,2),'y-');
set(gca, 'YScale', 'log')
ylim([1e-2 1e6])
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);
%% NATURAL GAS heater experiment
%clear; close all;
clearvars -except DATA
load('Data_processed_Heaters/Data_processed_NaturalGas_Heaters.mat')

figure(1); fig=gcf;
hold on
plot(CPC(:,1),CPC(:,2),'b.');
plot(Ptrak(:,1),Ptrak(:,2),'c.');
%plot(DustTrak(:,1),DustTrak(:,2),'c.');
%plot(SidePak(:,1),SidePak(:,2),'c*');
plot(ISEE_LCS_G1(:,1),ISEE_LCS_G1(:,2),'g-'); 
plot(ISEE_LCS_G201(:,1),ISEE_LCS_G201(:,2),'r-');
plot(ISEE_LCS_G202(:,1),ISEE_LCS_G202(:,2),'y-');
set(gca, 'YScale', 'log')
ylim([1e-2 1e6])
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

%% QUESTIONS AND ANSWERS BY TAREQ:
% Q1 Ptrak has one type of aerosol size in your data, what is the particle size?
% A1: This one is monitored total particle number concentration for
% everything larger than 25 nm, CPC is the same as Ptrak, but 10 nm instead of 25 nm
%
% Q2 SidePak has also one type of aerosol size in your data, 
% what is the particle size? Is it PM.5 like the first smoking experiment?
% A2: Yes PM2.5
%
% Q3 DustTrak has five type of aerosol size in your data, what is the particle size?
% A3: This one is the same as sidepak but provides 5 classes: 
% PM1, PM2.5, PM4, PM10, TOTAL
%
% Q4: What about PND and PMSD, what instruments did you use for obtaining 
% these concentrations? and what are they for each column?
% A4: These were a combination from CPC (1 col), PTRAK (1), AEROTRAK (6)
%
% Q5: OK, did you place them in this order?
% AeroTrak (6 channels: 0.3-0.5, 0.5-1, 1- 2.5, 2.5 -5, 5-10, 10-45 )
% CPC (above 10nm)
% Ptrak (above 25nm)
% the total should be 8 columns, instead for PMD, PMSD, PND and PNSD, 
% you have 9 columns for the concentrations
%
% A5: Yes, 
% CPC - PTRAK == 10 - 25 nm
% Ptrak - sum(aerotrak) == 25 - 300 nm
% These are the first two bins
% Then comes the aerotrak bins
%
% Q6: what are the differences between three four variables: PMD, PMSD, PND
% and PNSD?
% A6: PND is just number concentration in each size bin
% PNSD is the normalized concentration to the width of each size bin
% dN/dlog(Dp)
% Same for PMD PMSD
%
% Q7: so, the raw data from all instruments are either PND or PMD, right?
% Q7: the instruments do not give normalized numbers?
%
% A7: The raw data is PN
% I generate PND and PNSD
% then calculate PMD PMSD

%% UNDERSTANDING THE DATA:
clc
n = 2020

sum(PND(n,2:end-1))
PND(n,end)

PND(n,2)
CPC(n,2)

%mean(PND(n,3:4))
PND(n,4)
Ptrak(n,2)

PND(n,2:10)
[CPC(n,2),Ptrak(n,2),AeroTrak_PN(n,2:7)]

%%
clc
n=2000;%10000;%2002;
disp(['CPC+Ptrak+AeroTrak: ',num2str([CPC(n,2),Ptrak(n,2),AeroTrak_PN(n,2:7)])])
disp(['PND: ',num2str(PND(n,2:10))])

%Ptrak(n,2) - PND(n,3)

disp(['PNSD: ',num2str(PNSD(n,2:10))])
PND(n,5:end) == AeroTrak_PN(n,2:7)
%SidePak(n,2)
%DustTrak(n,2:6)