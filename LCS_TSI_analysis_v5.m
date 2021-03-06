%% This code is written in MATLAB
% The purpose aims to develop the estimator for ultra-fine particles based
% on low-cost sensors
% Written by:
% Martha Arbayani Bin Zaidan, Ph.D., Docent, CEng.
% Research Associate Professor, Nanjing University, China
% Senior Scientist, Helsinki University, Finland


%% UNDERSTANDING AEROSOL DATA: 
% To understand how PND was formed. 
% In practice, PND was formed using data from CPC, Ptrak and AeroTrak.
% We also need to understand how PMD was calculated.

clc;clear;close all;
% select any index in the data
n=10002;2002;2150;%10000;%2002;
r=0;%-2;
load('Data_clean_processed.mat'); D=1;
%load('Data_processed_Heaters/Data_processed_Kerosene_Heaters.mat'); D=2;
%load('Data_processed_Heaters/Data_processed_NaturalGas_Heaters.mat'); D=3;

if D == 1
    disp('Smoking Data')
    disp(['CPC+Ptrak+AeroTrak: ',num2str(roundn([CPC(n,8),nan,Ptrak(n,8),AeroTrak_PN(n,8:13)],r))])
    disp(['PND               : ',num2str(roundn(PND(n,8:16),r))])
    disp(['PNSD              : ',num2str(roundn(PNSD(n,8:16),r))])
    disp(['PMD               : ',num2str(roundn(PMD(n,8:16),r))])
    disp(['PMSD              : ',num2str(roundn(PMSD(n,8:16),r))])
    
    bin2 = CPC(n,8) - Ptrak(n,8);
    bin3 = Ptrak(n,8) - sum(AeroTrak_PN(n,8:13));
    disp(bin2)
    disp(bin3)
    
    figure(1); fig=gcf;
    X = CPC(:,8) - Ptrak(:,8); % CPC - PTRAK == 10 - 25 nm
    Y = PND(2:end,9);          
    subplot(211);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('CPC - Ptrak');title('Bin 1: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off   
    X = Ptrak(:,8) - sum(AeroTrak_PN(:,8:13),2); % Ptrak - sum(aerotrak) == 25 - 300 nm
    Y = PND(2:end,10);
    subplot(212);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('Ptrak - sum(aerotrak)');title('Bin 2: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off
    set(findall(fig,'-property','FontSize'),'FontSize',22);
    
    figure(2); fig=gcf;
    X = CPC(:,8); % CPC 
    Y = PND(2:end,8);          
    subplot(311);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('CPC');title('Bin 1: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off
    
    X = CPC(:,8) - Ptrak(:,8); % CPC - PTRAK == 10 - 25 nm
    Y = PND(2:end,9);
    subplot(312);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('CPC - Ptrak');title('Bin 2: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off
      
    X = Ptrak(:,8) - sum(AeroTrak_PN(:,8:13),2); % Ptrak - sum(aerotrak) == 25 - 300 nm
    Y = PND(2:end,10);
    subplot(313);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('Ptrak - sum(aerotrak)');title('Bin 3: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off
    set(findall(fig,'-property','FontSize'),'FontSize',22);
else 
    disp('Kerosene/Natural Gas Data')
    disp(['CPC+Ptrak+AeroTrak: ',num2str(roundn([CPC(n,2),nan,Ptrak(n,2),AeroTrak_PN(n,2:7)],r))])
    disp(['PND               : ',num2str(roundn(PND(n,2:10),r))])
    disp(['PNSD              : ',num2str(roundn(PNSD(n,2:10),r))])
    disp(['PMD               : ',num2str(roundn(PMD(n,2:10),r))])
    disp(['PMSD              : ',num2str(roundn(PMSD(n,2:10),r))])
    
    bin1 = CPC(n,2) - Ptrak(n,2);
    bin2 = Ptrak(n,2) - sum(AeroTrak_PN(n,2:7));
    disp(bin1)
    disp(bin2)
    
%     figure(1); fig=gcf;
%     subplot(211);scatter(CPC(:,2) - Ptrak(:,2),PND(:,3));
%     xlabel('PND');ylabel('CPC - Ptrak');title('Bin 1: PND')
%     subplot(212);scatter(Ptrak(:,2) - sum(AeroTrak_PN(:,2:7),2),PND(:,4));
%     xlabel('PND');ylabel('Ptrak - sum(aerotrak)');title('Bin 2: PND')
%     set(findall(fig,'-property','FontSize'),'FontSize',22);
    
    figure(2); fig=gcf;
    X = CPC(:,2) - Ptrak(:,2);
    Y = PND(:,3);
    subplot(211);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('CPC - Ptrak');title('Bin 1: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off   
    X = Ptrak(:,2) - sum(AeroTrak_PN(:,2:7),2);
    Y = PND(:,4);
    subplot(212);scatter(X,Y); hold on; grid on
    xlabel('PND');ylabel('Ptrak - sum(aerotrak)');title('Bin 2: PND')
    x = linspace(0,max(Y));
    y = linspace(0,max(Y));
    plot(x,y,'r');hold off
    set(findall(fig,'-property','FontSize'),'FontSize',22);
    
end



%% In this part, we pre-process the data
clear;clc;close all;

% SMOKING ACTIVITIES

load('Data_clean_processed.mat')


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
AT_met = timetable(T1,AT_T,AT_RH);
%clear AT_DoY1 AT_T AT_RH


% ClasOhlson_met
T1 = datetime(ClasOhlson_met(:,1:6));
CO_DoY = ClasOhlson_met(:,7);
CO_T   = ClasOhlson_met(:,8);
CO_RH  = ClasOhlson_met(:,9);
CO_P   = ClasOhlson_met(:,10);
CO_met = timetable(T1,CO_T,CO_RH,CO_P);
%clear T1 CO_DoY CO_T CO_RH CO_P


% DustTrak
T1 = datetime(DustTrak(:,1:6));
DustTrak_DoY = DustTrak(:,7);
DustTrak_c   = DustTrak(:,8:12);
DustTrak1 = timetable(T1,DustTrak_c);
%clear T1 DustTrak_DoY DustTrak_c


% SidePak
T1 = datetime(SidePak(:,1:6));
SidePak_DoY = SidePak(:,7);
SidePak_c   = SidePak(:,8);
SidePak1 = timetable(T1,SidePak_c);
%clear T1 DustTrak_DoY DustTrak_c


% PND
T1=datetime(PND(2:end,1:6));
PND_DoY = PND(2:end,7);
PND_c   = PND(2:end,8:end);
PND1 = timetable(T1,PND_c);
PND2 = sortrows(PND1);
PND3 = PND2(1:end-1,:);
%PND2 = unique(PND1);
%clear T1 PND_DoY PND_c


% PMD
T1=datetime(PMD(2:end,1:6));
PMD_DoY = PMD(2:end,7);
PMD_c   = PMD(2:end,8:end);
PMD1 = timetable(T1,PMD_c);
PMD2 = sortrows(PMD1);
PMD3 = PMD2(1:end-1,:);
%clear T1 PMD_DoY PMD_c

% LCS_G1:
T1 = datetime(ISEE_LCS_G1(:,1:6));
LCS_G1 = ISEE_LCS_G1(:,8);
LCS_G1_T = timetable(T1,LCS_G1);
%clear T1 LCS_G1


% LCS_G2_01 (PM2.5 and MET)
T1 = datetime(ISEE_LCS_G201(:,1:6));
LCS_G2_01 = ISEE_LCS_G201(:,8); % PM2.5
LCS_G2_01_met = ISEE_LCS_G201_met(:,8:end); % RH, T and P
LCS_G2_01_T = timetable(T1,LCS_G2_01,LCS_G2_01_met);
%clear T1 LCS_G2_01 LCS_G2_01_met


% LCS_G2_02 (PM2.5 and MET)
T1 = datetime(ISEE_LCS_G202(:,1:6));
LCS_G2_02 = ISEE_LCS_G202(:,8); % PM2.5
LCS_G2_02_met = ISEE_LCS_G202_met(:,8:end); % RH, T and P
LCS_G2_02_T = timetable(T1,LCS_G2_02,LCS_G2_02_met);

% https://www.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html

% DATAs = DATAsmoking

    DATAs1 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));

    DATAs2 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
        PND3, PMD3, ... % PND1, PMSD1,...
        LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');

    DATAs3 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
        PND3, PMD3, ... % PND1, PMSD1,...
        LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely');

    DATAs4 =[AeroTrak_met(:,1:6),AT_T,AT_RH,CO_T,CO_RH,CO_P, ...
        DustTrak_c,SidePak_c,PND_c,PMD_c, ...
        LCS_G1,LCS_G2_01,LCS_G2_01_met, ...
        LCS_G2_02,LCS_G2_02_met];

% KEROSENE HEATER

clearvars -except DATAs1 DATAs2 DATAs3 DATAs4% syn_s syn_k syn_n
load('Data_processed_Heaters/Data_processed_Kerosene_Heaters.mat');


% AeroTrak_Met

Date  = datevec(AeroTrak_met(:,1));
Datey = [repmat(2021,size(Date,1),1) ,Date(:,2:end)];
T1    = datetime(Datey);

AT_T = AeroTrak_met(:,2);
AT_RH = AeroTrak_met(:,3);
AT_met = timetable(T1,AT_T,AT_RH);

% ClasOhlson_met
CO_T   = ClasOhlson_met(:,2);
CO_RH  = ClasOhlson_met(:,3);
CO_P   = ClasOhlson_met(:,4);
CO_met = timetable(T1,CO_T,CO_RH,CO_P);

% DustTrak
DustTrak_c   = DustTrak(:,2:end);
DustTrak1 = timetable(T1,DustTrak_c);

% SidePak
SidePak_c   = SidePak(:,2);
SidePak1 = timetable(T1,SidePak_c);

% PND
PND_c   = PND(:,2:end);
PND1 = timetable(T1,PND_c);
PND2 = sortrows(PND1);
PND3 = PND2(1:end-1,:);

% PMD
PMD_c   = PMD(:,2:end);
PMD1 = timetable(T1,PMD_c);
PMD2 = sortrows(PMD1);
PMD3 = PMD2(1:end-1,:);

% LCS_G1:
LCS_G1 = ISEE_LCS_G1(:,2);
LCS_G1_T = timetable(T1,LCS_G1);

% LCS_G2_01 (PM2.5 and MET)
LCS_G2_01 = ISEE_LCS_G201(:,2); % PM2.5
LCS_G2_01_met = nan(size(ISEE_LCS_G201_met,1),3); %ISEE_LCS_G201_met(:,8:end); % RH, T and P
LCS_G2_01_T = timetable(T1,LCS_G2_01,LCS_G2_01_met);

% LCS_G2_02 (PM2.5 and MET)
LCS_G2_02 = ISEE_LCS_G202(:,2); % PM2.5
LCS_G2_02_met = nan(size(ISEE_LCS_G201_met,1),3);  % ISEE_LCS_G202_met(:,8:end); % RH, T and P
LCS_G2_02_T = timetable(T1,LCS_G2_02,LCS_G2_02_met);

    DATAk1 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));

    DATAk2 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');

    DATAk3 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely');

    DATAk4 =[Datey,AT_T,AT_RH,CO_T,CO_RH,CO_P, ...
        DustTrak_c,SidePak_c,PND_c,PMD_c, ...
        LCS_G1,LCS_G2_01,LCS_G2_01_met, ...
        LCS_G2_02,LCS_G2_02_met];

% NATURAL GAS HEATER
    
clearvars -except DATAs1 DATAs2 DATAs3 DATAs4 DATAk1 DATAk2 DATAk3 DATAk4 % syn_s syn_k syn_n
load('Data_processed_Heaters/Data_processed_NaturalGas_Heaters.mat')

% AeroTrak_Met

Date  = datevec(AeroTrak_met(:,1));
Datey = [repmat(2021,size(Date,1),1) ,Date(:,2:end)];
T1    = datetime(Datey);

AT_T = AeroTrak_met(:,2);
AT_RH = AeroTrak_met(:,3);
AT_met = timetable(T1,AT_T,AT_RH);

% ClasOhlson_met
CO_T   = nan(size(T1)); %ClasOhlson_met(:,2);
CO_RH  = nan(size(T1)); %ClasOhlson_met(:,3);
CO_P   = nan(size(T1)); %ClasOhlson_met(:,4);
CO_met = timetable(T1,CO_T,CO_RH,CO_P);

% DustTrak
DustTrak_c   = nan(size(T1,1),5); %DustTrak(:,2:end);
DustTrak1 = timetable(T1,DustTrak_c);

% SidePak
SidePak_c   = nan(size(T1)); %SidePak(:,2);
SidePak1 = timetable(T1,SidePak_c);

% PND
PND_c   = PND(:,2:end);
PND1 = timetable(T1,PND_c);
PND2 = sortrows(PND1);
PND3 = PND2(1:end-1,:);

% PMD
PMD_c   = PMD(:,2:end);
PMD1 = timetable(T1,PMD_c);
PMD2 = sortrows(PMD1);
PMD3 = PMD2(1:end-1,:);

% LCS_G1:
LCS_G1 = ISEE_LCS_G1(:,2);
LCS_G1_T = timetable(T1,LCS_G1);

% LCS_G2_01 (PM2.5 and MET)
LCS_G2_01 = ISEE_LCS_G201(:,2); % PM2.5
LCS_G2_01_met = ISEE_LCS_G201_met(:,2:end); % RH, T and P
LCS_G2_01_T = timetable(T1,LCS_G2_01,LCS_G2_01_met);

% LCS_G2_02 (PM2.5 and MET)
LCS_G2_02 = ISEE_LCS_G202(:,2); % PM2.5
LCS_G2_02_met = ISEE_LCS_G202_met(:,2:end); % RH, T and P
LCS_G2_02_T = timetable(T1,LCS_G2_02,LCS_G2_02_met);

DATAn1 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));

DATAn2 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');

DATAn3 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely');

DATAn4 =[Datey,AT_T,AT_RH,CO_T,CO_RH,CO_P, ...
        DustTrak_c,SidePak_c,PND_c,PMD_c, ...
        LCS_G1,LCS_G2_01,LCS_G2_01_met, ...
        LCS_G2_02,LCS_G2_02_met];
    
DATA4_label = {'year','month','day','hour','minute','second',...
    'AT_T','AT_RH','CO_T','CO_RH','CO_P', ...
        'DustTrak_c1','DustTrak_c2','DustTrak_c3','DustTrak_c4','DustTrak_c5', ...
        'SidePak_c', ...
        'PND_c1','PND_c2','PND_c3','PND_c4','PND_c5','PND_c6','PND_c7', ...
        'PND_c8','PND_c9', ...
        'PMD_c1','PMD_c2','PMD_c3','PMD_c4','PMD_c5','PMD_c6','PMD_c7', ...
        'PMD_c8','PMD_c9', ...
        'LCS_G1', ...
        'LCS_G2_01','LCS_G2_01_RH', 'LCS_G2_01_T', 'LCS_G2_01_P', ...
        'LCS_G2_02','LCS_G2_02_RH', 'LCS_G2_02_T', 'LCS_G2_02_P'};

% There are 4 different types of DATA for smoking, kerosine and natural gas
% The most reliable one is DATA2, because it syncronizes well
% DATA4 is the original data from Tareq

DATA1 = [DATAs1;DATAk1;DATAn1];
DATA2 = [DATAs2;DATAk2;DATAn2];
DATA3 = [DATAs3;DATAk3;DATAn3];
DATA4 = [DATAs4;DATAk4;DATAn4];

clearvars -except DATA1 DATA2 DATA3 DATA4 DATA4_label %DATAs DATAk DATAn syn_s syn_k syn_n
 

%% 
% PLOTTING LCS measurements 
% PLOTTING PM2.5 and UFP (PND) using DATA2 (after syncronization)
% PLOTTING PM2.5 and UFP (PND) using DATA4 (before syncronization)


DATA = DATA2;

clc
figure(1); fig =gcf; 
fig.Position = [100 100 540 400].*2.5;
ms=10;
subplot(511)
plot(DATA.AT_T,'b.','MarkerSize',ms);hold on
plot(DATA.CO_T,'r.','MarkerSize',ms);
plot(DATA.LCS_G2_01_met(:,2),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02_met(:,2),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
%title('Time-Series: Temperature','interpreter','latex')
legend('AeroTrak','ClasOhson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([10 40]); grid on
xlabel('Time Index','interpreter','latex'); ylabel('Temperature ($\circ$C)','interpreter','latex')
hold off

subplot(512)
plot(DATA.AT_RH,'b.','MarkerSize',ms);hold on
plot(DATA.CO_RH,'r.','MarkerSize',ms);
plot(DATA.LCS_G2_01_met(:,1),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02_met(:,1),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
%title('Time-Series: Relative Humidity','interpreter','latex')
legend('AeroTrak','ClasOhson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([10 70]); grid on
xlabel('Time Index','interpreter','latex'); ylabel('Relative Humidity (\%)','interpreter','latex')
hold off

subplot(513)
plot(DATA.CO_P,'r.','MarkerSize',ms);hold on
plot(DATA.LCS_G2_01_met(:,3),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02_met(:,3),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
%title('Time-Series: Pressure','interpreter','latex')
legend('ClasOhson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([850 950]); grid on
xlabel('Time Index','interpreter','latex'); ylabel('Pressure (mbar)','interpreter','latex')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

subplot(514)
plot(DATA.PMD_c(:,7),'b.','MarkerSize',ms); hold on
plot(DATA.DustTrak_c(:,2),'r.','MarkerSize',ms);
plot(DATA.SidePak_c(:,1),'m.','MarkerSize',ms);
plot(DATA.LCS_G1(:,1),'y.','MarkerSize',ms);
plot(DATA.LCS_G2_01(:,1),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02(:,1),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
%title('Time-Series: PM$_{2.5}$','interpreter','latex')
set(gca, 'YScale', 'log')
ylim([1e-2 5e4])
xlabel('Time Index','interpreter','latex'); ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
legend('PMD','DustTrak','SidePak','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
hold off


subplot(515)
plot(DATA.PND_c(:,1),'b.','MarkerSize',ms); hold on
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
set(gca, 'YScale', 'log')
ylim([0 7e5])
xlabel('Time Index','interpreter','latex'); ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);


% The plans
% Smoking experiment: we perform sensors validation for met variables 
% First : We compare between ref instruments (AT and CO) for Temp, RH
%         Conclusion: both results are consistent and become ref instruments
% Second: We compare between ref instruments and LCSs for Temp and RH
%         Conclusion: there are offsets, and we calibrate those
% We compare between CO_P and LCSs_P
%         Conclusion: there is a small offset and we calibrate it.
%
% Kerosene experiment: We use only CO, because LCSs were not measured
%
% Natural gas experiment: We use only LCSs (after calibrated)
% 
% How about PM2.5 vars
% We will validate all PM2.5 on the smoking experiment
% Results
% PM2.5 between PMD(:,3) and LCSs are good, 
% however PMD(:,3) and (SidePak and DustTrak) are bad 
% We will calibrate PM2.5 using PMD 
% 
% ABout UFP: We estimate small PND for CPC and Ptrak and AeroTrak (1 size)
% We then develop calibrators for met variables here
% Kerosine experiment: we use only CO_met data
% Natural gas experiment: we use only LCS_met data
% we used the CO met data for the kerosine experiment and 
% to train UFP models using PM2.5 and met vars
% 
% PROBLEMS:
% 1: the interpolation may not work well

% Aerosol size distribution

Exp_smoking = [1:1:11521]';
Exp_kerosine = [11521:1:30242]';
Exp_gas = [30243:1:54721]';

inst_res_PND = [0.025,0.3,0.5,1,2.5,5,10,45].*1e3;
Time_index   = [1:1:size(DATA.T1,1)]';

% AFTER SYNCRONIZATION (DATA2) 
figure(2);fig=gcf; 
fig.Position = [100 100 540 400].*2.5;
FS=26; PMx=6;
lc =1e0; hc=10e5;
tiledlayout(2,3);
nexttile % subplot(231)
plot(Time_index(Exp_smoking,1),DATA.PMD_c(Exp_smoking,PMx),'b.'); hold on
plot(Time_index(Exp_smoking,1),DATA.LCS_G1(Exp_smoking,1),'g.');
plot(Time_index(Exp_smoking,1),DATA.LCS_G2_01(Exp_smoking,1),'r.');
plot(Time_index(Exp_smoking,1),DATA.LCS_G2_02(Exp_smoking,1),'m.');
hold off
legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Smoking','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(232)
plot(Time_index(Exp_kerosine,1),DATA.PMD_c(Exp_kerosine,PMx),'b.'); hold on
plot(Time_index(Exp_kerosine,1),DATA.LCS_G1(Exp_kerosine,1),'g.');
plot(Time_index(Exp_kerosine,1),DATA.LCS_G2_01(Exp_kerosine,1),'r.');
plot(Time_index(Exp_kerosine,1),DATA.LCS_G2_02(Exp_kerosine,1),'m.');
hold off
legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Kerosene','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(233)
plot(Time_index(Exp_gas,1),DATA.PMD_c(Exp_gas,PMx),'b.'); hold on
plot(Time_index(Exp_gas,1),DATA.LCS_G1(Exp_gas,1),'g.');
plot(Time_index(Exp_gas,1),DATA.LCS_G2_01(Exp_gas,1),'r.');
plot(Time_index(Exp_gas,1),DATA.LCS_G2_02(Exp_gas,1),'m.');
hold off
legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Natureal Gas','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(234)
h1 = pcolor(Time_index(Exp_smoking,1)',inst_res_PND',[DATA1.PND_c(Exp_smoking,1),DATA1.PND_c(Exp_smoking,3:end)]');
set(h1, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC Smoking','interpreter','latex')
colormap jet
caxis([lc,hc])
set(gca, 'YScale', 'log','colorscale','log')

nexttile % subplot(235)
h2 = pcolor(Time_index(Exp_kerosine,1)',inst_res_PND',[DATA1.PND_c(Exp_kerosine,1),DATA1.PND_c(Exp_kerosine,3:end)]');
set(h2, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC kerosene','interpreter','latex')
colormap jet
caxis([lc,hc]) 
set(gca, 'YScale', 'log','colorscale','log')

nexttile % subplot(236)
h3 = pcolor(Time_index(Exp_gas,1)',inst_res_PND',[DATA1.PND_c(Exp_gas,1),DATA1.PND_c(Exp_gas,3:end)]');
set(h3, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC natural gas','interpreter','latex')
colormap jet
caxis([lc, hc])
set(gca, 'YScale', 'log','colorscale','log')

cb = colorbar;
%cb.Limits= [1e0, 1e5];
cb.Layout.Tile = 'east';
%cb.Location = 'southoutside';

% - Build title axes and title.
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text( 0.5, 0, 'After syncronize', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

set(findall(fig,'-property','FontSize'),'FontSize',FS);

% BEFORE SYNCRONIZATION (DATA4)
figure(3);fig=gcf; 
fig.Position = [100 100 540 400].*2.5;
FS=26; PMx4=33;
lc =1e0; hc=10e5;
tiledlayout(2,3);
nexttile % subplot(231)
plot(Time_index(Exp_smoking,1),DATA4(Exp_smoking,PMx4),'b.'); hold on
plot(Time_index(Exp_smoking,1),DATA4(Exp_smoking,36),'g.');
plot(Time_index(Exp_smoking,1),DATA4(Exp_smoking,37),'r.');
plot(Time_index(Exp_smoking,1),DATA4(Exp_smoking,41),'m.');
hold off
legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Smoking','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(232)
plot(Time_index(Exp_kerosine,1),DATA4(Exp_kerosine,PMx4),'b.'); hold on
plot(Time_index(Exp_kerosine,1),DATA4(Exp_kerosine,36),'g.');
plot(Time_index(Exp_kerosine,1),DATA4(Exp_kerosine,37),'r.');
plot(Time_index(Exp_kerosine,1),DATA4(Exp_kerosine,41),'m.');

hold off
legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Kerosene','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(233)
plot(Time_index(Exp_gas,1),DATA4(Exp_gas,PMx4),'b.'); hold on
plot(Time_index(Exp_gas,1),DATA4(Exp_gas,36),'g.');
plot(Time_index(Exp_gas,1),DATA4(Exp_gas,37),'r.');
plot(Time_index(Exp_gas,1),DATA4(Exp_gas,41),'m.');
hold off
legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Natureal Gas','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(234)
h1 = pcolor(Time_index(Exp_smoking,1)',inst_res_PND',[DATA4(Exp_smoking,18),DATA4(Exp_smoking,20:26)]');
set(h1, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC Smoking','interpreter','latex')
colormap jet
caxis([lc,hc])
set(gca, 'YScale', 'log','colorscale','log')

nexttile % subplot(235)
h2 = pcolor(Time_index(Exp_kerosine,1)',inst_res_PND',[DATA4(Exp_kerosine,18),DATA4(Exp_kerosine,20:26)]');
set(h2, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC kerosene','interpreter','latex')
colormap jet
caxis([lc,hc]) 
set(gca, 'YScale', 'log','colorscale','log')

nexttile % subplot(236)
h3 = pcolor(Time_index(Exp_gas,1)',inst_res_PND',[DATA4(Exp_gas,18),DATA4(Exp_gas,20:26)]');
set(h3, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC natural gas','interpreter','latex')
colormap jet
caxis([lc, hc])
set(gca, 'YScale', 'log','colorscale','log')

cb = colorbar;
cb.Layout.Tile = 'east';

% - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
 text( 0.5, 0, 'Before syncronize', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

set(findall(fig,'-property','FontSize'),'FontSize',FS);


%% Scatter plots between PMD and LCSs
clc
PMx=6;
figure(4); fig=gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('PM_{2.5}')
subplot(221)
scatter(DATA.PMD_c(:,PMx),DATA.DustTrak_c(:,2))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('DustTrak','interpreter','latex')
xlim([1e-2 1e4]);ylim([1e-2 1e4]);grid on
subplot(222)
scatter(DATA.PMD_c(:,PMx),DATA.SidePak_c(:,1))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('SidePak','interpreter','latex')
xlim([1e-2 1e4]);ylim([1e-2 1e4]);grid on
subplot(223)
scatter(DATA.PMD_c(Exp_smoking,PMx),DATA.LCS_G1(Exp_smoking,1),'b');hold on
scatter(DATA.PMD_c(Exp_kerosine,PMx),DATA.LCS_G1(Exp_kerosine,1),'r')
scatter(DATA.PMD_c(Exp_gas,PMx),DATA.LCS_G1(Exp_gas,1),'g'); hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('$\mathcal{L}_1$','interpreter','latex')
xlim([1e-2 1e4]);ylim([1e-2 1e4]);grid on
subplot(224)
scatter(DATA.PMD_c(Exp_smoking,PMx),DATA.LCS_G2_01(Exp_smoking,1),'b');hold on
scatter(DATA.PMD_c(Exp_kerosine,PMx),DATA.LCS_G2_01(Exp_kerosine,1),'r')
scatter(DATA.PMD_c(Exp_gas,PMx),DATA.LCS_G2_01(Exp_gas,1),'g'); hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([1e-2 1e4]);ylim([1e-2 1e4]);grid on
set(findall(fig,'-property','FontSize'),'FontSize',22);

figure(5); fig = gcf;
fig.Position = [100 100 540 400].*2.5;
subplot(241)
scatter(DATA.PMD_c(:,PMx),DATA.DustTrak_c(:,2))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('DustTrak','interpreter','latex')
subplot(242)
scatter(DATA.PMD_c(:,PMx),DATA.SidePak_c(:,1))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('SidePak','interpreter','latex')
subplot(243)
scatter(DATA.PMD_c(:,PMx),DATA.LCS_G1(:,1))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('$\mathcal{L}_{1}$','interpreter','latex')
subplot(244)
scatter(DATA.PMD_c(:,PMx),DATA.LCS_G2_01(:,1))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('AeroTrak','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
subplot(245)
scatter(DATA.DustTrak_c(:,2),DATA.LCS_G1(:,1))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('DustTrak','interpreter','latex')
ylabel('$\mathcal{L}_{1}$','interpreter','latex')
subplot(246)
scatter(DATA.DustTrak_c(:,2),DATA.LCS_G2_01(:,1))
xlabel('DustTrak','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
subplot(247)
scatter(DATA.SidePak_c(:,1),DATA.LCS_G1(:,1))
xlabel('SidePak','interpreter','latex')
ylabel('$\mathcal{L}_{1}$','interpreter','latex')
subplot(248)
scatter(DATA.SidePak_c(:,1),DATA.LCS_G2_01(:,1))
xlabel('SidePak','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
set(findall(fig,'-property','FontSize'),'FontSize',22);

%% SCATTER PLOTS BETWEEN MET VARS

figure(6);fig=gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('Temperature')
subplot(221)
scatter(DATA.AT_T(:,1),DATA.CO_T(:,1)); hold on
xlabel('AeroTrak','interpreter','latex')
ylabel('Clas Ohlson','interpreter','latex')
xlim([1e1 4e1]);ylim([1e1 4e1]);grid on
x = linspace(1e1,4e1);
y = linspace(1e1,4e1);
plot(x,y,'r');hold off
subplot(222)
scatter(DATA.AT_T(:,1),DATA.LCS_G2_01_met(:,2));hold on
xlabel('AeroTrak','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([1e1 4e1]);ylim([1e1 4e1]);grid on
x = linspace(1e1,4e1);
y = linspace(1e1,4e1);
plot(x,y,'r');hold off
subplot(223)
scatter(DATA.CO_T(:,1),DATA.AT_T(:,1));hold on
xlabel('Clas Ohlson','interpreter','latex')
ylabel('AeroTrak','interpreter','latex')
xlim([1e1 4e1]);ylim([1e1 4e1]);grid on
x = linspace(1e1,4e1);
y = linspace(1e1,4e1);
plot(x,y,'r');hold off
subplot(224)
scatter(DATA.CO_T(:,1),DATA.LCS_G2_01_met(:,2)); hold on
xlabel('Clas Ohlson','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([1e1 4e1]);ylim([1e1 4e1]);grid on
x = linspace(1e1,4e1);
y = linspace(1e1,4e1);
plot(x,y,'r');hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

figure(7); fig=gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('Relative Humidity')
subplot(221)
scatter(DATA.AT_RH(:,1),DATA.CO_RH(:,1));hold on
xlabel('AeroTrak','interpreter','latex')
ylabel('Clas Ohlson','interpreter','latex')
xlim([1e1 5e1]);ylim([1e1 5e1]);grid on
x = linspace(1e1,5e1);
y = linspace(1e1,5e1);
plot(x,y,'r');hold off
subplot(222)
scatter(DATA.AT_RH(:,1),DATA.LCS_G2_01_met(:,1));hold on
xlabel('AeroTrak','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([1e1 5e1]);ylim([1e1 5e1]);grid on
x = linspace(1e1,5e1);
y = linspace(1e1,5e1);
plot(x,y,'r');hold off
subplot(223)
scatter(DATA.CO_RH(:,1),DATA.AT_RH(:,1)); hold on
xlabel('Clas Ohlson','interpreter','latex')
ylabel('AeroTrak','interpreter','latex')
xlim([1e1 5e1]);ylim([1e1 5e1]);grid on
x = linspace(1e1,5e1);
y = linspace(1e1,5e1);
plot(x,y,'r');hold off
subplot(224)
scatter(DATA.CO_RH(:,1),DATA.LCS_G2_01_met(:,1)); hold on
xlabel('Clas Ohlson','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([1e1 5e1]);ylim([1e1 5e1]);grid on
x = linspace(1e1,5e1);
y = linspace(1e1,5e1);
plot(x,y,'r');hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);

figure(8);fig=gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('Pressure')
subplot(121)
scatter(DATA.CO_P(:,1),DATA.LCS_G2_01_met(:,3)); hold on
xlabel('Clas Ohlson','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([8.95e2 9.05e2]);ylim([8.95e2 9.05e2]);grid on
x = linspace(8.5e2,9.5e2);
y = linspace(8.5e2,9.5e2);
plot(x,y,'r');hold off
subplot(122)
scatter(DATA.CO_P(:,1),DATA.LCS_G2_02_met(:,3)); hold on
xlabel('Clas Ohlson','interpreter','latex')
ylabel('$\mathcal{L}_{2a}$','interpreter','latex')
xlim([8.95e2 9.05e2]);ylim([8.95e2 9.05e2]);grid on
x = linspace(8.5e2,9.5e2);
y = linspace(8.5e2,9.5e2);
plot(x,y,'r');hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);


%% SCATTER PLOTS between CPC and other vars
% WE SHOULD NOT USE DATA FROM Exp_kerosine and Exp_gas, because
% the reference devices do not work all in those experiments

Ds = Exp_smoking;
Dk = Exp_kerosine;
Dg = Exp_gas;


D = Exp_smoking;
%D = Exp_kerosine;
%D = Exp_gas;
%D=[Exp_smoking;Exp_kerosine];
%D=[Exp_smoking;Exp_gas];
%D=[Exp_kerosine;Exp_gas];
%D=[Exp_smoking;Exp_kerosine;Exp_gas];

figure(9); fig=gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('CPC to other vars')
subplot(241)
scatter(DATA.PND_c(D,1),DATA.AT_T(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('Temp (AT)','interpreter','latex')
xlim([1e3 7e5]);ylim([20 35]);grid on
set(gca, 'XScale', 'log')
subplot(242)
scatter(DATA.PND_c(D,1),DATA.AT_RH(D,1));
%scatter(DATA.PND_c(D,1),DATA.CO_RH(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('RH (AT)','interpreter','latex')
xlim([1e3 7e5]);ylim([10 50]);grid on
set(gca, 'XScale', 'log')

subplot(243)
%scatter(DATA.PND_c(D,1),DATA.CO_P(D,1));
scatter(DATA.PND_c(D,1),DATA.LCS_G2_01_met(D,3));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('P (CO)','interpreter','latex')
xlim([1e3 7e5]);ylim([895 911]);grid on
set(gca, 'XScale', 'log')

subplot(244)
scatter(DATA.PND_c(D,1),DATA.PMD_c(D,PMx));
%scatter(DATA.PND_c(D,1),DATA.LCS_G1(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('PM$_{2.5}$ (AT)','interpreter','latex')
xlim([1e3 7e5]); 
ylim([1e0 1e3]);
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(findall(fig,'-property','FontSize'),'FontSize',22);

subplot(245)
scatter(DATA.PND_c(D,1),DATA.LCS_G1(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('PM$_{2.5}$ (LCS)','interpreter','latex')
xlim([1e3 7e5]); 
ylim([1e0 1e3]);
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(findall(fig,'-property','FontSize'),'FontSize',22);

subplot(246)
scatter(DATA.PND_c(Ds,1),DATA.LCS_G1(Ds,1),'b.'); hold on
scatter(DATA.PND_c(Dk,1),DATA.LCS_G1(Dk,1),'r.');
scatter(DATA.PND_c(Dg,1),DATA.LCS_G1(Dg,1),'g.');
xlabel('PNC (CPC)','interpreter','latex')
ylabel('PM$_{2.5}$ (LCS)','interpreter','latex')
xlim([1e3 7e5]); 
ylim([1e0 1e3]);
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(findall(fig,'-property','FontSize'),'FontSize',22);


subplot(247)
scatter(DATA.PND_c(D,1),DATA.LCS_G1(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('PM$_{2.5}$ (LCS)','interpreter','latex')
xlim([1e3 1e6]); 
ylim([1e0 1e3]);
grid on
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
set(findall(fig,'-property','FontSize'),'FontSize',22);

subplot(248)
Dat_temp = [DATA.PND_c(D,1),DATA.LCS_G1(D,1)];
%Dat_temp(isnan(Dat_temp)) = 0;
Dat_temp = Dat_temp( ~any( isnan( Dat_temp ) | isinf( Dat_temp ), 2 ),: );
feature1 = abs(cwt(Dat_temp(:,1)))';
feature2 = abs(cwt(Dat_temp(:,2)))'; 
scatter(feature1(:,2),feature2(:,2));
%scatter(abs(cwt(Dat_temp(:,1))),abs(cwt(Dat_temp(:,2))));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('PM$_{2.5}$ (LCS)','interpreter','latex')
%xlim([1e3 7e5]); 
%ylim([1e0 1e3]);
grid on
scatter(feature1(:,2),feature2(:,2));
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(findall(fig,'-property','FontSize'),'FontSize',22);



%%

D=[Exp_smoking;Exp_kerosine;Exp_gas];

figure(10)
subplot(211)
plot(DATA.PND_c(D,1));
subplot(212)
plot(DATA.LCS_G1(D,1));


figure(11)
Dat_temp = [DATA.PND_c(D,1),DATA.LCS_G1(D,1)];
%Dat_temp(isnan(Dat_temp)) = 0;
Dat_temp = Dat_temp( ~any( isnan( Dat_temp ) | isinf( Dat_temp ), 2 ),: );
feature1 = abs(wdenoise(Dat_temp(:,1)));
feature2 = abs(wdenoise(Dat_temp(:,2))); 
%feature1 = abs(cwt(Dat_temp(:,1)))';
%feature2 = abs(cwt(Dat_temp(:,2)))'; 
scatter(feature1(:,1),feature2(:,1));
%scatter(abs(cwt(Dat_temp(:,1))),abs(cwt(Dat_temp(:,2))));
%xlabel('PNC (CPC)','interpreter','latex')
%ylabel('PM$_{2.5}$ (LCS)','interpreter','latex')
%xlim([0 2]);
%xlim([1e3 7e5]); 
%ylim([1e0 1e3]);
%grid on
%scatter(feature1(:,2),feature2(:,2));
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
set(findall(fig,'-property','FontSize'),'FontSize',22);


%C = xcorr(DATA.PND_c(D,1),DATA.LCS_G1(D,1));

%% MATRIX PLOT between all variables

DATAx = [DATA.PMD_c(:,2),DATA.PMD_c(:,6),DATA.DustTrak_c(:,2),DATA.SidePak_c(:,1), ...
    DATA.LCS_G1(:,1),DATA.LCS_G2_01(:,1),DATA.LCS_G2_02(:,1) ...
    DATA.AT_T,DATA.CO_T,DATA.LCS_G2_01_met(:,2),DATA.LCS_G2_01_met(:,2) ...
    DATA.AT_RH,DATA.CO_RH,DATA.LCS_G2_01_met(:,1),DATA.LCS_G2_01_met(:,1) ...
    DATA.CO_P,DATA.LCS_G2_01_met(:,3),DATA.LCS_G2_01_met(:,3)
    ];

labelX = {'CPC','PMD', ...
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

corrP = corr(DATAx,'Type',CT,'Rows','pairwise');

figure(10);fig=gcf; 
fig.Position = [100 100 540 400].*2.5;
FS=26;

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

%% HISTOGRAMS of the measurements
FS = 16;
figure(11); fig = gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('Standarization Algorithm col.3')
subplot(5,3,1);
histogram(DATA.PND_c(:,1),'BinWidth',.05,'FaceColor','b');
set(gca, 'XScale', 'log')
title('CPC')
xlabel('PNC [cm^{-3}]')
subplot(5,3,2);
histogram(log10(DATA.PND_c(:,1)),'BinWidth',.05,'FaceColor','b');
title('CPC (log)')
xlabel('log PNC [cm^{-3}]')
xlim([0 6])
subplot(5,3,3);
histogram(log10(DATA.PND_c(Exp_smoking,1)),'BinWidth',.05,'FaceColor','b');
hold on
histogram(log10(DATA.PND_c(Exp_kerosine,1)),'BinWidth',.05,'FaceColor','r');
histogram(log10(DATA.PND_c(Exp_gas,1)),'BinWidth',.05,'FaceColor','g');
hold off
title('CPC (log)')
xlabel('log PNC [cm^{-3}]')
xlim([0 6])

subplot(5,3,4);
histogram(DATA.PMD_c(:,PMx),'BinWidth',.05,'FaceColor','b');
set(gca, 'XScale', 'log')
title('PM_{2.5}')
xlabel('PM_{2.5} [\mug/m^{3}]')
subplot(5,3,5);
histogram(log10(DATA.PMD_c(:,PMx)),'BinWidth',.05,'FaceColor','b');
title('PM_{2.5} (log)')
xlabel('log PM_{2.5} [\mug/m^{3}]')
xlim([-1 2])
subplot(5,3,6);
histogram(log10(DATA.PMD_c(Exp_smoking,PMx)),'BinWidth',.05,'FaceColor','b');
hold on
histogram(log10(DATA.PMD_c(Exp_kerosine,PMx)),'BinWidth',.05,'FaceColor','r');
histogram(log10(DATA.PMD_c(Exp_gas,PMx)),'BinWidth',.05,'FaceColor','g');
hold off
title('PM_{2.5} (log)')
xlabel('log PM_{2.5} [\mug/m^{3}]')
xlim([-1 2])

% https://www.mathworks.com/help/matlab/ref/double.normalize.html
subplot(5,3,7);
histogram(DATA.AT_T(:,1),'BinWidth',.05,'FaceColor','b');
title('Temp')
xlabel('Temp [\circC]')
subplot(5,3,8);
Temp = normalize(DATA.AT_T(:,1));
histogram(Temp,'BinWidth',.05,'FaceColor','b');
title('Temp (norm)')
xlabel('Temp (norm)')
subplot(5,3,9);
mu = nanmean(DATA.AT_T(Exp_smoking,1));
sigma = nanstd(DATA.AT_T(Exp_smoking,1));
Temp_s = (DATA.AT_T(Exp_smoking,1) - mu)/sigma;
Temp_k = (DATA.AT_T(Exp_kerosine,1) - mu)/sigma;
Temp_n = (DATA.AT_T(Exp_gas,1) - mu)/sigma;
histogram(Temp_s,'BinWidth',.05,'FaceColor','b'); hold on
histogram(Temp_k,'BinWidth',.05,'FaceColor','r');
histogram(Temp_n,'BinWidth',.05,'FaceColor','g');
title('Temp (norm)')
xlabel('Temp (norm)')
hold off

subplot(5,3,10);
histogram(DATA.AT_RH(:,1),'BinWidth',.05,'FaceColor','b');
title('RH')
xlabel('RH [%]')
subplot(5,3,11);
Temp = normalize(DATA.AT_RH(:,1));
histogram(Temp,'BinWidth',.05,'FaceColor','b');
title('RH (norm)')
xlabel('RH (norm)')
subplot(5,3,12);
mu = nanmean(DATA.AT_RH(Exp_smoking,1));
sigma = nanstd(DATA.AT_RH(Exp_smoking,1));
RH_s = (DATA.AT_RH(Exp_smoking,1) - mu)/sigma;
RH_k = (DATA.AT_RH(Exp_kerosine,1) - mu)/sigma;
RH_n = (DATA.AT_RH(Exp_gas,1) - mu)/sigma;
histogram(RH_s,'BinWidth',.05,'FaceColor','b'); hold on
histogram(RH_k,'BinWidth',.05,'FaceColor','r');
histogram(RH_n,'BinWidth',.05,'FaceColor','g');
title('RH (norm)')
xlabel('RH (norm)')
hold off

subplot(5,3,13);
histogram(DATA.CO_P(:,1),'BinWidth',.05,'FaceColor','b');
title('P')
xlabel('P [mbar]')
subplot(5,3,14);
Temp = normalize(DATA.CO_P(:,1));
histogram(Temp,'BinWidth',.05,'FaceColor','b');
title('P (norm)')
xlabel('P (norm)')
subplot(5,3,15);
mu = nanmean(DATA.CO_P(Exp_smoking,1));
sigma = nanstd(DATA.CO_P(Exp_smoking,1));
P_s = (DATA.CO_P(Exp_smoking,1) - mu)/sigma;
P_k = (DATA.CO_P(Exp_kerosine,1) - mu)/sigma;
P_n = (DATA.CO_P(Exp_gas,1) - mu)/sigma;
histogram(P_s,'BinWidth',.05,'FaceColor','b'); hold on
histogram(P_k,'BinWidth',.05,'FaceColor','r');
histogram(P_n,'BinWidth',.05,'FaceColor','g');
title('P (norm)')
xlabel('P (norm)')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);

disp('For normalization of temp, RH and P, perhaps the best is to use max and min number')

% Using min-max algorithm
% https://www.mathworks.com/help/matlab/ref/rescale.html
% B = rescale(A,l,u,'InputMin',inmin,'InputMax',inmax) uses the formula
% B = l + [(A-inmin)./(inmax-inmin)].*(u-l)

figure(12); fig = gcf;
fig.Position = [100 100 540 400].*2.5;
sgtitle('Min-Max Algorithm col.3')
subplot(5,3,1);
histogram(DATA.PND_c(:,1),'BinWidth',.05,'FaceColor','b');
set(gca, 'XScale', 'log')
title('CPC')
xlabel('PNC [cm^{-3}]')
subplot(5,3,2);
histogram(log10(DATA.PND_c(:,1)),'BinWidth',.05,'FaceColor','b');
title('CPC (log)')
xlabel('log PNC [cm^{-3}]')
xlim([0 6])
subplot(5,3,3);
histogram(log10(DATA.PND_c(Exp_smoking,1)),'BinWidth',.05,'FaceColor','b');
hold on
histogram(log10(DATA.PND_c(Exp_kerosine,1)),'BinWidth',.05,'FaceColor','r');
histogram(log10(DATA.PND_c(Exp_gas,1)),'BinWidth',.05,'FaceColor','g');
hold off
title('CPC (log)')
xlabel('log PNC [cm^{-3}]')
xlim([0 6])

subplot(5,3,4);
histogram(DATA.PMD_c(:,PMx),'BinWidth',.05,'FaceColor','b');
set(gca, 'XScale', 'log')
title('PM_{2.5}')
xlabel('PM_{2.5} [\mug/m^{3}]')
subplot(5,3,5);
histogram(log10(DATA.PMD_c(:,PMx)),'BinWidth',.05,'FaceColor','b');
title('PM_{2.5} (log)')
xlabel('log PM_{2.5} [\mug/m^{3}]')
xlim([-1 2])
subplot(5,3,6);
histogram(log10(DATA.PMD_c(Exp_smoking,PMx)),'BinWidth',.05,'FaceColor','b');
hold on
histogram(log10(DATA.PMD_c(Exp_kerosine,PMx)),'BinWidth',.05,'FaceColor','r');
histogram(log10(DATA.PMD_c(Exp_gas,PMx)),'BinWidth',.05,'FaceColor','g');
hold off
title('PM_{2.5} (log)')
xlabel('log PM_{2.5} [\mug/m^{3}]')
xlim([-1 2])

% https://www.mathworks.com/help/matlab/ref/double.normalize.html
subplot(5,3,7);
histogram(DATA.AT_T(:,1),'BinWidth',.05,'FaceColor','b');
title('Temp')
xlabel('Temp [\circC]')
subplot(5,3,8);
inmin = 10; inmax = 40; l = 0; u = 1;
Temp = l + [ (DATA.AT_T(:,1) -inmin)./(inmax-inmin)].*(u-l) ;
histogram(Temp,'BinWidth',.05,'FaceColor','b');
title('Temp (norm)')
xlabel('Temp (norm)')
subplot(5,3,9);
Temp_s = l + [ (DATA.AT_T(Exp_smoking,1) -inmin)./(inmax-inmin)].*(u-l) ;
Temp_k = l + [ (DATA.AT_T(Exp_kerosine,1) -inmin)./(inmax-inmin)].*(u-l) ;
Temp_n = l + [ (DATA.AT_T(Exp_gas,1) -inmin)./(inmax-inmin)].*(u-l) ;
histogram(Temp_s,'BinWidth',.05,'FaceColor','b'); hold on
histogram(Temp_k,'BinWidth',.05,'FaceColor','r');
histogram(Temp_n,'BinWidth',.05,'FaceColor','g');
title('Temp (norm)')
xlabel('Temp (norm)')
hold off

% https://www.mathworks.com/help/matlab/ref/double.normalize.html
subplot(5,3,10);
histogram(DATA.AT_RH(:,1),'BinWidth',.05,'FaceColor','b');
title('RH')
xlabel('RH [%]')
subplot(5,3,11);
inmin = 10; inmax = 50; l = 0; u = 1;
Temp = l + [ (DATA.AT_RH(:,1) -inmin)./(inmax-inmin)].*(u-l) ;
histogram(Temp,'BinWidth',.05,'FaceColor','b');
title('RH (norm)')
xlabel('RH (norm)')
subplot(5,3,12);
RH_s = l + [ (DATA.AT_RH(Exp_smoking,1) -inmin)./(inmax-inmin)].*(u-l) ;
RH_k = l + [ (DATA.AT_RH(Exp_kerosine,1) -inmin)./(inmax-inmin)].*(u-l) ;
RH_n = l + [ (DATA.AT_RH(Exp_gas,1) -inmin)./(inmax-inmin)].*(u-l) ;
histogram(RH_s,'BinWidth',.05,'FaceColor','b'); hold on
histogram(RH_k,'BinWidth',.05,'FaceColor','r');
histogram(RH_n,'BinWidth',.05,'FaceColor','g');
title('RH (norm)')
xlabel('RH (norm)')
hold off

% https://www.mathworks.com/help/matlab/ref/double.normalize.html
subplot(5,3,13);
histogram(DATA.CO_P(:,1),'BinWidth',.05,'FaceColor','b');
title('P')
xlabel('P [mbar]')
subplot(5,3,14);
inmin = 890; inmax = 910; l = 0; u = 1;
Temp = l + [ (DATA.CO_P(:,1) -inmin)./(inmax-inmin)].*(u-l) ;
histogram(Temp,'BinWidth',.05,'FaceColor','b');
title('P (norm)')
xlabel('P (norm)')
subplot(5,3,15);
P_s = l + [ (DATA.CO_P(Exp_smoking,1) -inmin)./(inmax-inmin)].*(u-l) ;
P_k = l + [ (DATA.CO_P(Exp_kerosine,1) -inmin)./(inmax-inmin)].*(u-l) ;
P_n = l + [ (DATA.CO_P(Exp_gas,1) -inmin)./(inmax-inmin)].*(u-l) ;
histogram(P_s,'BinWidth',.05,'FaceColor','b'); hold on
histogram(P_k,'BinWidth',.05,'FaceColor','r');
histogram(P_n,'BinWidth',.05,'FaceColor','g');
title('P (norm)')
xlabel('P (norm)')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);


%% MODELLING Linear Models and Shallow Neural Networks (version 1)

% In this case, we use the same training and testing data,
% The idea is to ensure if the modelling concept works in general

% If we train in smoking and test in kerosine, 
% the below fail completely because there may be a lot of NaN data scatter
% around the data.



D = Exp_smoking; 
%D = Exp_kerosine;
Dt = Exp_gas;
%Dt =[Exp_smoking;Exp_kerosine;Exp_gas];

% For met vars
% B = rescale(A,l,u,'InputMin',inmin,'InputMax',inmax) uses the formula
% B = l + [(A-inmin)./(inmax-inmin)].*(u-l)
% For aerosol, simply take log10

DATAm1 = [DATA.AT_T(D,1),DATA.AT_RH(D,1),DATA.CO_P(D,1),DATA.LCS_G1(D,1),DATA.PND_c(D,1)];
%DATAm1 = [DATA.AT_T(D,1),DATA.LCS_G2_01_met(D,1),DATA.LCS_G2_01_met(D,3),DATA.PMD_c(D,6),DATA.PND_c(D,1)];
%DATAm1 = [DATA.AT_T(D,1),DATA.LCS_G2_01_met(D,1),DATA.CO_P(D,1),DATA.LCS_G1(D,1),DATA.PND_c(D,1)];
%DATAm1 = [DATA.AT_T(D,1),DATA.AT_RH(D,1),DATA.CO_P(D,1),DATA.LCS_G2_02(D,1),DATA.PND_c(D,1)];
DATAm2 = zeros(size(DATAm1));

DATAt1 = [DATA.AT_T(Dt,1),DATA.AT_RH(Dt,1),DATA.CO_P(Dt,1),DATA.LCS_G1(Dt,1),DATA.PND_c(Dt,1)];
%DATAt1 = [DATA.AT_T(Dt,1),DATA.AT_RH(Dt,1),DATA.CO_P(Dt,1),DATA.LCS_G1(Dt,1),DATA.PND_c(Dt,1)];
%DATAt1 = [DATA.CO_T(Dt,1),DATA.CO_RH(Dt,1),DATA.CO_P(Dt,1),DATA.LCS_G1(Dt,1),DATA.PND_c(Dt,1)];
DATAt2 = zeros(size(DATAt1));

for n=1:size(DATAm1,2)
    if n == 1
        disp('Temp')
        inmin = 10; inmax = 40; l = 0; u = 1;
    elseif n==2
        disp('RH')
        inmin = 10; inmax = 50; l = 0; u = 1;
    elseif n==3 
        disp('P')
        inmin = 890; inmax = 910; l = 0; u = 1;
    else
        disp('PM2.5')
    end
    
    if n==1 || n==2 || n==3
        disp('met var')
        %DATAm2(:,n) = (DATAm1(:,n) - nanmean(DATAm1(:,n))) ./ nanstd(DATAm1(:,n));
        %DATAt2(:,n) = (DATAt1(:,n) - nanmean(DATAt1(:,n))) ./ nanstd(DATAt1(:,n));
        DATAm2(:,n) = l + [(DATAm1(:,n)-inmin)./(inmax-inmin)].*(u-l);
        DATAt2(:,n) = l + [(DATAt1(:,n)-inmin)./(inmax-inmin)].*(u-l);
    elseif n == 4
        disp('Aerosol')
        %DATAm2(:,n) = (DATAm1(:,n) - nanmean(DATAm1(:,n))) ./ nanstd(DATAm1(:,n));
        %DATAt2(:,n) = (DATAt1(:,n) - nanmean(DATAt1(:,n))) ./ nanstd(DATAt1(:,n));
        DATAm2(:,n) = log10(DATAm1(:,n));
        DATAt2(:,n) = log10(DATAt1(:,n));
        %DATAm2(:,n) = DATAm1(:,n);
        %DATAt2(:,n) = DATAt1(:,n);
    elseif n == 5
        disp('Aerosol')
        DATAm2(:,n) = log10(DATAm1(:,n));
        DATAt2(:,n) = log10(DATAt1(:,n));
    end
end

% Linear model
X = DATAm2(:,1:4);
Y = DATAm2(:,5);

Xt = DATAt2(:,1:4);
Yt = DATAt2(:,5);

Xt = X;
Yt = Y;


% LINEAR MODEL:
mdl = fitlm(X,Y);
Ypred_lm = predict(mdl,Xt);

model = 2
if model == 1
    disp('ANN')
% ANN model
inputs = X';
targets = Y';
inputs_t = Xt';
% Create a Fitting Network
hiddenLayerSize = 25;%15;
net = fitnet(hiddenLayerSize);
% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 90;%70/100;
net.divideParam.valRatio = 10/100;%15/100;
net.divideParam.testRatio = 0;%15/100;
 % Train the Network
[net,tr] = train(net,inputs,targets);
 % Test the Network
outputs = net(inputs_t);
Ypred = outputs;

elseif model == 2
    disp('GLM')
    %b = glmfit(X,Y,'gamma')
    %mdl = fitglm(Xt,[Y X],'gamma')
    %mdl = fitglm(X,y,'y ~ x1 + x2','Distribution','poisson');
    mdl = fitglm(X,Y,'Distribution','gamma');
    %mdl = fitrobust(X,Y);
    Ypred = predict(mdl,Xt);
end
figure(13); fig = gcf;
fig.Position = [100 100 540 400].*2.5;
subplot(221);
scatter(Yt,Ypred_lm);hold on
Xlim1 = 3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log PND (CPC)'); ylabel('log Est. PND (CPC)')

subplot(222);
scatter(10.^Yt,10.^Ypred_lm); hold on
Xlim1 = 1e3;
Ylim1 = 5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r'); hold off
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('PND (CPC)');ylabel('Est. PND (CPC)')

subplot(223);
scatter(Yt,Ypred);hold on
Xlim1 = 3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log PND (CPC)');ylabel('log Est. PND (CPC)')

subplot(224);
scatter(10.^Yt,10.^Ypred); hold on
Xlim1 = 1e3;
Ylim1 = 5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r'); hold off
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('PND (CPC)');ylabel('Est. PND (CPC)')

set(findall(fig,'-property','FontSize'),'FontSize',FS);



%% MODELLING Linear Models and Shallow Neural Networks (version 2)

% 1) use all data, randominze, and predict CPC
% 2) use Exp_smoking, estimate Exp_kerosine and Exp_gas

Ds = Exp_smoking; 
Dk = Exp_kerosine;
Dg = Exp_gas;
Da =[Exp_smoking;Exp_kerosine;Exp_gas];

% CHOOSE THE DATA with the number of inputs and output
Di = 2;

if Di == 2
    disp('Temp and PM2.5')
    DATAi = [ [DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
              [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
    %disp('RH and PM2.5')
    %DATAi = [[DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
    %         [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 3
    disp('Temp, RH and PM2.5')
    DATAi = [[DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
             [DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
         
    %disp('RH, P and PM2.5')
    %DATAi = [[DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
    %         [DATA.LCS_G2_01_met(Ds,3);DATA.CO_P(Dk,1);DATA.LCS_G2_01_met(Dg,3)], ...
    %         [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]]; 
     
elseif Di == 4
    disp('Temp, RH, P and PM2.5')
    %DATAi = [[DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
    %         [DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
    %         [DATA.LCS_G2_01_met(Ds,3);DATA.CO_P(Dk,1);DATA.LCS_G2_01_met(Dg,3)], ...
    %         [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
    DATAi = [[DATA.AT_T(Ds,1);DATA.AT_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
             [DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G2_01_met(Ds,3);DATA.CO_P(Dk,1);DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
end

% NORMALIZATION
DATAi1 = zeros(size(DATAi));
l = 0; u = 1;
for n = 1: Di
    if n==1
        disp('Temp'); inmin = 10; inmax = 40;
    elseif n==2
        disp('RH'); inmin = 10; inmax = 50;
    elseif n==3
        disp('P'); inmin = 890; inmax = 910;
    else
        disp('PM2.5')
    end
    if n == Di
        disp('PM2.5')
        DATAi1(:,n) = log10(DATAi(:,n));
    else
        %DATAi1(:,n) = l + [(DATAi(:,n)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,n) = DATAi(:,n);
    end
end


%if Di == 4
       
    CPC = DATA.PND_c([Ds;Dk;Dg],1);
    %figure(100);plot(DATA.PND_c(Ds,1));hold on;plot(DATA.PND_c(Dk,1),'r');plot(DATA.PND_c(Dg,1),'g'); hold off
    CPClog = log10(CPC);
    CPCgradient = gradient(CPC);
    CPCdiff = diff(CPC);
    CPCloggradient = gradient(CPClog);
    idx = find(CPCgradient <= nanmedian(CPCgradient)+10);
    CPCclean =CPC;
    CPCclean(idx,:)=nan;
    %CPCclean = CPC(idx,:);
    
    %%%DATAo  = CPC; 
    DATAo  = CPCclean;
    %DATAo  = [DATA.PND_c(Da,1)];
    DATAo1 = log10(DATAo);
    % figure(100);plot(DATAo1(Ds,1));hold on;plot(DATAo1(Dk,1),'r');plot(DATAo1(Dg,1),'g'); hold off
    
    
    figure(13)
    for n=1:size(DATAi,2)
        subplot(2,5,n);plot(DATAi(:,n),'.');hold on
        subplot(2,5,n+5);plot(DATAi1(:,n),'.');
    end
    subplot(2,5,5);plot(DATAo(:,1),'.');
    subplot(2,5,10);plot(DATAo1(:,1),'.');
    hold off
    
    
    
    figure(14); 
    subplot(511);plot(CPC,'b.');
    hold on; plot(CPCgradient,'r.'); hold off
    subplot(512);plot(CPClog,'b.');
    hold on; plot(CPCloggradient,'r.'); hold off
    subplot(513);plot(CPCgradient,'b.');hold on
    subplot(513);plot(CPCloggradient,'r.');hold off
    subplot(514);plot(CPC,'b.');hold on;plot(CPCclean,'r.');hold off
    subplot(515);plot(log10(CPC),'b.');hold on;plot(log10(CPCclean),'r.');hold off
    %
%end

% SELECT TRAINING AND TESTING DATA
R     = zeros(1,12);
MAPE  = zeros(1,12);

for test_no=1:12
    if test_no == 1; TRAIN = Ds; TEST = Dk; end
    if test_no == 2; TRAIN = Ds; TEST = Dg; end
    if test_no == 3; TRAIN = Ds; TEST = [Dk;Dg]; end
    if test_no == 4; TRAIN = Dk; TEST = Ds; end
    if test_no == 5; TRAIN = Dk; TEST = Dg; end
    if test_no == 6; TRAIN = Dk; TEST = [Ds;Dg]; end
    if test_no == 7; TRAIN = Dg; TEST = Ds; end
    if test_no == 8; TRAIN = Dg; TEST = Dk; end
    if test_no == 9; TRAIN = Dg; TEST = [Ds;Dk]; end
    if test_no == 10; TRAIN = [Ds;Dk]; TEST = Dg; end
    if test_no == 11; TRAIN = [Ds;Dg]; TEST = Dk; end
    if test_no == 12; TRAIN = [Dk;Dg]; TEST = Ds; end

clc
method = 7;
if method == 1
    disp('Training and Testing data is the same')
    DATAt  = [DATAi1,DATAo1];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    Xtr = DATAt1(:,1:end-1);
    Ytr = DATAt1(:,end);
    Xte = DATAt1(:,1:end-1);
    Yte = DATAt1(:,end); 
elseif method == 2
    disp('Training and Testing data only for smoking experiment')
    
    DATAt  = [DATAi1(Ds,:),DATAo1(Ds,:)];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    Xtr = DATAt1(:,1:end-1);
    Ytr = DATAt1(:,end);
    Xte = DATAt1(:,1:end-1);
    Yte = DATAt1(:,end);
    
    %Xtr = DATAi1(Ds,:);
    %Ytr = DATAo1(Ds,:);
    %Xte = DATAi1(Ds,:);
    %Yte = DATAo1(Ds,:); 
    
elseif method == 3
    disp('Training and Testing data only for smoking and kerosene experiments')
    
    DATAt  = [DATAi1([Ds;Dk],:),DATAo1([Ds;Dk],:)];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    Xtr = DATAt1(:,1:end-1);
    Ytr = DATAt1(:,end);
    Xte = DATAt1(:,1:end-1);
    Yte = DATAt1(:,end);
    
    %Xtr = DATAi1([Ds:Dk],:);
    %Ytr = DATAo1([Ds:Dk],:);
    %Xte = DATAi1([Ds:Dk],:);
    %Yte = DATAo1([Ds:Dk],:); 

elseif method == 4
    
    disp('Training and Testing data only for smoking, kerosene and gas experiments')  
    
    DATAt  = [DATAi1([Ds;Dk;Dg],:), DATAo1([Ds;Dk;Dg],:)];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    Xtr = DATAt1(:,1:end-1);
    Ytr = DATAt1(:,end);
    Xte = DATAt1(:,1:end-1);
    Yte = DATAt1(:,end);
    
    %Xtr = DATAi1([Ds:Dk:Dg],:);
    %Ytr = DATAo1([Ds:Dk:Dg],:);
    %Xte = DATAi1([Ds:Dk:Dg],:);
    %Yte = DATAo1([Ds:Dk:Dg],:); 
    
elseif method == 5
    disp('Training and Testing data are randomized')
    
    DATAt  = [DATAi1,DATAo1];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    percent = 70/100;
    p = size(DATAt1,1);
    rng(1986);s = rng;
    c = randperm(p)';
    tr = c( 1 : roundn(percent * p,0)     , 1);
    te = c( roundn(percent * p,0)+1 : end , 1);
    
    Xtr = DATAt1(tr,1:end-1);
    Ytr = DATAt1(tr,end);
    Xte = DATAt1(te,1:end-1);
    Yte = DATAt1(te,end);
    
elseif method ==6
    disp('Wavelet feature')
    
    DATAt  = [DATAi1,DATAo1];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    for n =1:5
        if n==4
            %d = 10.^(DATAt1(:,n))';
            d = DATAt1(:,n)';
        elseif n == 5
            %d = 10.^(DATAt1(:,n))';
            d = DATAt1(:,n)';
        else
            d = DATAt1(:,n)';
        end
        wv = 'db2';
        [c,l] = wavedec(d,3,wv);
        % approx = appcoef(c,l,'db2');
        xs = waverec(c,l,wv);
        %xs = wrcoef('a',c,l,'sym4',2);
        err = norm(d-xs)
        % figure(100);plot(d,'b.');hold on;plot(xs,'r--');hold off
        
        [cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
        
        cfs = cwt(d,'amor');
        %cfs = waveletScattering(d');
        
        d1 = abs(cfs(1,:)) ;
        
        
        xden = wdenoise(d,4);
        % xden = wdenoise(d);
        %figure(100);plot(d);hold on;plot(xden,'r--');hold off
        
        DATAt1(:,n) = xden ;
        %DATAt1(:,n) = abs(cfs(1,:))' ;
        %DATAt1(:,n) = xs' ;
        %DATAt1(:,n) = c' ;
        %DATAt1(:,n) = d ;      
    end
    
    percent = 70/100;
    p = size(DATAt1,1);
    rng(1986);s = rng;
    c = randperm(p)';
    tr = c( 1 : roundn(percent * p,0)     , 1);
    te = c( roundn(percent * p,0)+1 : end , 1);
    
    Xtr = DATAt1(tr,1:end-1);
    Ytr = DATAt1(tr,end);
    Xte = DATAt1(te,1:end-1);
    Yte = DATAt1(te,end);
    
    
    elseif method == 7
    disp('Wavelet feature 2')
    
    %[Ds;Dk;Dg]
    %DATAt  = [DATAi1,DATAo1];
    no = linspace(1, size(DATAi1,1), size(DATAi1,1))' ;
    DATAt  = [DATAi1,DATAo1,no];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    for n =1:size(DATAt1 ,2)-1
        if n==4
            %d = 10.^(DATAt1(:,n))';
            d = DATAt1(:,n)';
        elseif n == 5
            %d = 10.^(DATAt1(:,n))';
            d = DATAt1(:,n)';
        else
            d = DATAt1(:,n)';
        end
        xden = d;
        xden = wdenoise(d,4);
        %figure(100);plot(d);hold on;plot(xden,'r--');hold off        
        DATAt1(:,n) = xden ;
    end
    
    %percent = 70/100;
    %p = size(DATAt1,1);
    %rng(1986);s = rng;
    %c = randperm(p)';
    %tr = c( 1 : roundn(percent * p,0)     , 1);
    %te = c( roundn(percent * p,0)+1 : end , 1);
    %%%tr = ismember(DATAt1(:,end),[Ds;Dg]);
    %%%te = ismember(DATAt1(:,end),Dk);
    %%%tr = ismember(DATAt1(:,end),[Ds]);
    %%%te = ismember(DATAt1(:,end),Dg);
    %%%tr = ismember(DATAt1(:,end),[Ds]);
    %%%te = ismember(DATAt1(:,end),[Dg]);
    
    tr = ismember(DATAt1(:,end),TRAIN);
    te = ismember(DATAt1(:,end),TEST);
    
    Xtr = DATAt1(tr,1:end-2);
    Ytr = DATAt1(tr,end-1);
    Xte = DATAt1(te,1:end-2);
    Yte = DATAt1(te,end-1);
    
end


% Linear model
X = Xtr;
Y = Ytr;

Xt = Xte;
Yt = Yte;

LM = 1;
% LINEAR MODEL:
if LM == 1
    disp('Standard Linear model')
    mdl = fitlm(X,Y);
    Ypred_lm = predict(mdl,Xt);
elseif LM == 2
end


inputs = X';
targets = Y';
inputs_t = Xt';
ANN = 5;
if ANN ==1
    disp('standard ANN')
% ANN model
% Create a Fitting Network
hiddenLayerSize = 50;100;25;20;%15;
net = fitnet(hiddenLayerSize);
net.trainfcn = 'trainlm';%'trainbr'
% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 95/100;%70/100;
net.divideParam.valRatio = 5/100;
net.divideParam.testRatio = 0/100;
 % Train the Network
[net,tr] = train(net,inputs,targets);
 % Test the Network
outputs = net(inputs_t);
Ypred_snn = outputs;
elseif ANN ==2
    disp('TDNN')
     % https://se.mathworks.com/help/deeplearning/ref/timedelaynet.html
    X = con2seq(inputs); % u =  con2seq(u0);
    T = num2cell(targets);
    Xnew = con2seq(inputs_t);
    %net = timedelaynet(1:2,20);
    net = timedelaynet(1:3,50);
    [Xs,Xi,Ai,Ts] = preparets(net,X,T);
    net = train(net,Xs,Ts,Xi,Ai);
    [Y,Xf,Af] = net(Xs,Xi,Ai);
    perf = perform(net,Ts,Y);
    [netc,Xic,Aic] = closeloop(net,Xf,Af);
    %view(netc)
    y2 = netc(Xnew,Xic,Aic);
    est_pm25 = cell2mat(y2');
    Ypred_snn = est_pm25;
elseif ANN == 3
    disp('LSTM')
   % Some very useful links for LSTM:
    % https://se.mathworks.com/matlabcentral/answers/393034-i-am-unable-to-resolve-this-error-invalid-training-data-predictors-must-be-a-n-by-1-cell-array-of
    % https://stats.stackexchange.com/questions/352036/what-should-i-do-when-my-neural-network-doesnt-learn/352037#352037
    % https://se.mathworks.com/help/deeplearning/ref/trainingoptions.html
    
    inputs = inputs';
    inputs_t = inputs_t';
    outputs = targets';
    
    XTrain=cell(size(inputs,1),1);
    XTest=cell(size(inputs_t,1),1);
    
    for n=1:size(inputs,1)
        XTrain{n,1}=inputs(n,:)';%num2cell(inputs(n,:));
    end
    for n=1:size(inputs_t,1)
        XTest{n,1}=inputs_t(n,:)';%num2cell(inputs_test1(n,:));
    end
    
    YTrain =  num2cell(outputs);
    
    % Normalise training predictors
    mu = mean([XTrain{:}],2);
    sig = std([XTrain{:}],0,2);
    
    for i = 1:numel(XTrain)
        XTrain{i} = (XTrain{i} - mu) ./ sig;
    end
    
    % Define Network Architecture
    numResponses = size(YTrain{1},1);
    featureDimension = size(XTrain{1},1);
    numHiddenUnits = 25;%200;%100;%200;
    
    layers = [ ...
        sequenceInputLayer(featureDimension)
        lstmLayer(numHiddenUnits,'OutputMode','sequence')
        fullyConnectedLayer(30)%fullyConnectedLayer(30)%fullyConnectedLayer(50)
        dropoutLayer(0.5)
        fullyConnectedLayer(numResponses)
        regressionLayer];
    
    maxEpochs = 100;%500; % 60;
    %miniBatchSize = 20;
    miniBatchSize = 300; %500;%1000;
    
    options = trainingOptions('sgdm', ...%'rmsprop', ...%'adam', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'InitialLearnRate',0.01, ...%0.01, ...
        'GradientThreshold',1, ...
        'Shuffle','never', ...
        'Plots','training-progress',...
        'Verbose',0);
    
    % Train the network
    net = trainNetwork(XTrain,YTrain,layers,options);
        
    for i = 1:numel(XTest)
        XTest{i} = (XTest{i} - mu) ./ sig;
    end
    
    YPred = predict(net,XTest,'MiniBatchSize',1);
    est_pm25 = cell2mat(YPred);
    Ypred_snn = est_pm25 ; 
elseif ANN == 4
    disp('Bayesian LASSO regression')
    
    % Run standard LASSO first
    [LassoBetaEstimates,FitInfo] = lasso(Xtr,Ytr','Standardize',false);
    
    
    
    PriorMdl = bayeslm(p,'ModelType','lasso');
    table(PriorMdl.Lambda,'RowNames',PriorMdl.VarNames);
    
    ismissing = any(isnan(Xtr'),2);
    n = sum(~ismissing); % Effective sample size
    lambda = FitInfo.Lambda*n./sqrt(FitInfo.MSE);
    
    numlambda = numel(lambda);

% Preallocate
BayesLassoCoefficients = zeros(p+1,numlambda);
BayesLassoCI95 = zeros(p+1,2,numlambda);
fmseBayesLasso = zeros(numlambda,1);
BLCPlot = zeros(p+1,numlambda);

% Estimate and forecast
rng(10); % For reproducibility
for j = 1:numlambda
    PriorMdl.Lambda = lambda(j);
    [EstMdl,Summary] = estimate(PriorMdl,Xtr,Ytr','Display',false);
    BayesLassoCoefficients(:,j) = Summary.Mean(1:(end - 1));
    BLCPlot(:,j) = Summary.Mean(1:(end - 1));
    BayesLassoCI95(:,:,j) = Summary.CI95(1:(end - 1),:);
    idx = BayesLassoCI95(:,2,j) > 0 & BayesLassoCI95(:,1,j) <= 0;
    BLCPlot(idx,j) = 0;
    yFBayesLasso = forecast(EstMdl,XF);
    fmseBayesLasso(j) = sqrt(mean((yF - yFBayesLasso).^2));
end
    
L1Vals = sum(abs(BLCPlot(2:end,:)),1)/max(sum(abs(BLCPlot(2:end,:)),1));

figure;
plot(L1Vals,BLCPlot(2:end,:))
xlabel('L1');
ylabel('Coefficient Estimates');
yyaxis right
h = plot(L1Vals,fmseBayesLasso,'LineWidth',2,'LineStyle','--');
legend(h,'FMSE','Location','SW');
ylabel('FMSE');
title('Bayesian Lasso')
    
    
elseif ANN == 5
    Ypred_snn = Ypred_lm;
    disp('We do not work on ANN')
end


    R0 =  corrcoef(Yt,Ypred_lm);
    R(1,test_no) = R0(2,1);

    %MAPE(1,test_no) = nanmean(abs((Yt-Ypred_lm)./Yt));
    MAPE(1,test_no)=errperf(Yt,Ypred_lm,'mape');
end


figure(15);
subplot(221);
scatter(Yt,Ypred_lm);hold on
Xlim1 = 0;3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')

subplot(222);
scatter(10.^Yt,10.^Ypred_lm);hold on
Xlim1 = 1e3;
Ylim1 = 1e6;7.5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Real PNC (CPC)');ylabel('Est PNC (CPC)')

subplot(223);
scatter(Yt,Ypred_snn);hold on
Xlim1 = 0;3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')

subplot(224);
scatter(10.^Yt,10.^Ypred_snn);hold on
Xlim1 = 1e3;
Ylim1 = 1e6;7.5e5;
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
plot(x,y,'r');hold off
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Real PNC (CPC)');ylabel('Est PNC (CPC)')

set(findall(fig,'-property','FontSize'),'FontSize',22);


figure(16); fig = gcf;
plot(10.^Yt,'b.','MarkerSize',12);hold on;grid on
plot(10.^Ypred_lm,'g.','MarkerSize',12);
plot(10.^Ypred_snn,'r.','MarkerSize',12);
ylabel('PND [/cm$^{-3}$]','interpreter','latex')
xlabel('Time Index')
ylim([0 1e6])
set(gca, 'YScale', 'log')
legend('Real','LM','ANN','interpreter','latex')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);


%% THE BEST SOLUTIONS SO FAR ARE:
% BY USING Smoking experiment only
% (i) it seems that the upper regime is difficult to predict,
% (ii) there may be a problem with gas experiments, check inputs and
% outputs again for gas experiment because it give warning:
% "Warning: Regression design matrix is rank deficient to within machine
% precision."
% solution: https://www.mathworks.com/matlabcentral/answers/637610-regression-design-matrix-is-rank-deficient-what-to-do-next
% SOLUTION (i)
% Time-frequency domains features, study wavelet, find correlation within
% wavelet transform
%
% SOLUTION (ii)
% We may need to establish a model, such as MoE, where aerosol models are 
% divided into different regimes
%




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


