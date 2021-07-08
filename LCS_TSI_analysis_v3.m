
%% UNDERSTANDING AEROSOL DATA:

clc;clear
n=2150;%10000;%2002;
r=-2;
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

else 
    disp('Kerosene/Natural Gas Data')
    disp(['CPC+Ptrak+AeroTrak: ',num2str(roundn([CPC(n,2),nan,Ptrak(n,2),AeroTrak_PN(n,2:7)],r))])
    disp(['PND               : ',num2str(roundn(PND(n,2:10),r))])
    disp(['PNSD              : ',num2str(roundn(PNSD(n,2:10),r))])
    disp(['PMD               : ',num2str(roundn(PMD(n,2:10),r))])
    disp(['PMSD              : ',num2str(roundn(PMSD(n,2:10),r))])
end


%%
clear;clc;close all

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
%clear T1 LCS_G2_02 LCS_G2_02_met

%DATA_ts = synchronize(AT_met,AT_PN,CO_met,CPC1,Ptrak1,PND1,PMSD1, ...
%    LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));

% https://www.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html

% DATAs = DATAsmoking

%if syn_s == 1
    DATAs1 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));
%elseif syn_s ==2
    DATAs2 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
        PND3, PMD3, ... % PND1, PMSD1,...
        LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');
%else
    DATAs3 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
        PND3, PMD3, ... % PND1, PMSD1,...
        LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely');
%end

    DATAs4 =[AeroTrak_met(:,1:6),AT_T,AT_RH,CO_T,CO_RH,CO_P, ...
        DustTrak_c,SidePak_c,PND_c,PMD_c, ...
        LCS_G1,LCS_G2_01,LCS_G2_01_met, ...
        LCS_G2_02,LCS_G2_02_met];


clearvars -except DATAs1 DATAs2 DATAs3 DATAs4% syn_s syn_k syn_n
load('Data_processed_Heaters/Data_processed_Kerosene_Heaters.mat');


% AeroTrak_Met

Date  = datevec(AeroTrak_met(:,1));
Datey = [repmat(2021,size(Date,1),1) ,Date(:,2:end)];
T1    = datetime(Datey);

AT_T = AeroTrak_met(:,2);
AT_RH = AeroTrak_met(:,3);
AT_met = timetable(T1,AT_T,AT_RH);
%clear AT_T AT_RH

% ClasOhlson_met
CO_T   = ClasOhlson_met(:,2);
CO_RH  = ClasOhlson_met(:,3);
CO_P   = ClasOhlson_met(:,4);
CO_met = timetable(T1,CO_T,CO_RH,CO_P);
%clear CO_T CO_RH CO_P

% DustTrak
DustTrak_c   = DustTrak(:,2:end);
DustTrak1 = timetable(T1,DustTrak_c);
%clear DustTrak_c

% SidePak
SidePak_c   = SidePak(:,2);
SidePak1 = timetable(T1,SidePak_c);
%clear DustTrak_c

% PND
PND_c   = PND(:,2:end);
PND1 = timetable(T1,PND_c);
PND2 = sortrows(PND1);
PND3 = PND2(1:end-1,:);
%clear PND_c

% PMD
PMD_c   = PMD(:,2:end);
PMD1 = timetable(T1,PMD_c);
PMD2 = sortrows(PMD1);
PMD3 = PMD2(1:end-1,:);
%clear PMD_c

% LCS_G1:
LCS_G1 = ISEE_LCS_G1(:,2);
LCS_G1_T = timetable(T1,LCS_G1);
%clear LCS_G1

% LCS_G2_01 (PM2.5 and MET)
LCS_G2_01 = ISEE_LCS_G201(:,2); % PM2.5
LCS_G2_01_met = nan(size(ISEE_LCS_G201_met,1),3); %ISEE_LCS_G201_met(:,8:end); % RH, T and P
LCS_G2_01_T = timetable(T1,LCS_G2_01,LCS_G2_01_met);
%clear LCS_G2_01 LCS_G2_01_met

% LCS_G2_02 (PM2.5 and MET)
%T1 = datetime(ISEE_LCS_G202(:,1:6));
LCS_G2_02 = ISEE_LCS_G202(:,2); % PM2.5
LCS_G2_02_met = nan(size(ISEE_LCS_G201_met,1),3);  % ISEE_LCS_G202_met(:,8:end); % RH, T and P
LCS_G2_02_T = timetable(T1,LCS_G2_02,LCS_G2_02_met);
%clear LCS_G2_02 LCS_G2_02_met

%if syn_k == 1
    DATAk1 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));
%elseif syn_k ==2 
    DATAk2 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');
%else
    DATAk3 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely');
%end

    DATAk4 =[Datey,AT_T,AT_RH,CO_T,CO_RH,CO_P, ...
        DustTrak_c,SidePak_c,PND_c,PMD_c, ...
        LCS_G1,LCS_G2_01,LCS_G2_01_met, ...
        LCS_G2_02,LCS_G2_02_met];


%DATAk = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
%    PND3, PMD3, ... % PND1, PMSD1,...
%    LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');


clearvars -except DATAs1 DATAs2 DATAs3 DATAs4 DATAk1 DATAk2 DATAk3 DATAk4 % syn_s syn_k syn_n
load('Data_processed_Heaters/Data_processed_NaturalGas_Heaters.mat')

% AeroTrak_Met

Date  = datevec(AeroTrak_met(:,1));
Datey = [repmat(2021,size(Date,1),1) ,Date(:,2:end)];
T1    = datetime(Datey);

AT_T = AeroTrak_met(:,2);
AT_RH = AeroTrak_met(:,3);
AT_met = timetable(T1,AT_T,AT_RH);
%clear AT_T AT_RH

% ClasOhlson_met
CO_T   = nan(size(T1)); %ClasOhlson_met(:,2);
CO_RH  = nan(size(T1)); %ClasOhlson_met(:,3);
CO_P   = nan(size(T1)); %ClasOhlson_met(:,4);
CO_met = timetable(T1,CO_T,CO_RH,CO_P);
%clear CO_T CO_RH CO_P

% DustTrak
DustTrak_c   = nan(size(T1,1),5); %DustTrak(:,2:end);
DustTrak1 = timetable(T1,DustTrak_c);
%clear DustTrak_c

% SidePak
SidePak_c   = nan(size(T1)); %SidePak(:,2);
SidePak1 = timetable(T1,SidePak_c);
%clear DustTrak_c

% PND
PND_c   = PND(:,2:end);
PND1 = timetable(T1,PND_c);
PND2 = sortrows(PND1);
PND3 = PND2(1:end-1,:);
%clear PND_c

% PMD
PMD_c   = PMD(:,2:end);
PMD1 = timetable(T1,PMD_c);
PMD2 = sortrows(PMD1);
PMD3 = PMD2(1:end-1,:);
%clear PMD_c

% LCS_G1:
LCS_G1 = ISEE_LCS_G1(:,2);
LCS_G1_T = timetable(T1,LCS_G1);
%clear LCS_G1

% LCS_G2_01 (PM2.5 and MET)
LCS_G2_01 = ISEE_LCS_G201(:,2); % PM2.5
LCS_G2_01_met = ISEE_LCS_G201_met(:,2:end); % RH, T and P
LCS_G2_01_T = timetable(T1,LCS_G2_01,LCS_G2_01_met);
%clear LCS_G2_01 LCS_G2_01_met

% LCS_G2_02 (PM2.5 and MET)
%T1 = datetime(ISEE_LCS_G202(:,1:6));
LCS_G2_02 = ISEE_LCS_G202(:,2); % PM2.5
LCS_G2_02_met = ISEE_LCS_G202_met(:,2:end); % RH, T and P
LCS_G2_02_T = timetable(T1,LCS_G2_02,LCS_G2_02_met);
%clear LCS_G2_02 LCS_G2_02_met

%if syn_n == 1
    DATAn1 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'regular','linear','TimeStep',minutes(1));
%elseif syn_n==2
    DATAn2 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely','mean');
%else
    DATAn3 = synchronize(AT_met, CO_met, DustTrak1, SidePak1, ...
            PND3, PMD3, ... % PND1, PMSD1,...
            LCS_G1_T,LCS_G2_01_T,LCS_G2_02_T,'minutely');
%end

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


 
DATA1 = [DATAs1;DATAk1;DATAn1];
DATA2 = [DATAs2;DATAk2;DATAn2];
DATA3 = [DATAs3;DATAk3;DATAn3];
DATA4 = [DATAs4;DATAk4;DATAn4];

clearvars -except DATA1 DATA2 DATA3 DATA4 DATA4_label %DATAs DATAk DATAn syn_s syn_k syn_n


%%
DATA = DATA2;

clc
figure(1); fig =gcf; ms=10;
subplot(411)
plot(DATA.AT_T,'b.','MarkerSize',ms);hold on
plot(DATA.CO_T,'r.','MarkerSize',ms);
plot(DATA.LCS_G2_01_met(:,2),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02_met(:,2),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
title('Time-Series: Temperature','interpreter','latex')
legend('AeroTrak','ClasOhson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([10 40]); grid on
xlabel('Time Index','interpreter','latex'); ylabel('Temperature ($\circ$C)','interpreter','latex')
hold off
%set(findall(fig,'-property','FontSize'),'FontSize',22);

%figure(2); fig =gcf;
subplot(412)
plot(DATA.AT_RH,'b.','MarkerSize',ms);hold on
plot(DATA.CO_RH,'r.','MarkerSize',ms);
plot(DATA.LCS_G2_01_met(:,1),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02_met(:,1),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
title('Time-Series: Relative Humidity','interpreter','latex')
legend('AeroTrak','ClasOhson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([10 70]); grid on
xlabel('Time Index','interpreter','latex'); ylabel('Relative Humidity (\%)','interpreter','latex')
hold off
%set(findall(fig,'-property','FontSize'),'FontSize',22);

%figure(3); fig =gcf;
subplot(413)
plot(DATA.CO_P,'r.','MarkerSize',ms);hold on
plot(DATA.LCS_G2_01_met(:,3),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02_met(:,3),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
title('Time-Series: Pressure','interpreter','latex')
legend('ClasOhson','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
ylim([850 950]); grid on
xlabel('Time Index','interpreter','latex'); ylabel('Pressure (mbar)','interpreter','latex')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);


clc
%figure(4); fig =gcf; ms = 10;
subplot(414)
plot(DATA.PMD_c(:,7),'b.','MarkerSize',ms); hold on
plot(DATA.DustTrak_c(:,2),'r.','MarkerSize',ms);
plot(DATA.SidePak_c(:,1),'m.','MarkerSize',ms);
plot(DATA.LCS_G1(:,1),'y.','MarkerSize',ms);
plot(DATA.LCS_G2_01(:,1),'g.','MarkerSize',ms);
plot(DATA.LCS_G2_02(:,1),'c.','MarkerSize',ms);
xline(0,'-',{'Smoking','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(11521,'-',{'Kerosene','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
xline(30243,'-',{'Gas','Experiment'},'interpreter','latex','LabelOrientation','horizontal');
title('Time-Series: PM$_{2.5}$','interpreter','latex')
set(gca, 'YScale', 'log')
ylim([1e-2 5e4])
xlabel('Time Index','interpreter','latex'); ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
legend('PMD','DustTrak','SidePak','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);


% figure(11);plot(PMD(:,10));hold on;plot(DATA.PMD_c(Exp_smoking,3),'r--');hold off

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


clc
figure(2);fig=gcf; FS=26; PMx=6;
lc =1e0; hc=10e5;
tiledlayout(2,3);
nexttile
%subplot(231)
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

nexttile
%subplot(232)
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

nexttile
%subplot(233)
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

nexttile
%subplot(234)
h1 = pcolor(Time_index(Exp_smoking,1)',inst_res_PND',[DATA1.PND_c(Exp_smoking,1),DATA1.PND_c(Exp_smoking,3:end)]');
%h = pcolor(Time_index(Exp_smoking,1)',inst_res_PND',[DATA.PND_c(Exp_smoking,1),DATA.PND_c(Exp_smoking,3:end)]');
set(h1, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
%ylabel('PNC (cm$^{-3}$)','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC Smoking','interpreter','latex')
%yticks(gca,inst_res_AT)
%yticks(gca,[1 2 3 4 5 6 7])
%yticklabels(gca,{'0.3 \mum','0.5 \mum','1 \mum','2.5 \mum','5 \mum','10 \mum','45 \mum','interpreter','latex'})
%colorbar
colormap jet
caxis([lc,hc])
set(gca, 'YScale', 'log','colorscale','log')

nexttile
%subplot(235)
h2 = pcolor(Time_index(Exp_kerosine,1)',inst_res_PND',[DATA1.PND_c(Exp_kerosine,1),DATA1.PND_c(Exp_kerosine,3:end)]');
set(h2, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
%ylabel('PNC (cm$^{-3}$)','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC kerosene','interpreter','latex')
%colorbar
colormap jet
caxis([lc,hc]) 
set(gca, 'YScale', 'log','colorscale','log')

nexttile
%subplot(236)
h3 = pcolor(Time_index(Exp_gas,1)',inst_res_PND',[DATA1.PND_c(Exp_gas,1),DATA1.PND_c(Exp_gas,3:end)]');
set(h3, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
%ylabel('PNC (cm$^{-3}$)','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNC natural gas','interpreter','latex')
%colorbar
colormap jet
caxis([lc, hc])
set(gca, 'YScale', 'log','colorscale','log')

cb = colorbar;
%cb.Limits= [1e0, 1e5];
cb.Layout.Tile = 'east';
%cb.Location = 'southoutside';
%hp4 = get(subplot(2,3,6),'Position');
%h=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.0455  hp4(2)+hp4(3)*2.1]);

%colorbar(h1,jet)
%suptitle('After syncronize the data')


% - Build title axes and title.
 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
 set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
 text( 0.5, 0, 'After syncronize', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;


set(findall(fig,'-property','FontSize'),'FontSize',FS);

figure(21);fig=gcf; FS=26; PMx4=33;
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

%%
%% Scatter plots between PMD and LCSs
clc
PMx=6;
figure(3); fig=gcf;
suptitle('PM_{2.5}')
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




figure(4)
subplot(241)
scatter(DATA.PMD_c(:,PMx),DATA.DustTrak_c(:,2))
subplot(242)
scatter(DATA.PMD_c(:,PMx),DATA.SidePak_c(:,1))
subplot(243)
scatter(DATA.PMD_c(:,PMx),DATA.LCS_G1(:,1))
subplot(244)
scatter(DATA.PMD_c(:,PMx),DATA.LCS_G2_01(:,1))

subplot(245)
scatter(DATA.DustTrak_c(:,2),DATA.LCS_G1(:,1))
subplot(246)
scatter(DATA.DustTrak_c(:,2),DATA.LCS_G2_01(:,1))
subplot(247)
scatter(DATA.SidePak_c(:,1),DATA.LCS_G1(:,1))
subplot(248)
scatter(DATA.SidePak_c(:,1),DATA.LCS_G2_01(:,1))


%%
figure(5);fig=gcf;
suptitle('Temperature')
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

figure(6); fig=gcf;
suptitle('Relative Humidity')
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



figure(7);fig=gcf;
suptitle('Pressure')
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


%%
D = Exp_smoking;
%D = Exp_kerosine;
%D = Exp_gas;


figure(8); fig=gcf;
suptitle('CPC to other vars')
subplot(221)
scatter(DATA.PND_c(D,1),DATA.AT_T(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('Temp (AT)','interpreter','latex')
xlim([1e3 7e5]);ylim([20 35]);grid on
set(gca, 'XScale', 'log')
subplot(222)
scatter(DATA.PND_c(D,1),DATA.AT_RH(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('RH (AT)','interpreter','latex')
xlim([1e3 7e5]);ylim([10 50]);grid on
set(gca, 'XScale', 'log')

subplot(223)
scatter(DATA.PND_c(D,1),DATA.CO_P(D,1));
xlabel('PNC (CPC)','interpreter','latex')
ylabel('P (CO)','interpreter','latex')
xlim([1e3 7e5]);ylim([895 901]);grid on
set(gca, 'XScale', 'log')

subplot(224)
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




%% MATRIX PLOT
%close all;
clc

DATAx = [DATA.PMD_c(:,1),DATA.PMD_c(:,6),DATA.DustTrak_c(:,2),DATA.SidePak_c(:,1), ...
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

clc
figure(9);fig=gcf; FS=26;

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




%% %% HISTOGRAMS
%close all; 
clc
figure(10); fig = gcf;
fig.Position = [100 100 540 400].*2.5;
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

disp('For normalization of temp, RH and P, perhaps the best is to use max and min number')

%% Using min-max algorithm
% https://www.mathworks.com/help/matlab/ref/rescale.html
% B = rescale(A,l,u,'InputMin',inmin,'InputMax',inmax) uses the formula
% B = l + [(A-inmin)./(inmax-inmin)].*(u-l)



%close all; clc
figure(11); fig = gcf;
fig.Position = [100 100 540 400].*2.5;
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



%% MODELLING: Linear regression and ANN
% 1) use all data, randominze, and predict CPC
% 2) use Exp_smoking, estimate Exp_kerosine and Exp_gas

D = Exp_smoking; 
Dt = Exp_kerosine;
%Dt = Exp_gas;
%D =[Exp_smokimg;Exp_kerosine;Exp_gas];

% For met vars
% B = rescale(A,l,u,'InputMin',inmin,'InputMax',inmax) uses the formula
% B = l + [(A-inmin)./(inmax-inmin)].*(u-l)
% For aerosol, simply take log10

DATAm1 = [DATA.AT_T(D,1),DATA.LCS_G2_01_met(D,1),DATA.LCS_G2_01_met(D,3),DATA.PMD_c(D,6),DATA.PND_c(D,1)];
%DATAm1 = [DATA.AT_T(D,1),DATA.LCS_G2_01_met(D,1),DATA.CO_P(D,1),DATA.LCS_G1(D,1),DATA.PND_c(D,1)];
%DATAm1 = [DATA.AT_T(D,1),DATA.AT_RH(D,1),DATA.CO_P(D,1),DATA.LCS_G2_02(D,1),DATA.PND_c(D,1)];
DATAm2 = zeros(size(DATAm1));

DATAt1 = [DATA.AT_T(Dt,1),DATA.AT_RH(Dt,1),DATA.CO_P(Dt,1),DATA.LCS_G1(Dt,1),DATA.PND_c(Dt,1)];
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
        DATAm2(:,n) = l + [(DATAm1(:,n)-inmin)./(inmax-inmin)].*(u-l);
        DATAt2(:,n) = l + [(DATAt1(:,n)-inmin)./(inmax-inmin)].*(u-l);
    else
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

% ANN model
inputs = X';
targets = Y';
inputs_t = Xt';
% Create a Fitting Network
hiddenLayerSize = 15;
net = fitnet(hiddenLayerSize);
% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
 % Train the Network
[net,tr] = train(net,inputs,targets);
 % Test the Network
outputs = net(inputs_t);
Ypred = outputs;


figure(12);
subplot(221);
scatter(Yt,Ypred_lm);hold on
Xlim1 = 3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')

subplot(222);
scatter(10.^Yt,10.^Ypred_lm);hold on
Xlim1 = 1e1;
Ylim1 = 5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('Real PNC (CPC)');ylabel('Est PNC (CPC)')

subplot(223);
scatter(Yt,Ypred);hold on
Xlim1 = 3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')

subplot(224);
scatter(10.^Yt,10.^Ypred);hold on
Xlim1 = 1e1;
Ylim1 = 5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('Real PNC (CPC)');ylabel('Est PNC (CPC)')


%% MODELLING Linear Models and Shallow Neural Networks

% 1) use all data, randominze, and predict CPC
% 2) use Exp_smoking, estimate Exp_kerosine and Exp_gas

Ds = Exp_smoking; 
Dk = Exp_kerosine;
Dg = Exp_gas;
Da =[Exp_smoking;Exp_kerosine;Exp_gas];

% CHOOSE THE DATA with the number of inputs and output
Di = 4 ;
if Di == 2
    disp('Temp and PM2.5')
    DATAi = [ [DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
              [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 3
    disp('Temp, RH and PM2.5')
    DATAi = [[DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
             [DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 4
    disp('Temp, RH, P and PM2.5')
    DATAi = [[DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
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

%
if Di == 4
    figure(13)
    for n=1:4
        subplot(2,5,n);plot(DATAi(:,n),'.');hold on
        subplot(2,5,n+5);plot(DATAi1(:,n),'.');
    end
    subplot(2,5,5);plot(DATAo(:,1),'.');
    subplot(2,5,10);plot(DATAo1(:,1),'.');
    hold off
    
    %
    %%CPC = DATA.PND_c([Ds;Dk;Dg],1);
    CPC = DATA.PND_c([Ds;Dk;Dg],1);
    CPClog = log10(CPC);
    CPCgradient = gradient(CPC);
    CPCdiff = diff(CPC);
    CPCloggradient = gradient(CPClog);
    idx = find(CPCgradient <= nanmedian(CPCgradient)+10);
    CPCclean =CPC;
    CPCclean(idx,:)=nan;
    %CPCclean = CPC(idx,:);
    
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
end

DATAo  = CPCclean;
%DATAo  = [DATA.PND_c(Da,1)];
DATAo1 = log10(DATAo);

% SELECT TRAINING AND TESTING DATA

clc
method = 4;
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
    Xtr = DATAi1(Ds,:);
    Ytr = DATAo1(Ds,:);
    Xte = DATAi1(Ds,:);
    Yte = DATAo1(Ds,:); 
    
elseif method == 3
    disp('Training and Testing data only for smoking and kerosene experiments')
    Xtr = DATAi1([Ds:Dk],:);
    Ytr = DATAo1([Ds:Dk],:);
    Xte = DATAi1([Ds:Dk],:);
    Yte = DATAo1([Ds:Dk],:); 

elseif method == 4
    disp('Training and Testing data only for smoking, kerosene and gas experiments')  
    Xtr = DATAi1([Ds:Dk:Dg],:);
    Ytr = DATAo1([Ds:Dk:Dg],:);
    Xte = DATAi1([Ds:Dk:Dg],:);
    Yte = DATAo1([Ds:Dk:Dg],:); 
    
elseif method == 5
    disp('Training and Testing data are randomized')
    
    DATAt  = [DATAi1,DATAo1];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    percent = 70/100;
    p = size(DATAt1,1);
    c = randperm(p)';
    tr = c( 1 : roundn(percent * p,0)     , 1);
    te = c( roundn(percent * p,0)+1 : end , 1);
    
    Xtr = DATAt1(tr,1:end-1);
    Ytr = DATAt1(tr,end);
    Xte = DATAt1(te,1:end-1);
    Yte = DATAt1(te,end); 
end


% Linear model
X = Xtr;
Y = Ytr;

Xt = Xte;
Yt = Yte;


% LINEAR MODEL:
mdl = fitlm(X,Y);
Ypred_lm = predict(mdl,Xt);

% ANN model
inputs = X';
targets = Y';
inputs_t = Xt';
% Create a Fitting Network
hiddenLayerSize = 20;%15;
net = fitnet(hiddenLayerSize);
% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
 % Train the Network
[net,tr] = train(net,inputs,targets);
 % Test the Network
outputs = net(inputs_t);
Ypred_snn = outputs;

figure(12);
subplot(221);
scatter(Yt,Ypred_lm);hold on
Xlim1 = 3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')

subplot(222);
scatter(10.^Yt,10.^Ypred_lm);hold on
Xlim1 = 1e0;
Ylim1 = 7.5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('Real PNC (CPC)');ylabel('Est PNC (CPC)')

subplot(223);
scatter(Yt,Ypred_snn);hold on
Xlim1 = 3;
Ylim1 = 6;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')

subplot(224);
scatter(10.^Yt,10.^Ypred_snn);hold on
Xlim1 = 1e0;
Ylim1 = 7.5e5;
xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
x = linspace(Xlim1,Ylim1);
y = linspace(Xlim1,Ylim1);
plot(x,y,'r');hold off
xlabel('Real PNC (CPC)');ylabel('Est PNC (CPC)')



%% THE BEST SOLUTIONS SO FAR ARE:
% BY USING Smoking experiment only
% (i) it seems that the upper regime is difficult to predict,
% (ii) there may be a problem with gas experiments, check inputs and
% outputs again for gas experiment because it give warning:
% "Warning: Regression design matrix is rank deficient to within machine
% precision."
% solution: https://www.mathworks.com/matlabcentral/answers/637610-regression-design-matrix-is-rank-deficient-what-to-do-next
% ANOTHER SOLUTION:
% We may need to establish a model, such as MoE, where aerosol models are 
% divided into different regimes


%%
clc
figure(2);
subplot(221)
scatter(DATA.AT_T,DATA.CO_T)
subplot(222)
scatter(DATA.CO_T,DATA.LCS_G2_01_met(:,2))
subplot(223)
scatter(DATA.AT_T,DATA.LCS_G2_01_met(:,2))
subplot(224)
scatter(DATA.AT_T,DATA.LCS_G2_02_met(:,2))

figure(3);
subplot(221)
scatter(DATA.AT_RH,DATA.CO_RH)
subplot(222)
scatter(DATA.CO_RH,DATA.LCS_G2_01_met(:,1))
subplot(223)
scatter(DATA.AT_RH,DATA.LCS_G2_01_met(:,1))
subplot(224)
scatter(DATA.AT_RH,DATA.LCS_G2_02_met(:,1))

%%
subplot(222)
scatter(DATA.AT_RH,DATA.CO_RH)


%%
subplot(212)
plot(DATA.LCS_G1(:,1),'b','LineWidth',2); hold on
plot(DATA.LCS_G2_01(:,1),'y','LineWidth',2); 
plot(DATA.LCS_G2_02(:,1),'c','LineWidth',2);
xlabel('Day of the year','interpreter','latex')
ylabel('PM$_{2.5}$ concentration','interpreter','latex')
title('Low-Cost Sensors PM$_{2.5}$','interpreter','latex');grid on
legend('$\mathcal{L}_1$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
set(gca, 'YScale', 'log')
ylim([1e-2 1e4])
colorbar
hold off
set(findall(fig,'-property','FontSize'),'FontSize',22);



%%
figure(2);
plot(DATA.T1,DATA.AT_T,'b.');hold on
plot(DATA.T1,DATA.CO_T,'r.');
plot(DATA.T1,DATA.LCS_G2_01_met(:,2),'g.');
plot(DATA.T1,DATA.LCS_G2_02_met(:,2),'c.');
ylim([10 40])
hold off
%plot(DATA.AT_RH)
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


