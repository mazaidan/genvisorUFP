DATA = DATA2;

% Aerosol size distribution

Exp_smoking = [1:1:11521]';
Exp_kerosine = [11521:1:30242]';
Exp_gas = [30243:1:54721]';

inst_res_PND = [0.025,0.3,0.5,1,2.5,5,10,45].*1e3;
Time_index   = [1:1:size(DATA.T1,1)]';

% AFTER SYNCRONIZATION (DATA2)
figure(6);fig=gcf; 
fig.Position = [100 100 540 400].*2.5;
FS=26; PMx=6;
lc =1e0; hc=10e5;
tiledlayout(2,2);
nexttile % subplot(231)
plot(Time_index(Exp_smoking,1),DATA.LCS_G1(Exp_smoking,1),'b.');
%plot(Time_index(Exp_smoking,1),DATA.PMD_c(Exp_smoking,PMx),'b.'); hold on
%plot(Time_index(Exp_smoking,1),DATA.LCS_G1(Exp_smoking,1),'g.');
%plot(Time_index(Exp_smoking,1),DATA.LCS_G2_01(Exp_smoking,1),'r.');
%plot(Time_index(Exp_smoking,1),DATA.LCS_G2_02(Exp_smoking,1),'m.');
hold off
%legend('PMD','$\mathcal{L}_{1}$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$','interpreter','latex')
legend('PM$_{2.5}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Smoking','interpreter','latex')
set(gca, 'YScale', 'log')


nexttile % subplot(222)
plot(Time_index(Exp_gas,1),DATA.LCS_G1(Exp_gas,1),'b.');
%plot(Time_index(Exp_gas,1),DATA.PMD_c(Exp_gas,PMx),'b.'); hold on
%plot(Time_index(Exp_gas,1),DATA.LCS_G1(Exp_gas,1),'g.');
%plot(Time_index(Exp_gas,1),DATA.LCS_G2_01(Exp_gas,1),'r.');
%plot(Time_index(Exp_gas,1),DATA.LCS_G2_02(Exp_gas,1),'m.');
hold off
legend('PM$_{2.5}$','interpreter','latex')
ylim([1e-1 1e3]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Natureal Gas','interpreter','latex')
set(gca, 'YScale', 'log')

nexttile % subplot(223)
h1 = pcolor(Time_index(Exp_smoking,1)',inst_res_PND',[DATA1.PND_c(Exp_smoking,1),DATA1.PND_c(Exp_smoking,3:end)]');
set(h1, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PND Smoking','interpreter','latex')
colormap jet
caxis([lc,hc])
set(gca, 'YScale', 'log','colorscale','log')

nexttile % subplot(224)
h3 = pcolor(Time_index(Exp_gas,1)',inst_res_PND',[DATA1.PND_c(Exp_gas,1),DATA1.PND_c(Exp_gas,3:end)]');
set(h3, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PND natural gas','interpreter','latex')
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
%text( 0.5, 0, 'After syncronize', 'FontSize', 14', 'FontWeight', 'Bold', ...
%      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

set(findall(fig,'-property','FontSize'),'FontSize',FS);

%% Tareq plot's data:

clear;clc;close all;

% SMOKING ACTIVITIES

load('Data_clean_processed.mat');

FS =26;
figure(1); fig = gcf;
%Y = PM2p5(:,8); 
Y = SidePak(:,8);
X = ISEE_LCS_G202(:,end); 
%X = PMSD(1:end-1,10);

MINx = 1e-2; MAXx = 1e5; MINy = 1e-2; MAXy = 1e5;  
%MINx = 0;0.01*min(X); MAXx = 100*max(X); MINy = 0;0.01*min(Y); MAXy = 100*max(Y);
scatter(X,Y);hold on
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('LCSs','interpreter','latex')
ylabel('SidePak','interpreter','latex')
xlim([MINx MAXx]);ylim([MINy MAXy]);
x = linspace(MINx,MAXx,1000);
y = linspace(MINy,MAXy,1000);
plot(x,y,'r');
set(findall(fig,'-property','FontSize'),'FontSize',FS);


%% Plot PNSD

clear;clc;close all;

% SMOKING ACTIVITIES

load('Data_clean_processed.mat');

%Exp_smoking = [1:1:11521]';
%Exp_kerosine = [11521:1:30242]';
%Exp_gas = [30243:1:54721]';

%inst_res_PND = [0.025,0.3,0.5,1,2.5,5,10,45].*1e3;
inst_res_PND = [Dp_min(1,1),Dp_max];
Time_index   = [1:1:size(PNSD,1)]';

% AFTER SYNCRONIZATION (DATA2)
figure(6);fig=gcf; 
fig.Position = [120 120 690 450].*2.5;
FS=26; PMx=6;
lc =1e0; hc=10e5;
tiledlayout(2,1);
nexttile % subplot(231)
plot(Time_index(2:end,1),DustTrak(:,9),'b.');
hold off
legend('PM$_{2.5}$','interpreter','latex')
ylim([1e-1 1e5]); grid on
xlabel('Time Index','interpreter','latex')
ylabel('PM$_{2.5}$ ($\mu$g/m$^3$)','interpreter','latex')
title('PM$_{2.5}$ Smoking','interpreter','latex')
set(gca, 'YScale', 'log')


nexttile % subplot(224)
%h3 = pcolor(Time_index(:,1)',inst_res_PND',[DATA1.PND_c(Exp_gas,1),DATA1.PND_c(Exp_gas,3:end)]');
h3 = pcolor(Time_index(:,1)',inst_res_PND',[PNSD(:,8:end)]');
set(h3, 'EdgeColor', 'none')
xlabel('Time Index','interpreter','latex')
ylabel('Particle size (nm)','interpreter','latex')
title('PNSD smoking','interpreter','latex')
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

set(findall(fig,'-property','FontSize'),'FontSize',FS);

%% Tareq plot's data:

clear;clc;close all;

% SMOKING ACTIVITIES

load('Data_clean_processed.mat');

FS =26;
figure(1); fig = gcf;
%Y = PM2p5(:,8); 
Y = SidePak(:,8);
X = ISEE_LCS_G202(:,end); 
%X = PMSD(1:end-1,10);

MINx = 1e-2; MAXx = 1e5; MINy = 1e-2; MAXy = 1e5;  
%MINx = 0;0.01*min(X); MAXx = 100*max(X); MINy = 0;0.01*min(Y); MAXy = 100*max(Y);
scatter(X,Y);hold on
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('LCSs','interpreter','latex')
ylabel('SidePak','interpreter','latex')
xlim([MINx MAXx]);ylim([MINy MAXy]);
x = linspace(MINx,MAXx,1000);
y = linspace(MINy,MAXy,1000);
plot(x,y,'r');
set(findall(fig,'-property','FontSize'),'FontSize',FS);


%% Tareq's computations for derived variables 

clear;clc;close all;

% SMOKING ACTIVITIES

load('Data_clean_processed.mat');

% How to calculate PND

figure(1); fig=gcf;
X = CPC(:,8); % CPC
Y = PND(2:end,8);
subplot(421);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('CPC');
title(['Bin 1: ',num2str(Dp_min(1,1)),' - ',num2str(Dp_max(1,1)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = CPC(:,8) - Ptrak(:,8); % CPC - PTRAK == 10 - 25 nm
Y = PND(2:end,9);
subplot(422);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('CPC - Ptrak'); %title('Bin 2: PND')
title(['Bin 2: ',num2str(Dp_min(1,2)),' - ',num2str(Dp_max(1,2)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = Ptrak(:,8) - sum(AeroTrak_PN(:,8:13),2); % Ptrak - sum(aerotrak) == 25 - 300 nm
Y = PND(2:end,10);
subplot(423);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('Ptrak - sum(aerotrak)');%title('Bin 3: PND')
title(['Bin 3: ',num2str(Dp_min(1,3)),' - ',num2str(Dp_max(1,3)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = AeroTrak_PN(:,8); 
Y = PND(2:end,11);
subplot(424);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('AeoTrak');%title('Bin 4: PND')
title(['Bin 4: ',num2str(Dp_min(1,4)),' - ',num2str(Dp_max(1,4)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = AeroTrak_PN(:,9); 
Y = PND(2:end,12);
subplot(425);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('AeoTrak');%title('Bin 5: PND')
title(['Bin 5: ',num2str(Dp_min(1,5)),' - ',num2str(Dp_max(1,5)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = AeroTrak_PN(:,10); 
Y = PND(2:end,13);
subplot(426);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('AeoTrak');%title('Bin 6: PND')
title(['Bin 6: ',num2str(Dp_min(1,6)),' - ',num2str(Dp_max(1,6)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = AeroTrak_PN(:,11); 
Y = PND(2:end,14);
subplot(427);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('AeoTrak');%title('Bin 7: PND')
title(['Bin 7: ',num2str(Dp_min(1,7)),' - ',num2str(Dp_max(1,7)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

X = AeroTrak_PN(:,12); 
Y = PND(2:end,15);
subplot(428);scatter(X,Y); hold on; grid on
xlabel('PND');ylabel('AeoTrak');%title('Bin 8: PND')
title(['Bin 8: ',num2str(Dp_min(1,8)),' - ',num2str(Dp_max(1,8)),' m'])
x = linspace(0,max(Y));
y = linspace(0,max(Y));
plot(x,y,'r');hold off

set(findall(fig,'-property','FontSize'),'FontSize',22);

% From PND, how to calculate PNSD:
figure(2); fig = gcf;
for n= 1:8
    subplot(4,2,n);
    PNSD1 = (PND(:,n+8))./log10(Dp_max(1,n)/Dp_min(1,n));
    %PNSD1 = abs(PND(:,n+7) - PND(:,n+8))./log10(Dp_max(1,n)/Dp_min(1,n));
    scatter(PNSD(:,n+8),PNSD1); hold on; grid on
    x = linspace(0,max(PNSD1));
    y = linspace(0,max(PNSD1));
    plot(x,y,'r');hold off
    xlabel('PNSD');ylabel('PNSD1');
    title(['Bin ',num2str(n), ' : ',num2str(Dp_min(1,n)*1e6),' - ',num2str(Dp_max(1,n)*1e6),'$\mu$m'], 'interpreter','latex')
end
set(findall(fig,'-property','FontSize'),'FontSize',22);

%% Computing PMD

PN = pi * (PND(:,8).*(Dp_gmd(1,1)^3)/6)*rho_unit;

% From PND, how to calculate PNSD:

for n= 1:8
    
    figure(3); fig = gcf;
    fig.Position = [120 120 690 690].*1.5;
    subplot(4,2,n);
    PMD1 = 1e6*PND(:,n+8) .* (pi/6)*((Dp_gmd(1,n))^3)*(1e9.*rho_unit);
    scatter(PMD(:,n+8),PMD1); hold on
    x = linspace(0,max(PMD1));
    y = linspace(0,max(PMD1));
    plot(x,y,'r');hold off
    xlabel('PMSD');ylabel('PMSD1');
    title(['Bin ',num2str(n), ' : ',num2str(Dp_min(1,n)*1e6),' - ',num2str(Dp_max(1,n)*1e6),' $\mu$m'], 'interpreter','latex')
    set(findall(fig,'-property','FontSize'),'FontSize',22);
    
    figure(4); fig = gcf;
    fig.Position = [120 120 690 690].*1.5;
    subplot(4,2,n);
    PMSD1 = 1e6*PNSD(:,n+8) .* (pi/6)*((Dp_gmd(1,n))^3)*(1e9.*rho_unit);
    scatter(PMSD(:,n+8),PMSD1); hold on
    x = linspace(0,max(PMSD1));
    y = linspace(0,max(PMSD1));
    plot(x,y,'r');hold off
    xlabel('PMSD');ylabel('PMSD1');
    title(['Bin ',num2str(n), ' : ',num2str(Dp_min(1,n)*1e6),' - ',num2str(Dp_max(1,n)*1e6),' $\mu$m'], 'interpreter','latex')
    set(findall(fig,'-property','FontSize'),'FontSize',22);
    
end

%% Calculating PM2.5 from PMD
clc
Dp = [Dp_min',Dp_gmd',Dp_max']

figure(5); fig = gcf;
title('Calculating PM$_{2.5}$','interpreter','latex');
PM25b  =[];
for n = 1:5
    PM25a = PMSD(2:end,n+8) * log10(Dp_max(1,n)/Dp_min(1,n));
    PM25b = [PM25b,PM25a];
end
PM25 = sum(PM25b,2);
scatter(PM2p5(:,8),PM25); grid on
hold on;
x = linspace(0,max(PM25));
y = linspace(0,max(PM25));
plot(x,y,'r');hold off
xlabel('PM$_{2.5}$ [Data]','interpreter','latex');
ylabel('PM$_{2.5}$ [Recalculate]','interpreter','latex');
set(findall(fig,'-property','FontSize'),'FontSize',22);


%%
Dp
figure(6);
PM25b  =[];
for n = 5%6:7
    %PM25a = PMSD(2:end,n+1) * log10(Dp_max(1,n)/Dp_min(1,n));
    PM25a = PMSD(2:end,n+8)*log10(Dp_max(1,n)/Dp_min(1,n));
    PM25b = [PM25b,PM25a];
end
%PM25 = sum(PM25b,2);
%PM25 = abs(PM25b(:,2) - PM25b(:,1));
%PM25 = sum(PMD(2:end,9:end),2);
PM25 = PM25b;
scatter(PM25,PM2p5(:,8));





%% CALIBRATIONS LCS PM2.5

addpath(genpath('Functions'));
addpath(genpath('Functions_special/BayesianNeuralNetworks'));
rmpath('Functions_special/BayesianNeuralNetworks/netlab/Garbages');

clear; close all; clc;

load('DATA2.mat') ;

Exp_smoking = [1:1:11521]';
Exp_kerosine = [11522:1:30242]';
Exp_gas = [30243:1:54721]';

%DATA = DATA2([Exp_smoking;Exp_gas],:);
DATAs2 = DATA2(Exp_smoking,:);
%DATAn2 = DATA2(Exp_gas,:);


FS =20;
figure(2); fig = gcf;
PMx = 6;
%Y = SidePak(:,8);
Y = DATAs2.DustTrak_c(:,2); %DATAs2.SidePak_c(:,1);%DATAs2.PMD_c(:,PMx); %DATAs2.DustTrak_c(:,2);  % 
%X = ISEE_LCS_G202(:,end); 
X = DATAs2.LCS_G2_02(:,1); Temp =DATAs2.LCS_G2_02_met(:,2);% 
index = linspace(1,size(X,1),size(X,1))';

O = rmmissing([X,Y,Temp,index]);
X = O(:,1);
Y = O(:,2);
Temp = O(:,3);
index = O(:,4);


MINx = 1e-2; MAXx = 1e5; MINy = 1e-2; MAXy = 1e5;  
%MINx = 0;0.01*min(X); MAXx = 100*max(X); MINy = 0;0.01*min(Y); MAXy = 100*max(Y);
scatter(X,Y,'r');hold on; grid on;
set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');
xlabel('PM$_{2.5}$ [$\mu$g/m$^3$] (DustTrak)','interpreter','latex');
ylabel('PM$_{2.5}$ [$\mu$g/m$^3$] ($\mathcal{L}_{2a}$)','interpreter','latex');
xlim([MINx MAXx]);ylim([MINy MAXy]);

Xlog = abs(log10(X));
Ylog = abs(log10(Y));

idx = isinf(Xlog);
Xlog(idx) = [];
Ylog(idx) = [];
Temp(idx) = [];
index(idx) = [];

idx = isinf(Ylog);
Xlog(idx) = [];
Ylog(idx) = [];
Temp(idx) = [];
index(idx) = [];

Te =1;
if Te ==0
    mdl = fitlm(Xlog,Ylog);
    Ylm = predict(mdl,Xlog);
else
    mdl = fitlm([Xlog, Temp],Ylog);
    Ylm = predict(mdl,[Xlog, Temp]);
end
scatter(10.^Ylog,10.^Ylm,'g.');

x = linspace(MINx,MAXx,1000);
y = linspace(MINy,MAXy,1000);
plot(x,y,'r');

legend('Before calibration','After calibration','interpreter','latex')

hold off

set(findall(fig,'-property','FontSize'),'FontSize',FS);

anova(mdl,'summary')

figure(3); fig = gcf;
plot(DATAs2.T1(index),10.^Ylog,'b.');
hold on; grid on;
plot(DATAs2.T1(index),10.^Xlog,'r.');
plot(DATAs2.T1(index),10.^Ylm,'g'); 
ylabel('PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
legend('Reference instrument','$\mathcal{L}_{2a}$ before calibration', ...
    '$\mathcal{L}_{2a}$ after calibration','interpreter','latex');
set(gca, 'YScale', 'log');
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);


%% CALIBRATIONS LCS PM2.5 (lopping on different calibrators)
addpath(genpath('Functions'));

clear; close all; clc;

load('DATA2.mat') ;

Exp_smoking = [1:1:11521]';
Exp_kerosine = [11522:1:30242]';
Exp_gas = [30243:1:54721]';

%DATA = DATA2([Exp_smoking;Exp_gas],:);
DATAs2 = DATA2(Exp_smoking,:);
%DATAn2 = DATA2(Exp_gas,:);


FS =20;
figure(2); fig = gcf;
PMx = 6;
X = DATAs2.LCS_G2_02(:,1); 
Y = DATAs2.DustTrak_c(:,2); %DATAs2.SidePak_c(:,1);%DATAs2.PMD_c(:,PMx); %DATAs2.DustTrak_c(:,2);  % 
Temp = DATAs2.LCS_G2_02_met(:,2);
RH   = DATAs2.LCS_G2_02_met(:,1);
P    = DATAs2.LCS_G2_02_met(:,3);
index = linspace(1,size(X,1),size(X,1))';


%%%%%%%%%%%%%%%%%%%

O = rmmissing([X,Y,Temp,RH,P,index]);
X     = O(:,1);
Y     = O(:,2);
Temp  = O(:,3);
RH    = O(:,4);
P     = O(:,5);
index = O(:,6);

Xlog = abs(log10(X));
Ylog = abs(log10(Y));

idx = isinf(Xlog);
Xlog(idx)  = [];
Ylog(idx)  = [];
Temp(idx)  = [];
RH(idx)    = [];
P(idx)     = [];
index(idx) = [];

idx = isinf(Ylog);
Xlog(idx) = [];
Ylog(idx) = [];
Temp(idx) = [];
RH(idx)    = [];
P(idx)     = [];
index(idx) = [];

for Cal = 0:4
    
%Cal =2;
if Cal == 0
    Ylm = Xlog;
elseif Cal == 1
    mdl = fitlm(Xlog,Ylog);
    Ylm = predict(mdl,Xlog);
elseif Cal == 2
    mdl = fitlm([Xlog, Temp],Ylog);
    Ylm = predict(mdl,[Xlog, Temp]);
elseif Cal == 3
    mdl = fitlm([Xlog, Temp, RH],Ylog);
    Ylm = predict(mdl,[Xlog, Temp, RH]);
elseif Cal == 4
    mdl = fitlm([Xlog, Temp, RH, P],Ylog);
    Ylm = predict(mdl,[Xlog, Temp, RH, P]);
end

Rp    = corr(Ylog,Ylm,'Type','Pearson','Rows','complete');
MAPE  = errperf(Ylog,Ylm,'mape');
MAE   = errperf(10.^Ylog,10.^Ylm,'mae');

Rp1(1,Cal+1) = Rp;
MAPE1(1,Cal+1) = MAPE;
MAE1(1,Cal+1) = MAE;

end
%%%%%%

%

MINx = 1e-2; MAXx = 1e5; MINy = 1e-2; MAXy = 1e5;  
scatter(X,Y,'r');hold on; grid on;
set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');
xlabel('PM$_{2.5}$ [$\mu$g/m$^3$] (DustTrak)','interpreter','latex');
ylabel('PM$_{2.5}$ [$\mu$g/m$^3$] ($\mathcal{L}_{2a}$)','interpreter','latex');
xlim([MINx MAXx]);ylim([MINy MAXy]);

scatter(10.^Ylog,10.^Ylm,'g.');

x = linspace(MINx,MAXx,1000);
y = linspace(MINy,MAXy,1000);
plot(x,y,'r');

legend('Before calibration','After calibration','interpreter','latex')

hold off

set(findall(fig,'-property','FontSize'),'FontSize',FS);

anova(mdl,'summary')

figure(3); fig = gcf;
plot(DATAs2.T1(index),10.^Ylog,'b.');
hold on; grid on;
plot(DATAs2.T1(index),10.^Xlog,'r.');
plot(DATAs2.T1(index),10.^Ylm,'g'); 
ylabel('PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
legend('Reference instrument','$\mathcal{L}_{2a}$ before calibration', ...
    '$\mathcal{L}_{2a}$ after calibration','interpreter','latex');
set(gca, 'YScale', 'log');
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);



%%
%% CALIBRATIONS LCS PM2.5 (lopping on different scenarios and 
%% different calibrators)
addpath(genpath('Functions'));

clear; close all; clc;

load('DATA2.mat') ;

Exp_smoking = [1:1:11521]';
Exp_kerosine = [11522:1:30242]';
Exp_gas = [30243:1:54721]';

%DATA = DATA2([Exp_smoking;Exp_gas],:);
DATAs2 = DATA2(Exp_smoking,:);
%DATAn2 = DATA2(Exp_gas,:);


FS =20;

%PMx = 6;
R      = DATAs2.DustTrak_c(:,2);
Rl     = abs(log10(R));

A1     = DATAs2.LCS_G2_01(:,1);
A1l    = abs(log10(A1));
T1     = DATAs2.LCS_G2_01_met(:,2);
RH1    = DATAs2.LCS_G2_01_met(:,1);
P1     = DATAs2.LCS_G2_01_met(:,3);
index1 = linspace(1,size(A1,1),size(A1,1))';

A2     = DATAs2.LCS_G2_01(:,1);
A2l    = abs(log10(A2));
T2     = DATAs2.LCS_G2_02_met(:,2);
RH2    = DATAs2.LCS_G2_02_met(:,1);
P2     = DATAs2.LCS_G2_02_met(:,3);
index2 = linspace(1,size(A2,1),size(A2,1))';

%%%%%%%%%%%%%%%%%%%

%D  = [R,Rl,A1,A1l,A2,A2l,T1,T2,RH1,RH2,P1,P2,index1,index2];
D  = [R,Rl, A1,A1l,T1,RH1,P1,index1, A2,A2l,T2,RH2,P2,index2]; 
D1 = rmmissing(D);

idx = isinf(D1);
idx1 = sum(idx,2);
idx2 = idx1>0;
D1(idx2,:) = [];

%S = 1;
for S = 1:4
if S == 1
    Y  = D1(:,2);
    Y1 = Y;
    X  = D1(:,[3:8]);  % LCS2a
    X1 = D1(:,[9:14]); % LCS2b
elseif S == 2
    Y  = D1(:,2);
    Y1 = Y;
    X  = D1(:,[9:14]);  % LCS2b
    X1 = D1(:,[3:8]); % LCS2a
elseif S == 3
    th = roundn(0.7 .* size(D1,1),0);
    RandIdx = randperm(size(D1,1))';
    Tr = RandIdx(1:th,1);
    Te = RandIdx(th+1:end,1);
    Y  = D1(Tr,2);
    Y1 = D1(Te,2);
    X  = D1(Tr,[3:8]);   % LCS2a
    X1 = D1(Te,[9:14]);  % LCS2b
elseif S == 4
    th = roundn(0.7 .* size(D1,1),0);
    RandIdx = randperm(size(D1,1))';
    Tr = RandIdx(1:th,1);
    Te = RandIdx(th+1:end,1);
    Y  = D1(Tr,2);
    Y1 = D1(Te,2);
    X   = D1(Tr,[9:14]);  % LCS2b
    X1  = D1(Te,[3:8]);   % LCS2a
else
    error('Please select S Scenario')
end

for Cal = 0:2%0:4
    
if Cal == 0
    Ylm = X1(:,2);
elseif Cal == 1
    mdl = fitlm(X(:,2),Y);
    Ylm = predict(mdl,X1(:,2));
elseif Cal == 2
    mdl = fitlm(X(:,[2:3]),Y);
    Ylm = predict(mdl,X1(:,[2:3]));
elseif Cal == 3
    mdl = fitlm(X(:,[2:4]),Y);
    Ylm = predict(mdl,X1(:,[2:4]));
elseif Cal == 4
    mdl = fitlm(X(:,[2:5]),Y);
    Ylm = predict(mdl,X1(:,[2:5]));
end

Rp    = corr(Y1,Ylm,'Type','Pearson','Rows','complete');
MAPE  = errperf(Y1,Ylm,'mape');
MAE   = errperf(10.^Y1,10.^Ylm,'mae');

Rp1(S,Cal+1)   = roundn(Rp,-2);
MAPE1(S,Cal+1) = roundn(MAPE,-2);
MAE1(S,Cal+1)  = roundn(MAE,-2);

end
end
%%%%%%

index = X1(:,6);

figure(4); fig = gcf;
%plot(DATAs2.T1(index),10.^Ylog,'b.');
plot(DATAs2.T1(index),DATAs2.LCS_G2_01(index,1),'b.');
hold on; grid on;
%%plot(DATAs2.T1(index),10.^Xlog,'r.');
plot(DATAs2.T1(index),10.^Ylm,'g.'); hold off
ylabel('PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
%legend('Reference instrument','$\mathcal{L}_{2a}$ before calibration', ...
%    '$\mathcal{L}_{2a}$ after calibration','interpreter','latex');
legend('$\mathcal{L}_{2a}$ before calibration', ...
    '$\mathcal{L}_{2a}$ after calibration','interpreter','latex');
set(gca, 'YScale', 'log');
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);

figure(5); fig = gcf;
scatter(DATAs2.LCS_G2_01(index,1),10.^Ylm);
xlabel('Uncalibrated PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
ylabel('Calibrated PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([1 1e3]);ylim([1 1e4])
set(findall(fig,'-property','FontSize'),'FontSize',FS);


%% To test the developed calibration on natural gas heater data

DATAs2 = DATA2(Exp_gas,:);

%PMx = 6;
R      = DATAs2.DustTrak_c(:,2);
Rl     = abs(log10(R));

A1     = DATAs2.LCS_G2_01(:,1);
A1l    = abs(log10(A1));
T1     = DATAs2.LCS_G2_01_met(:,2);
RH1    = DATAs2.LCS_G2_01_met(:,1);
P1     = DATAs2.LCS_G2_01_met(:,3);
index1 = linspace(1,size(A1,1),size(A1,1))';

A2     = DATAs2.LCS_G2_01(:,1);
A2l    = abs(log10(A2));
T2     = DATAs2.LCS_G2_02_met(:,2);
RH2    = DATAs2.LCS_G2_02_met(:,1);
P2     = DATAs2.LCS_G2_02_met(:,3);
index2 = linspace(1,size(A2,1),size(A2,1))';

%%%%%%%%%%%%%%%%%%%

%D  = [R,Rl, A1,A1l,T1,RH1,P1,index1, A2,A2l,T2,RH2,P2,index2];
D0  = [A1,A1l,T1,RH1,P1,index1, A2,A2l,T2,RH2,P2,index2];
D00 = rmmissing(D0);
D1 = [nan(size(D00,1),2),D00];

idx = isinf(D1);
idx1 = sum(idx,2);
idx2 = idx1>0;
D1(idx2,:) = [];

X1  = D1(:,[3:8]);   % LCS2a
X1  = D1(:,[9:14]);  % LCS2b

Ylm = predict(mdl,X1(:,[2:3])); % only PM2.5 and Temp 


index = X1(:,6);

figure(6); fig = gcf;
%plot(DATAs2.T1,DATAs2.LCS_G1(:,1),'b.');
plot(DATAs2.T1,DATAs2.LCS_G2_02(:,1),'b.');
%%%plot(DATAs2.T1(index),DATAs2.LCS_G2_01(index,1),'b.');
hold on; grid on;
%%plot(DATAs2.T1(index),10.^Xlog,'r.');
plot(DATAs2.T1(index),10.^Ylm,'g.'); hold off
ylabel('PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
%legend('Reference instrument','$\mathcal{L}_{2a}$ before calibration', ...
%    '$\mathcal{L}_{2a}$ after calibration','interpreter','latex');
legend('$\mathcal{L}_{2a}$ before calibration', ...
    '$\mathcal{L}_{2a}$ after calibration','interpreter','latex');
set(gca, 'YScale', 'log');
hold off
set(findall(fig,'-property','FontSize'),'FontSize',FS);

figure(7); fig = gcf;
%scatter(DATAs2.LCS_G2_01(index,1),10.^Ylm);
scatter(DATAs2.LCS_G2_02(index,1),10.^Ylm);
xlabel('Uncalibrated PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
ylabel('Calibrated PM$_{2.5}$ [$\mu$g/m$^3$]','interpreter','latex');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([1 1e3]);ylim([1 1e4]);
set(findall(fig,'-property','FontSize'),'FontSize',FS);

%% FIG.2: MATRIX PLOT BETWEEN AEROSOL SENSORS
% We do not use this.



labelX = {'$\mathcal{R}_1$','$\mathcal{R}_2$', ...
            '$\mathcal{L}_1$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$'};
        
%labelX = {'$\mathcal{R}_1$','$\mathcal{R}_2$', ...
%            '$\mathcal{S}_1$','$\mathcal{S}_2$','$\mathcal{S}_3$'};

labelY = fliplr({'$\mathcal{R}_1$','$\mathcal{R}_2$', ...
          '$\mathcal{L}_1$','$\mathcal{L}_{2a}$','$\mathcal{L}_{2b}$'});

corrP = corr(DATA1,'Type','Pearson','Rows','pairwise');
corrP = abs(corrP);
corrP = tril(corrP,-1);
for i=1:5
for j=1:5
mape(i,j) = nanmean(abs(DATA1(:,i)-DATA1(:,j))./(abs(DATA1(:,i))+abs(DATA1(:,j)))/2);
end
end
mape = triu(mape,1);
clc
figure(2);fig=gcf;
set(fig,'Position',[10.3333 41.6667 1280 599.3333])

Rms= corrP+mape;
Rms(Rms==0)=NaN;
Rms = fliplr(Rms);
Rms(:,end+1)=Rms(:,end);
Rms(end+1,:)=Rms(end,:);
p = pcolor(Rms);
p.LineStyle='None';
hold on
fill([1 2 2 3 3 4 4 5 5 6 6 5 5 4 4 3 3 2 2 1 1],...
     [6 6 5 5 4 4 3 3 2 2 1 1 2 2 3 3 4 4 5 5 6],[1 1 1],'LineStyle','None')
hold on
fill([1 1 2 2 3 3 4 4 5 5 1],[1 5 5 4 4 3 3 2 2 1 1],'k','FaceAlpha',0,'LineWidth',1)
hold on
fill([2 2 6 6 5 5 4 4 3 3 2],[5 6 6 2 2 3 3 4 4 5 5],'k','FaceAlpha',0,'LineWidth',1)
xlabel({'Mean absolute percentage error ($MAPE$)','Pearson coefficient ($R$)'}, 'Interpreter', 'latex')
h = colorbar;
h.Label.String = '$R$/$MAPE$';
h.Label.Interpreter = 'latex';
colormap(jet)
ax = gca;
caxis([0 1])
set(ax,'Position',[0.15 0.2 0.6 0.65])
set(ax, 'XTick', 1.5:length(labelX)+0.5, 'XTickLabel', labelX, 'TickLabelInterpreter', 'latex')
set(ax, 'YTick', 1.5:length(labelY)+0.5, 'YTickLabel', labelY, 'TickLabelInterpreter', 'latex')
set(h,'Ticks',0:0.2:1,'TickLabelInterpreter', 'latex')
%colorbar off
box off
set(findall(fig,'-property','FontSize'),'FontSize',FS);
print(gcf,'Fig2.png','-dpng','-r1000')

