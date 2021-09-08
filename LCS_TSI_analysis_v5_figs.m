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
Y = SidePak(:,8);
X = ISEE_LCS_G202(:,end); 

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