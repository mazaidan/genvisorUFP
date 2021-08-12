%% This code is written in MATLAB
% The purpose aims to develop the estimator for ultra-fine particles based
% on low-cost sensors
% Written by:
% Martha Arbayani Bin Zaidan, Ph.D., Docent, CEng.
% Research Associate Professor, Nanjing University, China
% Senior Scientist, Helsinki University, Finland

addpath(genpath('Functions'));

clear; close all;clc;

load('DATA2.mat') 
DATA = DATA2;

Exp_smoking = [1:1:11521]';
Exp_kerosine = [11522:1:30242]';
Exp_gas = [30243:1:54721]';

%% QUICK CORRELATION ANALYSIS including time-delayed features
Ds = Exp_smoking; 
Dg = Exp_gas;
%Da =[Exp_smoking;Exp_gas];
Da = Ds; [Ds;Dg];

% % time-delayed 0 (no time delayed)
% Da_o = Da(1:end,1);
% Da_i = Da(1:end,1);
% % time-delayed 1
% Da_o = Da(2:end,1);
% Da_i = Da(1:end-1,1);
% % time-delayed 2
% Da_o = Da(3:end,1);
% Da_i = Da(1:end-2,1);

for T = 0:5
    % time-delayed T
    Da_o = Da(T+1:end,1);
    Da_i = Da(1:end-T,1);
    disp(['Time Delayed: ', num2str(T)])
    
    %CPC = DATA.PND_c([Ds;Dk;Dg],1);
    CPC = DATA.PND_c(Da_o,1);
    CPClog = log10(CPC);
    CPCgradient = gradient(CPC);
    CPCdiff = diff(CPC);
    CPCloggradient = gradient(CPClog);
    idx = find(CPCgradient <= nanmedian(CPCgradient)+10);
    CPCclean =CPC;
    CPCclean(idx,:)=nan;
    DATAo  = CPCclean;
    DATAo1 = log10(DATAo);
            
    X0 = [DATA.AT_T(Da_i,1),DATA.AT_RH(Da_i,1),DATA.LCS_G2_01_met(Da_i,3),DATA.LCS_G1(Da_i,1)];
    
    for n = 1:4
        Y = DATAo;DATAo1;
        X = X0(:,n);
        Rp = corr(X,Y,'Type','Pearson','Rows','complete');
        Rs = corr(X,Y,'Type','Spearman','Rows','complete');
        disp(['Rp: ',num2str(Rp), ' Rs: ',num2str(Rs), ' '])
    end
    disp('  ')
    
end

%% Cross-correlation

T = 0;
Da_o = Da(T+1:end,1);
Da_i = Da(1:end-T,1);

% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [DATA.AT_T(Da_i,1),DATA.AT_RH(Da_i,1),DATA.LCS_G2_01_met(Da_i,3),DATA.LCS_G1(Da_i,1)];

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CPC = DATA.PND_c(Da_o,1);
CPClog = log10(CPC);
CPCgradient = gradient(CPC);
CPCdiff = diff(CPC);
CPCloggradient = gradient(CPClog);
idx = find(CPCgradient <= nanmedian(CPCgradient)+10);
CPCclean =CPC;
CPCclean(idx,:)=nan;
DATAo  = CPCclean;
DATAo1 = log10(DATAo);
Y = DATAo;DATAo1;

% Compute XCORR and PLOT %%%%

Vars_names = {'Temp','RH','Pressure','PM_{2.5}'};
Corr_Type = 'Pearson';%'Spearman';
Lag =10;
Lags = [-Lag:1:Lag]';
for n = 1:4
    [R,L,pvalue] = crosscorrelation(X(:,n)',Y',Lag,Corr_Type);
    figure(1);
    subplot(2,2,n)
    stem(Lags,R)
    title(Vars_names{1,n})
end

%% MODELLING Linear Models and Shallow Neural Networks (version 2)

% 1) use all data, randominze, and predict CPC
% 2) use Exp_smoking, estimate Exp_kerosine and Exp_gas

Ds = Exp_smoking; 
Dg = Exp_gas;
Da =[Exp_smoking;Exp_gas];

%% Available Output (with/without cleaning) and Input Features

% PND data cleaning
CLEAN_PND = 1;
if CLEAN_PND == 0
    disp('OUTPUT: We do not remove PND data which is not clean')
    DATAo  = [DATA.PND_c([Ds;Dg],1)];
    DATAo1 = log10(DATAo);
elseif CLEAN_PND == 1
    disp('OUTPUT: We remove PND data which is not clean and normalize it')
    CPC = DATA.PND_c([Ds;Dg],1);
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
end

% PND data cleaning

x1 = [DATA.AT_T(Ds,1);  DATA.LCS_G2_01_met(Dg,2)]; % Temp
x1n = normalize_UFPsensors(x1,'Temp');

x2 = [DATA.AT_RH(Ds,1); DATA.LCS_G2_01_met(Dg,1)]; % RH
x2n = normalize_UFPsensors(x2,'RH');

x3 = [DATA.LCS_G2_01_met(Ds,3); DATA.LCS_G2_01_met(Dg,3)]; % P
x3n = normalize_UFPsensors(x3,'P');

x4 = [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]; % PM2.5
x4n = normalize_UFPsensors(x4,'PM25');

Xd  = []; % Input delayed
Xdn = []; % normalized input delayed
for T = 1 : 2
    %T =1; % time delayed = 1
    Ds_T = Ds(T+1:end,1);
    Dg_T = Dg(T+1:end,1);
    % Only Temp and PM2.5 included, necause the auto-correlation which we
    % obtain indicate that these two variables influence PND
    Xd0 = [[DATA.AT_T(Ds_T,1); nan(T,1); DATA.LCS_G2_01_met(Dg_T,2); nan(T,1)], ... 
        [DATA.LCS_G1(Ds_T,1);nan(T,1);DATA.LCS_G1(Dg_T,1); nan(T,1)]];
    %%%Xd0 = [DATA.LCS_G1(Ds_T,1);nan(T,1);DATA.LCS_G1(Dg_T,1); nan(T,1)];
    Xd = [Xd,Xd0];
    
    % Normalized values:
    xdn0a = normalize_UFPsensors(Xd0(:,1),'Temp');
    xdn0b = normalize_UFPsensors(Xd0(:,2),'PM25');
    Xdn0 = [xdn0a,xdn0b];
    Xdn = [Xdn,Xdn0]; 
end

DATAi  = [x1,x2,x3,x4,Xd];      %  DATA Input
%DATAi1 = [x1n,x2n,x3n,x4n,Xdn]; %  DATA input (with normalization)
%DATAi1 = [x1n,x2n,x4n,Xdn]; %  DATA input (with normalization)
%DATAi1 = [x1n.*x2n,x4n,Xdn]; %  DATA input (with normalization)
DATAi1 = [x1n.*x2n.*x4n,Xdn]; %  DATA input (with normalization)



% MODELLING
% INPUT   = DATAi, DATAi1
% OUTPUT  = DATAo, DATAo1

R     = zeros(2,2);%zeros(2,12);
MAPE  = zeros(2,2);%zeros(2,12);
for test_no=1:2%12
    if test_no == 1 
        Ds_half = roundn(size(Ds,1)/2,0);
        Dg_half = roundn(size(Dg,1)/2,0);
        TRAIN = [Ds(1:Ds_half);Dg(1:Dg_half)]; 
        TEST  = [Ds(Ds_half+1:end);Dg(Dg_half+1:end)]; 
    end
    if test_no == 2 
        TRAIN = [Ds(Ds_half+1:end);Dg(Dg_half+1:end)];  
        TEST  = [Ds(1:Ds_half);Dg(1:Dg_half)]; 
    end
    %if test_no == 1; TRAIN = Ds; TEST = Dg; end
    %if test_no == 2; TRAIN = Dg; TEST = Ds; end
    
    
    % Make the number for tracking
    no = [Ds;Dg];
    % Remove NaN data
    DATAt  = [DATAi1,DATAo1,no];
    DATAt1 = DATAt( ~any( isnan( DATAt ) | isinf( DATAt ), 2 ),: );
    
    Feature = 1;
    if Feature == 1
        disp('Wavelet filtering')
        d = DATAt1(:,end-1);
        xden = wdenoise(d,4);
        DATAt1(:,end-1) = xden ;
    else
        disp('No feature extraction')
    end
    
    
    tr = ismember(DATAt1(:,end),TRAIN);
    te = ismember(DATAt1(:,end),TEST);
    
    X = DATAt1(tr,1:end-2);
    Y = DATAt1(tr,end-1);
    Xt = DATAt1(te,1:end-2);
    Yt = DATAt1(te,end-1);
     
    Model1 = 'LM1'; 
    Model2 = 'ANN1';
    [Ypred_lm] = UFPmodelling(X,Y,Xt,Model1);
    [Ypred_snn] = UFPmodelling(X,Y,Xt,Model2);

   
    % RESULT PLOTS and METRICS
    figure(1); fig =gcf;
    subplot(2,2, test_no*2 - 1 );
    scatter(Yt,Ypred_lm);hold on
    Xlim1 = 0;3;
    Ylim1 = 6;
    xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
    x = linspace(Xlim1,Ylim1);
    y = linspace(Xlim1,Ylim1);
    plot(x,y,'r');hold off
    title(['M1: Test No: ',num2str(test_no)])
    xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')
    subplot(2,2, test_no*2 );
    scatter(Yt,Ypred_snn);hold on
    Xlim1 = 0;3;
    Ylim1 = 6;
    xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
    x = linspace(Xlim1,Ylim1);
    y = linspace(Xlim1,Ylim1);
    plot(x,y,'r');hold off
    title(['M2: Test No: ',num2str(test_no)])
    xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')
    set(findall(fig,'-property','FontSize'),'FontSize',22);
    
    figure(2); fig = gcf;
    subplot(1,2,test_no); %
    %subplot(4,3,test_no);
    plot(10.^Yt,'b.','MarkerSize',12);hold on;grid on
    plot(10.^Ypred_lm,'g.','MarkerSize',12);
    plot(10.^Ypred_snn,'r.','MarkerSize',12);
    ylabel('PND [/cm$^{-3}$]','interpreter','latex')
    xlabel('Time Index')
    ylim([0 1e6])
    set(gca, 'YScale', 'log')
    legend('Real','M$_1$','M$_2$','interpreter','latex')
    title(['Test No: ',num2str(test_no)])
    hold off
    set(findall(fig,'-property','FontSize'),'FontSize',22);

 
    R(1,test_no) = corr(Yt,Ypred_lm,'Type','Spearman','Rows','complete');
    R(2,test_no) = corr(Yt,Ypred_snn,'Type','Spearman','Rows','complete');

    MAPE(1,test_no)=errperf(Yt,Ypred_lm,'mape');
    MAPE(2,test_no)=errperf(Yt,Ypred_snn,'mape');
         
end

Rmean     = mean(R,2);
MAPEmean  = mean(MAPE,2);
disp(['Rmean of M1: ', num2str(Rmean(1)), ' and Rmean of M2: ', num2str(Rmean(2))])
disp(['MAPEmean of M1: ', num2str(MAPEmean(1)), ' and MAPEmean of M2: ', num2str(MAPEmean(2))])


%%

% CHOOSE THE DATA with the different types of inputs
Di = 9;

if Di == 1
    disp('Temp and PM2.5')
    DATAi = [ [DATA.AT_T(Ds,1);DATA.LCS_G2_01_met(Dg,2)], ...
              [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 2
    disp('RH and PM2.5')
    DATAi = [[DATA.AT_RH(Ds,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 3
    disp('P and PM2.5')
    DATAi = [[DATA.LCS_G2_01_met(Ds,3);DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 4
    disp('Temp, RH and PM2.5')
    DATAi = [[DATA.AT_T(Ds,1);DATA.LCS_G2_01_met(Dg,2)], ...
             [DATA.AT_RH(Ds,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 5         
    disp('Temp, P and PM2.5')
    DATAi = [[DATA.AT_T(Ds,1);DATA.LCS_G2_01_met(Dg,2)], ...
             [DATA.LCS_G2_01_met(Ds,3);DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 6         
    disp('RH, P and PM2.5')
    DATAi = [[DATA.AT_RH(Ds,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G2_01_met(Ds,3);DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 7
    disp('Temp, RH, P and PM2.5')
    %DATAi = [[DATA.AT_T(Ds,1);DATA.CO_T(Dk,1);DATA.LCS_G2_01_met(Dg,2)], ...
    %         [DATA.AT_RH(Ds,1);DATA.CO_RH(Dk,1);DATA.LCS_G2_01_met(Dg,1)], ...
    %         [DATA.LCS_G2_01_met(Ds,3);DATA.CO_P(Dk,1);DATA.LCS_G2_01_met(Dg,3)], ...
    %         [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dk,1);DATA.LCS_G1(Dg,1)]];
    DATAi = [[DATA.AT_T(Ds,1);DATA.LCS_G2_01_met(Dg,2)], ...
             [DATA.AT_RH(Ds,1);DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G2_01_met(Ds,3);DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 8
    disp('Temp, Temp(t-1), RH, P and PM2.5')
    % Da_o = Da(T+1:end,1); Da_i = Da(1:end-T,1);
    T = 1 ;
    Ds_T = Ds(T+1:end,1);
    Dg_T = Dg(T+1:end,1);
    DATAi = [[DATA.AT_T(Ds,1);  DATA.LCS_G2_01_met(Dg,2)], ... 
            [DATA.AT_T(Ds_T,1); nan(T,1); DATA.LCS_G2_01_met(Dg_T,2); nan(T,1)] ...
             [DATA.AT_RH(Ds,1); DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G2_01_met(Ds,3); DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)]];
elseif Di == 9
    disp('Temp, Temp(t-1), RH, P and PM2.5, PM2.5(t-1)')
    % Da_o = Da(T+1:end,1); Da_i = Da(1:end-T,1);
    T = 1 ;
    Ds_T = Ds(T+1:end,1);
    Dg_T = Dg(T+1:end,1);
    DATAi = [[DATA.AT_T(Ds,1);  DATA.LCS_G2_01_met(Dg,2)], ... 
            [DATA.AT_T(Ds_T,1); nan(T,1); DATA.LCS_G2_01_met(Dg_T,2); nan(T,1)] ...
             [DATA.AT_RH(Ds,1); DATA.LCS_G2_01_met(Dg,1)], ...
             [DATA.LCS_G2_01_met(Ds,3); DATA.LCS_G2_01_met(Dg,3)], ...
             [DATA.LCS_G1(Ds,1);DATA.LCS_G1(Dg,1)] ...
             [DATA.LCS_G1(Ds_T,1);nan(T,1);DATA.LCS_G1(Dg_T,1); nan(T,1)]];
end

% NORMALIZATION
NORM = 1; l = -1; u = 1;
DATAi1 = zeros(size(DATAi));
if NORM == 0
    disp('No NORMALIZATION for MET vars')
    DATAi1(:,1:end-1) = DATAi(:,1:end-1);
    DATAi1(:,end) = log10(DATAi(:,end));
elseif NORM == 1
    disp('NORMALIZATION for MET vars')
    if Di == 1
        disp('Temp and PM2.5'); inmin = 10; inmax = 40;
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
    elseif Di == 2
        disp('RH and PM2.5'); inmin = 10; inmax = 50;
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
    elseif Di == 3
        disp('P and PM2.5'); inmin = 890; inmax = 910;
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
    elseif Di == 4
        disp('Temp, RH and PM2.5'); 
        inmin = 10; inmax = 40; % Temp
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 10; inmax = 50; % RH
        DATAi1(:,2) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
    elseif Di == 5
        disp('Temp, RH and PM2.5'); 
        inmin = 10; inmax = 40; % Temp
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 890; inmax = 910; % P
        DATAi1(:,2) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
    elseif Di == 6
        disp('RH, P and PM2.5'); 
        inmin = 10; inmax = 50; % RH
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 890; inmax = 910; % P
        DATAi1(:,2) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
    elseif Di == 7
        disp('Temp, RH, P and PM2.5');
        inmin = 10; inmax = 40; % T
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 10; inmax = 50; % RH
        DATAi1(:,2) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 890; inmax = 910; % P
        DATAi1(:,3) = l + [(DATAi(:,3)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,end) = log10(DATAi(:,end));
        
     elseif Di == 8
        disp('Temp, Temp(t-1), RH, P and PM2.5');
        inmin = 10; inmax = 40; % T
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,2) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 10; inmax = 50; % RH
        DATAi1(:,3) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 890; inmax = 910; % P
        DATAi1(:,4) = l + [(DATAi(:,3)-inmin)./(inmax-inmin)].*(u-l);
        % PM2.5
        DATAi1(:,end) = log10(DATAi(:,end));
        
     elseif Di == 9
        disp('Temp, Temp(t-1), RH, P and PM2.5, PM2.5(t-1)');
        inmin = 10; inmax = 40; % T
        DATAi1(:,1) = l + [(DATAi(:,1)-inmin)./(inmax-inmin)].*(u-l);
        DATAi1(:,2) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 10; inmax = 50; % RH
        DATAi1(:,3) = l + [(DATAi(:,2)-inmin)./(inmax-inmin)].*(u-l);
        inmin = 890; inmax = 910; % P
        DATAi1(:,4) = l + [(DATAi(:,3)-inmin)./(inmax-inmin)].*(u-l);
        % PM2.5
        DATAi1(:,end) = log10(DATAi(:,end-1));
        DATAi1(:,end) = log10(DATAi(:,end));
    end
end

% PND data cleaning
CLEAN_PND = 1;
if CLEAN_PND == 0
    disp('We do not remove PND data which is not clean')
    DATAo  = [DATA.PND_c([Ds;Dg],1)];
    DATAo1 = log10(DATAo);
elseif CLEAN_PND == 1
    disp('We remove PND data which is not clean')
    CPC = DATA.PND_c([Ds;Dg],1);
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
    
% % %     figure(100)
% % %     for n=1:size(DATAi,2)
% % %         subplot(2,5,n);plot(DATAi(:,n),'.');hold on
% % %         subplot(2,5,n+5);plot(DATAi1(:,n),'.');
% % %     end
% % %     subplot(2,5,5);plot(DATAo(:,1),'.');
% % %     subplot(2,5,10);plot(DATAo1(:,1),'.');
% % %     hold off
    
    figure(101);
    subplot(511);plot(CPC,'b.');
    hold on; plot(CPCgradient,'r.'); hold off
    subplot(512);plot(CPClog,'b.');
    hold on; plot(CPCloggradient,'r.'); hold off
    subplot(513);plot(CPCgradient,'b.');hold on
    subplot(513);plot(CPCloggradient,'r.');hold off
    subplot(514);plot(CPC,'b.');hold on;plot(CPCclean,'r.');hold off
    subplot(515);plot(log10(CPC),'b.');hold on;plot(log10(CPCclean),'r.');hold off
else
    disp('We need to choose CLEAN_PND either 0 or 1')
end
Dk = Dg;
% SELECT TRAINING AND TESTING DATA
R     = zeros(2,2);%zeros(2,12);
MAPE  = zeros(2,2);%zeros(2,12);
%R     = zeros(2,12);%zeros(2,12);
%MAPE  = zeros(2,12);%zeros(2,12);


for test_no=1:2%12
    if test_no == 1; TRAIN = Ds; TEST = Dg; end
    if test_no == 2; TRAIN = Dg; TEST = Ds; end
    
%     if test_no == 1; TRAIN = Ds; TEST = Dk; end
%     if test_no == 2; TRAIN = Ds; TEST = Dg; end
%     if test_no == 3; TRAIN = Ds; TEST = [Dk;Dg]; end
%     if test_no == 4; TRAIN = Dk; TEST = Ds; end
%     if test_no == 5; TRAIN = Dk; TEST = Dg; end
%     if test_no == 6; TRAIN = Dk; TEST = [Ds;Dg]; end
%     if test_no == 7; TRAIN = Dg; TEST = Ds; end
%     if test_no == 8; TRAIN = Dg; TEST = Dk; end
%     if test_no == 9; TRAIN = Dg; TEST = [Ds;Dk]; end
%     if test_no == 10; TRAIN = [Ds;Dk]; TEST = Dg; end
%     if test_no == 11; TRAIN = [Ds;Dg]; TEST = Dk; end
%     if test_no == 12; TRAIN = [Dk;Dg]; TEST = Ds; end

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
    %%%no = linspace(1, size(DATAi1,1), size(DATAi1,1))' ;
    no = [Ds;Dg];
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


    %R0 =  corrcoef(Yt,Ypred_lm);
    %R(1,test_no) = R0(2,1);
    R(1,test_no) = corr(Yt,Ypred_lm,'Type','Spearman','Rows','complete');
    
    %R0 =  corrcoef(Yt,Ypred_snn);
    %R(2,test_no) = R0(2,1);
    R(2,test_no) = corr(Yt,Ypred_snn,'Type','Spearman','Rows','complete');

    %MAPE(1,test_no) = nanmean(abs((Yt-Ypred_lm)./Yt));
    MAPE(1,test_no)=errperf(Yt,Ypred_lm,'mape');
    MAPE(2,test_no)=errperf(Yt,Ypred_snn,'mape');
    
    figure(1);
    subplot(1,2,test_no);
    %subplot(4,3,test_no);
    scatter(Yt,Ypred_lm);hold on
    Xlim1 = 0;3;
    Ylim1 = 6;
    xlim([Xlim1 Ylim1]);ylim([Xlim1 Ylim1]);grid on
    x = linspace(Xlim1,Ylim1);
    y = linspace(Xlim1,Ylim1);
    plot(x,y,'r');hold off
    title(['Test No: ',num2str(test_no)])
    xlabel('log Real PNC (CPC)');ylabel('log Est PNC (CPC)')
    
    figure(2); %fig = gcf;
    subplot(1,2,test_no); %
    %subplot(4,3,test_no);
    plot(10.^Yt,'b.','MarkerSize',12);hold on;grid on
    plot(10.^Ypred_lm,'g.','MarkerSize',12);
    plot(10.^Ypred_snn,'r.','MarkerSize',12);
    ylabel('PND [/cm$^{-3}$]','interpreter','latex')
    xlabel('Time Index')
    ylim([0 1e6])
    set(gca, 'YScale', 'log')
    legend('Real','M$_1$','M$_2$','interpreter','latex')
    title(['Test No: ',num2str(test_no)])
    hold off
    %set(findall(fig,'-property','FontSize'),'FontSize',22);

    
    
end

Rmean     = mean(R,2);
MAPEmean  = mean(MAPE,2); 
disp(['Rmean of M1: ', num2str(Rmean(1)), ' and Rmean of M2: ', num2str(Rmean(2))])
disp(['MAPEmean of M1: ', num2str(MAPEmean(1)), ' and MAPEmean of M2: ', num2str(MAPEmean(2))])


figure(15); fig = gcf;
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



