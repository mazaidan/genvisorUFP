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