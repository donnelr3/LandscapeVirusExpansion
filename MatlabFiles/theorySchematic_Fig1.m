% "The role of pathogen mediated insect superabundance in the east-African emergence of a plant virus" Fig.1
%                                                                                                R. Donnelly
% this script simply produces illustrative curves for schematic 
% only calculation is to take ratio of vector abundance and infected plant incidence 

Fxx=-200:200;

kapp=0.4;

% abundance per plant
A0=2;
A1=20;

% max number of infected plants per field
Lmax=100;

% delay between insect and plant disease waves
invasive_tau=15;
modif_tau=15;

xx0_1=50;
xx0_2=xx0_1-modif_tau;
xx0_3=xx0_1+invasive_tau;

I_yy1=Lmax-(Lmax./(1+exp(-kapp*(Fxx-xx0_1))));

A_yy2=A0+(A1-A0)*(Lmax-(Lmax./(1+exp(-kapp*(Fxx-xx0_2)))))/Lmax;

A_yy3=A0+(A1-A0)*(Lmax-(Lmax./(1+exp(-kapp*(Fxx-xx0_3)))))/Lmax;

A_yy4=A1*ones(size(Fxx));
con=1;
xticmarks=0:20:100;

figure;
hold on;

subplot(5,1,1)
plot([1,2],[1,2]);
ax1=gca;

subplot(5,1,2)
yyaxis left
plot(Fxx,I_yy1);
ylim([0 Lmax*1.05])
ylabel('#Infected plants (per field')
yyaxis right
hold on;
plot(Fxx,A_yy2,'-') %,'m'
plot(Fxx,A_yy3,'-') %,'b'
plot(Fxx,A_yy4,'-') %,'--g'
set(gca,'ytick',0:10:20,'xtick',xticmarks,'xticklabel',[]);
%ylim([0 A1*1.15])
ylabel('#Insect vectors (per plant)')
xlim([0 100])
text(10,13,'A');
ax2=gca;

subplot(5,1,3)
plot(Fxx,(A_yy2./(I_yy1+con)),'m')
set(gca,'ytick',0:5,'xtick',xticmarks,'xticklabel',[]);
xlim([0 100])
%ylim([-0.1 2.5])
text(10,1.5,'B');
ax3=gca;

subplot(5,1,4)
plot(Fxx,(A_yy3./(I_yy1+con)),'b')
set(gca,'ytick',0:5:10,'xtick',xticmarks);
xlim([0 100])
%ylim([0 11])
text(10,9,'C');
xlabel('Field position')
ax4=gca;

subplot(5,1,5)
plot(Fxx,(A_yy4./(I_yy1+con)),'g')
set(gca,'ytick',0:10:20,'xtick',xticmarks);
xlim([0 100])
%ylim([0 25])
text(10,9,'D');
xlabel('Field position')
ax5=gca;

set(gcf, 'Units','centimeters', 'Position',1.5*[0 0 5*(1/1.275) (5/3)*13.5])

set(gcf,'color', 'w');
set([ax1 ax2 ax3 ax4 ax5], ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'FontSize', 12, ...
  'YColor'      , [.3 .3 .3], ...
  'FontName', 'Times New Roman', ...
  'LineWidth'   , 0.2         );

h=gcf;
h.Renderer='Painters';
orient(h,'landscape')
