% "The role of pathogen mediated insect superabundance in the east-African emergence of a plant virus" Fig.3
%                                                                                                R. Donnelly
% this script simply produces graphs from digitised data
% only calculation is to use least scquares to fit epidemic curve from infected plant
% incidence data - note that this is fitting is to produce reference curves
% only in Fig. 3

% digistised from Colvin et al '04 Plant Pathology ("DigitiseFiles_Experiment" folder)
% (see folder "rGithub_v2\DigitiseFiles_Experiment" for source of the following numbers)
inc96=[74.215247 90.807175 53.587444 28.699552 14.573991  4.260090  6.278027  4.932735];
nymph96=[19.917012 53.112033 17.842324  7.883817  6.639004  6.224066 14.522822  7.883817];
adult96=[4.832347 13.017751  7.790927  2.761341  1.380671  1.676529  3.254438  3.846154];

inc97=[99.55157 98.65471 97.08520 44.39462 60.76233 90.35874 73.54260 97.75785];
nymph97=[79.253112 136.929461  48.132780   4.564315  29.045643 103.734440 118.257261 155.601660];
adult97=[16.469428 35.207101 15.779093  2.366864  6.213018 28.106509 30.276134 38.067061];

fieldkm=[0,9.202666667,18.40533333,27.53066667,36.73333333,45.78133333,51.81333333,58];

figure;
subplot(3,1,1);
subplot(3,1,2);
hold on;
yyaxis left
plot(fieldkm,adult96./(inc96+1),'color',[0.6 0.6 0.6])
scale1=max(adult96./(inc96+1));
set(gca,'ytick',[],'xtick',fieldkm);
xlim([-2 58])
ylim([0 0.75])
xlabel('#Insect vectors (per plant)')
ylabel('#Infected plants (per field')
yyaxis right
plot(fieldkm,inc96/100,'b.')
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',[0,10,20,30,40,50],'xticklabel',[0,10,20,30,40,50]);
xlim([-2 58])
ax1=gca;




subplot(3,1,3);
hold on;
yyaxis left
scale2=max(nymph96./(inc96+1));
plot(fieldkm,(scale1/scale2)*nymph96./(inc96+1),'--','color',[0.6 0.6 0.6])
set(gca,'ytick',[],'xtick',fieldkm);
xlim([-2 58])
ylim([0 0.75])
xlabel('#Insect vectors (per plant)')
ylabel('#Infected plants (per field')
yyaxis right
plot(fieldkm,inc96/100,'b.')
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',[0,10,20,30,40,50],'xticklabel',[0,10,20,30,40,50]);
xlim([-2 58])
ax2=gca;


%%%%%%%%% with fits %%%%%%%%%%%%%%%%%%
inc962=100-inc96;
x096=[1 30 100 5];

fun=@(xin,xdom)(xin(4)+(xin(3)-xin(4))./(1 + exp(-xin(1)*(xdom-xin(2)))));

xout96 = lsqcurvefit(fun,x096,fieldkm,inc962);
tmp=(fieldkm(1)-50):(fieldkm(end)+50);
y=(xout96(4)+(xout96(3)-xout96(4))./(1 + exp(-xout96(1)*(tmp-xout96(2)))));

subplot(3,1,3);
yyaxis right
hold on;
plot1=plot(tmp,(100-y)/100,'b');
plot1.Color(4) = 0.4;
subplot(3,1,2);
yyaxis right
hold on;
plot1=plot(tmp,(100-y)/100,'b');
plot1.Color(4) = 0.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


set(gcf, 'Units','centimeters', 'Position',1.5*[0 0 (7.25/8)*5*0.91 1*13.5])

set(gcf,'color', 'w');
set([ax1 ax2], ...
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

