% "The role of pathogen mediated insect superabundance in the east-African emergence of a plant virus" Fig.4
%                                                                                                R. Donnelly
% this script simply produces graphs from digitised data
% only calculation is to use least scquares to fit epidemic curve from infected plant
% incidence data - note that this is fitting is to produce reference curves
% only in Fig. 4

% digistised from Legg and Ogwal '98 Journal of Applied Entomology ("DigitiseFiles_Survey" folder)
% (see folder "rGithub_v2\DigitiseFiles_Survey" for source of the following numbers)

% SurveyDataConversions.csv column V
mnAdultC=[1.446886447,2.003663004,2.164835165,1.619047619,3.655677656,2.450549451,1.545787546,0.611721612,0.648351648,0.805860806];   % mean of 5 'seasons' Fig 2
% SurveyDataConversions.csv column Z
mnNymphC=[2.65,2.91,2,3.54,6.06,3.96,2.47,0.8,0.9,1.07]*160/6.57;   % from Fig 2C (c(NakitomaLength,MageeraLength...)*yextent/yaxissize)
% SurveyDataConversions.csv column W
mnIncC=[92.60504202,77.71908764,85.882353,88.49939976,82.97719088,50.18007203,29.81992797,7.154861945,25.81032413,44.53781513]; % mean of 5 'seasons' Fig 2

% SurveyDataConversions.csv column X
mnAdultE=[1.051788376,0.89418778,1.171385991,0.369970194,0.489567809,0.746646796,0.376676602,0.34314456,0.292846498,0.384500745];   % mean of 4 'seasons' Fig 3
% SurveyDataConversions.csv column AC
mnNymphE=[3.6,6.33,4.52,2.7,4.1,6.56,3.23,1.76,2.24,3.74]*35/6.72;   % from Fig 3C (c(KumiLength,AttuturLength...)*yextent/yaxissize)
% SurveyDataConversions.csv column Y
mnIncE=[70.7510 79.7431 74.7036 66.9631 73.88011 62.4835 21.0145 16.5679 16.69960 32.3123]; % mean of 4 'seasons' Fig 3

mnAdultC=mnAdultC(1:8);
mnNymphC=mnNymphC(1:8);
mnIncC=mnIncC(1:8);
mnAdultE=mnAdultE(1:8);
mnNymphE=mnNymphE(1:8);
mnIncE=mnIncE(1:8);

wpC=mnAdultC./(mnIncC+1);
NwpC=mnNymphC./(mnIncC+1);
wpE=mnAdultE./(mnIncE+1);
NwpE=mnNymphE./(mnIncE+1);


fieldkmE=[0,9.613636364,14.11363636,28.43181818,39.06818182,48.06818182,59.11363636,69.13636364,79.36363636,90];
fieldkmE=fieldkmE(1:8);

fieldkmC=[0,9.404934688,28.47605225,54.07837446,81.50943396,93.26560232,113.3817126,131.9303338,153.3526851,180];
fieldkmC=fieldkmC(1:8);

Fxx=1:8;




















%%%%%%%%% with fits %%%%%%%%%%%%%%%%%%
% For reference only, fitting a logistic curve to the incidence data points
incCT=100-mnIncC;
incET=100-mnIncE;
x0=[1 100 50 90 10];

fun=@(xin,xdom)([xin(5)+(xin(4)-xin(5))./(1 + exp(-xin(1)*(xdom(1:end/2)-xin(2)))) xin(5)+(xin(4)-xin(5))./(1 + exp(-xin(1)*(xdom((end/2)+1:end)-xin(3))))]);
xout = lsqcurvefit(fun,x0,[fieldkmC fieldkmE],[incCT incET]);
tmpC=(fieldkmC(1)-50):(fieldkmC(end)+50);
tmpE=(fieldkmE(1)-50):(fieldkmE(end)+50);
yC=xout(5)+(xout(4)-xout(5))./(1 + exp(-xout(1)*(tmpC-xout(2))));
yE=xout(5)+(xout(4)-xout(5))./(1 + exp(-xout(1)*(tmpE-xout(3))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
subplot(2,3,1);
subplot(2,3,2); % CENTRAL ADULT
hold on;
yyaxis left
plot(fieldkmC,wpC,'color',[0.6 0.6 0.6])
text(10,1.5,'A');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
xlabel('#Insect vectors (per plant)')
ylabel('#Infected plants (per field')
scale1=wpC(end);
ylim([scale1*0.15 scale1*1.1])
yyaxis right
plot(fieldkmC,mnIncC/100,'b.')
xlim([fieldkmC(1)-4 fieldkmC(end)+5])
ax1=gca;
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:20:140);
plot1=plot(tmpC,(100-yC)/100,'b');
plot1.Color(4) = 0.4;

subplot(2,3,5); % CENTRAL NYMPH
hold on;
yyaxis left
text(10,1.5,'A');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
xlabel('#Insect vectors (per plant)')
ylabel('#Infected plants (per field')
scale2=NwpC(end);
plot(fieldkmC,(scale1/scale2)*NwpC,'--','color',[0.6 0.6 0.6])
ylim([scale1*0.15 scale1*1.1])
yyaxis right
plot(fieldkmC,mnIncC/100,'b.')
xlim([fieldkmC(1)-4 fieldkmC(end)+5])
ax2=gca;
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:20:140);
plot1=plot(tmpC,(100-yC)/100,'b');
plot1.Color(4) = 0.4;


subplot(2,3,3);
hold on;
yyaxis left
plot(fieldkmE,wpE,'color',[0.6 0.6 0.6])
text(10,1.5,'B');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
xlabel('#Insect vectors (per plant)')
scale1=max(wpE);
ylim([scale1*0.15 scale1*1.1])

yyaxis right
plot(fieldkmE,mnIncE/100,'b.')
xlim([fieldkmE(1)-2 fieldkmE(end)+2])
ax3=gca;
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:10:70);
plot1=plot(tmpE,(100-yE)/100,'b');
plot1.Color(4) = 0.4;

subplot(2,3,6);
hold on;
yyaxis left
text(10,1.5,'B');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
xlabel('#Insect vectors (per plant)')
scale2=max(NwpE);
plot(fieldkmE,(scale1/scale2)*NwpE,'--','color',[0.6 0.6 0.6])
ylim([scale1*0.15 scale1*1.1])

yyaxis right
plot(fieldkmE,mnIncE/100,'b.')
xlim([fieldkmE(1)-2 fieldkmE(end)+2])
ax4=gca;
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:10:70);
plot1=plot(tmpE,(100-yE)/100,'b');
plot1.Color(4) = 0.4;

set(gcf, 'Units','centimeters', 'Position',1.5*[0 0 0.99*13.5*1.125 2*(7.25/8)*5*0.95])

set(gcf,'color', 'w');
set([ax1 ax2 ax3 ax4], ...
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
