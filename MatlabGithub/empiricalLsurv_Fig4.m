% digitised data from Legg and Ogwal 98 Journal Applied Entomology
% landscape invasion manuscript (contrasting scenarios for the role of insect superabundance)

% NEW
mnAdultC=[1.446886447,2.003663004,2.164835165,1.619047619,3.655677656,2.450549451,1.545787546,0.611721612,0.648351648,0.805860806];   % mean of 5 'seasons' Fig 2
mnNymphC=[2.65,2.91,2,3.54,6.06,3.96,2.47,0.8,0.9,1.07]*160/6.57;   % from Fig 2C (c(NakitomaLength,MageeraLength...)*yextent/yaxissize)
mnIncC=[92.60504202,77.71908764,85.882353,88.49939976,82.97719088,50.18007203,29.81992797,7.154861945,25.81032413,44.53781513]; % mean of 5 'seasons' Fig 2

mnAdultE=[1.051788376,0.89418778,1.171385991,0.369970194,0.489567809,0.746646796,0.376676602,0.34314456,0.292846498,0.384500745];   % mean of 4 'seasons' Fig 3
mnNymphE=[3.6,6.33,4.52,2.7,4.1,6.56,3.23,1.76,2.24,3.74]*35/6.72;   % from Fig 3C (c(KumiLength,AttuturLength...)*yextent/yaxissize)
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

fieldkmE=[0,10.73752711,16.78958785,31.82212581,42.55965293,51.73535792,61.4967462,70.86767896,79.65292842,90];
fieldkmE=fieldkmE(1:8);
fieldkmC=[0,18.5046729,44.69959947,66.80907877,91.08144192,101.6555407,118.7182911,136.5020027,153.564753,180];
fieldkmC=fieldkmC(1:8);


Fxx=1:8;





%%%%%%%%% with fits %%%%%%%%%%%%%%%%%%
incCT=100-mnIncC;
incET=100-mnIncE;
x096=[1 100 50 90 10];
% r is x(1)
% thalf is x(2)
% p0 is x(2)

fun=@(xin,xdom)([xin(5)+(xin(4)-xin(5))./(1 + exp(-xin(1)*(xdom(1:end/2)-xin(2)))) xin(5)+(xin(4)-xin(5))./(1 + exp(-xin(1)*(xdom((end/2)+1:end)-xin(3))))]);
xout96 = lsqcurvefit(fun,x096,[fieldkmC fieldkmE],[incCT incET]);
tmpC=(fieldkmC(1)-50):(fieldkmC(end)+50);
tmpE=(fieldkmE(1)-50):(fieldkmE(end)+50);
yC=xout96(5)+(xout96(4)-xout96(5))./(1 + exp(-xout96(1)*(tmpC-xout96(2))));
yE=xout96(5)+(xout96(4)-xout96(5))./(1 + exp(-xout96(1)*(tmpE-xout96(3))));
%alt_xoutCEyr24=[0.3591   21.8969   97.8652   57.4474   60.2731   90.1880   26.6238];
%y=(alt_xoutCEyr24(7)+(alt_xoutCEyr24(6)-alt_xoutCEyr24(7))./(1 + exp(-alt_xoutCEyr24(1)*(tmp-alt_xoutCEyr24(2)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
subplot(2,3,1);
subplot(2,3,2); % CENTRAL ADULT
hold on;
yyaxis left
%plot(fieldkmC,wpC,'color',[0.6 0.6 0.6])
plot(fieldkmC,wpC,'color',[0.6 0.6 0.6])

text(10,1.5,'A');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
%xlim([0.75 8.25])
%ylim([0 0.75])
xlabel('#Insect vectors (per plant)')
ylabel('#Infected plants (per field')
scale1=wpC(end);
ylim([scale1*0.15 scale1*1.1])
yyaxis right
plot(fieldkmC,mnIncC/100,'b.')
%set(gca,'ytick',0:10:100,'yticklabel',[],'xtick',1:8);
xlim([fieldkmC(1)-4 fieldkmC(end)+5])
%ylim([0 Lmax*1.05])
ax1=gca;
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:20:140);

plot1=plot(tmpC,(100-yC)/100,'b');
plot1.Color(4) = 0.4;




subplot(2,3,5); % CENTRAL NYMPH
hold on;
yyaxis left

text(10,1.5,'A');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
%xlim([0.75 8.25])
%ylim([0 0.75])
xlabel('#Insect vectors (per plant)')
ylabel('#Infected plants (per field')
scale2=NwpC(end);
plot(fieldkmC,(scale1/scale2)*NwpC,'--','color',[0.6 0.6 0.6])
ylim([scale1*0.15 scale1*1.1])
yyaxis right
plot(fieldkmC,mnIncC/100,'b.')
%set(gca,'ytick',0:10:100,'yticklabel',[],'xtick',1:8);
xlim([fieldkmC(1)-4 fieldkmC(end)+5])
%ylim([0 Lmax*1.05])
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
%xlim([0.75 8.25])
%ylim([0 0.75])
xlabel('#Insect vectors (per plant)')
%ylabel('#Infected plants (per field')
scale1=max(wpE);
ylim([scale1*0.15 scale1*1.1])

yyaxis right
plot(fieldkmE,mnIncE/100,'b.')
%set(gca,'ytick',0:0.5:1,'yticklabel',[],'xtick',1:8);
xlim([fieldkmE(1)-2 fieldkmE(end)+2])
%ylim([0 Lmax*1.05])
ax3=gca;
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:10:70);

plot1=plot(tmpE,(100-yE)/100,'b');
plot1.Color(4) = 0.4;



subplot(2,3,6);
hold on;
yyaxis left
text(10,1.5,'B');
set(gca,'ytick',0:0.5:1,'xtick',1:8);
%xlim([0.75 8.25])
%ylim([0 0.75])
xlabel('#Insect vectors (per plant)')
%ylabel('#Infected plants (per field')
scale2=max(NwpE);
plot(fieldkmE,(scale1/scale2)*NwpE,'--','color',[0.6 0.6 0.6])
ylim([scale1*0.15 scale1*1.1])

yyaxis right
plot(fieldkmE,mnIncE/100,'b.')
%set(gca,'ytick',0:0.5:1,'yticklabel',[],'xtick',1:8);
xlim([fieldkmE(1)-2 fieldkmE(end)+2])
%ylim([0 Lmax*1.05])
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
