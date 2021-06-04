clear all
global a K thet H muu Pacq Pinoc plD sig pMigrate pDisappear eps1 epsSP2 muuN plE epslscape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Switched to coding format based on midscape stopping condition %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note SP2 switch must be on for invasive wf scenario
muu=1/50; % vector mortality
a=1*(100*muu/2); % vector repro
K=40; % carrying capacity
H=100; % number plants in field

Pacq=0.032;
Pinoc=0.032; % rate inoculate

thet=1; % dispersal rate
sig=0; % assuming life long retention period
plE=1/360; % plant mortality
plD=plE; % infected, plant mortality roguing 
muuN=2*muu; % vector mortality for nymphs

% default migration parameter values
pMigrate=0.5;
pDisappear=0.5;
prSwitch=0;

icbnMain=1/30;% disease latency
dNmain=1/25; %vector development
extraFig=1;
propCleanSeed=1;
cutFromWithin=0.5;
range=119;

numFieldsIn=120;
denC=1;  % wave-profile: A/(I+denC)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% emergence through PATHOGEN MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps1=20;    
epsSP2=1; 
sp2On=0;
epslscape=1;
tMaxSS=1000;
tMax=10000;
startField=25;
%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%


%%% two-step approach, - run to equil. in first step 
invIndex=1; % field location of invader
initVals=[zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) K*(1-(muu/a))*ones(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RUN first TO insect STEADY-STATE (no pathogen!) %%%%%%%%%%%%%%%%%%%%
%%% Note incubation period does not impact outcome %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tnewSS,ynewSS] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,icbnMain,dNmain,range,propCleanSeed,cutFromWithin),[0 tMaxSS],initVals);
nymphsSS=ynewSS(:,0*numFieldsIn+(1:numFieldsIn));
insectsSS=ynewSS(:,3*numFieldsIn+(1:numFieldsIn));
nymphsSS2=ynewSS(:,11*numFieldsIn+(1:numFieldsIn));
insectsSS2=ynewSS(:,14*numFieldsIn+(1:numFieldsIn));
if prSwitch
    figure;
    plot(tnewSS,nymphsSS);
    figure;
    plot(tnewSS,insectsSS);
    figure;
    plot(tnewSS,nymphsSS2);
    figure;
    plot(tnewSS,insectsSS2);
end
initValsFromSS=ynewSS(end,:);
initValsFromSS(7*numFieldsIn+invIndex)=1;    % <- disease invasion
initValsFromSS(14*numFieldsIn+invIndex)=sp2On*1; % <- insect invasion (susceptible)
initValsFromSS(3*numFieldsIn+invIndex)=initValsFromSS(3*numFieldsIn+invIndex)-sp2On*1; % <- simply reflecting invasion in loss of indiv. from WT (actually this represents insect mutation rather than invasion)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PLOT incubation curves at diff. snapshots %%%%%%%%%%
options = odeset('Events',@midScapeEvents);
[tnew,ynew,te1,ye,ie] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,icbnMain,dNmain,range,propCleanSeed,cutFromWithin),[0 tMax],initValsFromSS,options);
infVec=ynew(:,7*numFieldsIn+(1:numFieldsIn));
exposVecF=ye(6*numFieldsIn+(1:numFieldsIn/2));
infVecF=ye(7*numFieldsIn+(1:numFieldsIn/2));
expVecF=exposVecF/H;
incVecF=infVecF/H;
if sum(sum((incVecF+expVecF)>1.001))
    disp('ERROR!!! unusual incidence');
    for qq=1:numFieldsIn
        tmpfield=incVecF(:,qq)+expVecF(:,qq);
        inds=find(tmpfield>1.001);
        if ~isempty(inds)
            disp(['Field is ' num2str(qq)]);
            disp(['At time ' num2str(inds')]);
            disp(incVecF(inds,qq)');
            disp(expVecF(inds,qq)');
            return;
        end
    end
end
superAbundAF=H*(1-incVecF-expVecF).*(ye(3*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(4*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(:,5*numFieldsIn+(1:numFieldsIn/2)));%sp#1
superAbundBF=H*(1-incVecF-expVecF).*(ye(14*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(15*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(:,16*numFieldsIn+(1:numFieldsIn/2)));%sp#2
n_superAbundAF=H*(1-incVecF-expVecF).*(ye(0*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(1*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(:,2*numFieldsIn+(1:numFieldsIn/2)));%sp#1nymph
n_superAbundBF=H*(1-incVecF-expVecF).*(ye(11*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(12*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(:,13*numFieldsIn+(1:numFieldsIn/2)));%sp#2nymph
superAbundF=superAbundAF+superAbundBF; % both sp
n_superAbundF=n_superAbundAF+n_superAbundBF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CLOSE-UP
xticmarks=0:10:100;
figure;
hold on;
yyaxis left
whichFields=startField:(startField+15);
incy=infVecF;
incy2=100-incy;
options = optimoptions('lsqcurvefit'); %,'Algorithm','levenberg-marquardt'

alt_lb = [0, 10, 50, 0];
alt_ub = [3, 50, 150, 50];
alt_x0y=(alt_lb+alt_ub)/2;
funy=@(xin,xdom1)(xin(4)+((xin(3)-xin(4))./(1 + exp(-xin(1)*(xdom1-xin(2))))));
alt_xouty = lsqcurvefit(funy,alt_x0y, 1:numFieldsIn/2,incy2,alt_lb,alt_ub,options);

wProf=(1/H)*superAbundF./(infVecF+denC);
newX=1:0.001:(numFieldsIn/2);

yy = spline(1:(numFieldsIn/2),wProf,newX);
xMin=newX(find(yy==min(yy)));

plot(1:(numFieldsIn/2),wProf);
set(gca,'ytick',0:1:10,'xtick',0:50);
plot([xMin xMin],[0 10])
ylim([0 2.5])
hold on;
yyaxis right
plot(1:(numFieldsIn/2),incVecF);
set(gca,'ytick',0:0.2:1,'yticklabel',[],'xtick',0:50);
title(['Close-up A\infty= ' num2str(round(superAbundF(1)/H))])
xlim([whichFields(1) whichFields(end)])

tmpVec=find(incVecF<max(incVecF)/2);
if ~isempty(tmpVec)
    halfSatLoc=tmpVec(1);
    plot([alt_xouty(2) alt_xouty(2)],[0 1])
    text(alt_xouty(2)*1.01,0.1,['\delta=  ' num2str(alt_xouty(2)-xMin)])
    text(alt_xouty(2)*1.01,0.2,['r=  ' num2str(alt_xouty(1))])
end

figure
hold on;
subplot(4,2,3)
yyaxis left
p1=plot(1:(numFieldsIn/2),(1/H)*superAbundF./(infVecF+denC),'k');
normTo1=max((1/H)*superAbundF./(infVecF+denC));
hold on;
normTo2=max((1/H)*n_superAbundF./(infVecF+denC));
p2=plot(1:(numFieldsIn/2),(normTo1/(normTo2*H))*n_superAbundF./(infVecF+denC),'g--');
xlim([20 55])
set(gca,'ytick',0:10:200,'xtick',xticmarks,'xticklabel',[]);
ylabel('Wave-profile (Y/(I+1))');
xlabel('Field position');
lgd=legend([p1 p2],'adult','nymph','location','west');
lgd.FontSize = 9;
yyaxis right
p3=plot(1:(numFieldsIn/2),incVecF,'color',[0.68 0.85 0.9]);
ylim([0 1]);
set(gca,'ytick',[],'xtick',xticmarks,'xticklabel',[]);
ax1=gca;
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend boxoff;  
ax1.YAxis(1).Color = 'k';

subplot(4,2,4)
yyaxis left
plot(1:(numFieldsIn/2),incVecF);
%ylabel('#Insect vectors (per plant)')
ylim([0 1])
set(gca,'ytick',0:0.2:1,'xtick',xticmarks,'xticklabel',[]);
yyaxis right
plot(1:(numFieldsIn/2),(1/H)*superAbundF);
%plot(1:(numFieldsIn/2),superAbundF);
set(gca,'ytick',0:100:1000,'xtick',xticmarks,'xticklabel',[]);
ylim([0 350])
xlim([20 55])
ax1b=gca;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%% emergence through INVASIVE VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps1=1;
epsSP2=5; 
sp2On=1;
epslscape=1;
%%%% INVASIVE VECTOR %%%%%%%%%%%%%
%%%% INVASIVE VECTOR %%%%%%%%%%%%%
%%%% INVASIVE VECTOR %%%%%%%%%%%%%
%%%% INVASIVE VECTOR %%%%%%%%%%%%%
%%%% INVASIVE VECTOR %%%%%%%%%%%%%

%%% two-step approach, - run to equil. in first step
invIndex=1; % field location of invader
initVals=[zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) K*(1-(muu/a))*ones(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RUN first TO insect STEADY-STATE (no pathogen!) %%%%%%%%%%%%%%%%%%%%
%%% Note incubation period does not impact outcome %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
propCleanSeed=1;
[tnewSS,ynewSS] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,icbnMain,dNmain,range,propCleanSeed,cutFromWithin),[0 tMaxSS],initVals);
nymphsSS=ynewSS(:,0*numFieldsIn+(1:numFieldsIn));
insectsSS=ynewSS(:,3*numFieldsIn+(1:numFieldsIn));
nymphsSS2=ynewSS(:,11*numFieldsIn+(1:numFieldsIn));
insectsSS2=ynewSS(:,14*numFieldsIn+(1:numFieldsIn));
if prSwitch
    figure;
    plot(tnewSS,nymphsSS);
    figure;
    plot(tnewSS,insectsSS);
    figure;
    plot(tnewSS,nymphsSS2);
    figure;
    plot(tnewSS,insectsSS2);
end
initValsFromSS=ynewSS(end,:);
initValsFromSS(7*numFieldsIn+invIndex)=1;    % <- disease invasion
initValsFromSS(14*numFieldsIn+invIndex)=sp2On*1; % <- insect invasion (susceptible)
initValsFromSS(3*numFieldsIn+invIndex)=initValsFromSS(3*numFieldsIn+invIndex)-sp2On*1; % <- simply reflecting invasion in loss of indiv. from WT (actually this represents insect mutation rather than invasion)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PLOT incubation curves at diff. snapshots %%%%%%%%%%
options = odeset('Events',@midScapeEvents);
[tnew,ynew,te2,ye,ie] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,icbnMain,dNmain,range,propCleanSeed,cutFromWithin),[0 tMax],initValsFromSS,options);
infVecF=ye(7*numFieldsIn+(1:numFieldsIn/2));
exposVecF=ye(6*numFieldsIn+(1:numFieldsIn/2));
expVecF=exposVecF/H;
incVecF=infVecF/H;
if sum(sum((incVecF+expVecF)>1.001))
    disp('ERROR!!! unusual incidence');
    for qq=1:numFieldsIn/2
        tmpfield=incVecF(qq)+expVecF(qq);
        inds=find(tmpfield>1.001);
        if ~isempty(inds)
            disp(['Field is ' num2str(qq)]);
            disp(['At time ' num2str(inds')]);
            disp(incVecF(qq)');
            disp(expVecF(qq)');
            return;
        end
    end
end
superAbundAF=H*(1-incVecF-expVecF).*(ye(3*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(4*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(5*numFieldsIn+(1:numFieldsIn/2)));%sp#1
superAbundBF=H*(1-incVecF-expVecF).*(ye(14*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(15*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(16*numFieldsIn+(1:numFieldsIn/2)));%sp#2
n_superAbundAF=H*(1-incVecF-expVecF).*(ye(0*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(1*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(2*numFieldsIn+(1:numFieldsIn/2)));%sp#1nymph
n_superAbundBF=H*(1-incVecF-expVecF).*(ye(11*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(12*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(13*numFieldsIn+(1:numFieldsIn/2)));%sp#2nymph
superAbundF=superAbundAF+superAbundBF;
n_superAbundF=n_superAbundAF+n_superAbundBF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xticmarks=0:10:100;
subplot(4,2,5)
yyaxis left
p1=plot(1:(numFieldsIn/2),(1/H)*superAbundF./(infVecF+denC),'k');
normTo1=max((1/H)*superAbundF./(infVecF+denC));
hold on;
normTo2=max((1/H)*n_superAbundF./(infVecF+denC));
p2=plot(1:(numFieldsIn/2),(normTo1/(normTo2*H))*n_superAbundF./(infVecF+denC),'g--');
set(gca,'ytick',0:20:200,'xtick',xticmarks,'xticklabel',[]);
ylabel('Wave-profile (Y/(I+1))')
xlabel('Field position')
lgd=legend([p1 p2],'adult','nymph','location','west');
lgd.FontSize = 9;
%ylim([0 90])
yyaxis right
p3=plot(1:(numFieldsIn/2),incVecF,'color',[0.68 0.85 0.9]);
set(gca,'ytick',[],'xtick',xticmarks,'xticklabel',[]);
xlim([20 55])
ylim([0 1])
ax2=gca;
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend boxoff  
ax2.YAxis(1).Color = 'k';

subplot(4,2,6)
yyaxis left
plot(1:(numFieldsIn/2),incVecF);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks,'xticklabel',[]);
ylim([0 1])
yyaxis right
plot(1:(numFieldsIn/2),superAbundF/H);
%plot(1:(numFieldsIn/2),superAbundF);
set(gca,'ytick',0:100:600,'xtick',xticmarks,'xticklabel',[]);
xlim([20 55])
ylim([0 350])
ax2b=gca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%% emergence through PATHOGEN MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps1=1;     
epsSP2=1; 
sp2On=0;
epslscape=3; 
%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
%%% two-step approach, - run to equil. in first step
invIndex=1; % field location of invader
initVals=[zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) K*(1-(muu/a))*ones(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% RUN first TO insect STEADY-STATE (no pathogen!) %%%%%%%%%%%%%%%%%%%%
%%% Note incubation period does not impact outcome %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
propCleanSeed=1;
[tnewSS,ynewSS] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,icbnMain,dNmain,range,propCleanSeed,cutFromWithin),[0 tMaxSS],initVals);
nymphsSS=ynewSS(:,0*numFieldsIn+(1:numFieldsIn));
insectsSS=ynewSS(:,3*numFieldsIn+(1:numFieldsIn));
nymphsSS2=ynewSS(:,11*numFieldsIn+(1:numFieldsIn));
insectsSS2=ynewSS(:,14*numFieldsIn+(1:numFieldsIn));
initValsFromSS=ynewSS(end,:);
initValsFromSS(7*numFieldsIn+invIndex)=1;    % <- disease invasion
initValsFromSS(14*numFieldsIn+invIndex)=sp2On*1; % <- insect invasion (susceptible)
initValsFromSS(3*numFieldsIn+invIndex)=initValsFromSS(3*numFieldsIn+invIndex)-sp2On*1; % <- simply reflecting invasion in loss of indiv. from WT (actually this represents insect mutation rather than invasion)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PLOT incubation curves at diff. snapshots %%%%%%%%%%
options = odeset('Events',@midScapeEvents);
[tnew,ynew,te3,ye,ie] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,icbnMain,dNmain,range,propCleanSeed,cutFromWithin),[0 tMax],initValsFromSS,options);
infVecF=ye(7*numFieldsIn+(1:numFieldsIn/2));
exposVecF=ye(6*numFieldsIn+(1:numFieldsIn/2));
expVecF=exposVecF/H;
incVecF=infVecF/H;
if sum(sum((incVecF+expVecF)>1.001))
    disp('ERROR!!! unusual incidence');
    for qq=1:numFieldsIn/2
        tmpfield=incVecF(qq)+expVecF(qq);
        inds=find(tmpfield>1.001);
        if ~isempty(inds)
            disp(['Field is ' num2str(qq)]);
            disp(['At time ' num2str(inds')]);
            disp(incVecF(qq)');
            disp(expVecF(qq)');
            return;
        end
    end
end
superAbundAF=H*(1-incVecF-expVecF).*(ye(3*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(4*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(5*numFieldsIn+(1:numFieldsIn/2)));%sp#1
superAbundBF=H*(1-incVecF-expVecF).*(ye(14*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(15*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(16*numFieldsIn+(1:numFieldsIn/2)));%sp#2
n_superAbundAF=H*(1-incVecF-expVecF).*(ye(0*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(1*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(2*numFieldsIn+(1:numFieldsIn/2)));%sp#1nymph
n_superAbundBF=H*(1-incVecF-expVecF).*(ye(11*numFieldsIn+(1:numFieldsIn/2)))+H*expVecF.*(ye(12*numFieldsIn+(1:numFieldsIn/2)))+H*incVecF.*(ye(13*numFieldsIn+(1:numFieldsIn/2)));%sp#2nymph
superAbundF=superAbundAF+superAbundBF;
n_superAbundF=n_superAbundAF+n_superAbundBF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xticmarks=0:10:100;
subplot(4,2,7)
yyaxis left
p1=plot(1:(numFieldsIn/2),(1/H)*superAbundF./(infVecF+denC),'k');
normTo1=max((1/H)*superAbundF./(infVecF+denC));
hold on;
normTo2=max((1/H)*n_superAbundF./(infVecF+denC));
p2=plot(1:(numFieldsIn/2),(normTo1/(normTo2*H))*n_superAbundF./(infVecF+denC),'g--');
set(gca,'ytick',0:20:200,'xtick',xticmarks);
ylabel('Wave-profile (Y/(I+1))')
xlabel('Field position')
lgd=legend([p1 p2],'adult','nymph','location','west');
lgd.FontSize = 9;
%ylim([0 80])
yyaxis right
p3=plot(1:(numFieldsIn/2),incVecF,'color',[0.68 0.85 0.9]);
set(gca,'ytick',[],'xtick',xticmarks);
xlim([20 55])
ylim([0 1])
ax3=gca;
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend boxoff  
ax3.YAxis(1).Color = 'k';


subplot(4,2,8)
yyaxis left
plot(1:(numFieldsIn/2),incVecF);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks);
ylim([0 1])
yyaxis right
plot(1:(numFieldsIn/2),superAbundF/H);
%plot(1:(numFieldsIn/2),superAbundF);
set(gca,'ytick',0:100:600,'xtick',xticmarks);
xlim([20 55])
ylim([0 350])
ax3b=gca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%% dummy plots
subplot(4,2,1)
plot([0,1],[0,1]);
ax0=gca;
subplot(4,2,2)
plot([0,1],[0,1]);
ax0b=gca;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'Units','centimeters', 'Position',[0 0 1.225*10.5 0.8*30.85])

set(gcf,'color', 'w');
set([ax0 ax0b ax1 ax2 ax3 ax1b ax2b ax3b], ...
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
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% orient(h,'portrait')
h.Renderer='Painters';
orient(h,'landscape')